import multiprocessing as mp
import queue
import codecs
import gzip
import sys
import binascii
import signal
import zlib
import re
import pickle

from .model import _Call, _Record, make_calldata_tuple
from .model import _Substitution, _Breakend, _SingleBreakend, _SV


# print(__name__+" => "+__file__)

class LineWorkerProcess(mp.Process):
	def __init__(self, id, queue_in, queue_out, lock_finished, file_info, options):
		super(LineWorkerProcess, self).__init__()

		self.options = options
		self.header_info = file_info["headerInfo"]


		# Conversion between value in file and Python value
		self.field_counts = {
			'.': None,  # Unknown number of values
			'A': -1,  # Equal to the number of alternate alleles in a given record
			'G': -2,  # Equal to the number of genotypes in a given record
			'R': -3,  # Equal to the number of alleles including reference in a given record
		}
		self.id = id
		self.q_in = queue_in
		self.q_out = queue_out
		self.l_finished = lock_finished
		self._row_pattern = options["row_pattern"]
		self._format_cache = {}
		self._sample_indexes = options["_sample_indexes"]
		self.samples = options["samples"]
		self.formats = file_info["headerInfo"]["header_data"]["FORMAT"]

	def _abort_sig_handler(self, signum, frame):
#		self.done = True
		pass

	def run(self):
		signal.signal(signal.SIGINT, self._abort_sig_handler)
		signal.signal(signal.SIGTERM, self._abort_sig_handler)

		self.done = False
		while not self.done:
			try:
				# get a line from the input queue
				work_buffer = self.q_in.get(False)
#				work_buffer = zlib.decompress(bytes(work_buffer)).decode().split("\n")
				work_buffer = work_buffer.split("\n")
				for work_line in work_buffer:
					work_line = work_line.strip()

					# split the line into major elements
					work_elements = work_line.split()

	#					print(work_elements)

					# From https://github.com/jamescasbon/PyVCF/blob/master/vcf/parser.py#L555
					parsed_line = {}

					# SITE's CHROMOSOME
					parsed_line["chrom"] = work_elements[0]
					if self.options["prepend_chr"] == True:
						parsed_line["chrom"] = "chr" + parsed_line["chrom"]

					# SITE's POSITION
					parsed_line["pos"] = int(work_elements[1])

					# SITE's IDs
					if (work_elements[2]) != '.':
						IDs = work_elements[2].split(",")
					else:
						IDs = None
					parsed_line["id"] = IDs

					#SITE's REFERENCE ALLELES
					parsed_line["ref"] = work_elements[3]

					#SITE's Alternate Alleles
					parsed_line["alt"] = self._map(self._parse_alt, work_elements[4].split(','))

					# SITE'S QUAL
					try:
						qual = int(work_elements[5])
					except ValueError:
						try:
							qual = float(work_elements[5])
						except ValueError:
							qual = None
					parsed_line["qual"] = qual

					# SITE'S FILTER
					filt = work_elements[6]
					if filt == '.':
						filt = None
					elif filt == 'PASS':
						filt = []
					else:
						filt = filt.split(';')
					parsed_line["filter"] = filt

					parsed_line["info"] = self._parse_info(work_elements[7])

					try:
						fmt = work_elements[8]
					except IndexError:
						fmt = None
					else:
						if fmt == '.':
							fmt = None
					parsed_line["format"] = fmt

					record = _Record(parsed_line["chrom"], parsed_line["pos"], parsed_line["id"], parsed_line["ref"], parsed_line["alt"], parsed_line["qual"], parsed_line["filter"],
									 parsed_line["info"], parsed_line["format"], self._sample_indexes)

					if fmt is not None:
						samples = self._parse_samples(work_elements[9:], parsed_line["format"], record)
						record.samples = samples

                    #
					# # SITE'S INFO (entered under the corresponding "alt" entries, or ref+alt entries)
					# cnt_A = len(parsed_line["alt"])
					# cnt_R = len(parsed_line["ref"]) + cnt_A
					# cnt_G = len(self.header_info["SAMPLES"])
					# parsed_line["info"] = {}
					# for rec in work_elements[7].split(';'):
					# 	temp = rec.split('=')
					# 	id = temp[0]
					# 	if len(temp) == 1:
					# 		parsed_line["info"][id] = True
					# 	else:
					# 		temp = temp[1].split(',')
					# 		if len(temp) == 1:
					# 			parsed_line["info"][id] = temp[0]
					# 		else:
					# 			if len(temp) == cnt_A:
					# 				for idx in range(len(temp)):
					# 					parsed_line["alt"][alt_index[idx]][id] = temp[idx]
					# 			elif len(temp) == cnt_R:
					# 				splt_idx = len(parsed_line["ref"])
					# 				# for idx in range(splt_idx):
					# 				# 	parsed_line["ref"][alt_index[idx]][id] = temp[idx]
                     #                #
					# 				# parsed_line["ref"][alt_index[idx]][id] = temp[idx]
					# 				# for idx in range(len(temp)):
					# 				# 	parsed_line["alt"][alt_index[idx]][id] = temp[idx]
                    #
					# # FORMAT FOR EACH SAMPLE OF THE SITE
					# if "FORMAT" in self.header_info:
					# 	try:
					# 		fmt = work_elements[8]
					# 	except IndexError:
					# 		fmt = None
					# 	else:
					# 		if fmt == '.':
					# 			fmt = None
                    #
					# 	if fmt is not None:
					# 		parsed_line["format"] = work_elements[8]
					# 		format_spec = work_elements[8].split(':')
                    #
					# # PROCESS EACH SAMPLE FOR THE CURRENT SITE
					# for sample_idx in range(len(work_elements) - 9):
					# 	# TODO: DO MORE WORK HERE
					# 	sample_data = work_elements[sample_idx + 9].split(':')
					# 	gt_record = {"phased": False, "alleles": []}
					# 	for sample_format_idx in range(len(sample_data)):
					# 		if format_spec[sample_format_idx] == "GT":
					# 			pass

					# push the processed line into our output buffer
					self.q_out.put(pickle.dumps(record))
#					self.q_out.put("test")

			except queue.Empty:
				if self.l_finished.acquire(False):
					self.l_finished.release()
					print("worker "+str(self.id)+" FINISHED")
					self.done = True

	def _map(self, func, iterable, bad='.'):
		'''``map``, but make bad values None.'''
		return [func(x) if x != bad else None
				for x in iterable]

	def _parse_alt(self, str):
		if self.options['alt_pattern'].search(str) is not None:
			# Paired breakend
			items = self.options['alt_pattern'].split(str)
			remoteCoords = items[1].split(':')
			chr = remoteCoords[0]
			if chr[0] == '<':
				chr = chr[1:-1]
				withinMainAssembly = False
			else:
				withinMainAssembly = True
			pos = remoteCoords[1]
			orientation = (str[0] == '[' or str[0] == ']')
			remoteOrientation = (re.search('\[', str) is not None)
			if orientation:
				connectingSequence = items[2]
			else:
				connectingSequence = items[0]
			return _Breakend(chr, pos, orientation, remoteOrientation, connectingSequence, withinMainAssembly)
		elif str[0] == '.' and len(str) > 1:
			return _SingleBreakend(True, str[1:])
		elif str[-1] == '.' and len(str) > 1:
			return _SingleBreakend(False, str[:-1])
		elif str[0] == "<" and str[-1] == ">":
			return _SV(str[1:-1])
		else:
			return _Substitution(str)


	def _parse_info(self, info_str):
		'''Parse the INFO field of a VCF entry into a dictionary of Python
		types.

		'''
		if info_str == '.':
			return {}

		entries = info_str.split(';')
		retdict = {}

		for entry in entries:
			entry = entry.split('=', 1)
			ID = entry[0]
			try:
				entry_type = self.header_info["header_data"]["INFO"][ID]["Type"]
#				entry_type = self.infos[ID].type
			except KeyError:
				if entry[1:]:
					entry_type = 'String'
				else:
					entry_type = 'Flag'

			if entry_type == 'Integer':
				vals = entry[1].split(',')
				try:
					val = self._map(int, vals)
				# Allow specified integers to be flexibly parsed as floats.
				# Handles cases with incorrectly specified header types.
				except ValueError:
					val = self._map(float, vals)
			elif entry_type == 'Float':
				vals = entry[1].split(',')
				val = self._map(float, vals)
			elif entry_type == 'Flag':
				val = True
			elif entry_type in ('String', 'Character'):
				try:
					vals = entry[1].split(',')  # commas are reserved characters indicating multiple values
					val = self._map(str, vals)
				except IndexError:
					entry_type = 'Flag'
					val = True

			try:
				if self.header_info["header_data"]["INFO"][ID]["Number"] == 1 and entry_type not in ('Flag',):
#				if self.infos[ID].num == 1 and entry_type not in ('Flag',):
					val = val[0]
			except KeyError:
				pass

			retdict[ID] = val

		return retdict

	def _parse_samples(self, samples, samp_fmt, site):
		'''Parse a sample entry according to the format specified in the FORMAT
        column.

        NOTE: this method has a cython equivalent and care must be taken
        to keep the two methods equivalent
        '''

		# check whether we already know how to parse this format
		if samp_fmt not in self._format_cache:
			self._format_cache[samp_fmt] = self._parse_sample_format(samp_fmt)
		samp_fmt = self._format_cache[samp_fmt]

		samp_data = []
		_map = self._map

		nfields = len(samp_fmt._fields)

		for name, sample in zip(self.samples, samples):

			# parse the data for this sample
			sampdat = [None] * nfields

			for i, vals in enumerate(sample.split(':')):

				# short circuit the most common
				if samp_fmt._fields[i] == 'GT':
					sampdat[i] = vals
					continue
				elif not vals or vals == ".":
					sampdat[i] = None
					continue

				entry_num = samp_fmt._nums[i]
				entry_type = samp_fmt._types[i]

				# we don't need to split single entries
				if entry_num == 1 or ',' not in vals:

					if entry_type == 'Integer':
						try:
							sampdat[i] = int(vals)
						except ValueError:
							sampdat[i] = float(vals)
					elif entry_type == 'Float':
						sampdat[i] = float(vals)
					else:
						sampdat[i] = vals

					if entry_num != 1:
						sampdat[i] = (sampdat[i])

					continue

				vals = vals.split(',')

				if entry_type == 'Integer':
					try:
						sampdat[i] = _map(int, vals)
					except ValueError:
						sampdat[i] = _map(float, vals)
				elif entry_type == 'Float' or entry_type == 'Numeric':
					sampdat[i] = _map(float, vals)
				else:
					sampdat[i] = vals

			# create a call object
			call = _Call(site, name, samp_fmt(*sampdat))
			samp_data.append(call)

		return samp_data

	def _parse_sample_format(self, samp_fmt):
		""" Parse the format of the calls in this _Record """
		samp_fmt = make_calldata_tuple(samp_fmt.split(':'))

		for fmt in samp_fmt._fields:
			try:
				entry_type = self.formats[fmt]["Type"]
				entry_num = self.formats[fmt]["Number"]
			except KeyError:
				entry_num = None
				entry_type = 'String'
			samp_fmt._types.append(entry_type)
			samp_fmt._nums.append(entry_num)
		return samp_fmt


class FileLoaderProcess(mp.Process):
	def __init__(self, filename, queue_in, lock_finished, lock_abort):
		super(FileLoaderProcess, self).__init__()
		self.q_in = queue_in
		self.l_finished = lock_finished
		self.l_abort = lock_abort
		self.filename = filename

	def _abort_sig_handler(self, signum, frame):
		self.done = True

	def run(self):
		signal.signal(signal.SIGINT, self._abort_sig_handler)
		signal.signal(signal.SIGTERM, self._abort_sig_handler)

		# open file
		fsock = open(self.filename, 'rt')
		if self.filename.endswith('.gz'):
			fsock = gzip.GzipFile(fileobj=fsock)
			if sys.version > '3':
				fsock = codecs.getreader('ascii')(fsock)

		# stream lines into the processing queue
		print("populate start")
		self.done = False
		self.l_finished.acquire()
		print("populate locked")
		while not self.done:
			if self.l_abort.acquire(False):
				print("FileReader Abort")
				self.l_abort.release()
				self.done = True
			else:
				# read line
				line = fsock.readline()
				if not line:
					self.done = True;
				else:
					line = line.strip()
					if not line.startswith("#"):
						self.q_in.put(line)
		self.q_in.close()
		print("populate unlocked")
		self.l_finished.release()

class FileChunkLoaderProcess(mp.Process):
	def __init__(self, filename, queue_in, lock_finished, lock_abort):
		super(FileChunkLoaderProcess, self).__init__()
		self.q_in = queue_in
		self.l_finished = lock_finished
		self.l_abort = lock_abort
		self.filename = filename

	def _abort_sig_handler(self, signum, frame):
		self.l_finished.release()
		self.done = True

	def run(self, chunk_size = 32768 ): # 65536 or 32768
		signal.signal(signal.SIGINT, self._abort_sig_handler)
		signal.signal(signal.SIGTERM, self._abort_sig_handler)

		# open file
		fsock = open(self.filename, 'rt')
		if self.filename.endswith('.gz'):
			fsock = gzip.GzipFile(fileobj=fsock)
			if sys.version > '3':
				fsock = codecs.getreader('ascii')(fsock)

		# process past all the VCF header lines
		print("priming past VCF header")
		line = fsock.readline().strip()
		while line.startswith("#"):
			line = fsock.readline().strip()

		# stream lines into the processing queue
		print("populate start")
		self.done = False
		self.chunk_sz = chunk_size
		self.l_finished.acquire()
		chunk_buf = ""
		chunk_remains = ""
		print("populate locked")
		while not self.done:
			if self.l_abort.acquire(False):
				print("FileReader Abort")
				self.l_abort.release()
				self.done = True
			else:
				# read chunk
				chunk_buf = fsock.read(self.chunk_sz)

				# see if we are at the end of the file
				if len(chunk_buf) == 0:
					self.done = True
				else:
					# find the end of the last whole line captured
					chunk_end = chunk_buf.rindex("\n")

					# assemble the chunk to send and queue it up
					chunk_send = ''.join([chunk_remains, chunk_buf[:chunk_end]])
					self.q_in.put(chunk_send)
#					chunk_send = bytearray(chunk_remains.decode()).append(chunk_buf[:chunk_end].decode())
#					self.q_in.put(zlib.compress(chunk_send.encode(), level=1))

					# extract the chunk's trailing chars
					chunk_remains = chunk_buf[chunk_end+1:]

		self.q_in.close()
		print("populate unlocked")
		self.l_finished.release()