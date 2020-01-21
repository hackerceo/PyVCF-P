import multiprocessing as mp
import queue
import codecs
import gzip
import sys
import binascii
import signal
import zlib


# print(__name__+" => "+__file__)

class LineWorkerProcess(mp.Process):
	def __init__(self, id, queue_in, queue_out, lock_finished, header_info):
		super(LineWorkerProcess, self).__init__()
		self.id = id
		self.q_in = queue_in
		self.q_out = queue_out
		self.l_finished = lock_finished
		self.header_info = header_info
		from .model import _Call, _Record, make_calldata_tuple
		from .model import _Substitution, _Breakend, _SingleBreakend, _SV

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

					# SITE's POSITION
					parsed_line["pos"] = int(work_elements[1])

					# SITE's IDs
					if (work_elements[2]) != '.':
						IDs = work_elements[2].split(",")
					else:
						IDs = []
					parsed_line["id"] = IDs

					#SITE's REFERENCE ALLELES
					parsed_line["ref"] = work_elements[3].split(',')

					# SITE'S ALT ALLELES
					alt = {}
					idx_pos = 0
					alt_index = list()
					for temp in work_elements[4].split(','):
						alt[temp] = {'index': idx_pos}
						alt_index.append(temp)
						idx_pos = idx_pos + 1
					parsed_line["alt"] = alt

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
					parsed_line["filter"] = work_elements[6].split(';')

					# SITE'S INFO (entered under the corresponding "alt" entries, or ref+alt entries)
					cnt_A = len(parsed_line["alt"])
					cnt_R = len(parsed_line["ref"]) + cnt_A
					cnt_G = len(self.header_info["SAMPLES"])
					parsed_line["info"] = {}
					for rec in work_elements[7].split(';'):
						temp = rec.split('=')
						id = temp[0]
						if len(temp) == 1:
							parsed_line["info"][id] = True
						else:
							temp = temp[1].split(',')
							if len(temp) == 1:
								parsed_line["info"][id] = temp[0]
							else:
								if len(temp) == cnt_A:
									for idx in range(len(temp)):
										parsed_line["alt"][alt_index[idx]][id] = temp[idx]
								elif len(temp) == cnt_R:
									splt_idx = len(parsed_line["ref"])
									# for idx in range(splt_idx):
									# 	parsed_line["ref"][alt_index[idx]][id] = temp[idx]
                                    #
									# parsed_line["ref"][alt_index[idx]][id] = temp[idx]
									# for idx in range(len(temp)):
									# 	parsed_line["alt"][alt_index[idx]][id] = temp[idx]

					# FORMAT FOR EACH SAMPLE OF THE SITE
					if "FORMAT" in self.header_info:
						try:
							fmt = work_elements[8]
						except IndexError:
							fmt = None
						else:
							if fmt == '.':
								fmt = None

						if fmt is not None:
							parsed_line["format"] = work_elements[8]
							format_spec = work_elements[8].split(':')

					# PROCESS EACH SAMPLE FOR THE CURRENT SITE
					for sample_idx in range(len(work_elements) - 9):
						# TODO: DO MORE WORK HERE
						sample_data = work_elements[sample_idx + 9].split(':')
						gt_record = {"phased": False, "alleles": []}
						for sample_format_idx in range(len(sample_data)):
							if format_spec[sample_format_idx] == "GT":
								pass

					# push the processed line into our output buffer
					self.q_out.put(parsed_line)

			except queue.Empty:
				if self.l_finished.acquire(False):
					self.l_finished.release()
					print("worker "+str(self.id)+" FINISHED")
					self.done = True






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