import multiprocessing as mp
import queue
import codecs
import gzip
import sys
import binascii
import signal

# print(__name__+" => "+__file__)

class LineWorkerProcess(mp.Process):
	def __init__(self, id, queue_in, queue_out, lock_finished):
		super(LineWorkerProcess, self).__init__()
		self.id = id
		self.q_in = queue_in
		self.q_out = queue_out
		self.l_finished = lock_finished

	def _abort_sig_handler(self, signum, frame):
#		self.done = True
		pass

	def run(self):
		signal.signal(signal.SIGINT, self._abort_sig_handler)
		signal.signal(signal.SIGTERM, self._abort_sig_handler)

		self.done = False
		while not self.done:
			try:
				work_line = self.q_in.get(False)
				# split the line into major elements
				work_elements = work_line.split()

#					print(work_elements)

				# From https://github.com/jamescasbon/PyVCF/blob/master/vcf/parser.py#L555
				chrom = work_elements[0]
				pos = int(work_elements[1])
				if (work_elements[2]) != '.':
					IDs = work_elements[2].split(",")
				else:
					IDs = []
				ref = work_elements[3]

				alt = {}
				for temp in work_elements[4].split(','):
					alt[temp] = {}


				try:
					qual = int(work_elements[5])
				except ValueError:
					try:
						qual = float(work_elements[5])
					except ValueError:
						qual = None


				# TODO: DO WORK HERE
#					ignored_var = math.factorial(random.randint(9999,99999))
				self.q_out.put(binascii.crc32(work_line.encode('ascii')))
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