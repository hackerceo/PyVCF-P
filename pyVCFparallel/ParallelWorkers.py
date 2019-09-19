import multiprocessing as mp
import queue
import math
import random
import codecs
import gzip
import sys
import binascii

print(__name__+" => "+__file__)

class LineWorkerProcess(mp.Process):
	def __init__(self, id, queue_in, queue_out, lock_finished):
		super(LineWorkerProcess, self).__init__()
		self.id = id
		self.q_in = queue_in
		self.q_out = queue_out
		self.l_finished = lock_finished

	def run(self):
			finished = False
			while not finished:
				try:
					work_line = self.q_in.get(False)
					# TODO: DO WORK HERE
#					ignored_var = math.factorial(random.randint(9999,99999))
					self.q_out.put(binascii.crc32(work_line.encode('ascii')))
				except queue.Empty:
					if self.l_finished.acquire(False):
						self.l_finished.release()
						print("worker "+str(self.id)+" FINISHED")
						finished = True

class FileLoaderProcess(mp.Process):
	def __init__(self, filename, queue_in, lock_finished):
		super(FileLoaderProcess, self).__init__()
		self.q_in = queue_in
		self.l_finished = lock_finished
		self.filename = filename
	def run(self):
		# open file
		fsock = open(self.filename, 'rt')
		if self.filename.endswith('.gz'):
			fsock = gzip.GzipFile(fileobj=fsock)
			if sys.version > '3':
				fsock = codecs.getreader('ascii')(fsock)

		# stream lines into the processing queue
		print("populate start")
		self.l_finished.acquire()
		print("populate locked")
		while True:
			line = fsock.readline()
			if not line:
				break;
			self.q_in.put(line.strip())
		self.q_in.close()
		print("populate unlocked")
		self.l_finished.release()