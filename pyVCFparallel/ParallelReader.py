import multiprocessing as mp
import queue
import gzip
import codecs
import sys
import signal
import collections

# print(__name__+" => "+__file__)


class Reader(object):
    def __init__(self, fsock=None, filename=None, compressed=None, prepend_chr=False, strict_whitespace=False, encoding='ascii'):
        """ Create a new Reader for a VCF file.

            You must specify either fsock (stream) or filename.  Gzipped streams
            or files are attempted to be recogized by the file extension, or gzipped
            can be forced with ``compressed=True``

            'prepend_chr=True' will put 'chr' before all the CHROM values, useful
            for different sources.

            'strict_whitespace=True' will split records on tabs only (as with VCF
            spec) which allows you to parse files with spaces in the sample names.
        """
        from pyVCFparallel.ParallelWorkers import LineWorkerProcess, FileChunkLoaderProcess
        self.LineWorkerProcess = LineWorkerProcess
        self.FileLoaderProcess = FileChunkLoaderProcess
        self.running = False

        # Define the processing queues
        self.queue_in = None
        self.queue_out = None

        # Define the locks
        self.done_lock = mp.Lock()
        self.flush_lock = mp.Lock()
        self.abort_lock = mp.Lock()


        # super(Reader, self).__init__() <====== is this neeeded?

        self.ctxFile = {}

        if not (fsock or filename):
            raise Exception('You must provide at least fsock or filename')

        if fsock:
            self.ctxFile["reader"] = fsock
            if filename is None and hasattr(fsock, 'name'):
                filename = fsock.name
                if compressed is None:
                    compressed = filename.endswith('.gz')
        elif filename:
            if compressed is None:
                compressed = filename.endswith('.gz')
                self.ctxFile["reader"] = open(filename, 'rb' if compressed else 'rt')
        self.ctxFile["filename"] = filename
        self.ctxFile["compressed"] = False
        if compressed:
            self.ctxFile["reader"] = gzip.GzipFile(fileobj=self.ctxFile["reader"])
            self.ctxFile["compressed"] = True
            if sys.version > '3':
                self.ctxFile["reader"] = codecs.getreader(encoding)(self.ctxFile["reader"])

        # parse the header information if we are only passed a filepathname
        if filename:
            from pyVCFparallel.Parsers import Parsers
            self.ctxFile["headerInfo"], self._header_lines, self._column_headers = Parsers.ParseHeader(fsock=self.ctxFile["reader"], compressed=self.ctxFile["compressed"])


            # ___match PyVCF for compatibility____________________________
            self.metadata = None
            self.infos = None
            self.filters = None
            self.alts = None
            self.formats = None
            self.contigs = None
            self.samples = None


            # build the samples list being sure to keep the correct index order
            self.samples = []
            self._sample_indexes = {}
            for sampleID in self.ctxFile["headerInfo"]["SAMPLES"]:
                self._sample_indexes[sampleID] = self.ctxFile["headerInfo"]["SAMPLES"][sampleID]["index"]
                self.samples.append(sampleID)
            del self.ctxFile["headerInfo"]["SAMPLES"]


            setup_info = {"infos":{"file":"INFO", "obj":"Info"},
                          "filters":{"file": "FILTER", "obj": "Filter"},
                          "alts":{"file":"ALT", "obj":"Alt"},
                          "contigs":{"file":"contig", "obj":"Contig"},
                          "formats":{"file":"FORMAT", "obj":"Format"}
                          }

            for attr in setup_info:
                if setup_info[attr]["file"] in self.ctxFile["headerInfo"]:
                    if type(self.ctxFile["headerInfo"][setup_info[attr]["file"]]) == dict:
                        # extract a collection of attributes used by the records
                        attr_counts = collections.Counter()
                        # generate a list of all attributes used by this type of record
                        for rec_key in self.ctxFile["headerInfo"][setup_info[attr]["file"]]:
                            attr_counts.update(self.ctxFile["headerInfo"][setup_info[attr]["file"]][rec_key].keys())
                        generator = collections.namedtuple(setup_info[attr]["obj"], attr_counts.keys())
                        # create a defaults dict to fall back on
                        default_data = {}
                        for name in attr_counts.keys():
                            default_data[name] = None
                        # use generator to convert all the records
                        setup_info[attr]["data"] = collections.OrderedDict()
                        for rec_key in self.ctxFile["headerInfo"][setup_info[attr]["file"]]:
                            safe_data = collections.ChainMap(self.ctxFile["headerInfo"][setup_info[attr]["file"]][rec_key], default_data)
                            setup_info[attr]["data"][rec_key] = generator(**safe_data)
                        setattr(self, attr, setup_info[attr]["data"])
                        del setup_info[attr]["data"]
                        del self.ctxFile["headerInfo"][setup_info[attr]["file"]]
                    else:
                        pass # this is where we need to save this to metadata object

            # assume the remaining information is metadata
            metadata = collections.OrderedDict()
            source_list = []
            for rec_key in self.ctxFile["headerInfo"]:
                if rec_key == "source":
                    source_list.append(self.ctxFile["headerInfo"][rec_key])
                else:
                    metadata[rec_key] = self.ctxFile["headerInfo"][rec_key]
                metadata["source"] = source_list
            del self.ctxFile["headerInfo"]

            # _Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc', 'source', 'version'])
            # _Filter = collections.namedtuple('Filter', ['id', 'desc'])
            # _Alt = collections.namedtuple('Alt', ['id', 'desc'])
            # _Format = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])
            # _Contig = collections.namedtuple('Contig', ['id', 'length'])

            _SampleInfo = collections.namedtuple('SampleInfo', ['samples', 'gt_bases', 'gt_types', 'gt_phases'])


                # if strict_whitespace:
        #     self.ctxFile._separator = '\t'
        # else:
        #     self.ctxFile._separator = '\t| +'

#        self.ctxFile._row_pattern = re.compile(self.ctxFile._separator)
#        self.ctxFile._alt_pattern = re.compile('[\[\]]')

# TODO: THIS NEEDS TO CHANGE!
#        self.ctxFile.reader = (line.strip() for line in self.ctxFile._reader if line.strip())




    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        if self.queue_out == None:
            raise Exception('You must start execute ParallelReader.run(...) before you try to read from it!')
        try:
            record = self.queue_out.get(True, 3)
            return record
        except queue.Empty:
            if self.done_lock.acquire(False):
                self.done_lock.release()
                # the loader process has released the "working" lock
                # once it is released it means that no more items will
                # be added to queue_in

                print(self.queue_in.empty())
                print(self.queue_out.empty())
                aborting = False
                if self.abort_lock.acquire(False):
                    aborting = True
                    self.abort_lock.release()

                if (self.queue_in.empty() and self.queue_out.empty()) or aborting:
                    self.running = False
                    for p in self.line_processes:
                        if p.exitcode == None:
                            print("Shutting down worker...")
                            p.terminate()
                    raise StopIteration
        except InterruptedError:
            pass
        return None



    def _abort_sig_handler(self, signum, frame):
        print("abort!!!")
        self.abort_lock.release()
        for worker in self.line_processes:
            worker.terminate()
#        self.done_lock.release()


    def run(self, threadCnt, queueSize):
        # Define the processing queues
        self.queue_in = mp.Queue(queueSize)
        self.queue_out = mp.Queue(queueSize)

        self.running = True

        # handling for clean abort process
        self.abort_lock.acquire()
        signal.signal(signal.SIGINT, self._abort_sig_handler)
        signal.signal(signal.SIGTERM, self._abort_sig_handler)

        # Start filling the input queue
        # process_fill = mp.Process(target=populate_queue, args=(done_lock, queue_in))
        self.file_process = self.FileLoaderProcess(self.ctxFile["filename"], self.queue_in, self.done_lock, self.abort_lock)
        self.file_process.daemon = True
        self.file_process.start()

        # make sure it has control of the finish lock before starting the workers
        while self.done_lock.acquire(False):
            self.done_lock.release()

        # Setup a list of processing workers
        #	processes = [mp.Process(target=processing_worker, args=(5, x, queue_in, queue_out, done_lock), daemon=True) for x in range(4)]
        self.line_processes = [self.LineWorkerProcess(x, self.queue_in, self.queue_out, self.done_lock, self.ctxFile["headerInfo"]) for x in range(threadCnt)]
        # Run worker processes
        for p in self.line_processes:
            print("starting worker...")
            p.daemon = True
            p.start()

        # # Get the results from the output queue
        # finished = False
        # while not finished:
        #     try:
        #         result = queue_out.get(True, 3)
        #         print(result)
        #     except queue.Empty:
        #         print("Output Queue is EMPTY")
        #         if done_lock.acquire(False):
        #             done_lock.release()
        #             # the loader process has released the "working" lock
        #             # once it is released it means that no more items will
        #             # be added to queue_in
        #             if queue_in.empty() and queue_out.empty():
        #                 finished = True


# class Reader(object):
#     """ Reader for a VCF v 4.0 file, an iterator returning ``_Record objects`` """
#
#     def __init__(self, fsock=None, filename=None, compressed=None, prepend_chr=False,
#                  strict_whitespace=False, encoding='ascii'):
#         """ Create a new Reader for a VCF file.
#
#             You must specify either fsock (stream) or filename.  Gzipped streams
#             or files are attempted to be recogized by the file extension, or gzipped
#             can be forced with ``compressed=True``
#
#             'prepend_chr=True' will put 'chr' before all the CHROM values, useful
#             for different sources.
#
#             'strict_whitespace=True' will split records on tabs only (as with VCF
#             spec) which allows you to parse files with spaces in the sample names.
#         """
#         super(Reader, self).__init__()
#
#         if not (fsock or filename):
#             raise Exception('You must provide at least fsock or filename')
#
#         if fsock:
#             self._reader = fsock
#             if filename is None and hasattr(fsock, 'name'):
#                 filename = fsock.name
#                 if compressed is None:
#                     compressed = filename.endswith('.gz')
#         elif filename:
#             if compressed is None:
#                 compressed = filename.endswith('.gz')
#             self._reader = open(filename, 'rb' if compressed else 'rt')
#         self.filename = filename
#         if compressed:
#             self._reader = gzip.GzipFile(fileobj=self._reader)
#             if sys.version > '3':
#                 self._reader = codecs.getreader(encoding)(self._reader)
#
#         if strict_whitespace:
#             self._separator = '\t'
#         else:
#             self._separator = '\t| +'
#
#         self._row_pattern = re.compile(self._separator)
#         self._alt_pattern = re.compile('[\[\]]')
#
#         self.reader = (line.strip() for line in self._reader if line.strip())
#
#         #: metadata fields from header (string or hash, depending)
#         self.metadata = None
#         #: INFO fields from header
#         self.infos = None
#         #: FILTER fields from header
#         self.filters = None
#         #: ALT fields from header
#         self.alts = None
#         #: FORMAT fields from header
#         self.formats = None
#         #: contig fields from header
#         self.contigs = None
#         self.samples = None
#         self._sample_indexes = None
#         self._header_lines = []
#         self._column_headers = []
#         self._tabix = None
#         self._prepend_chr = prepend_chr
#         self._parse_metainfo()
#         self._format_cache = {}
#         self.encoding = encoding
#
#     def __iter__(self):
#         return self
#
#     def _parse_metainfo(self):
#         '''Parse the information stored in the metainfo of the VCF.
#
#         The end user shouldn't have to use this.  She can access the metainfo
#         directly with ``self.metadata``.'''
#         for attr in ('metadata', 'infos', 'filters', 'alts', 'contigs', 'formats'):
#             setattr(self, attr, OrderedDict())
#
#         parser = _vcf_metadata_parser()
#
#         line = next(self.reader)
#         while line.startswith('##'):
#             self._header_lines.append(line)
#
#             if line.startswith('##INFO'):
#                 key, val = parser.read_info(line)
#                 self.infos[key] = val
#
#             elif line.startswith('##FILTER'):
#                 key, val = parser.read_filter(line)
#                 self.filters[key] = val
#
#             elif line.startswith('##ALT'):
#                 key, val = parser.read_alt(line)
#                 self.alts[key] = val
#
#             elif line.startswith('##FORMAT'):
#                 key, val = parser.read_format(line)
#                 self.formats[key] = val
#
#             elif line.startswith('##contig'):
#                 key, val = parser.read_contig(line)
#                 self.contigs[key] = val
#
#             else:
#                 key, val = parser.read_meta(line)
#                 if key in SINGULAR_METADATA:
#                     self.metadata[key] = val
#                 else:
#                     if key not in self.metadata:
#                         self.metadata[key] = []
#                     self.metadata[key].append(val)
#
#             line = next(self.reader)
#
#         fields = self._row_pattern.split(line[1:])
#         self._column_headers = fields[:9]
#         self.samples = fields[9:]
#         self._sample_indexes = dict([(x,i) for (i,x) in enumerate(self.samples)])
#
#     def _map(self, func, iterable, bad=['.', '']):
#         '''``map``, but make bad values None.'''
#         return [func(x) if x not in bad else None
#                 for x in iterable]
#
#     def _parse_filter(self, filt_str):
#         '''Parse the FILTER field of a VCF entry into a Python list
#
#         NOTE: this method has a cython equivalent and care must be taken
#         to keep the two methods equivalent
#         '''
#         if filt_str == '.':
#             return None
#         elif filt_str == 'PASS':
#             return []
#         else:
#             return filt_str.split(';')
#
#     def _parse_info(self, info_str):
#         '''Parse the INFO field of a VCF entry into a dictionary of Python
#         types.
#
#         '''
#         if info_str == '.':
#             return {}
#
#         entries = info_str.split(';')
#         retdict = {}
#
#         for entry in entries:
#             entry = entry.split('=', 1)
#             ID = entry[0]
#             try:
#                 entry_type = self.infos[ID].type
#             except KeyError:
#                 try:
#                     entry_type = RESERVED_INFO[ID]
#                 except KeyError:
#                     if entry[1:]:
#                         entry_type = 'String'
#                     else:
#                         entry_type = 'Flag'
#
#             if entry_type == 'Integer':
#                 vals = entry[1].split(',')
#                 try:
#                     val = self._map(int, vals)
#                 # Allow specified integers to be flexibly parsed as floats.
#                 # Handles cases with incorrectly specified header types.
#                 except ValueError:
#                     val = self._map(float, vals)
#             elif entry_type == 'Float':
#                 vals = entry[1].split(',')
#                 val = self._map(float, vals)
#             elif entry_type == 'Flag':
#                 val = True
#             elif entry_type in ('String', 'Character'):
#                 try:
#                     vals = entry[1].split(',') # commas are reserved characters indicating multiple values
#                     val = self._map(str, vals)
#                 except IndexError:
#                     entry_type = 'Flag'
#                     val = True
#
#             try:
#                 if self.infos[ID].num == 1 and entry_type not in ( 'Flag', ):
#                     val = val[0]
#             except KeyError:
#                 pass
#
#             retdict[ID] = val
#
#         return retdict
#
#     def _parse_sample_format(self, samp_fmt):
#         """ Parse the format of the calls in this _Record """
#         samp_fmt = make_calldata_tuple(samp_fmt.split(':'))
#
#         for fmt in samp_fmt._fields:
#             try:
#                 entry_type = self.formats[fmt].type
#                 entry_num = self.formats[fmt].num
#             except KeyError:
#                 entry_num = None
#                 try:
#                     entry_type = RESERVED_FORMAT[fmt]
#                 except KeyError:
#                     entry_type = 'String'
#             samp_fmt._types.append(entry_type)
#             samp_fmt._nums.append(entry_num)
#         return samp_fmt
#
#     def _parse_samples(self, samples, samp_fmt, site):
#         '''Parse a sample entry according to the format specified in the FORMAT
#         column.
#
#         NOTE: this method has a cython equivalent and care must be taken
#         to keep the two methods equivalent
#         '''
#
#         # check whether we already know how to parse this format
#         if samp_fmt not in self._format_cache:
#             self._format_cache[samp_fmt] = self._parse_sample_format(samp_fmt)
#         samp_fmt = self._format_cache[samp_fmt]
#
#         if cparse:
#             return cparse.parse_samples(
#                 self.samples, samples, samp_fmt, samp_fmt._types, samp_fmt._nums, site)
#
#         samp_data = []
#         _map = self._map
#
#         nfields = len(samp_fmt._fields)
#
#         for name, sample in itertools.izip(self.samples, samples):
#
#             # parse the data for this sample
#             sampdat = [None] * nfields
#
#             for i, vals in enumerate(sample.split(':')):
#
#                 # short circuit the most common
#                 if samp_fmt._fields[i] == 'GT':
#                     sampdat[i] = vals
#                     continue
#                 # genotype filters are a special case
#                 elif samp_fmt._fields[i] == 'FT':
#                     sampdat[i] = self._parse_filter(vals)
#                     continue
#                 elif not vals or vals == ".":
#                     sampdat[i] = None
#                     continue
#
#                 entry_num = samp_fmt._nums[i]
#                 entry_type = samp_fmt._types[i]
#
#                 # we don't need to split single entries
#                 if entry_num == 1:
#                     if entry_type == 'Integer':
#                         try:
#                             sampdat[i] = int(vals)
#                         except ValueError:
#                             sampdat[i] = float(vals)
#                     elif entry_type == 'Float' or entry_type == 'Numeric':
#                         sampdat[i] = float(vals)
#                     else:
#                         sampdat[i] = vals
#                     continue
#
#                 vals = vals.split(',')
#                 if entry_type == 'Integer':
#                     try:
#                         sampdat[i] = _map(int, vals)
#                     except ValueError:
#                         sampdat[i] = _map(float, vals)
#                 elif entry_type == 'Float' or entry_type == 'Numeric':
#                     sampdat[i] = _map(float, vals)
#                 else:
#                     sampdat[i] = vals
#
#             # create a call object
#             call = _Call(site, name, samp_fmt(*sampdat))
#             samp_data.append(call)
#
#         return samp_data
#
#     def _parse_alt(self, str):
#         if self._alt_pattern.search(str) is not None:
#             # Paired breakend
#             items = self._alt_pattern.split(str)
#             remoteCoords = items[1].split(':')
#             chr = remoteCoords[0]
#             if chr[0] == '<':
#                 chr = chr[1:-1]
#                 withinMainAssembly = False
#             else:
#                 withinMainAssembly = True
#             pos = remoteCoords[1]
#             orientation = (str[0] == '[' or str[0] == ']')
#             remoteOrientation = (re.search('\[', str) is not None)
#             if orientation:
#                 connectingSequence = items[2]
#             else:
#                 connectingSequence = items[0]
#             return _Breakend(chr, pos, orientation, remoteOrientation, connectingSequence, withinMainAssembly)
#         elif str[0] == '.' and len(str) > 1:
#             return _SingleBreakend(True, str[1:])
#         elif str[-1] == '.' and len(str) > 1:
#             return _SingleBreakend(False, str[:-1])
#         elif str[0] == "<" and str[-1] == ">":
#             return _SV(str[1:-1])
#         else:
#             return _Substitution(str)
#
#     def next(self):
#         '''Return the next record in the file.'''
#         line = next(self.reader)
#         row = self._row_pattern.split(line.rstrip())
#         chrom = row[0]
#         if self._prepend_chr:
#             chrom = 'chr' + chrom
#         pos = int(row[1])
#
#         if row[2] != '.':
#             ID = row[2]
#         else:
#             ID = None
#
#         ref = row[3]
#         alt = self._map(self._parse_alt, row[4].split(','))
#
#         try:
#             qual = int(row[5])
#         except ValueError:
#             try:
#                 qual = float(row[5])
#             except ValueError:
#                 qual = None
#
#         filt = self._parse_filter(row[6])
#         info = self._parse_info(row[7])
#
#         try:
#             fmt = row[8]
#         except IndexError:
#             fmt = None
#         else:
#             if fmt == '.':
#                 fmt = None
#
#         record = _Record(chrom, pos, ID, ref, alt, qual, filt,
#                 info, fmt, self._sample_indexes)
#
#         if fmt is not None:
#             samples = self._parse_samples(row[9:], fmt, record)
#             record.samples = samples
#
#         return record