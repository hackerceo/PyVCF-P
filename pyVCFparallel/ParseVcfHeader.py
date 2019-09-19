import codecs
import gzip
import re
import sys


class Parsers(object):
    '''Parse the metadata in the header of a VCF file.'''
    @staticmethod
    def ParseHeader(fsock=None, filename=None, compressed=None, encoding='ascii'):
        """ Create a new Reader for a VCF file.

            You must specify either fsock (stream) or filename.  Gzipped streams
            or files are attempted to be recogized by the file extension, or gzipped
            can be forced with ``compressed=True``
        """

        if not (fsock or filename):
            raise Exception('You must provide at least fsock or filename')

        if fsock:
            _reader = fsock
            if filename is None and hasattr(fsock, 'name'):
                filename = fsock.name
                if compressed is None:
                    compressed = filename.endswith('.gz')

        elif filename:
            if compressed is None:
                compressed = filename.endswith('.gz')
            _reader = open(filename, 'rb' if compressed else 'rt')
        filename = filename
        if compressed:
            _reader = gzip.GzipFile(fileobj=_reader)
            if sys.version > '3':
                _reader = codecs.getreader(encoding)(_reader)


        # rewind and read from the begining of the file.
        _reader.seek(0)
        file_header_info = {}
        _line = _reader.readline().strip()
        while (_line.startswith("#")):
            # header data starts with "#"
            if _line.startswith("##"):
                if re.match(".*<.*>", _line):
                    #--------------------------------------- Line type A
                    items = _line.split('=', 1)
                    # Removing initial hash marks
                    key = items[0].lstrip('#')
                    # N.B., items can have quoted values, so cannot just split on comma
                    val = {}
                    state = 0
                    k = ''
                    v = ''
                    for c in items[1].strip('[<>]'):

                        if state == 0:  # reading item key
                            if c == '=':
                                state = 1  # end of key, start reading value
                            else:
                                k += c  # extend key
                        elif state == 1:  # reading item value
                            if v == '' and c == '"':
                                v += c  # include quote mark in value
                                state = 2  # start reading quoted value
                            elif c == ',':
                                val[k] = v  # store parsed item
                                state = 0  # read next key
                                k = ''
                                v = ''
                            else:
                                v += c
                        elif state == 2:  # reading quoted item value
                            if c == '"':
                                v += c  # include quote mark in value
                                state = 1  # end quoting
                            else:
                                v += c
                    if k != '':
                        val[k] = v
                    #---------------------------------------
                    # save the data to return data structure
                    if key not in file_header_info:
                        file_header_info[key] = {}
                    file_header_info[key][val["ID"]] = val
                else:
                    #--------------------------------------- Line type B
                    items = _line.split('=', 1)
                    # Removing initial hash marks
                    key = items[0].lstrip('#')
                    # N.B., items can have quoted values, so cannot just split on comma
                    val = {}
                    state = 0
                    k = ''
                    v = ''
                    for c in items[1].strip('[<>]'):

                        if state == 0:  # reading item key
                            if c == '=':
                                state = 1  # end of key, start reading value
                            else:
                                k += c  # extend key
                        elif state == 1:  # reading item value
                            if v == '' and c == '"':
                                v += c  # include quote mark in value
                                state = 2  # start reading quoted value
                            elif c == ',':
                                val[k] = v  # store parsed item
                                state = 0  # read next key
                                k = ''
                                v = ''
                            else:
                                v += c
                        elif state == 2:  # reading quoted item value
                            if c == '"':
                                v += c  # include quote mark in value
                                state = 1  # end quoting
                            else:
                                v += c
                    if k != '':
                        val[k] = v
                    #---------------------------------------
                    # save the data to return data structure
                    if key not in file_header_info:
                        file_header_info[key] = []
                    for k, v in val.items():
                        file_header_info[key].append(k)
            else:
                # the final header line provides a list of sample identifiers
                file_header_info["SAMPLES"] = {}
                samples = _line.split()[9::]
                for id in samples:
                    file_header_info["SAMPLES"][id] = {}

            # get next line
            _line = _reader.readline().strip()

        # finished parsing the header information
        return file_header_info

