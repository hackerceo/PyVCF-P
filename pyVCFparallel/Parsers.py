import codecs
import gzip
import re
import sys


class Parsers(object):
    @staticmethod
    def getInfoDistributions(self):
        return {
            "AA": "1",
            "AC":"A",
            "AD":"R",
            "ADF":"R",
            "ADR":"R",
            "AF":"A",
            "AN":"1",
            "BQ":"1",
            "CIGAR":"A",
            "DB":"0",
            "DP":"1",
            "END":"1",
            "H2":"0",
            "H3":"0",
            "MQ":"1",
            "MQ0":"1",
            "NS":"1",
            "SB":"4",
            "SOMATIC":"0",
            "VALIDATED":"0",
            "1000G":"0"
        }

    @staticmethod
    def ParseHeader(fsock=None, filename=None, compressed=None, encoding='ascii', options=None):
        '''Parse the metadata in the header of a VCF file.'''
        """ 
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
        return_data = {}
        return_data["header_lines"] = []
        file_header_info = {}
        _line = _reader.readline().strip()
        while (_line.startswith("#")):
            # header data starts with "#"
            if _line.startswith("##"):
                return_data["header_lines"].append(_line)
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
                # the final header line provides a list of sample identifiers and field identifiers
                fields = options["row_pattern"].split(_line[1:])
                fields_cnt = 0
                fields_map = {"CHROM":None, "POS":None, "ID":None, "REF":None, "ALT":None, "QUAL":None, "FILTER":None, "INFO":None, "FORMAT":None, "END":None}
                for fieldname in fields:
                    if fieldname.upper() in fields_map:
                        fields_map[fieldname.upper()] = fields_cnt
                        fields_cnt += 1
                return_data["column_names"] = fields[0:fields_cnt]
                return_data["column_index"] = fields_map

                # extract the samples list
                return_data["samples"] = fields[fields_cnt:]
                return_data["sample_indexes"] = dict([(x, i) for (i, x) in enumerate(return_data["samples"])])

            # get next line
            _line = _reader.readline().strip()

        # finished parsing the header information (FILL IN DEFAULTS AND MISSING COMMON DEFINITIONS)
        return_data["header_data"] = Parsers._inject_common_header_data(file_header_info)

        return return_data

    @staticmethod
    def _parse_sample_format(self, samp_fmt):
        """ Parse the format of the calls in this _Record """
        samp_fmt = make_calldata_tuple(samp_fmt.split(':'))

        for fmt in samp_fmt._fields:
            try:
                entry_type = self.formats[fmt].type
                entry_num = self.formats[fmt].num
            except KeyError:
                entry_num = None
                try:
                    entry_type = RESERVED_FORMAT[fmt]
                except KeyError:
                    entry_type = 'String'
            samp_fmt._types.append(entry_type)
            samp_fmt._nums.append(entry_num)
        return samp_fmt

    @staticmethod
    def _inject_common_header_data(header_tag_array):
        # ========[ FORMAT ELEMENTS ]===================================================================================
        if not "FORMAT" in header_tag_array:
            header_tag_array["FORMAT"] = {}
        defaults = {
            "AD": {"Number":"R","Type":"Integer", "Description":"Read depth for each allele"},
            "ADF": {"Number":"R","Type":"Integer", "Description":"Read depth for each allele on the forward strand"},
            "ADR": {"Number":"R","Type":"Integer", "Description":"Read depth for each allele on the reverse strand"},
            "AHAP": {"Number":"1","Type":"Integer", "Description":"Unique identifier of ancestral haplotype"},
            "CN": {"Number": "1", "Type": "Integer", "Description": "Copy number genotype for imprecise events"},
            "CNQ": {"Number": "1", "Type": "Float", "Description": "Copy number genotype quality for imprecise events"},
            "CNL": {"Number": "G", "Type": "Float", "Description": "Copy number genotype likelihood for imprecise events"},
            "CNP": {"Number": "G", "Type": "Float", "Description": "Copy number posterior probabilities"},
            "DP": {"Number": "1", "Type": "Integer", "Description": "Read depth"},
            "EC": {"Number": "A", "Type": "Integer", "Description": "Expected alternate allele counts"},
            "FT": {"Number": "1", "Type": "String", "Description": "Filter indicating if this genotype was \"\called\""},
            "GL": {"Number": "G", "Type": "Float", "Description": "Genotype likelihoods"},
            "GP": {"Number": "G", "Type": "Float", "Description": "Genotype posterior probabilities"},
            "GQ": {"Number":"1","Type":"Integer", "Description":"Genotype quality"},
            "GT": {"Number":"1","Type":"String", "Description":"Genotype"},
            "HAP": {"Number": "1", "Type": "Integer", "Description": "Unique haplotype identifier"},
            "HQ": {"Number": "2", "Type": "Integer", "Description": "Haplotype quality"},
            "MQ": {"Number": "1", "Type": "Integer", "Description": "RMS mapping quality"},
            "NQ": {"Number": "1", "Type": "Integer", "Description": "Phred style probability score that the variant is novel"},
            "PL": {"Number": "G", "Type": "Integer", "Description": "Phred-scaled genotype likelihoods rounded to the closest integer"},
            "PQ": {"Number": "1", "Type": "Integer", "Description": "Phasing quality"},
            "PS": {"Number": "1", "Type": "Integer", "Description": "Phase set"},
        }
        for key in header_tag_array["FORMAT"]:
            if key in defaults:
                header_tag_array["FORMAT"][key] = {**defaults[key], **header_tag_array["FORMAT"][key]}
        for key in defaults:
            if key not in header_tag_array["FORMAT"]:
                header_tag_array["FORMAT"][key] = defaults[key]

        # ========[ INFO ELEMENTS ]===================================================================================
        if not "INFO" in header_tag_array:
            header_tag_array["INFO"] = {}
        defaults = {
            "AA": {"Number":"1", "Type":"String", "Description":"Ancestral allele"},
            "AC": {"Number":"A", "Type":"Integer", "Description":"Allele count in genotypes, for each ALT allele, in the same order as listed"},
            "AD": {"Number":"R", "Type":"Integer", "Description":"Total read depth for each allele"},
            "ADF": {"Number":"R", "Type":"Integer", "Description":"Read depth for each allele on the forward strand"},
            "ADR": {"Number":"R", "Type":"Integer", "Description":"Read depth for each allele on the reverse strand"},
            "AF": {"Number":"A", "Type":"Float", "Description":"Allele frequency for each ALT allele in the same order as listed (estimated from primary data, not called genotypes)"},
            "AN": {"Number":"1", "Type":"Integer", "Description":"Total number of alleles in called genotypes"},
            "BQ": {"Number":"1", "Type":"Float", "Description":"RMS base quality"},
            "CIGAR": {"Number":"A", "Type":"String", "Description":"Cigar string describing how to align an alternate allele to the reference allele"},
            "DB": {"Number":"0", "Type":"Flag", "Description":"dbSNP membership"},
            "DP": {"Number":"1", "Type":"Integer", "Description":"Combined depth across samples"},
            "END": {"Number":"1", "Type":"Integer", "Description":"End position on CHROM(used with symbolic alleles; see below)"},
            "H2:": {"Number":"0", "Type":"Flag", "Description":"HapMap2 membership"},
            "H3:": {"Number":"0", "Type":"Flag", "Description":"HapMap3 membership"},
            "MQ:": {"Number":"1", "Type":"Float", "Description":"RMS mapping quality"},
            "MQ0:": {"Number":"1", "Type":"Integer", "Description":"Number of MAPQ == 0 reads"},
            "NS:": {"Number":"1", "Type":"Integer", "Description":"Number of samples with data"},
            "SB": {"Number":"4", "Type":"Integer", "Description":"Strand bias"},
            "SOMATIC": {"Number":"0", "Type":"Flag", "Description":"Somatic mutation (for cancer genomics)"},
            "VALIDATED": {"Number":"0", "Type":"Flag", "Description":"Validated by follow-up experiment"},
            "1000G": {"Number":"0", "Type":"Flag", "Description":"1000 Genomes membership"}
        }
        for key in header_tag_array["INFO"]:
            if key in defaults:
                header_tag_array["INFO"][key] = {**defaults[key], **header_tag_array["INFO"][key]}
        for key in defaults:
            if key not in header_tag_array["INFO"]:
                header_tag_array["INFO"][key] = defaults[key]
        return header_tag_array


    @staticmethod
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

        if cparse:
            return cparse.parse_samples(
                self.samples, samples, samp_fmt, samp_fmt._types, samp_fmt._nums, site)

        samp_data = []
        _map = self._map

        nfields = len(samp_fmt._fields)

        for name, sample in itertools.izip(self.samples, samples):

            # parse the data for this sample
            sampdat = [None] * nfields

            for i, vals in enumerate(sample.split(':')):

                # short circuit the most common
                if samp_fmt._fields[i] == 'GT':
                    sampdat[i] = vals
                    continue
                # genotype filters are a special case
                elif samp_fmt._fields[i] == 'FT':
                    sampdat[i] = self._parse_filter(vals)
                    continue
                elif not vals or vals == ".":
                    sampdat[i] = None
                    continue

                entry_num = samp_fmt._nums[i]
                entry_type = samp_fmt._types[i]

                # we don't need to split single entries
                if entry_num == 1:
                    if entry_type == 'Integer':
                        try:
                            sampdat[i] = int(vals)
                        except ValueError:
                            sampdat[i] = float(vals)
                    elif entry_type == 'Float' or entry_type == 'Numeric':
                        sampdat[i] = float(vals)
                    else:
                        sampdat[i] = vals
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

    @staticmethod
    def _parse_alt(self, str):
        if self._alt_pattern.search(str) is not None:
            # Paired breakend
            items = self._alt_pattern.split(str)
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

    @staticmethod
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
                entry_type = self.infos[ID].type
            except KeyError:
                try:
                    entry_type = RESERVED_INFO[ID]
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
                if self.infos[ID].num == 1 and entry_type not in ('Flag',):
                    val = val[0]
            except KeyError:
                pass

            retdict[ID] = val

        return retdict