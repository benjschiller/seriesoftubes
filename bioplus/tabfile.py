'''tools for dealing with files that have tabular data'''
import os
import itertools
import re
import gzip
import bz2

def merge_files(left, right, output, comments='left'):
    '''
    merge_files merges the tab-delimited files named left and right,
    which may have commmented lines. the output is directed to the
    file named output.

    There are few modes. If comments='left', comments in left are preserved.
    If comments='right', comments in right are preserved. If comments='none', 
    no comments are preserved. If comments='all', all comments in left and right
    are appended to the beginning of output, although they may previously have
    been contained within the data in left or right.

    Use merge_tab_files if you need to pass custom parameters to TabFile.
    '''
    file1 = TabFile(left)
    file2 = TabFile(right)
    return merge_tab_files(file1, file2, output, comments)

def merge_tab_files(file1, file2, output_filename, comments='left'):
    '''merge_tab_files merges the tab-delimited files represented by TabFile
       objects left and right, which may have commmented lines. the output is
       directed to the file named output.

        There are few modes. If comments='left', comments in left are preserved.
        If comments='right', comments in right are preserved.
        If comments='none', no comments are preserved.
        If comments='all', all comments in left and right are appended to the
            beginning of output, although they may previously have been
            contained within the data in left or right.

    See also merge_files
    '''
    if comments == 'all':
        file1.process_table(file1, lambda x: x + file2.read_row())
    elif comments == 'right':
        file2.process_table(file2, lambda x: x + file1.read_row())
    elif comments == 'none':
        outfile = TabFile(output_filename, write=True)
        for x in file1: outfile.write_row(x + file2.read_row())
    elif comments == 'all':
        # write comments first
        for x in file1.comment_line_contents(): outfile.write_row(x)
        for x in file2.comment_line_contents(): outfile.write_row(x)
        for x in file1: outfile.write_row( x+ file2.read_row())
    else:
        raise ValueError("comments must be one of 'left', 'right', 'none', 'all'")

class TabFileError(Exception):
    def __init__(self, *args):
        self.msg = ' '.join(args)
        return
    
    def __str__(self):
        return repr(self.msg)

class DetectCommentsError(TabFileError):
    def __init__(self, *args):
        super(DetectCommentsError, self).__init__(*args)
        return
                
class TabFile(object):
    '''Usage: f = TabFile('filename', convert_spaces=True,
                          comments=[],
                          column_names = False)

       TabFile is a class for handling tab-delimited files.

       Use convert_spaces=False if you're file is tab-delimited and you wish to
       preserve other whitespace.
       
       TabFile suports commented lines. Commented lines are not recognzied as
       part of the table

       By default, only lines beginning with '#' will be recognized as comments
       (not part of the table). You may specify a list of additional keywords
       using comments=['keyword1','keyword2',etc.]. All lines containing that
       keyword will be recognized as a comment. keywords may be regular
       expressions.
['(?i)track','(?i)browser']
       if column_names = True, the first properly formatted row will be treated as
       column names (i.e. ignored as a comment)
    '''

    def __init__(self, filename, mode='r', convert_spaces=True,
                 compression=None, comments=[], column_names=False):
        self._filename = filename
        self._file_extension = os.path.splitext(filename)[1].lstrip(os.extsep)
        self.mode = mode
        self.open(mode)
        self._previous_line = 0
        # _detect_comments sets _comment_line_numbers,
        #                       _comment_line_contents
        if 'r' in mode :
            try:
                self._detect_comments(comments, column_names)
            except:
                raise DetectCommentsError("Failed with comments=", *comments)
            
    def previous_line(self):
        '''returns the line number of the last line read'''
        return self._previous_line

    def _detect_comments(self, comments, column_names=False):
        self.close()
        self.open()
        if not os.path.isfile(self._filename):
            self._column_names = None
            self._comment_line_numbers = []
            self._comment_line_contents = []
            return

        number_keywords = len(comments)
        if number_keywords is not 0:
            searchables = [re.compile(keyword) for keyword in comments]

        if column_names:
            first_valid_line = True
        else:
            first_valid_line = False
            self._column_names = None

        line_number = 0
        comment_line_numbers = []
        comment_line_contents = []
        for line in self.__rawiter__():
            line_number += 1
            # blank lines are comments
            if line.lstrip() == '':
                comment_line_numbers.append(line_number)
                comment_line_contents.append(line)
            # lines starting with # are comments
            elif line.lstrip()[0] == '#':
                comment_line_numbers.append(line_number)
                comment_line_contents.append(line)
            elif number_keywords is not 0:
                for searchable in searchables:
                    if searchable.search(line) is not None:
                        comment_line_numbers.append(line_number)
                        comment_line_contents.append(line)
                        break
            elif first_valid_line:
                first_valid_line = False
                self._column_names = self._parse_line(line)
                comment_line_numbers.append(line_number)
                comment_line_contents.append(line)

        self._comment_line_numbers = comment_line_numbers
        self._comment_line_contents = comment_line_contents
        self.close()
        self.open()

    def _parse_line(self, input_string, convert_spaces=True):
        '''
        parses tab-separated elements in a string and returns a list.
        if convert_spaces=True, _parse_line will treat contiguous whitespace as
        a tab
        '''
        if convert_spaces:
            return input_string.split()
        else:
            return input_string.split('\t')
    
    def _make_line(self, input_array):
        '''takes a list or array and returns a tab-delimited string'''
        return '\t'.join([str(item) for item in input_array]) + os.linesep 

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._file_pointer.close()

    def open(self, mode=None):
        '''
        mode can be overriden here but defaults to TabFile.mode
        
        acts just like the built-in open method in the file class. use
        write=True to write to a file, otherwise it will be opened in read-only
        mode
        '''
        if mode is None: mode = self.mode
        self._file_pointer = open(self._filename, mode)

    def zap(self):
        '''forces status to not open. use with caution. this may destroy data'''
        self._file_pointer = None
        self._previous_line = 0

    def close(self):
        '''works just like the built-in close method in the file class'''
        if self._file_pointer is None:
            raise IOError('{!s} not open'.format(self._filename))
        else:
            self._file_pointer.close()
        self._previous_line = 0

    def read_table(self, override=False):
        '''returns the contents of a file as a list of rows (with each row as a list). will ignore any lines that begin with a "#" symbol and truncate any lines that contain a "#" symbol'''
        if self._file_pointer is None: self.open()
        elif not override:
            raise UserWarning('File already opened and must be closed and re-opened to read whole table.')
        return [ row for row in self ]

    def read_cols(self, L):
        '''read_cols behaves like read_col but instead of taking a single column number (numbering starts at 0), it takes a list of column numbers and returns a list of partial rows, where each partial row is a list with entries from the appropriate columns IN THE ORDER SPECIFIED.

tip: use range() to create lists of ordered integers. e.g., range(2,6)=[2,3,4,5]
'''
        if self._file_pointer is None: self.open()
        elif type(L)==list:
            raise UserWarning('File already opened and must be closed and \
re-opened to read all rows of the column.')
        else: raise TypeError('read_cols requires a list of column numbers (0-based) as input')

        my_array = []
        for row in self:
            my_array.append([])
            for n in L: my_array[-1].append(row[n])
        return my_array

    def read_col(self, n):
        '''
        returns a list of items in column n (numbering starts at 0) as items
        rather than lists. Unlike read_cols and read_table, elements of the
        read_col list are not lists, but strings
        '''
        if self._file_pointer is None: self.open()
        elif type(n) is int:
            raise UserWarning('File already opened and must be close and \
re-open to read all rows of the column')
        else:
            raise TypeError('read_col requires an integer \
(column number, 0-based) as input')

        return [ row[n] for row in self ]

    def read_first_col(self):
        '''return a list of items in the first column'''
        return self.read_col(0)

    def read_last_col(self):
        '''returns a list of items in the last column'''
        return self.read_col(-1)

    def readline(self):
        '''
        reads one line and returns it. uses a generator, and will raise
        StopIteration if it reaches the EOF.

        readline is deprecated. use x = self.__rawiter__() and x.next()
        '''
        lines = self.__rawiter__()
        return lines.next()

    def __rawiter__(self):
        '''like iter, but does not process rows'''
        if self._file_pointer is None: self.open()
        while True:
            self._previous_line += 1
            thisline = self._file_pointer.readline()
            if thisline=='':
                self.close()
                raise StopIteration
            else:
                yield thisline

    def __iter__(self):
        '''returns the next (or first) line that is not a comment, parsed'''
        for line in self.__rawiter__():
            if not self.previous_line() in self.comment_line_numbers():
                yield self._parse_line( line )

    def comment_line_contents(self):
        '''returns the list of lines that are comments'''
        return self._comment_line_contents

    def comment_line_numbers(self):
        '''returns the list of lines that are comments'''
        return self._comment_line_numbers

    def read_row(self):
        '''returns the next (or first) line that is not a comment, parsed. uses __iter__ as a generator, and simply returns the next value

        read_row is deprecated. Use x = self.__iter__() and x.next()
        '''
        gen = self.__iter__()
        return gen.next()
    
    def write(self, s):
        '''
        writes a string directly to a file,
        without modification (user must supply \n if desired
        '''
        if self._file_pointer is None:
            raise IOError('File not open for writing.')
        else:
            self._file_pointer.write(s)

    def write_rows(self, iterable):
        for row in iterable: self.write_row(row)
        return
    
    def write_row(self, row, separator='\t'):
        '''Writes a list to the file as a line (Tab-delimited). A different separator may also be specified with separator='x'. (Note: uses file writelines method)'''
        if self._file_pointer is None:
            raise IOError('File not open for writing.')
        else:
            self._file_pointer.writelines([self._make_line(row)])

    def write_table(self, table, separator='\t', override=False,
                    column_names=True):
        '''
        Writes a table to a file (tab-delimited).
        An alternative separator may be specified with separator='x'.
        if column_names = True, column_names will be included as the first line
        unless they do not exist.
        '''
        if self._file_pointer is None:
            raise IOError('File not open for writing.')
        elif not override:
            raise UserWarning("File already opened and must be closed and \
re-opened to write the table. You may override this behavior by passing \
override=True to the method')")
        try:
            # write column_names if applicable
            if column_names and not self._column_names==None:
                self._file_pointer.wiritelines(
                    self._make_line(self._column_names))
            write_list = itertools.imap(self._make_line, table)
            self._file_pointer.writelines(write_list)
        finally:
            self.close()

    def process_table(self, output_filename, fnc, column_names=None):
        '''
        process_tables2 writes a new file (name is specified with new_file), 
        which applies a user-defined function fnc to each row of data in the 
        original file. fnc should yield a row (i.e. a list, array, or something 
        else finitely iterable).

        process_tables2 preserves all commented lines and also the line which 
        column names, if applicable. The user may specify new column names 
        using column_names (a list or other finite iterable), or we will use 
        the old column_names, which might not preserve the column labels if 
        columns were inserted in the middle of the table
        '''
        if self._file_pointer is None: self.open()
        output_file = TabFile(output_filename, 'w')

        if not column_names == None:
            first_valid_line = True
            output_file.set_column_names(column_names)
        elif not self._column_names==None:
            first_valid_line = True
            output_file.set_column_names(column_names)
        else: first_valid_line = False

        # treat rawiter as a generator
        # we're going to thread it
        lines = self.__rawiter__()

        # make sure we set column_names correctly
        if first_valid_line:
            while True:
                line = lines.next()
                if self.previous_line() in self.comment_line_numbers():
                    output_file.write(line)
                else:
                    first_valid_line = False
                    output_file.write_row(column_names)
                    break

        # proceed with the remainder
        while True:
            line = lines.next()
            if self.previous_line() in self.comment_line_numbers():
                output_file.write(line)
            else:
                output_file.write_row(fnc(self._parse_line(line)))
        output_file.close()

    def get_column_names(self):
        '''returns the column names'''
        return self._column_names

    def column_dict(self):
        '''returns a dictionary which gives the index corresponding to a particular column name'''
        d = {}
        for i in range(len(self._column_names)):
            d[self._column_names[i]] = i
        return d

    def set_column_names(self, L):
        self._column_names = L

    def write_column_names(self):
        '''writes the stored column names'''
        self.write_row(self._column_names)

    def mergesort(self,f, n,numerical=False):
        raise NotImplementedError
#        if self._file_pointer is None:
#            raise UserWarning('File already opened and must be closed and re-opened to mergesort the file. This behavior may not be overridden.')
#        fcontents = self.read_table(f)
#        if numerical==True: fcontents = self._numerize(fcontents,n)
##        newcontents = merge.mergesort(fcontents,n)
##        writeFile(newcontents,f)

class BedRow(list):
    '''BEDrows are list, but you can access their chromStart(), chromEnd(), etc. use help for a full list. Uses the same conventions as http://genome.ucsc.edu/FAQ/FAQformat#format1. Note that only the first three entries (chrom, chromStart, chromEnd) are required, so the others may not be defined.'''
    def chrom(self):
        '''returns the name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671)'''
        return self[0]

    def chrom_start(self):
        '''returns the starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0'''
        return self[1]

    def chromStart(self):
        raise DeprecationWarning()
        return self.chrom_start()

    def chrom_end(self):
        '''returns the ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99'''
        return self[2]

    def chromEnd(self):
        raise DeprecationWarning()
        return self.chrom_end()

    def name(self):
        '''returns the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.'''
        return self[3]

    def score(self):
        '''returns the score, a number between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray)'''
        return self[4]

    def strand(self):
        '''returns the strand, either '+' or '-' '''
        return self[5]

    def thickStart(self):
        '''returns the starting position at which the feature is drawn thickly (for example, the start codon in gene displays)'''
        return self[6]


    def thickEnd(self):
        '''returns the ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).'''
        return self[7]

    def itemRgb(self):
        '''returns itemRgb, An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browse'''
        return self[8]


    def blockCount(self):
        '''return the number of blocks (exons) in the BED line'''
        return self[9]


    def blockSizes(self):
        '''returns a comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.'''
        return self[10]

    def blockStart(self):
        '''returns a comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.'''
        return self[11]

class GzipTabFile(TabFile):
    '''
    For gzip-compressed tab-delimited files
    
    See Tabfile for usage info
    '''
    def __init__(self, *args, **kwargs):
        super(TabFile, self).__init__(*args, **kwargs)
        
    def _open(self, mode=None):
        '''
        mode can be overriden here but defaults to TabFile.mode
        
        acts just like the built-in open method in the file class. use
        write=True to write to a file, otherwise it will be opened in read-only
        mode
        '''
        if mode is None: mode = self.mode
        if not 'b' in mode: mode += 'b'
        self._file_pointer = gzip.open(self._filename, mode)
    
class Bzip2TabFile(TabFile):
    '''
    For bzip2-compressed tab-delimited files
    
    See Tabfile for usage info
    '''
    def __init__(self, *args, **kwargs):
        super(TabFile, self).__init__(*args, **kwargs)
        
    def _open(self, mode=None):
        '''
        mode can be overriden here but defaults to TabFile.mode
        
        acts just like the built-in open method in the file class. use
        write=True to write to a file, otherwise it will be opened in read-only
        mode
        '''
        if mode is None: mode = self.mode
        if not 'b' in mode: mode += 'b'
        self._file_pointer = bz2.BZ2File(self._filename, mode)
                
class BedFile(TabFile):
    '''
    A BED file is a type of TabFile, but also defines a method for working with
    rows. rows are given as instances of BedRow, instead of lists. BEDrows 
    inerhit all list methods and therefore are compatible with write_row. 
    BedRow has additional methods for chrom, chromStart, chromEnd, etc. 
    For more info, see BedRow

    track, browser lines are treated as comments

    Assumes track row is a comment. Use getTrackLine to see the track info
    '''
    DEFAULT_BED_COMMENTS = ['(?i)track','(?i)browser']

    def __init__(self, f, additional_comments=[], **kwargs):
        comments = self.DEFAULT_BED_COMMENTS + additional_comments
        TabFile.__init__(self, f, comments=comments,
                         **kwargs)
        self._track_line = None
        # check again for track line, and save it #
        for x in self._comment_line_contents:
            expr = re.compile('(?i)track')
            if expr.search(x) is not None:
                self._track_line = x
                break

    def get_track_line(self):
        '''returns the current track line, if any'''
        return self._track_line

    def __iter__(self):
        '''returns the next (or first) line that is not a comment, parsed'''
        for line in self.__rawiter__():
            if not self.previous_line() in self.comment_line_numbers():
                yield BedRow(self._parse_line( line ))

class MacsRow(list):
    '''MACSrows are list, but you can access their features as follows:        

    chr() or chrom() -- chromosome name
    start() or chromStart() -- start position, start() is 1-based, chromStart is 0-based (BED)
    end() or chromEnd() -- end position, equivalent but chromEnd (BED) is defined as 0-based, exclusive
    length() -- length
    summit() -- position of summit
    tags() -- number of unique tags in the peak region
    pvalue() -- returns the -10*log10(pvalue)
    fold_enrichment -- returns the fold enrichment
    FDR -- returns the FDR in %
    '''
    def chr(self):
        '''returns the name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671)'''
        return self[0]

    def chrom(self):
        '''returns the name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671)'''
        return self[0]

    def start(self):
        '''returns the starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 1'''
        return int(self[1])

    def chrom_start(self):
        '''returns the starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0'''
        return int(self[1]) - 1

    def chromStart(self):
        raise DeprecationWarning()
        return self.chrom_start()

    def end(self):
        '''returns the end position, 1-based, inclusive'''
        return int(self[2])

    def chrom_end(self):
        '''returns the ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99'''
        return int(self[2])

    def chromEnd(self):
        raise DeprecationWarning()
        return self.chrom_end()

    def length(self):
        '''returns the length'''
        return int(self[3])
    
    def summit(self):
        '''returns the position of the summit'''
        return int(self[4]) - 1

    def tags(self, type_=int):
        '''returns the number of unique tags in the peak region'''
        return type_(self[5])
    
    def tagsv1(self):
        return self.tags(int)
    
    def tagsv2(self):
        return self.tags(lambda x: int(float(x)))

    def pvalue(self,type_=str):
        '''returns the -10*log10(pvalue). preserves the str to eliminate rounding error. use type=float to get a decimal value'''
        return type_(self[6])

    def fold_enrichment(self,type=str):
        '''returns the fold_enrichment vs control. preserves the str to eliminate rounding error. use type=float to get a decimal value'''
        return type(self[7])

    def FDR(self,type_=str):
        '''returns the FDR (%). preserves the str to eliminate rounding error. use type=float to get a decimal value'''
        return type_(self[8])    

class MacsFile(TabFile):
    '''A MACS file is a type of TabFile, but also defines a method for working with rows. rows are given as instances of MACSRow, instead of lists. MACSrows inerhit all list methods and therefore are compatible with write_row. MacsRow has additional methods for chrom, chromStart, chromEnd, etc. For more info, see MacsRow'''

    def __init__(self, f, convert_spaces=True, **kwargs):
        super(MacsFile, self).__init__(f, column_names=True, **kwargs)

    def __iter__(self):
        '''returns the next (or first) line that is not a comment, parsed'''
        for line in self.__rawiter__():
            if not self.previous_line() in self.comment_line_numbers():
                yield MacsRow(self._parse_line(line))

def shift_peaks(f, peak_lengths=2):
    '''
    shift_peaks takes a file f (foo.bed)
    and produces a new file (foo_shifted.bed)
    with all the sequences shifted (left) by peak_lengths times their length.
    If peak_lengths, is negative they are shifted to the right.
    comments are stripped
    '''
    x = BedFile(f)
    # CHANGE TO USE FILENAME CORRECTION SCHEME
    y = BedFile(f[0:-4]+'_shifted.bed')
    y.open(write=True)
    for peak in x:
        peak_start, peak_end = int(peak[1]), int(peak[2])
        new_peak = peak.copy()
        peak_shift = peak_lengths * (peak_end - peak_start)
        # update start, end
        new_peak[1] = peak_start - peak_shift
        new_peak[2] = peak_end - peak_shift
        if int(new_peak[1])<0: new_peak[1] = 0
        y.write_row(new_peak)
    return

def _quote(s):
    return ''.join(["'", s ,"'"])
