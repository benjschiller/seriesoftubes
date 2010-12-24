#!/usr/bin/env python
'''
--references=ref1,ref2  use reference genomes ref1, ref2

# specify uniqueness of alignments
default is --unique (-m 1) and --random (-M 1)
--no-unique             equivalent to -m 1
--no-random             equivalent to -M 1
'''
import sys
import os
import multiprocessing
import subprocess
import scripter
from scripter import assert_path, print_debug
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.2"
BOOLEAN_OPTS = ["long-reads"]
scripter.SCRIPT_LONG_OPTS = ["no-random", "no-unique",
                             "references="] + BOOLEAN_OPTS
scripter.SOURCE_DIR = 'sequences.FASTQ'
scripter.TARGET_DIR = 'alignments.SAM'

PATH_TO_BOWTIE = '/usr/local/bowtie-0.12.5/bowtie'
NUM_CPUS = 1 # we'll change 
COMMON_FLAGS = ['-y', '-a', '--time', '--best', '--chunkmbs', '1024',
                    '--strata', '--sam', '--phred64-quals']

def check_script_options(options):
    specific_options = {}
    if options.has_key('references'):
        specific_options['references'] = options['references'].split(',')
    else: # default to hg18,hg19
        specific_options['references'] = ['hg18', 'hg19']

    specific_options['unique'] = not options.has_key('no-unique')
    specific_options['random'] = not options.has_key('no-random')

    for option in BOOLEAN_OPTS:
        pyoption = "_".join(option.split("-"))
        specific_options[pyoption] = options.has_key(option)

    return specific_options

def action(filename, references=[], random=True, unique=True,
           long_reads = False, **kwargs):
    if unique and random:
        uniqueness= {'unique': ['-m','1'], 'random': ['-M','1']}
    elif unique: uniqueness= {'unique': ['-m','1']}
    elif random: uniqueness= {'random': ['-M','1']}

    for match_type, match_flag in uniqueness.items():
        flags = [item for item in COMMON_FLAGS]
        flags.extend(match_flag)

        # In case we're from that bad day...
        source = filename.fastq_source
        if source=='081124_HWI-EAS355_0001_Meghan' and match_type=='unique':
            flags.extend(['-v','3'])
        elif long_reads:
            flags.extend(['-n','2'])
        else:
            flags.extend(['-v','2'])

        if filename.paired_end: flags.extend(['-X','600'])

        for ref in references:
            path_to_output = filename.output_filename(ref, match_type)
            filename.check_output_dir(os.path.split(path_to_output)[0])

            if filename.paired_end:
                bowtie_args = [PATH_TO_BOWTIE] + flags + [ref, '-1',
                                                          filename.input_file,
                                                          '-2',
                                                          filename.second_file,
                                                          path_to_output]
            else:
                bowtie_args = [PATH_TO_BOWTIE] + flags + [ref, 
                                                          filename.input_file, 
                                                          path_to_output]
            sys.stdout.write(os.linesep)
            sys.stdout.write(' '.join(bowtie_args))
            sys.stdout.write(os.linesep)
            sys.stdout.flush()
            bowtie_job = subprocess.Popen(bowtie_args, stdout=sys.stdout,
                                            stderr=sys.stdout)
            bowtie_job.wait()
    return

def get_pair_info(illumina_name):
    name_parts = illumina_name.split('.')
    for i in range(len(name_parts)):
        part = name_parts[i]
        if part.count('_') is not 0:
            if part.startswith('s_'):
                if part.count('_') is 1: continue
            pair_index = part.split('_')[-1]
            if pair_index == '1' or pair_index == '2':
                new_name = '.'.join(name_parts[0:i] +
                                     ['_'.join(part.split('_')[0:-1])])
                second_name = '.'.join(name_parts[0:i] +
                                     ['_'.join(part.split('_')[0:-1] + ['2'])])
                if i < len(name_parts):
                    new_name = '.'.join([new_name] + name_parts[i+1:])
                    second_name = '.'.join([second_name] + name_parts[i+1:])
                return (pair_index, second_name, new_name)
    return None

class BowtieFilenameParser(scripter.FilenameParser):
    def __init__(self, filename, verbose=False, *args, **kwargs):
        super(BowtieFilenameParser, self).__init__(filename, *args, **kwargs)
        if len(self.protoname.split('.')) > 6:
            self.fastq_source = self.protoname.split('.')[6]
        else:
            self.fastq_source = 'Unknown'

        # check if this is a paired-end file
        # if so, grab its partner
        illumina_name = os.path.basename(self.input_file)
        pair_info = get_pair_info(illumina_name)
        if pair_info is not None:
            pair_index = pair_info[0]
            second_name = pair_info[1]
            new_name = pair_info[2]
            if verbose: print_debug('NOTICE: Detected paired read file.')
            if pair_index == '1':
                if verbose: print_debug('Attempting to find second file.')

                self.second_file = os.sep.join([self.input_dir, second_name])
                self.protoname = os.path.splitext(new_name)[0]
                if verbose: print_debug('Found', self.second_file)
                try:
                    assert_path(self.second_file)
                    self.paired_end = True
                except IOError:
                    if verbose: print_debug('Failed to find paired end file')
                    self.paired_end = False
            elif pair_index == '2':
                if verbose: print_debug('This is the second file, ignoring it.')
                self.is_invalid = True
            else:
                if verbose: print_debug('Failed to find paired end')
                self.paired_end = False
        else: self.paired_end = False

    def output_filename(self, ref, match_type, mode='sam'):
        path_to_output = os.path.join(self.output_dir, ref, match_type,
                                       self.with_extension('sam'))
        return path_to_output

if __name__=="__main__":
    NUM_CPUS = scripter.NUM_CPUS # get the number of cpus
    scripter.NUM_CPUS = 1 # OVERRIDE NUM_CPUS
    COMMON_FLAGS.extend(['-p', str(NUM_CPUS)]) # Apply to bowtie

    scripter.check_script_options = check_script_options

    FilenameParser = BowtieFilenameParser
    scripter.perform(action, FilenameParser)
