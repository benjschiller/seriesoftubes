#!/usr/bin/env python
"""
--references=ref1,ref2  use reference genomes ref1, ref2
--mismatches=n          allow n mismatches (if use_quality, allow n mismatches
                                            in the seed)
                                           (default: 2)
--use-quality           use quality scores (equivalent to -n mode)
                        (default is to use -n mode)
--quals-type=           integer, solexa1.3, solexa, phred33,
                          or phred64 (see bowtie)
                        (default is --quals-type=solexa1.3)
--max-quality=          specify maximum quality scores of all mismatched
                        positions (default: 70)

--long-reads            equivalent to --use-quality

# specify uniqueness of alignments
default is --unique (-m 1) and --random (-M 1)
--no-unique             equivalent to -m 1
--no-random             equivalent to -M 1
"""
import sys
import os
import multiprocessing
import subprocess
import scripter
from scripter import assert_path, print_debug, path_to_executable

SOURCE_DIR = 'sequences.FASTQ'
TARGET_DIR = 'alignments.SAM'
VERSION = "2.4"

def main():
    boolean_opts = ["long-reads", "use-quality"]
    long_opts = ["no-random", "no-unique", "mismatches=",
                 "quals-type=", "references=", "max-quality="] + boolean_opts
    e = scripter.Environment(long_opts=long_opts, version=VERSION, doc=__doc__)
    e.set_filename_parser(BowtieFilenameParser)
    num_cpus = e.get_num_cpus()
    e.set_num_cpus(1) # we'll let bowtie do the multiprocessing
    common_flags = ['-y', '-a', '--time', '--best', '--chunkmbs', '1024',
                    '--strata', '--sam', '-p', str(num_cpus)]
    path_to_bowtie = path_to_executable('bowtie', '/usr/local/bowtie-*')
    e.update_script_kwargs({'common_flags': common_flags,
                        'path_to_bowtie': path_to_bowtie})
    e.parse_boolean_opts(boolean_opts)
    e.update_script_kwargs(check_script_options(e.get_options()))
    e.source_dir = SOURCE_DIR
    e.target_dir = TARGET_DIR
    e.do_action(action)

def check_script_options(options):
    specific_options = {}
    if options.has_key('references'):
        specific_options['references'] = options['references'].split(',')
    else: # default to hg18,hg19
        specific_options['references'] = ['hg18', 'hg19']

    specific_options['unique'] = not options.has_key('no-unique')
    specific_options['random'] = not options.has_key('no-random')

    if options.has_key('mismatches'):
        specific_options['mismatches'] = options['mismatches']
    else:
        specific_options['mismatches'] = '2'

    if options.has_key('quals-type'):
        if options['quals-type'] in ['solexa', 'solexa1.3', 'phred64', 
                                       'phred33','integer']:
            specific_options['quals_type'] = options['quals-type']
        else: raise scripter.Usage('invalid quals-type')
    else:
        specific_options['quals_type'] = 'solexa1.3' # same as bowtie

    if options.has_key('max-quality'):
        specific_options['max_quality'] = options['max-quality']
    else:
        specific_options['max_quality'] = '70' # same as bowtie

    return specific_options

def action(filename, references=[], random=True, unique=True,
           long_reads=False, use_quality=False, max_quality='70',
           mismatches='2',  quals_type='solexa1.3', verbose=False,
           common_flags = [],
           **kwargs):
    if unique and random:
        uniqueness= {'unique': ['-m','1'], 'random': ['-M','1']}
    elif unique: uniqueness= {'unique': ['-m','1']}
    elif random: uniqueness= {'random': ['-M','1']}

    for match_type, match_flag in uniqueness.items():
        flags = [item for item in common_flags]
        flags.extend(match_flag)

        # In case we're from that bad day...
        source = filename.fastq_source
        if source=='081124_HWI-EAS355_0001_Meghan' and match_type=='unique':
            flags.extend(['-v', '3'])
        elif use_quality or long_reads:
            flags.extend(['-e', max_quality])
            flags.extend(['-n', mismatches])
            flags.append(''.join(['--', quals_type, '-quals']))
        else:
            flags.extend(['-v', mismatches])

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
            bowtie_job = subprocess.call(bowtie_args,
                                         stdout=sys.stdout,
                                         stderr=subprocess.STDOUT)
            if verbose: print_debug("Finished")
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

if __name__=="__main__": main()
