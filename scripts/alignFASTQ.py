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
scripter.SCRIPT_DOC = __doc__
scripter.SCRIPT_VERSION = "2.1"
scripter.SCRIPT_LONG_OPTS = ["no-random", "no-unique", "references="]
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

    return specific_options

def action(filename, references=[], random=True, unique=True, **kwargs):
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
        else:
            flags.extend(['-v','2'])

        for ref in references:
            path_to_output = filename.output_filename(ref, match_type)
            filename.check_output_dir(os.path.split(path_to_output)[0])

            bowtie_args = [PATH_TO_BOWTIE] + flags + [ref, 
                                                      filename.input_file, 
                                                      path_to_output]
            sys.stdout.write(os.linesep)
            sys.stdout.write(' '.join(bowtie_args))
            sys.stdout.write(os.linesep)
            sys.stdout.flush()
            continue # delete this when done
            bowtie_job = subprocess.Popen(bowtie_args, stdout=sys.stdout,
                                            stderr=sys.stdout)
            bowtie_job.wait()
    return

class BowtieFilenameParser(scripter.FilenameParser):
    def __init__(self, filename, *args, **kwargs):
        super(BowtieFilenameParser, self).__init__(filename, *args, **kwargs)
        self.fastq_source = self.protoname.split('.')[6]

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
