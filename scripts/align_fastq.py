#!/usr/bin/env python
"""
align FASTQ files with bowtie, produces BAM files (sorted and indexed)
Output is two folders
    all.BAM/        sorted, indexed BAM files with mapped + unmapped reads
    aligned.BAM/    sorted, indexed BAM files with mapped reads only
In all/aligned.BAM, there will be a folder for each reference genome
    e.g. all.BAM/ref1 , all.BAM/ref2, aligned.BAM/ref1, aligned.BAM/ref2
In reference folder, there will be a folder for unique or random
    alignments (random is maq-like behavior)
    e.g. all.BAM/ref1/random, all.BAM/ref1/unique, etc.

--references=ref1,ref2  use reference genomes ref1, ref2
                        (default is hg19,hg18)
--mismatches=n          allow n mismatches in the seed (default is 2)
                        (allow n mismatches total if --ignore-quality)
--seed-length=m         use seed length of m (default is 28)
--ignore-quality        Use -v mode with bowtie, allows only n mismatches total
--quals-type=           integer, solexa1.3, solexa, phred33, or phred64 (see bowtie)
                        (default is --quals-type=solexa1.3)
--max-quality=          specify maximum quality scores of all mismatched
                        positions (default is 70), ignored in -v mode

--no-unique           do not produce unique/ folder 
--no-random           do not produce random/ folder 
--no-filtering        do not produce aligned.BAM/ folder
"""
import os
import platform
import re
import subprocess
import sys
import pysam
import scripter
from scripter import assert_path, print_debug, path_to_executable, Usage, \
                     exit_on_Usage, InvalidFileException

SOURCE_DIR = 'sequences.FASTQ'
TARGET_DIR = 'all.BAM'
FILTERED_DIR = 'aligned.BAM'
VERSION = "2.4"

def main():
    long_opts = ["no-random", "no-unique", "no-filtering", "mismatches=",
                 "seed-length=", "quals-type=", "references=", "max-quality="]
    e = scripter.Environment(long_opts=long_opts, version=VERSION, doc=__doc__)
    e.set_filename_parser(BowtieFilenameParser)
    if platform.system() == 'Windows':
        # then let python do the multprocessing
        num_cpus = 1
    else:
        # then let bowtie do the multiprocessing
        num_cpus = e.get_num_cpus()
        e.set_num_cpus(1) # we'll let bowtie do the multiprocessing
    common_flags = ['-y', '-a', '--time', '--best', '--chunkmbs', '1024',
                    '--strata', '--sam', '-p', str(num_cpus)]
    path_to_bowtie = path_to_executable('bowtie', '/usr/local/bowtie-*')
    e.update_script_kwargs({'common_flags': common_flags,
                        'path_to_bowtie': path_to_bowtie})
    e.update_script_kwargs(check_script_options(e.get_options()))
    e.set_source_dir(SOURCE_DIR)
    e.set_target_dir(TARGET_DIR)
    e.do_action(align, stay_open=True)
    e.set_num_cpus(num_cpus)
    
    e.do_action(produce_bam_files)

@exit_on_Usage
def check_script_options(options):
    specific_options = {}
    if options.has_key('references'):
        specific_options['references'] = options['references'].split(',')
    else: # default to hg18,hg19
        specific_options['references'] = ['hg19', 'hg18']

    unique = not options.has_key('no-unique')
    random = not options.has_key('no-random')
    if not unique and not random: raise Usage('Nothing to do')
    specific_options.update({'unique': unique, 'random': random})

    specific_options['filtering'] = not options.has_key('no-filtering')  
    specific_options['use_quality'] = not options.has_key('ignore-quality')  

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

    if options.has_key('seed-length'):
        specific_options['seed_len'] = options['seed-length']
    else:
        specific_options['seed_len'] = '28' # same as bowtie
        
    return specific_options

def align(fp_obj, references=[], random=True, unique=True, max_quality='70',
          quals_type='solexa1.3',  mismatches='2', seed_len='28',
          use_quality=True, verbose=False, debug=False,
          common_flags=[], **kwargs):
    uniqueness = {}
    if unique: uniqueness.update({'unique': ['-m','1']})
    if random: uniqueness.update({'random': ['-M','1']})
    stdout_buffer = []
    common_flags.extend(['-l', seed_len])
    for match_type, match_flag in uniqueness.items():
        flags = [item for item in common_flags]
        flags.extend(match_flag)
        
        source = fp_obj.fastq_source
        # In case we're from that bad day...
        if source=='081124_HWI-EAS355_0001_Meghan' and match_type=='unique':
            flags.extend(['-v', '3'])
        elif use_quality:
            flags.extend(['-e', max_quality])
            flags.extend(['-n', mismatches])
            flags.append(''.join(['--', quals_type, '-quals']))
        else:
            flags.extend(['-v', mismatches])

        if fp_obj.paired_end: flags.extend(['-X','600'])
        
        for ref in references:
            path_to_output = fp_obj.sam_filename(ref, match_type)
            output_dir = os.path.split(path_to_output)[0]
            fp_obj.check_output_dir(output_dir)
            fp_obj.sam_files.append(path_to_output)
            s = align_once(fp_obj.input_file, fp_obj.second_file, flags,
                           ref, path_to_output, debug=debug, **kwargs)
            stdout_buffer.append(s)
    return '\n'.join([s for s in stdout_buffer if s is not None])

def align_once(filename1, filename2, flags, ref, path_to_output,
          use_quality=False, path_to_bowtie= None, debug=True, **kwargs):
    if path_to_bowtie is None: path_to_executable('bowtie')
    
    if filename2 is not None:
        file_args = [ref, '-1', filename1, '-2', filename2, path_to_output]
    else:
        file_args = [ref, filename1, path_to_output]
    bowtie_args = [path_to_bowtie] + flags + file_args
    run_title = '\n{!s}\n'.format(' '.join(bowtie_args))
    if debug:
        print_debug('Launching botwie as ' + run_title)
        bowtie_job = subprocess.Popen(bowtie_args,
                                      stdout=sys.stdout,
                                      stderr=subprocess.STDOUT)
        bowtie_job.wait()
        return
    else:
        bowtie_job = subprocess.Popen(bowtie_args,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
        stdout, stderr = bowtie_job.communicate()
        return '{!s}{!s}{!s}'.format(run_title, stdout, stderr)

@exit_on_Usage
def produce_bam_files(fp_obj, verbose=False, debug=False,
                      filtering=True, *args, **kwargs):
    stdout_buffer = []
    bam_files = []
    # convert
    for sam_file in fp_obj.sam_files:
        name = os.path.splitext(sam_file)[0]
        temp_file = os.extsep.join([name,'tmp'])
        bam_file = os.extsep.join([name,'bam'])
        if debug: print_debug('Converting', sam_file, 
                              'to temporary BAM file', temp_file)
        pysam.view('-b', '-S', '-o' + temp_file, sam_file)
        
        if debug: print_debug('Sorting', temp_file, 'into', bam_file)
        pysam.sort(temp_file, name)
        
        if debug: print_debug('Removing temporary files', sam_file, temp_file)
        os.remove(sam_file)
        os.remove(temp_file)
        
        if debug: print_debug('Indexing', bam_file) 
        pysam.index(bam_file)
        
        bam_files.append(bam_file)
        stdout_buffer.append('Converted {!s} to {!s} (sorted and indexed)'.\
                             format(sam_file, bam_file))
    #filter
    if filtering:
        for bam_file in bam_files:
            new_bam_file = re.sub(TARGET_DIR, FILTERED_DIR, bam_file, 1)
            new_bam_dir = os.path.split(new_bam_file)[0]
            if not os.path.exists(new_bam_dir): os.makedirs(new_bam_dir)
            if debug: print_debug('Copying aligned reads from {!s} to {!s}'.\
                                  format(bam_file, new_bam_file))
            pysam.view('-F 0x4', '-b', '-o' + new_bam_file, bam_file)
        if debug: print_debug('Indexing', new_bam_file) 
        pysam.index(new_bam_file)
        stdout_buffer.append('Filtered aligned reads from {!s} to {!s} \
(sorted and indexed)'.format(bam_file, new_bam_file))
    
    return '\n'.join(stdout_buffer)

def get_pair_info(illumina_name):
    """
    take a filename from GERALD output and figure out whether it's the first
    or second read (or single-end read) and return a tuple
    (pair_index, second_file_name, new_output_name)
    """
    name_parts = illumina_name.split('.')
    for i in range(len(name_parts)):
        part = name_parts[i]
        subparts = part.split('_') 
        if len(subparts) > 0:
            if subparts[0] == 's':
                if len(subparts) == 1: continue
                pair_index = subparts[2]
                # lane = subparts[1]
                if pair_index == '1':
                    join_part = '_'.join(subparts[0:2] + subparts[3:])
                    new_output_name = '.'.join(name_parts[0:i] + [join_part] +\
                                                name_parts[i+1:])
                    second_part = '_'.join(subparts[0:2] + ['2'] + subparts[3:])
                    second_name = '.'.join(name_parts[0:i] + [second_part] + 
                                           name_parts[i+1:])
                    return (pair_index, second_name, new_output_name)
                elif pair_index == '2':
                    raise InvalidFileException
                    return
    return None

class BowtieFilenameParser(scripter.FilenameParser):
    def __init__(self, filename, debug=False,
                 *args, **kwargs):
        super(BowtieFilenameParser, self).__init__(filename, *args, **kwargs)
        self.sam_files = []
        self.second_file = None
        self.check_paired_end(debug=debug)
        if len(self.protoname.split('.')) > 6:
            self.fastq_source = self.protoname.split('.')[6]
        else:
            self.fastq_source = 'Unknown'
            
    def check_paired_end(self, debug=False):
        # check if this is a paired-end file
        # if so, grab its partner
        seqfile_name = os.path.basename(self.input_file)
        pair_info = get_pair_info(seqfile_name)
        if pair_info is not None:
            pair_index = pair_info[0]
            second_name = pair_info[1]
            new_name = pair_info[2]
            if debug: print_debug('NOTICE: Detected paired read file.')
            if pair_index == '1':
                if debug: print_debug('Attempting to find second file.')

                self.second_file = os.sep.join([self.input_dir, second_name])
                self.protoname = os.path.splitext(new_name)[0]
                if debug: print_debug('Found', self.second_file)
                try:
                    assert_path(self.second_file)
                    self.paired_end = True
                except IOError:
                    if debug: print_debug('Failed to find paired end file')
                    self.paired_end = False
            elif pair_index == '2':
                if debug: print_debug('This is the second file, ignoring it.')
                raise InvalidFileException
            else:
                if debug: print_debug('Failed to find paired end')
                self.paired_end = False
        else: self.paired_end = False

    def sam_filename(self, ref, match_type):
        return os.path.join(self.output_dir, ref, match_type,
                                       self.with_extension('sam'))

if __name__=="__main__": main()
