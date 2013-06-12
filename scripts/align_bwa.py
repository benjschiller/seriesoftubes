# !/usr/bin/env python
"""
align FASTQ, SAM, or BAM file (gzip and bzip2 supported)
with bwa, produces BAM files (sorted and indexed)
output will be in ./align_bwa
In align_bwa, there will be a folder for each reference genome
    e.g. align_bwa/ref1 , align2/ref2, align2/ref1, align2/ref2
"""
import logging
import os
import errno
from os import environ, curdir, devnull  # , makedirs
from os.path import join, exists
from textwrap import dedent
import sys
from subprocess import Popen, PIPE, STDOUT
from scripter import assert_path, path_to_executable, \
    exit_on_Usage, get_logger, Environment, critical, debug, info
from seriesoftubes.tubes.polledpipe import PolledPipe
from seriesoftubes.tubes import wait_for_job
from seriesoftubes.fnparsers import BowtieFilenameParser

from pkg_resources import get_distribution
__version__ = get_distribution('seriesoftubes').version
VERSION = __version__


def main():
    e = Environment(version=VERSION, doc=__doc__)
    e.set_filename_parser(BowtieFilenameParser)
    # let bwa do the multiprocessing
    parser = e.argument_parser
    parser.add_argument('--path-to-bwa', nargs='?',
                        default=path_to_executable('bwa', '/usr/local/bwa*',
                                                   environ='SOT_PATH_TO_BWA'),
                        help='The path to the bwa executable')
    parser.add_argument('--path-to-samtools', nargs='?',
                        default=path_to_executable('samtools',
                                                   '/usr/local/samtools*',
                                                   environ=
                                                   'SOT_PATH_TO_SAMTOOLS'),
                        help='The path to the samtools executable')
    # fix aliases, should be --ref too
    parser.add_argument('--reference', dest='references', action='append',
                        help=dedent('''\
    Reference genome to align against (should be a
    fasta file indexed by bwa). This flag may be called multiple times
    (which will cause each reference to be aligned to separately). If no
    references are specified, we'll look the for environment variable
    SOT_DEFAULT_REFERENCES, which should be given as a list,
    e.g. "foo foo2 foo3"'''),
                        )
    parser.add_argument('--passthru-args', nargs='*',
                        help='A list of arguments to be passed through to bwa '
                             'Substitute + '
                             'for - (e.g., --passthru-args +m 4 50')
    context = e.get_context()
    new_references = validate_references(**context)
    e.update_context({'references': new_references})
    sequence = e.get_sequence(**context)
    e._sequence = merge_pairs(sequence)
    e.do_action(align_bwa)


def fasta_to_bwa(fasta_file, path_to_bwa='bwa'):
    """given a filename, makes a bwa index
    if that file is a FASTA file
    """
    if exists(fasta_file):
        f = open(fasta_file, 'rU')
        for line in f:
            if line.startswith('#'):
                continue
            elif line.startswith('>'):
                args = [path_to_bwa, 'index', fasta_file]
                debug(" ".join(args))
                P = Popen(args, stdout=open(devnull, 'w'), stderr=PIPE)
                stderr = P.communicate()[1]
                if stderr.splitlines()[0].startswith('Error'):
                    return None
                else:
                    return fasta_file
    return None


def merge_pairs(old_sequence):
    """pair sequence files that differ only by 1/2 in one position
    """
    new_sequence = []
    sequence = []
    for fp_obj in old_sequence:
        if fp_obj.paired_end:
            new_sequence.append(fp_obj)
        else:
            sequence.append(fp_obj)
    input_files = [fp_obj.input_file for fp_obj in sequence]
    mates = []
    for i in xrange(len(input_files)):
        if i in mates:
            continue
        these_mates = []
        for j in xrange(i, len(input_files)):
            identity = [(not x == y) for x, y in
                        zip(input_files[i], input_files[j])]
            if sum(identity) == 1:
                index = identity.index(True)
                if (input_files[i][index] == '1' and
                        input_files[j][index] == '2')\
                   or \
                   (input_files[i][index] == '2' and
                        input_files[j][index] == '1'):
                        these_mates.append(j)
                        print i, j, input_files[i], input_files[j]
        if len(these_mates) == 0:
            continue
        elif len(these_mates) == 1:
            j = these_mates[0]
        else:
            # User input required
            print "Ambiguous filename pairing"
            print "Please select the correct mate pair for %s:" % \
                  input_files[i]
            while True:
                for j in these_mates:
                    print "[%d] %s" % (j, input_files[j])
                rawinput = raw_input("Enter your choice: ")
                try:
                    choice = int(rawinput.strip())
                except ValueError:
                    continue
                if not choice in these_mates:
                    print "%d is not a valid choice" % choice
                else:
                    j = choice
                    break
        identity = [(not x == y) for x, y in
                    zip(input_files[i], input_files[j])]
        index = identity.index(True)
        if input_files[i][index] == '1' and input_files[j][index] == '2':
            sequence[i].second_file = sequence[j].input_file
            sequence[i].paired_end = True
            mates.append(j)
        elif input_files[j][index] == '1' and input_files[i][index] == '2':
            sequence[j].second_file = sequence[i].input_file
            sequence[j].paired_end = True
            mates.append(i)

    new_sequence.extend([x for i, x in enumerate(sequence) if not i in mates])
    return new_sequence


def validate_references(references=None, path_to_bwa='bwa',
                        logger=None, environ_key='SOT_DEFAULT_REFERENCES',
                        target_dir=curdir,
                        **kwargs):
    ## Make the output directory, complain if we fail
    #if os.path.exists(target_dir):
    #    debug('Output directory %s already exists', target_dir)
    #else:
    #    debug('Creating directory "%s"', target_dir)
    #    makedirs(target_dir, mode=0755)
    #    if not os.path.exists(target_dir):
    #        raise IOError('Could not create directory %s' % target_dir)
    debug('Validating references')
    new_references = []
    if references is None:
        if environ_key in environ:
            references = environ[environ_key].split()
        else:
            critical('no reference genomes specified')
            return []

    for r in references:
        if exists(r):
            if not all(map(exists, [r + '.amb', r + '.ann', r + '.bwt',
                                    r + '.pac', r + '.sa'])):
                info('Attempting to build bwa index from %s' % r)
                new_index = fasta_to_bwa(r, target_dir=target_dir,
                                         path_to_bwa=path_to_bwa)
                if new_index is not None:
                    new_references.append(new_index)
                    continue
                else:
                    critical('Failed to build bwa index.')
            else:
                debug('Found bwa index for %s' % r)
                new_references.append(r)
        else:
            critical('bwa could not find the reference %s', r)
            critical('we will not align to %s', r)
    return new_references


def convert_to_fastq(fp_obj, logger=None):
    """Convert a SAM or BAM file to FASTQ file(s) for alignment
    """
    input_file = fp_obj.input_file
    output_dir = fp_obj.output_dir
    fastq_dir = join(output_dir, 'fastq_input')
    fp_obj.check_output_dir(fastq_dir)
    protoname = fp_obj.protoname
    if fp_obj.paired_end:
        fastq_filenames = (join(fastq_dir, '%s.1.txt.gz' % protoname),
                           join(fastq_dir, '%s.2.txt.gz' % protoname))
        logger.info('Converting file %s to FASTQ files %s, %s',
                    input_file, fastq_filenames[0], fastq_filenames[1])
        in_args = [sys.executable, '-m',
                   'seriesoftubes.converters.bamtofastq2',
                   '--gzip', input_file, fastq_filenames[0],
                   fastq_filenames[1]]
    else:
        fastq_filename = join(fastq_dir, '%s.txt.gz' % protoname)
        logger.info('Converting file %s to FASTQ file %s',
                    input_file, fastq_filename)
        in_args = [sys.executable, '-m',
                   'seriesoftubes.converters.bamtofastq2',
                   input_file, fastq_filename]
    logger.debug('Launching %s', ' '.join(in_args))
    polledpipe = PolledPipe(logger=logger, level=logging.ERROR)
    job = Popen(in_args, stdout=polledpipe.w, stderr=STDOUT)
    wait_for_job(job, [polledpipe], logger)
    if fp_obj.paired_end:
        logger.debug('Settings input_file to %s', fastq_filenames[0])
        fp_obj.input_file = fastq_filenames[0]
        logger.debug('Settings second_file to %s', fastq_filenames[1])
        fp_obj.second_file = fastq_filenames[1]
    else:
        logger.debug('Settings input_file to %s', fastq_filename)
        fp_obj.input_file = fastq_filenames[0]
    logger.debug('Setting use_pysam to False')
    fp_obj.use_pysam = False
    logger.debug('Setting format to FASTQ')
    fp_obj.format = 'FASTQ'
    logger.debug('Ignoring open_func, it will not be used')
    if not job.returncode == 0:
        logger.critical('Conversion FAILED!')
    else:
        logger.info('Conversion successful')
    return


@exit_on_Usage
def align_bwa(fp_obj, references=[],
              unique=True, seed_len='28',
              use_quality=True, logging_level=10, num_threads=1,
              passthru_args=None,
              **kwargs):
    common_flags = ['-t', str(num_threads), '-M']

    logger = get_logger(logging_level)
    if references is None or references == []:
        logger.critical('Nothing to do')
        return
    stdout_buffer = []

    if fp_obj.format == 'BAM' or fp_obj.format == 'SAM':
        convert_to_fastq(fp_obj, logger=logger)

    if not fp_obj.format == 'FASTQ':
        logger.critical('%s only supports FASTQ files' % __file__)
        return

    if passthru_args is not None:
        for i in range(len(passthru_args)):
            passthru_args[i] = passthru_args[i].replace('+', '-')
        logger.debug('Passing thru arguments %s', ' '.join(passthru_args))
    kwargs['passthru_args'] = passthru_args

    flags = [item for item in common_flags]

    for ref in references:
        s = align_once(fp_obj, flags, ref, logger=logger, **kwargs)
        stdout_buffer.append(s)
    return '\n'.join([s for s in stdout_buffer if s is not None])


def make_paired_name(name1, name2):
    identity = [(not x == y) for x, y in zip(name1, name2)]
    if not sum(identity) == 1:
        raise ValueError('Not valid names')
    else:
        index = identity.index(True)
        return identity[0:index] + '%' + identity[(index + 1):]


def align_once(fp_obj, flags, ref, use_quality=False,
               path_to_bwa=None, path_to_samtools=None, logger=None,
               passthru_args=None,
               **kwargs):
    refname = os.path.basename(ref)
    path_to_unsorted = fp_obj.tmp_filename(refname)
    output_dir = os.path.split(path_to_unsorted)[0]
    try:
        fp_obj.check_output_dir(output_dir)
    except OSError, e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise
    filename1 = os.path.abspath(fp_obj.input_file)
    second_file = fp_obj.second_file
    file_args = [ref, filename1]
    if second_file is not None:
        file_args.append(os.path.abspath(second_file))

    if passthru_args:
        bwa_args = [path_to_bwa, 'mem'] + flags + passthru_args + file_args
    else:
        bwa_args = [path_to_bwa, 'mem'] + flags + file_args

    # finish parsing input here
    bwa_stderr = PolledPipe(logger=logger, level=logging.ERROR)
    logger.info('Launching bwa '
                '(output will be piped to samtools for BAM encoding)')
    logger.info(' '.join(bwa_args))
    bwa_aligner = Popen(bwa_args, stdout=PIPE, stderr=bwa_stderr.w,
                        bufsize=-1)

    samtools_args = [path_to_samtools, 'view', '-b', '-S', '-o',
                     path_to_unsorted, '-']
    logger.info('Launching samtools to encode bwa output as BAM')
    logger.info(' '.join(samtools_args))
    samtools_stdout = PolledPipe(logger=logger, level=logging.WARN)
    samtools_stderr = PolledPipe(logger=logger, level=logging.ERROR)
    samtools_viewer = Popen(samtools_args, stdin=bwa_aligner.stdout,
                            stdout=samtools_stdout.w,
                            stderr=samtools_stderr.w, bufsize=-1)

    logger.debug('Waiting for bwa to finish')
    pollables = [bwa_stderr, samtools_stdout, samtools_stderr]
    wait_for_job(bwa_aligner, pollables, logger)

    if not bwa_aligner.returncode == 0:
        logger.critical("bwa did not run properly [%d]",
                        bwa_aligner.returncode)
        samtools_viewer.terminate()
        samtools_viewer.poll()
        logger.critical("samtools terminated")
        return

    logger.debug('Alignment successfully completed')
    logger.debug('Waiting for samtools to finish')
    wait_for_job(samtools_viewer, [samtools_stdout, samtools_stderr], logger)
    if not samtools_viewer.returncode == 0:
        logger.critical("samtools view did not run properly [%d]",
                        samtools_viewer.returncode)
        return

    logger.debug('Unsorted BAM file successfully written')

    logger.info('Launching samtools again to sort BAM output')
    output_dir, output_file = os.path.split(path_to_unsorted)
    bam_file = os.path.splitext(output_file)[0]
    sorter_args = [path_to_samtools, 'sort', output_file, bam_file]
    logger.info(' '.join(sorter_args))
    samtools_stdout = PolledPipe(logger=logger, level=logging.WARN)
    samtools_stderr = PolledPipe(logger=logger, level=logging.ERROR)
    samtools_sorter = Popen(sorter_args, stdout=samtools_stdout.w,
                            stderr=samtools_stderr.w, cwd=output_dir)
    wait_for_job(samtools_sorter, [samtools_stdout, samtools_stderr], logger)
    if not samtools_sorter.returncode == 0:
        logger.critical("samtools sort did not run properly [%d]",
                        samtools_sorter.returncode)
        return

    # don't destroy the files until we're sure we succeeded!
    assert_path(os.path.join(output_dir, bam_file + '.bam'))
    logger.debug('Removing unsorted file %s', path_to_unsorted)
    os.remove(path_to_unsorted)

    logger.debug('Launching samtools again to index sorted BAM output')
    samtools_stdout = PolledPipe(logger=logger, level=logging.WARN)
    samtools_stderr = PolledPipe(logger=logger, level=logging.ERROR)
    index_args = [path_to_samtools, 'index', bam_file + '.bam']
    samtools_indexer = Popen(index_args, stdout=samtools_stdout.w,
                             stderr=samtools_stderr.w, cwd=output_dir)
    wait_for_job(samtools_indexer, [samtools_stdout, samtools_stderr], logger)
    if not samtools_indexer.returncode == 0:
        logger.critical("samtools index did not run properly [%d]",
                        samtools_indexer.returncode)
        return

    # Make sure indexing succeeds
    assert_path(os.path.join(output_dir, bam_file + '.bam.bai'))
    return

if __name__ == "__main__":
    main()
