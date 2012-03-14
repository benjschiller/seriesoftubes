#!/usr/bin/env python
'''
make track hub from MACS output
'''
raise NotImplementedError
import sys
from os import getenv, listdir
from os.path import splitext
from subprocess import Popen
from scripter import path_to_executable, get_logger
from seriesoftubes.fnparsers import BAMFilenameParser
from bioplus.genometools import guess_bam_genome, genome, NoMatchFoundError#, TemporaryGenomeFile
from pkg_resources import get_distribution, VersionConflict
__version__ = get_distribution('seriesoftubes').version
VERSION = __version__
TARGET_DIR = 'trackhub'

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--version',
                        help='show version info and exit',
                        action='version',
                        version='%(prog)s {0!s}'.format(__version__))
    parser.add_argument('-g', '--genome', dest='genome_build',
                        required=True, help='UCSC genome build (ex. hg19)')
    parser.add_argument('--ssh', nargs='?',
                        default=getenv("SERIESOFTUBES_SERVER_PATH"),
                        help="a valid path to a webserver where files can be uploaded using the ssh2 protocol. Passwords are not supported, you must use secure key authentication. This may provided as an environment variable SERIESOFTUBES_SERVER_PATH")
    parser.add_argument('--http', nargs='?',
                        default=getenv("SERIESOFTUBES_SERVER_URL"),
                        help='A valid URL for the webserver directory where files are uploaded to.  This may provided as an environment variable SERIESOFTUBES_SERVER_URL"')
    parser.set_defaults(**{'target': 'peaks'})
    args = vars(parser.parse_args())
    if args['ssh'] is None:
        raise Usage, 'No server path for upload'
    if args['http'] is None:
        raise Usage, 'No http URL given. Track hub will not be produced'
    try: bg2bw = path_to_executable("bedGraphToBigWig")
    except Usage:
        raise Usage, 'bedGraphToBigWig must be in PATH in order to make a track hub'
    #e.set_filename_parser(BAMFilenameParser)
    sys.exit()
          
def f(x, **kwargs):
    print x
    return x

def make_track_hub(name, genome_build, macs_dir,
                   ssh=None, http=None, logger=None):
    """
    make a track hub out of the MACS files on your local webserver
    """
    try: chrom_sizes = TemporaryGenomeFile(genome_build)
    except NoMatchFoundError:
        raise Usage('Could not determine genome size for %s' % bam_file)
    chrom_file = chrom_sizes.name
    files = listdir(macs_dir)
    bdg_files = [f for f in files if splitext(f)[1]=='.bdg']
    for f in bdg_files:
        bwf = '%s.bw' % splitext(f)[0]
        logger.info('Converting bedgraph file %s to bigWig file %s (%s)',
                    f, bwf, genome_build)
        Popen([bg2bw, f, chrom_file, bwf], cwd=macs_dir).wait()
        

if __name__=="__main__": main()