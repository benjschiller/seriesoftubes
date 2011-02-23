#!/usr/bin/env python
# Creates a session for the genome browser
import sys
import os
import glob
URL_PREFIX = 'http://kryweb.ucsf.edu/data'
URL_SEP ='/'

def main():
    '''the main function'''

    if len(sys.argv) == 1:
        sys.stderr.write( 'Usage: ./convertSamToBam.py file1.bam file2.bam ...\n')
        sys.exit()
# else
    filenames = sys.argv[1:]
    for wildfilename in filenames:
        for filename in glob.glob(wildfilename):
            if not os.path.exists(os.getcwd() + os.sep + filename):
                sys.stderr.write('Error: file not found (' \
                    + os.curdir + os.sep + filename \
                    + '). Exiting...\n')
                sys.exit(1)

            track_line = get_track_line(filename)
            if track_line is not None: print track_line
	
def get_track_line(filename):
    url_parts = filename.split(os.sep)
    filename_parts = url_parts[-1].split(os.extsep)
    file_extension = filename_parts[-1]

    uniqueness = url_parts[-2]
    genome = url_parts[-3]
    parent_folder = url_parts[-4]
    sample_name = filename_parts[0]
    ChIP = filename_parts[1]
    display_name_parts = [sample_name]

    if filename_parts[-2]=='all' or filename_parts[-2].startswith('barcode'):
        display_name_parts.extend(['.', filename_parts[-2]])

    if ChIP=='Input': display_name_parts.extend([' ', ChIP])
    display_name_parts.extend([' ','(', parent_folder, '/', uniqueness, ')'])
    display_name = ''.join(display_name_parts)

    full_url = URL_SEP.join([URL_PREFIX]+url_parts[-4:])

    if file_extension=='bam':
        track_line = ' '.join(['track', 'type=bam',
                               'name="'+ display_name + '"',
                               'bigDataUrl="' + full_url + '"',
                               'description="' + url_parts[-1] + '"',
                               'visibility=squish', 'db=' + genome ,
                               'bamColorMode=strand'])

        return track_line

if __name__=='__main__': main()
