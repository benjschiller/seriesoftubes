#!/usr/bin/env python
# Creates a session for the genome browser
import sys,os


def main():
	'''the main function'''
	url_prefix = 'http://kryweb.ucsf.edu/data'
	url_sep='/'

	if len(sys.argv) == 1:
		sys.stderr.write( 'Usage: ./convertSamToBam.py file1.bam file2.bam ...' + os.linesep )
		sys.exit()
	# else
	filenames = sys.argv[1:]
	for filename in filenames:
		if not os.path.exists(os.getcwd() + os.sep + filename):
			sys.stderr.write('Error: file not found (' \
				+ os.curdir + os.sep + filename \
				+ '). Exiting...' + os.linesep)
			sys.exit(1)
		url_parts = filename.split(os.sep)
		filename_parts = url_parts[-1].split(os.extsep)
		if not filename_parts[-1]=='bam': continue

		uniqueness = url_parts[-2]
		genome = url_parts[-3]
		parent_folder = url_parts[-4]
		sample_name = filename_parts[0]
		display_name_parts = [sample_name]

		if filename_parts[-2]=='all' or filename_parts[-2].startswith('barcode'):
			display_name_parts.extend(['.', filename_parts[-2]])

		display_name_parts.extend([' ','(', parent_folder, '/', uniqueness, ')'])
		display_name = ''.join(display_name_parts)

		full_url = url_sep.join([url_prefix]+url_parts[-4:])
		print ' '.join([ \
		'track','type=bam','name="'+ display_name + '"', \
		'bigDataUrl="' + full_url + '"', \
		'description="' + url_parts[-1] + '"', \
		'visibility=squish','db=' + genome ,\
		'bamColorMode=strand' \
		])
	
if __name__=='__main__': main()
