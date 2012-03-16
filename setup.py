#!/usr/bin/env python
# Setup script

"""Description

Setup script for seriesoftubes 

Copyright (c) 2010 Benjamin Schiller <benjamin.schiller@ucsf.edu>

All rights reserved.
"""

import os
import sys
from distutils.core import setup
try: import py2exe
except ImportError: pass
try: import py2app
except ImportError: pass
cmdclass = {}
try:
    from sphinx.setup_command import BuildDoc
    cmdclass['build_sphinx'] = BuildDoc
except ImportError: pass
try:
    from sphinx_pypi_upload import UploadDoc
    cmdclass['upload_sphinx'] = UploadDoc
except ImportError: pass
name = 'seriesoftubes'
version = '0.9.0'

def main():
	if not float(sys.version[:3])>=2.7:
		sys.exit("CRITICAL: Python version must greater than or equal to 2.7! python 2.7.2 is recommended!\n")
	setup(name=name, version = version,
	      description='An extended pipeline for Solexa ChIP-seq data',
	      author='Benjamin Schiller',
	      author_email='benjamin.schiller@ucsf.edu',
	      url = 'https://github.com/benjschiller/seriesoftubes',
          requires = ['scripter (>=3.1)', 'biopython (>=1.56)',
					  'twobitreader (>=2.4)', 'pysam (>=0.4)',
					  'MACS (>=2.0.10)'],
	      packages = ['bioplus', 'seriesoftubes', 'seriesoftubes.converters',
					  'seriesoftubes.fnparsers', 'seriesoftubes.tubes'],
		  package_data= {'bioplus': ['data/genomes.db']},
          scripts = [os.path.join('scripts', x) for x in os.listdir('scripts') 
                     if not x.startswith('.')],
		  cmdclass=cmdclass,
		  command_options={
			  'project': ('setup.py', name),
			  'version': ('setup.py', version),
		  },
  	      classifiers = [
				'Development Status :: 3 - Alpha',
				'License :: OSI Approved :: Artistic License',
				'Intended Audience :: Developers',
				'Intended Audience :: End Users/Desktop',
				'Intended Audience :: Science/Research',
				'Operating System :: MacOS :: MacOS X',
				'Operating System :: POSIX',
				'Programming Language :: Python :: 2.7',
				'Topic :: Scientific/Engineering :: Bio-Informatics'
				]
	      )
if __name__ == '__main__':
	main()
