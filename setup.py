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

def main():
	if not float(sys.version[:3])>=2.5:
		sys.stderr.write("CRITICAL: Python version must greater than or equal to 2.5! python 2.6.4 is recommended!\n")
		sys.exit(1)
	setup(name='seriesoftubes',
          version = "2.4-pre1",
	      description='An extended pipeline for Solexa ChIP-seq data',
	      author='Benjamin Schiller',
	      author_email='benjamin.schiller@ucsf.edu',
          install_requires = ['scripter>=2.0.2', 'biopython>=1.56',
							  'twobitreader>=1.0.2', 'pysam>=0.4'],
	      packages = ['bioplus'],
          scripts = [os.path.join('scripts', x) for x in os.listdir('scripts') 
                     if not x.startswith('.')]
	      )
if __name__ == '__main__':
	main()
