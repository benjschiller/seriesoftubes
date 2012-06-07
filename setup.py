#!/usr/bin/env python
# Setup script

"""Description

Setup script for seriesoftubes 

Copyright (c) 2010-2012 Benjamin Schiller <benjamin.schiller@ucsf.edu>,
                        University of California, San Francisco

All rights reserved.
"""

import os
import sys
from distutils.core import setup
from distutils.extension import Extension

has_cython = True
try:
    import Cython.Distutils
    command_classes = {'build_ext': Cython.Distutils.build_ext}
except ImportError:
    has_cython = False
    command_classes =  {}

name = 'seriesoftubes'
version = '0.9.3.6'

def main():
    if not float(sys.version[:3])>=2.7:
    	sys.exit("CRITICAL: Python version must greater than or equal to 2.7! python 2.7.2 is recommended!\n")
          
    if has_cython:
        ext_modules=[Extension('seriesoftubes.cPreprocess', ["seriesoftubes/cPreprocess.pyx"])]
    else:
        ext_modules=[Extension('seriesoftubes.cPreprocess', ["seriesoftubes/cPreprocess.c"])]
          
    setup(name=name, version = version,
          description='An extended pipeline for Solexa ChIP-seq data',
          author='Benjamin Schiller',
          author_email='benjamin.schiller@ucsf.edu',
          url = 'https://github.com/benjschiller/seriesoftubes',
            cmdclass= command_classes,
            install_requires = ['scripter>=3.3.0',
    						  'biopython>=1.56',
    				          'pysam>=0.4',
    				          'twobitreader>=2.5',
    				          'MACS>=2.0.10'
    				         ],
          packages = ['bioplus', 'seriesoftubes', 'seriesoftubes.converters',
    				  'seriesoftubes.fnparsers', 'seriesoftubes.tubes'],
    	  package_data= {'bioplus': ['data/genomes.db']},
            ext_modules = ext_modules,                
            scripts = [os.path.join('scripts', x) for x in os.listdir('scripts') 
                       if not x.startswith('.')],
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
    			],
    	  dependency_links = ['https://github.com/downloads/benjschiller/MACS/benjschiller-MACS-v2.0.10pre2.tar.gz#egg=MACS-2.0.10'],
          )
if __name__ == '__main__':
	main()
