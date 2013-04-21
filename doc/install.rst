.. seriesoftubes documentation master file, created by
   sphinx-quickstart on Thu Mar 15 18:58:51 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Installation of seriesoftubes
==============================================================

.. toctree::
   :maxdepth: 2

Installing seriesoftubes is easy if you already have Python and setuptools
installed:

.. code-block:: bash

   easy_install seriesoftubes
   
However, in order to use the scripts, you will need to manually install the
:ref:`external dependencies <extern>` too.

Requirements
============ 

Python 2.7
----------
seriesoftubes is written in and requires `Python <http://www.python.org>`_.
Currently version >= 2.7 is required
(`Python 2.7.3 <http://www.python.org/download/releases/2.7.3/>`_ is recommended)

Setuptools
----------
`Setuptools <http://pypi.python.org/pypi/setuptools>`_
(or `distribute <http://pypi.python.org/pypi/distribute#downloads>`_) is now
required for installation. It is probably already installed, but if not,
you must download and install it. It will install all the other python packages
for you.

Operating System
----------------

Several parts of seriesoftubes will not function properly on Windows, if at all.
Mac OS X and Linux/\*nix should work, and the pipeline has been most recently tested on `Ubuntu 12.04 LTS <http://releases.ubuntu.com/precise/>`_.

Python Packages (installed with easy_install)
---------------------------------------------

If you install seriesoftubes using :program:`easy_install` from `setuptools
<http://pypi.python.org/pypi/setuptools>`_ (or an equivalent package), then it
will automatically install most of the dependencies. If not, then you must
install the following packages

    - `Biopython <http://pypi.python.org/pypi/biopython>`_
    - `pysam <http://pypi.python.org/pypi/pysam>`_
    - `scripter <http://pypi.python.org/pypi/scripter>`_    
    - `MACS >=2.0.10 beta <https://github.com/taoliu/MACS/>`_

.. _extern:

External software dependencies
------------------------------

Bowtie2 (recommended) or Bowtie
-------------------------------

`Bowtie2`_ is used by :program:`align2.py` (recommended) to align reads to the reference genome. Legacy support for :program:`bowtie` is provided by `align.py`. `align2.py` expects :program:`bowtie2` to be in your PATH or in /usr/local/bowtie2(-\*). Similarly, `align.py` expects :program:`bowtie` to be in your PATH or in /usr/local/bowtie(-\*).

.. _Bowtie: http://bowtie-bio.sourceforge.net/index.shtml
.. _Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

Samtools
--------

`Samtools <http://samtools.sourceforge.net/>`_ is required for 
:program:`preprocess_reads.py`, :program:`align.py`, and :program:`align2.py`.
seriesoftubes expects :program:`samtools` to be in your PATH or in
/usr/local/samtools(-\*). 
