"""
requires genome.db from bioplus, an sqlite3 database w/ the following structure:
tables are named according to genome name
each table has fields VARCHAR 'name' and INT UNSIGNED 'length' 
  which describes the chromosomes/scaffolds in a genome
"""
from pysam import Samfile
from pkg_resources import resource_filename, cleanup_resources
import sqlite3
_CONNECTION = sqlite3.connect(resource_filename('bioplus','data/genomes.db'))
import atexit
atexit.register(cleanup_resources)
from operator import itemgetter
from tempfile import TemporaryFile

def TemporaryGenomeFile(genome_name):
    """
    returns a file-like object pointing to a temporary file containing 
    the chromsome names and sizes
    
    the current file position will be 0
    
    it will be deleted when the object is garbage collected
    """
    f = NamedTemporaryFile()
    genome_rows = genome(genome_name).iteritems()
    f.writelines(('%s\t%d\n' for chrom, size in genome_rows))
    f.flush()
    f.seek(0)
    return f
    
def _populate_available_genomes():
    name_tuples = _CONNECTION.execute("select name from sqlite_master WHERE type='table'").fetchall()
    return map(itemgetter(0), name_tuples)

AVAILABLE_GENOMES = _populate_available_genomes()

def import_db(db, replace=False, ignore_warnings=False):
    '''
    imports a db like genome.db into the existing database
    
    note: this may not be persistent if your egg is zipped
    '''
    connection = sqlite3.connect(db)
    name_tuples = connection.execute("select name from sqlite_master WHERE type='table'").fetchall()
    available_genomes = map(itemgetter(0), name_tuples)
    for genome_name in available_genomes:
        if not ignore_warnings:
            assert genome_name in AVAILABLE_GENOMES and not replace, \
                "Genome %s already in database" % genome_name
        _CONNECTION.execute("create table %s (name VARCHAR, length INT UNSIGNED)" % 
                            genome_name)
        dict(connection.execute("select * from %s" % genome_name))
        _CONNECTION.executemany("insert into %s values (?,?)" % 
                                genome_name, dict(genome_dict).items())
    return

def add_genome(genome_name, genome_dict, replace=False):
    """
    WARNING: THIS WILL PERMANENTLY ADD A GENOME
    
    make a call to pybedtool and tries to add to genome registry
    then reloads module and repopulates globals    
    note: this may not be persistent if your egg is zipped
    """
    assert genome_name in AVAILABLE_GENOMES and not replace, "Genome %s already in database" % genome_name
    _CONNECTION.execute("create table %s (name VARCHAR, length INT UNSIGNED)" % 
                        genome_name)
    _CONNECTION.executemany("insert into %s values (?,?)" % 
                            genome_name, dict(genome_dict).items())
    return


def genome(genome_name):
    """
    if available,
    returns a dict of {'chr_name': length} for all chromosomes in that genome
    """
    if not genome_name in AVAILABLE_GENOMES: raise GenomeNotAvailableError
    else:
        try:
            return dict(_CONNECTION.execute("select * from %s" % genome_name))
        except sqlite3.OperationalError:
            raise GenomeNotAvailableError(genome_name)
    
class GenomeNotAvailableError(ValueError):
    pass

class NoMatchFoundError(LookupError):
    pass

def guess_bam_genome(bam_file_or_filename):
    return guess_genome(genome_from_bam(bam_file_or_filename))

def genome_from_bam(bam_file_or_filename):
    if isinstance(bam_file_or_filename, Samfile):
        bam_file = bam_file_or_filename
    else:
        bam_file = Samfile(bam_file_or_filename)
    return dict(zip(bam_file.references, bam_file.lengths))

def guess_genome(genome1):
    """
    expects a dictionary or iterable of ('name', length)
    
    checks for a match in the database
    """
    for name in AVAILABLE_GENOMES:
        genome2 = genome(name)
        if matches_genome(genome1, genome2): return name
    raise NoMatchFoundError
    
def matches_genome(genome1, genome2, symmetric=False):
    """
    expects two dictionaries or iterables of ('name', length)
    tells you if genome1 matches genome2
    (not all the entries in genome2 need to be matched, but all the entries in
    genome1 must be)
    
    symmetric=True means matches(genome1,genome2) and matches(genome2,genome1)
    """
    g1 = dict(genome1)
    g2 = dict(genome2)
    for k in g1.keys():
        if not g2.has_key(k): return False
        elif not g2[k] == g1[k]: return False
        #else: continue
    if symmetric:
        return matches_genome(genome2, genome1)
    else:
        return True

#def guess_genome_from_dict(dict):
    