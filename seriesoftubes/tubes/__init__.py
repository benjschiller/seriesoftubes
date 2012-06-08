"""
contains useful functions and classes for stringing together input and output
via pipes and converting data formats on the fly
"""
from polledpipe import PolledPipe
from cTubes import wait_for_job
__all__ = ['PolledPipe', 'wait_for_job']
#from time import time, sleep

#def wait_for_job(job, logs=[], logger=None):
#    """wait for a job (Popen instance) to complete
#    """
#    otime = time()
#    poll = job.poll
#    while True:
#        atime = time()
#        if otime - atime > 60 and logger is not None:
#            logger.debug('Still running')
#        otime = atime
#        sleep(3)
#        if poll() is not None: break
#        for log in logs: log.log()
#    for log in logs: log.log()
#    
#def run_job(job, readable_input, logs=[], logger=None):
#    """run a job (Popen instance) to completion
#    by piping in readable_input
#    """
#    w = job.stdin.write
#    i = 0
#    for line in readable_input:
#        i += 1
#        m = i % 10**5
#        if m == 0: logger.debug('Read in %d lines', i)
#        w(line)
#        for log in logs: log.log()
#    w('\x1a')
#    otime = time()
#    while True:
#        atime = time()
#        if otime - atime > 60 and logger is not None:
#            logger.debug('Still running')
#        otime = atime
#        sleep(3)
#        if job.poll() is not None: break
#        for log in logs: log.log()
#    for log in logs: log.log()