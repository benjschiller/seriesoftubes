"""
contains useful functions and classes for stringing together input and output
via pipes and converting data formats on the fly
"""
from polledpipe import PolledPipe
import select
from select import poll
from os import fdopen
__all__ = ['wait_for_job', 'PolledPipe']
from time import time, sleep

def wait_for_job(job, logs=[], logger=None):
    """wait for a job (Popen instance) to complete
    """
    otime = time()
    poll = job.poll
    while True:
        atime = time()
        if otime - atime > 60 and logger is not None:
            logger.debug('Still running')
            otime = atime
        sleep(3)
        if poll() is not None: break
        for log in logs: log.log()
    for log in logs: log.log()
#    
#def run_job(job, readable_input, input_pipe, logs=[], logger=None):
#    """run a job (Popen instance) to completion
#    by piping in readable_input
#    """
#    cdef:
#        long i
#        int m
#        double otime, atime
#        bytes line
#        list events
#    i = 0
#    input = fdopen(input_pipe.w, 'wb')
#    for line in readable_input:
#        i += 1
#        m = i % (5 * 10 ** 5)
#        if m == 0: logger.debug('Read in %d lines', i)
#        while True:
#            events = input_pipe.poll()
#            if len(events) > 0:
#                if not (events[0][1] & select.POLLIN) == select.POLLIN: break 
#            else: break
##        logger.debug(line)
#        input.write(line)
#        for log in logs: log.log()
#    logger.debug('Read in %d lines', i)
#    logger.debug('Closing stdin/stdout/stderr')
##    input.write('\x1a')
#    input.flush()
##    job.stdin.flush()
##    job.stdin.close()
#    if job.stdout is not None:
#        job.stdout.flush()
##        job.stdout.close()
#    if job.stderr is not None:
#        job.stderr.flush()
##        job.stderr.close()
#    otime = time()
#    while True:
#        job.poll()
#        if job.returncode is not None: break
#        sleep(10)
#        logger.debug('Waiting for job to finish %d seconds', int(time() - otime))
#    for log in logs: log.log()
