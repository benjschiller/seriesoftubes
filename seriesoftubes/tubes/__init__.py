"""
contains useful functions and classes for stringing together input and output
via pipes and converting data formats on the fly
"""
from polledpipe import PolledPipe
__all__ = ['wait_for_job', 'PolledPipe']
import time

def wait_for_job(job, logs=[], logger=None):
    otime = time.time()
    while True:
        atime = time.time()
        if otime - atime > 60 and logger is not None:
            logger.debug('Still running')
        otime = atime
        time.sleep(3)
        if job.poll() is not None: break
        for log in logs: log.log()
    for log in logs: log.log()