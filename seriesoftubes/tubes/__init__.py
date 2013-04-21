"""
contains useful functions and classes for stringing together input and output
via pipes and converting data formats on the fly
"""
from polledpipe import PolledPipe
from time import time, sleep
__all__ = ['wait_for_job', 'PolledPipe']


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
        if poll() is not None:
            break
        for log in logs:
            log.log()
    for log in logs:
        log.log()
