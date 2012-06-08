import select
import os
import logging
import platform
if platform.system() == 'Windows':
    raise RuntimeError("""Microsoft Windows is not and will not be supported.
This is due to an underlying design problem in the OS.
Use cygwin if you need to run this on Windows.""")

class PolledPipe(object):
    """
    A PolledPipe object has two attributes:
    
    r - pipe read file descriptor,
    w - pipe write file descriptor
    
    and three methods:
    
    poll() -- poll the read end of the pipe,
    readlines() -- read availble lines if the poll says the pipe is ready,
    log(level=logging.error) -- emit the result of readline (if not None)
    
    a select.poll() object has r registered to it
    """
    def __init__(self, logger=None, level=None):
        self._logger = logger
        self._level = level
        r, w = os.pipe()
        self.r = r
        self._r_file = os.fdopen(r, 'r', 0)
        self.w = w
        self._poll = select.poll()
        self._poll.register(r)
    
    def poll(self, timeout=0):
        return self._poll.poll(timeout)
        
    def readlines(self, timeout=0):
        results = self.poll(timeout)
        if results is None: raise StopIteration
        for result in results:
            event = result[1]
            if (event & select.POLLERR) == select.POLLERR:
                raise IOError('Something went wrong trying to read from a pipe')
            elif (event & select.POLLNVAL) == select.POLLNVAL:
                raise IOError('Invalid request: descriptor for pipe not open')
            elif (event & select.POLLHUP) == select.POLLHUP:
                raise IOError('Pipe file descriptor already hung up')
            elif (event & select.POLLIN) == select.POLLIN or \
                 (event & select.POLLPRI) == select.POLLPRI:
                yield self._r_file.readline()
        raise StopIteration
    
    def log(self, level=None):
        """
        emit the results of readlines
        """
        if self._logger is None:
            raise RuntimeWarning('PolledPipe did not have a logger attached')
        if level is None: level = self._level
        if level is None: level = logging.ERROR
        while True:
            for line in self.readlines():
                self._logger.log(level, line.rstrip())
            results = self.poll()
            if results is None: break
            if len(results) == 0: break
            event = results[0][1]
            if (event & select.POLLIN) == select.POLLIN or \
               (event & select.POLLPRI) == select.POLLPRI: continue
            else: break
        return