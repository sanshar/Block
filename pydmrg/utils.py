#
# File: utils.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import os, sys
import tempfile

class capture_stdout:
    '''redirect stdout to a string

    Examples
    --------
    with capture_stdout() as stdout:
        run_sth()
    print stdout.read()
    '''
    def __enter__(self):
        sys.stdout.flush()
        self._contents = None
        self.old_stdout_fileno = sys.stdout.fileno()
        self.bak_stdout_fd = os.dup(self.old_stdout_fileno)
        self.fd, self.ftmp = tempfile.mkstemp()
        os.dup2(self.fd, self.old_stdout_fileno)
        return self
    def __exit__(self, type, value, traceback):
        sys.stdout.flush()
        self._contents = open(self.ftmp, 'r').read()
        os.dup2(self.bak_stdout_fd, self.old_stdout_fileno)
        os.close(self.fd)
        os.remove(self.ftmp)

    def read(self):
        if self._contents:
            return self._contents
        else:
            sys.stdout.flush()
            return open(self.ftmp, 'r').read()
