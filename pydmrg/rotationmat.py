#
# File: quanta.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import os, sys
import _dmrg
import wavefunction

class RotationMatrix(object):
    def __init__(self, dmrg_env=None):
        self._env = dmrg_env
        self.size = None

    def load(self, start_id, end_id, state_id=0, prefix=None):
        if prefix is None:
            if self._env is None:
                prefix = os.environ['TMPDIR'] + '/'
            else:
                prefix = self._env.scratch_prefix + '/'
        matfile = '%sRotation-%d-%d.0.state%d.tmp' \
                % (prefix, start_id, end_id, state_id)
        if not os.path.isfile(matfile):
            raise OSError('file %s does not exist' % matfile)
        self._raw = _dmrg.NewRawRotationMatrix()
        self._raw.load(matfile)
        self._sync_raw2self()

    def save(self, start_id, end_id, state_id=0, prefix=None):
        if prefix is None:
            if self._env is None:
                prefix = os.environ['TMPDIR'] + '/'
            else:
                prefix = self._env.scratch_prefix + '/'
        matfile = '%sRotation-%d-%d.0.state%d.tmp' \
                % (prefix, start_id, end_id, state_id)
        self._raw.save(matfile)

    def refresh_by(self, rawmat):
        assert(isinstance(rawmat, _dmrg.NewRawRotationMatrix))
        self._raw = rawmat
        self._sync_raw2self()

    def _sync_raw2self(self):
        self.size = self._raw.get_size()

    def get_matrix_by_quanta_id(self, quanta_id):
        '''2D numpy array'''
        assert(quanta_id < self.size)
        return self._raw.get_matrix_by_quanta_id(quanta_id)


def update_rotmat(dmrg_env, wfn, sys, big, keep_states, noise):
    rmat = RotationMatrix(dmrg_env)
    rmat.refresh_by(_dmrg.Pyupdate_rotmat(wfn._raw, sys._raw, big._raw,
                                          keep_states, 0, noise))
    return rmat


if __name__ == '__main__':
    rotmats = RotationMatrix()
    rotmats.load(0, 5, prefix='/dev/shm/')
    for i in range(8):
        mat = rotmats.get_matrix_by_quanta_id(i)
        print mat.shape
