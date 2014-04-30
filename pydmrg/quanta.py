#
# File: quanta.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import _dmrg
import quanta

class SpinQuantum(object):
    def __init__(self):
        self._raw = _dmrg.NewRawSpinQuantum()
        self.particleNumber = None
        self.totalSpin = None
        self.irrep = None

    def __str__(self):
        return '%d:%d:sym=%d' % (self.particleNumber, self.totalSpin, self.irrep)

    def refresh_by(self, rawquanta):
        assert(isinstance(rawquanta, _dmrg.RawSpinQuantum))
        self._raw = rawquanta
        self._sync_raw2self()

    def _sync_raw2self(self):
        self.particleNumber = self._raw.particleNumber
        self.totalSpin = self._raw.totalSpin
        self.irrep = self._raw.irrep()

    def _sync_self2raw(self):
        pass

    def init(self, nparticle, spin, irrep_id):
        self._raw.init(nparticle, spin, irrep_id)
        self.particleNumber = nparticle
        self.totalSpin = spin
        self.irrep = irrep_id
