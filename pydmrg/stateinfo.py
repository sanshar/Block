#
# File: stateinfo.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import os
import numpy
import _dmrg
import quanta

class StateInfo(object):
    def __init__(self, dmrg_env=None):
        self._env = dmrg_env
        self._raw = _dmrg.NewRawStateInfo()
        self.totalStates = 0
        self.quantaStates = []
        self.quanta = None
        self.allowedQuanta = numpy.ones((0,0),dtype=bool)
        self.unBlockedIndex = []
        self.leftUnMapQuanta = []
        self.rightUnMapQuanta = []

        # point to another StateInfo, (not RawStateInfo)
        self.leftStateInfo = None
        self.rightStateInfo = None
        self.unCollectedStateInfo = None

    def load(self, start_id, end_id, forward=True, state_id=-1, prefix=None):
        if prefix is None:
            if self._env is None:
                prefix = os.environ['TMPDIR'] + '/'
            else:
                prefix = self._env.scratch_prefix + '/'
        if forward:
            sifile = '%sStateInfo-forward-%d-%d.%d.0.tmp' \
                    % (prefix, start_id, end_id, state_id)
        else:
            sifile = '%sStateInfo-backward-%d-%d.%d.0.tmp' \
                    % (prefix, start_id, end_id, state_id)
        self._raw.load(sifile)
        self._sync_raw2self()

    def save(self, start_id, end_id, forward=True, state_id=-1, prefix=None):
        if prefix is None:
            if self._env is None:
                prefix = os.environ['TMPDIR'] + '/'
            else:
                prefix = self._env.scratch_prefix + '/'
        if forward:
            sifile = '%sStateInfo-forward-%d-%d.%d.0.tmp' \
                    % (prefix, start_id, end_id, state_id)
        else:
            sifile = '%sStateInfo-backward-%d-%d.%d.0.tmp' \
                    % (prefix, start_id, end_id, state_id)
        #self._sync_self2raw()
        self._raw.save(sifile)

    def refresh_by(self, rawstateinfo, with_left_right=False):
        assert(isinstance(rawstateinfo, _dmrg.RawStateInfo))
        self._raw = rawstateinfo
        self._sync_raw2self()
        if with_left_right:
            self.leftStateInfo = StateInfo(self._env)
            self.leftStateInfo._raw = self._raw.leftStateInfo
            self.leftStateInfo._sync_raw2self()
            self.rightStateInfo = StateInfo(self._env)
            self.rightStateInfo._raw = self._raw.rightStateInfo
            self.rightStateInfo._sync_raw2self()
        #FIXME: carefully sync leftStateInfo etc.
        #self.allowedQuanta = self._raw.get_whole_allowedQuanta()
        #self.leftUnMapQuanta = self._raw.leftUnMapQuanta
        #self.rightUnMapQuanta = self._raw.leftUnMapQuanta
        # rightStateInfo.leftStateInfo, rightStateInfo.rightStateInfo

        #FIXME: how to initialize?
        # self.unCollectedStateInfo self.leftStateInfo self.rightStateInfo

        #FIXME: I didn't find anywhere to use previousStateInfo

    def init_by_spinquantum(self, n, sq, qs_lst):
        if n > 1 or len(qs_lst) > 1:
            raise ValueError('TODO: initialize StateInfo')
        self._raw = _dmrg.NewRawStateInfo()
        self._raw.init_by_a_spinquantum(sq._raw)

    def init_by_spinquantum1(self, n, sq, qs_lst):
        if n > 1 or len(qs_lst) > 1:
            raise ValueError('TODO: initialize StateInfo')
        self._raw = _dmrg.NewRawStateInfo()
        self._raw.init_by_a_spinquantum(sq._raw)

    def _sync_raw2self(self):
        self.totalStates = sum(self._raw.quantaStates)
        self.quantaStates = self._raw.quantaStates
        self.quanta = _SpinQuantumList(self)
        self.unBlockedIndex = [sum(self.quantaStates[:i]) \
                               for i,_ in enumerate(self.quantaStates)]
        #if self.leftStateInfo is not None \
        #   and self.rightStateInfo is not None:
        #    self.allowedQuanta = self._raw.get_whole_allowedQuanta()

    def _sync_self2raw(self):
        qlist = [self.get_quanta(i) for i in range(self.quanta.size)]
        self._raw.sync(self.totalStates, self.quantaStates, qlist, \
                       self.leftStateInfo._raw, self.rightStateInfo._raw, \
                       self.unCollectedStateInfo._raw)

    def get_quanta(self, quanta_id):
        assert(quanta_id < len(self.quantaStates))
        spinquanta = quanta.SpinQuantum()
        spinquanta.refresh_by(self._raw.get_quanta(quanta_id))
        return spinquanta

    def get_quantaMap(self, lquanta_id, rquanta_id):
        '''return a list'''
        # TODO: test me
        #assert(lquanta_id < len(self.leftStateInfo.quantaStates))
        #assert(rquanta_id < len(self.rightStateInfo.quantaStates))
        return self._raw.get_quantaMap(lquanta_id, rquanta_id)

    def get_allowedQuanta(self, lquanta_id, rquanta_id):
        '''return True/False'''
        # TODO: test me
        #assert(lquanta_id < len(self.leftStateInfo.quantaStates))
        #assert(rquanta_id < len(self.rightStateInfo.quantaStates))
        return self._raw.get_allowedQuanta(lquanta_id, rquanta_id)

    def __str__(self):
        s = ['totalStates: %d' % self.totalStates,
             'quantaStates: %s' % self.quantaStates,
             #'quanta.size: %d' % self.quanta.size,
            ]
        return '\n'.join(s)

class _SpinQuantumList(object):
    def __init__(self, super_stateinfo):
        self.size = len(super_stateinfo.quantaStates)
        self._super_stateinfo = super_stateinfo
    def __getitem__(self, quanta_id):
        return self._super_stateinfo.get_quanta(quanta_id)

def CollectQuanta(old_stateinfo):
    # make self as an unCollectedStateInfo for a new StateInfo
    newsi = StateInfo(old_stateinfo._env)
    _dmrg.Pyunion_StateInfo_quanta(newsi._raw, old_stateinfo._raw)
    newsi.unCollectedStateInfo = old_stateinfo
    newsi._raw.set_unCollectedStateInfo(old_stateinfo._raw)

    # NOTE avoid the following four being overwritten by _sync_raw2self
    newsi.leftStateInfo = old_stateinfo.leftStateInfo
    newsi.rightStateInfo = old_stateinfo.rightStateInfo
    newsi.leftUnMapQuanta = old_stateinfo.leftUnMapQuanta
    newsi.rightUnMapQuanta = old_stateinfo.rightUnMapQuanta
    newsi._sync_raw2self()
    return newsi


def TensorProduct(a, b, constraint=0):
    '''return StateInfo c, c = (StateInfo a) * (StateInfo b)

    constraint = 0 for NO_PARTICLE_SPIN_NUMBER_CONSTRAINT
    constraint = 1 for PARTICLE_SPIN_NUMBER_CONSTRAINT
    '''
    c = StateInfo(a._env)
    c.refresh_by(_dmrg.PyTensorProduct(a._raw, b._raw, constraint), True)
    return c

if __name__ == '__main__':
    si = StateInfo()
    si.load(0, 4, True, -1, prefix='/dev/shm/')
    print si.totalStates

