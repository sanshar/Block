#
# File: spinblock.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import os, sys
import copy
import _dmrg
import stateinfo
import param


class SpinBlock(object):
    def __init__(self, dmrg_env=None):
        self._env = dmrg_env
        self._raw = _dmrg.NewRawSpinBlock()
        self.braStateInfo = stateinfo.StateInfo(dmrg_env)
        self.ketStateInfo = stateinfo.StateInfo(dmrg_env)
        self.sites = []
        self.leftBlock = None
        self.rightBlock = None
        #self.loopblock = False
    @property
    def stateInfo(self):
        return self.braStateInfo

    def load(self, start_id, end_id, forward=True, left_id=-1, right_id=-1, \
             prefix=None):
        if prefix is None:
            if self._env is None:
                prefix = os.environ['TMPDIR'] + '/'
            else:
                prefix = self._env.scratch_prefix + '/'
        if forward:
            blockfile = '%sSpinBlock-forward-%d-%d.0.tmp' \
                    % (prefix, start_id, end_id)
        else:
            blockfile = '%sSpinBlock-backward-%d-%d.0.tmp' \
                    % (prefix, start_id, end_id)

        if not os.path.isfile(blockfile):
            raise OSError('file %s does not exist' % blockfile)
        rawbra, rawket = self._raw.load(blockfile)
        self.braStateInfo._raw = rawbra # use the non-mem _raw as it is one
                                        # part of class SbinBlock
        self.braStateInfo.load(start_id, end_id, forward, left_id, \
                               prefix=prefix)
        self.ketStateInfo._raw = rawket
        self.ketStateInfo.load(start_id, end_id, forward, right_id, \
                               prefix=prefix)
        self._sync_raw2self()

    def save(self, start_id, end_id, forward=True, left_id=-1, right_id=-1, \
             prefix=None):
        #TODO:self._raw.sync:
        #TODO: localstorage
        #TODO: name
        #TODO: complementary
        #TODO: hasMemoryAllocated
        #TODO: normal
        #TODO: direct
        #TODO: loopblock
        #TODO: sites
        #TODO: complementary_sites
        #TODO: braStateInfo
        #TODO: ketStateInfo
        #TODO: ops
        if prefix is None:
            if self._env is None:
                prefix = os.environ['TMPDIR'] + '/'
            else:
                prefix = self._env.scratch_prefix + '/'
        if forward:
            blockfile = '%sSpinBlock-forward-%d-%d.0.tmp' \
                    % (prefix, start_id, end_id)
        else:
            blockfile = '%sSpinBlock-backward-%d-%d.0.tmp' \
                    % (prefix, start_id, end_id)
        #self._sync_self2raw()
        self._raw.save(blockfile)
        self.braStateInfo.save(start_id, end_id, forward, left_id, \
                               prefix=prefix)
        self.ketStateInfo.save(start_id, end_id, forward, right_id, \
                               prefix=prefix)

    def _sync_raw2self(self):
        self.sites = self._raw.sites

    def _sync_self2raw(self):
        self._raw.sync(self.leftBlock._raw, self.rightBlock._raw, self.sites, \
                       self.braStateInfo._raw, self.ketStateInfo._raw)

    def init_dot(self, forward, start_id, dot_size=1, implicitTranspose=True, \
                 is_complement=False):
        # see e.g. sweep.C, BlockAndDecimate
        if forward:
            end_id = start_id + (dot_size-1)
            self.sites = range(start_id, end_id+1)
        else:
            end_id = start_id - (dot_size-1)
            self.sites = range(end_id, start_id+1)
        self._raw.init_by_dot_id(start_id, end_id, implicitTranspose, is_complement)
        self.braStateInfo.refresh_by(self._raw.get_braStateInfo())
        self.ketStateInfo.refresh_by(self._raw.get_ketStateInfo())

    def init_by_stateinfo(self, si):
        self._raw.init_by_stateinfo(si._raw)
        self.braStateInfo = si
        self.ketStateInfo = si
        self._sync_raw2self()

    def BuildSumBlock(self, constraint, lBlock, rBlock):
        self.leftBlock = lBlock
        self.rightBlock = rBlock
        #maybe initialize self._raw._this.twoInt here
        self.sites = sorted(self.leftBlock.sites + self.rightBlock.sites)
        c_sites = self.get_complementary_sites()
        self._raw.set_complementary_sites(c_sites)
        self.braStateInfo = stateinfo.TensorProduct(lBlock.braStateInfo, \
                                                    rBlock.braStateInfo, \
                                                    constraint)
        self.ketStateInfo = stateinfo.TensorProduct(lBlock.ketStateInfo, \
                                                    rBlock.ketStateInfo, \
                                                    constraint)
        if constraint != param.PARTICLE_SPIN_NUMBER_CONSTRAINT:
            self.braStateInfo = stateinfo.CollectQuanta(self.braStateInfo)
            self.ketStateInfo = stateinfo.CollectQuanta(self.ketStateInfo)

        self._sync_self2raw() # must synb before build_ops
        self._raw.set_twoInt()
        self._raw.build_ops()
        return self

    def BuildTensorProductBlock(self, sites):
        self._raw.BuildTensorProductBlock(sites)
        self._sync_raw2self()
        self.braStateInfo.refresh_by(self._raw.get_braStateInfo())
        self.ketStateInfo.refresh_by(self._raw.get_ketStateInfo())
        return self

    def transform_operators(self, rotmat):
        #old_si = self.stateInfo
        self._raw.transform_operators(rotmat._raw)
        self.braStateInfo.refresh_by(self._raw.get_braStateInfo())
        self.ketStateInfo.refresh_by(self._raw.get_ketStateInfo())
        #for i in self.stateInfo.quanta.size():
        #    assert(self.stateInfo.quanta[i] == old_si.quanta[newQuantaMap[i]])

        self.leftBlock = None
        self.rightBlock = None

    def default_op_components_compl(self, complementary=False, \
                                    implicitTranspose=True):
        self._raw.default_op_components_compl(complementary, implicitTranspose)

    def default_op_components(self, direct, sys, sysDot, haveNormops=False, \
                              haveCompops=True, implicitTranspose=True,
                              storagetype=param.LOCAL_STORAGE):
        self._raw.default_op_components(direct, sys._raw, sysDot._raw, \
                                        haveNormops, haveCompops, \
                                        implicitTranspose, storagetype)

    def print_summary(self):
        self._raw.printOperatorSummary()
        print 'self.sites', self.sites
        print 'self.braStateInfo', self.braStateInfo
        #print 'self.braStateInfo.braStateInfo', self.braStateInfo.braStateInfo
        #print 'self.braStateInfo.ketStateInfo', self.braStateInfo.ketStateInfo
        print 'self.ketStateInfo', self.ketStateInfo
        #print 'self.ketStateInfo.braStateInfo', self.ketStateInfo.braStateInfo
        #print 'self.ketStateInfo.ketStateInfo', self.ketStateInfo.ketStateInfo

    def get_complementary_sites(self):
        return [i for i in range(self._env.tot_sites) if i not in self.sites]

    def addAdditionalCompOps(self):
        self._raw.addAdditionalCompOps()

    def BuildSlaterBlock(self, si, env_sites, haveNormops):
        # lots of things initialized in BuildSlaterBlock
        _dmrg.PyBuildSlaterBlock_with_stateinfo(self._raw, si._raw, env_sites,
                                                haveNormops)
        self._sync_raw2self()
        self.braStateInfo.refresh_by(self._raw.get_braStateInfo())
        self.ketStateInfo.refresh_by(self._raw.get_ketStateInfo())
        return self

    def set_loopblock(self, tf):
        self._raw.set_loopblock(tf)


def InitStartingBlock(dmrg_env, forward=True):
    # usually molecule_quantum_tot_spin = 0
    molecule_quantum_tot_spin = dmrg_env.spin
    # usually forward_starting_size = backward_starting_size = 1
    forward_starting_size = dmrg_env.forward_starting_size
    backward_starting_size = dmrg_env.backward_starting_size

    startingBlock = SpinBlock(dmrg_env)
    if forward:
        startingBlock.init_dot(True, 0, forward_starting_size, is_complement=True)
        # dmrg_env.add_noninteracting_orbs is always True
        if dmrg_env.add_noninteracting_orbs and molecule_quantum_tot_spin != 0:
            s = quanta.SpinQuantum()
            s.init(nparticle, spin, irrep_id) # fixme, nparticle =?= spin, see initblocks.C
            addstate = stateinfo.StateInfo(dmrg_env)
            addstate.init_by_a_spinquantum(s)
            dummyblock = SpinBlock(dmrg_env)
            dummyblock.init_by_stateinfo(addstate)
            newblk = SpinBlock(dmrg_env)
            newblk.default_op_components(False, startingBlock, dummyblock, \
                                         True, True,dmrg_env.implicitTranspose)
            newblk.BuildSumBlock(param.NO_PARTICLE_SPIN_NUMBER_CONSTRAINT,
                                 startingBlock, dummyblock)
            startingBlock = newblk
    else:
        backwardSites = range(dmrg_env.tot_sites-backward_starting_size,
                              dmrg_env.tot_sites)
        startingBlock.default_op_components_compl(False)
        startingBlock.BuildTensorProductBlock(backwardSites)
    return startingBlock

# haveNormops whether construct CRE_DESCOMP
# haveCompops whether construct CRE_DESCOMP
def InitNewSystemBlock(dmrg_env, system, systemDot, haveNormops=False, haveCompops=True):
    direct = dmrg_env.direct # direct is obtained form dmrginp, true by default
    storagetype = 0 # mostly = DISTRIBUTED_STORAGE, but we use local for testing

    newsys = SpinBlock(dmrg_env)
    newsys.default_op_components(direct, system, systemDot, haveNormops,
                                 haveCompops, dmrg_env.implicitTranspose)
    newsys.BuildSumBlock(param.NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, system, systemDot)
    return newsys

def InitNewEnvironmentBlock(dmrg_env, isweep, environDot, system, systemDot, \
                            dot_with_sys, warmUp=False):
    forward = (system.sites[0] == 0)
    if forward:
        env_start_id = system.sites[-1]+dmrg_env.sys_add \
                + dmrg_env.env_add(isweep) + 1
        env_sites = range(env_start_id, dmrg_env.tot_sites)
    else:
        env_start_id = system.sites[0]-dmrg_env.sys_add \
                - dmrg_env.env_add(isweep) - 1
        env_sites = range(0, env_start_id+1)

    # usually forward_starting_size = backward_starting_size = 1
    forward_starting_size = dmrg_env.forward_starting_size
    backward_starting_size = dmrg_env.backward_starting_size
    nexact = forward_starting_size # or backward_starting_size

    environ = SpinBlock(dmrg_env)
    if dmrg_env.outputlevel > 0:
        print 'dot_with_sys, onedot, warmUp, env_sites = ', \
                dot_with_sys, dmrg_env.onedot(isweep), warmUp, env_sites
    if dot_with_sys and dmrg_env.onedot(isweep):
        newenviron = SpinBlock(dmrg_env)
        if warmUp:
            if len(env_sites) == nexact:
                newenviron.default_op_components_compl(not forward)
                newenviron.BuildTensorProductBlock(env_sites)
                newenviron.save(env_sites[0], env_sites[-1], forward=True)
            else:
                si = stateinfo.TensorProduct(system.stateInfo,
                                             systemDot.stateInfo,
                                             param.NO_PARTICLE_SPIN_NUMBER_CONSTRAINT)
                si = stateinfo.CollectQuanta(si)
                newenviron.BuildSlaterBlock(si, env_sites, False)
        else:
            newenviron.load(env_sites[0], env_sites[-1], not forward)
    else:
        haveNormops = not dot_with_sys # see initblocks.C
        if warmUp:
            if len(env_sites) == nexact:
                environ.default_op_components_compl(not forward)
                environ.BuildTensorProductBlock(env_sites)
                environ.save(env_sites[0], env_sites[-1], forward=True)
            else:
                si = stateinfo.TensorProduct(system.stateInfo,
                                             systemDot.stateInfo,
                                             param.NO_PARTICLE_SPIN_NUMBER_CONSTRAINT)
                si = stateinfo.CollectQuanta(si)
                if not dmrg_env.onedot(isweep):
                    si = stateinfo.TensorProduct(si, environDot.stateInfo,
                                                 param.NO_PARTICLE_SPIN_NUMBER_CONSTRAINT)
                    si = stateinfo.CollectQuanta(si)
                environ.BuildSlaterBlock(si, env_sites, haveNormops)
        else:
            environ.load(env_sites[0], env_sites[-1], not forward)
        environ.addAdditionalCompOps()
        # initialize newenv as that did in sys-block
        if dmrg_env.outputlevel > 0:
            print 'environ.print_summary()'
            environ.print_summary()
            print 'environDot.print_summary()'
            environDot.print_summary()
        newenviron = InitNewSystemBlock(dmrg_env, environ, environDot,
                                        haveNormops)
    return environ, newenviron

def InitBigBlock(dmrg_env, left, right):
    big = SpinBlock(dmrg_env)
    big._raw.set_big_components() # TODO: direct access self.ops
    #print 'start InitBigBlock', \
    #        left.braStateInfo.quanta.size,  left.braStateInfo.totalStates, \
    #        left.ketStateInfo.quanta.size,  left.ketStateInfo.totalStates, \
    #        right.braStateInfo.quanta.size, right.braStateInfo.totalStates,\
    #        right.ketStateInfo.quanta.size, right.ketStateInfo.totalStates
    big.BuildSumBlock(param.PARTICLE_SPIN_NUMBER_CONSTRAINT, left, right)
    return big



if __name__ == '__main__':
    block = SpinBlock()
    block.load(0, 4, True, prefix='/dev/shm/')
    print block.sites
    block.print_summary()

