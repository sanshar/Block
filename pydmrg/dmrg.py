#
# File: dmrg.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import os, sys
import tempfile
import sweep
import param
import _dmrg


class DMRGEnv(object):
    def __init__(self, **keywords):
        self.outputlevel = 0
        self.scratch_prefix = '/tmp'
        self.sys_add = 1
        self.tot_sites = 1
        self.sym = 'c1' # need interface to IrrepSpace
        self.nelec = -1
        self.spin = -1
        self.spinAdapted = True
        self.hf_occupancy = []
        self.hf_occ_user = 'integral'
        self.m_Bogoliubov = False
        self.sweep_iter_schedule   = [0   ]  # these lists
        self.sweep_state_schedule  = [2   ]  #  ...
        self.davidson_tol_schedule = [1e-5]  #  ...
        self.noise_schedule        = [1e-5]  # should have same lenth
        self.maxM = 10000
        self.m_lastM = 500
        self.m_startM = 500
#        self.integral_disk_storage_thresh = 100
#        self.do_diis = False
#        self.diis_error = 1e-2
#        self.start_diis_iter = 8
#        self.diis_keep_states = 6
#        self.diis_error_tol = 1e-8
        self.calc_type = param.DMRG
        self.noise_type = param.RANDOM
        self.ham_type = param.QUANTUM_CHEMISTRY
        self.nroots = 1
        self.weights = [1]
        self.solve_type = param.DAVIDSON
        self.do_fci = False
        self.do_npdm_ops = False
        self.do_npdm_in_core = False
        self.new_npdm_code = False
#        self.set_Sz = False
        self.maxiter = 10
        self.oneindex_screen_tol = 1e-15
        self.twoindex_screen_tol = 1e-15
#        self.no_transform = False
        self.add_noninteracting_orbs = True
        self.nquanta = 2
        self.deflation_min_size = 2
        self.deflation_max_size = 20
        self.algorithm_type = param.TWODOT_TO_ONEDOT
        self.onedot_start_cycle = 20 # effective only if algorithm_type == TWODOT_TO_ONEDOT
        self.direct = True
#        self.orbenergies = []
#        self.m_maxj = 15 for 9j CG coefficients
        self.max_lanczos_dimension = 5000
        self.sweep_tol = 1e-5
#        self.spin_vector = []
        self.implicitTranspose = True
        self.spin_orbs_symmetry = []
        self.spatial_to_spin = []
        self.spin_to_spatial = []
        self.core_energy = 0.
        self.orbformat = param.MOLPROFORM
        self.reorderType = param.FIEDLER
        self.reorderfile = ''
        self.reorder = []
        self.gaconffile = ''
#        self.molecule_quantum = None

        # they are from sweepParams
        self.forward_starting_size = 1
        self.backward_starting_size = 1
        self.max_blk_cyc = 0


    # call me after setting dmrg input parameters to update dmrginp in Block
    def update_dmrginp(self, fcidump):
        # TODO: parse and add more inputs

        line1 = open(fcidump, 'r').readline()[:-1].upper().split(',')
        for dat in line1:
            if 'NELEC' in dat:
                nelec = int(dat.split('=')[1])
            elif 'MS2' in dat:
                spin = int(dat.split('=')[1])
        if self.nelec == -1:
            self.nelec = nelec
        if self.spin == -1:
            self.spin = spin

        tmpinp = tempfile.mktemp()
        finp = open(tmpinp, 'w')
        finp.write('nelec %d\n' % self.nelec)
        finp.write('spin %d\n' % self.spin)
        #finp.write('irrep %d\n' % self.irrep)
        finp.write('hf_occ %s\n' % self.hf_occ_user)
        finp.write('outputlevel %d\n' % self.outputlevel)
        finp.write('schedule\n')
        for i,k in enumerate(self.sweep_iter_schedule):
            finp.write('%d %d  %g %g\n' % (k, self.sweep_state_schedule[i],
                                           self.davidson_tol_schedule[i],
                                           self.noise_schedule[i]))
        finp.write('end\n')
        finp.write('twodot_to_onedot %d\n' % self.onedot_start_cycle)
        finp.write('maxiter %d\n' % self.maxiter)
        finp.write('sweep_tol %g\n' % self.sweep_tol)
        finp.write('sym %s\n' % self.sym)
        finp.write('orbitals %s\n' % fcidump)
        finp.write('scratch %s\n' % self.scratch_prefix)
        finp.close()
        _dmrg.Pyinitialize_defaults(tmpinp)
        os.remove(tmpinp)

        _dmrg.Pysync_from_dmrginp(self)

    def build(self, fcidump):
        self.update_dmrginp(fcidump)


    def sweep_schedule(self, sweep_iter):
        i = 0
        i_ls = 0
        for k,v in enumerate(self.sweep_iter_schedule):
            if sweep_iter >= v:
                i = k
            if sweep_iter+1 >= v:
                i_ls = k
        return self.sweep_state_schedule[i], self.sweep_state_schedule[i_ls], \
                self.davidson_tol_schedule[i], self.noise_schedule[i],

    def max_block_cycle(self, sweep_iter=0):
        n_iters = (self.tot_sites - 2*self.forward_starting_size \
                   - self.sys_add - self.env_add(sweep_iter)) / self.sys_add \
                + 1
        return n_iters

    def env_add(self, sweep_iter=0):
        if self.algorithm_type == param.ONEDOT:
            return 0
        elif self.algorithm_type == param.TWODOT_TO_ONEDOT \
                and sweep_iter >= self.onedot_start_cycle:
            return 0
        else:
            return 1

    def onedot(self, sweep_iter=0):
        if self.algorithm_type == param.ONEDOT:
            return True
        elif self.algorithm_type == param.TWODOT_TO_ONEDOT \
                and sweep_iter >= self.onedot_start_cycle:
            return True
        else:
            return False

    def guesstype(self, warmUp, isweep, iblkcyc):
        if warmUp or (isweep < 2 and iblkcyc == 0):
            return param.BASIC
        else:
            if iblkcyc == 0:
                if self.algorithm_type != param.TWODOT_TO_ONEDOT:
                    return param.TRANSPOSE
                elif isweep != self.onedot_start_cycle:
                    return param.TRANSPOSE
                else:
                    return param.BASIC
            else:
                return param.TRANSFORM

    def sync2dmrginp(self):
        _dmrg.Pysync2dmrginp(self)

def parse(dmrg_conf_file):
    dmrg_env = DMRGEnv()
    #TODO: parse the dmrg.conf file and set relevant keys to the dmrg_env
    return dmrg_env

def _find_index(test, lst):
    for i,j in enumerate(lst):
        if test(j):
            return i
    return None


def dmrg_single(dmrg_env):
    eforward = sweep.do_one(dmrg_env, 0, forward=True, warmUp=True)
    ebackward = 0
    for isweep in range(1, dmrg_env.maxiter, 2):
        old_ef = eforward
        old_eb = ebackward
        ebackward = sweep.do_one(dmrg_env, isweep, forward=False, warmUp=False)
        #TODO: extapolate energy

        eforward = sweep.do_one(dmrg_env, isweep+1, forward=True, warmUp=False)
        if abs(eforward-old_ef) < dmrg_env.sweep_tol \
           and abs(ebackward-old_eb) < dmrg_env.sweep_tol \
           and ((dmrg_env.algorithm_type != param.TWODOT_TO_ONEDOT) \
                or (isweep+1 > dmrg_env.onedot_start_cycle)):
            break
    return eforward


if __name__ == '__main__':
    import os, shutil
    dmrg_tmpdir = './dmrg_tmp'
    if not os.path.exists(dmrg_tmpdir):
        os.mkdir(dmrg_tmpdir)

    dmrg_env = DMRGEnv()
    dmrg_env.nelec = 4
    dmrg_env.scratch_prefix = dmrg_tmpdir
    dmrg_env.sym = 'd2h'
    dmrg_env.outputlevel = 0
    dmrg_env.sweep_iter_schedule   = [0     ,8     ,18    ]
    dmrg_env.sweep_state_schedule  = [20    ,50    ,100   ]
    dmrg_env.sweep_qstate_schedule = [0     ,0     ,0     ]
    dmrg_env.davidson_tol_schedule = [1e-6  ,1e-8  ,1e-8  ]
    dmrg_env.noise_schedule        = [1e-6  ,1e-7  ,1e-8  ]
    dmrg_env.additional_noise      = [0     ,0     ,0     ]
    dmrg_env.maxiter = 30
    dmrg_env.update_dmrginp('tests/FCIDUMP-1')

    print 'E=', dmrg_single(dmrg_env)

    import wavefunction
    wfn = wavefunction.Wavefunction()
    wfn.load(2, 7, prefix=dmrg_tmpdir+'/')
    print wfn

    shutil.rmtree(dmrg_tmpdir)
