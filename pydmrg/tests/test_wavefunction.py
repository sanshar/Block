#!/usr/bin/env python

import os, shutil
import unittest
import numpy
import pydmrg


dmrg_tmpdir = 'dmrg_tmp/'
if not os.path.exists(dmrg_tmpdir):
    os.mkdir(dmrg_tmpdir)

with pydmrg.capture_stdout() as stdout:
    dmrg_env = pydmrg.DMRGEnv()
    dmrg_env.scratch_prefix = dmrg_tmpdir
    dmrg_env.sym = 'd2h'
    dmrg_env.outputlevel = 1
    dmrg_env.sweep_iter_schedule   = [0     ,8     ,18    ]
    dmrg_env.sweep_state_schedule  = [20    ,50    ,100   ]
    dmrg_env.sweep_qstate_schedule = [0     ,0     ,0     ]
    dmrg_env.davidson_tol_schedule = [1e-6  ,1e-8  ,1e-8  ]
    dmrg_env.noise_schedule        = [1e-6  ,1e-7  ,1e-8  ]
    dmrg_env.additional_noise      = [0     ,0     ,0     ]
    dmrg_env.maxiter = 30
    dmrg_env.update_dmrginp('FCIDUMP-1')

    pydmrg.dmrg_single(dmrg_env)


class KnowValues(unittest.TestCase):
    def test_load(self):
        wfn = pydmrg.Wavefunction()
        wfn.load(5, 7, 0, prefix=dmrg_tmpdir)
        self.assertEqual(wfn.deltaQuantum.particleNumber, 4)
        self.assertEqual(wfn.deltaQuantum.totalSpin, 0)
        self.assertEqual(wfn.stateInfo.totalStates, 3)
        #mat = wfn.stateInfo.allowedQuanta.copy()
        #mat[1,1] = mat[2,5] = mat[3,0] = mat[3,6] = False
        #self.assertTrue(not mat.any())

        diff = wfn.stateInfo.quantaStates - numpy.array([1, 1, 1],dtype=int)
        self.assertEqual(abs(diff).sum(), 0)
        diff = wfn.stateInfo.leftStateInfo.quantaStates \
                - numpy.array([1, 1, 1],dtype=int)
        self.assertEqual(abs(diff).sum(), 0)
        diff = wfn.stateInfo.rightStateInfo.quantaStates \
                - numpy.array([1, 1, 1],dtype=int)
        self.assertEqual(abs(diff).sum(), 0)

        spinquanta = wfn.stateInfo.get_quanta(0)
        self.assertEqual(spinquanta.particleNumber, 4)
        self.assertEqual(spinquanta.totalSpin, 0)
        self.assertEqual(spinquanta.irrep, 0)
        self.assertEqual(len(wfn.stateInfo.get_quantaMap(1,1)), 1)


if __name__ == "__main__":
    print "Full Tests for pydmrg load/read"
    unittest.main()

