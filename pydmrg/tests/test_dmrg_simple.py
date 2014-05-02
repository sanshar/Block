#
# File: dmrg.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import os, shutil
import unittest
import pydmrg


class KnowValues(unittest.TestCase):
    def test_dmrg_simple(self):
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
            e = pydmrg.dmrg_single(dmrg_env)

        shutil.rmtree(dmrg_tmpdir)

        self.assertAlmostEqual(e, -17.10516373008, 10)

if __name__ == "__main__":
    print "Full Tests for pydmrg scf"
    unittest.main()



