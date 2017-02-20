Overview
********

Features
========

New features in Block 1.5 (stackblock)
--------------------------------------

* Optimized memory usage and efficiency
* Supported (OpenMP) threads and shared memory

New features in Block 1.1
-------------------------

* Perturbation methods including NEVPT2 and MPSPT
* One-, two-, three- and four-particle density matrices
* One- and two-particle transition density matrices between two states

Features in Block 1.0
---------------------

* DMRG sweep algorithm for quantum chemistry, Hubbard and Heisenberg hamiltonians
* Full spin-adaptation (SU(2) symmetry) and Abelian point-group symmetries
* State-averaged and state-specific excited states
* DMRG-SCF and/or DMRG-NEVPT2 interfaces to the PySCF, Molpro, ORCA, Q-Chem and Molcas program packages


Downloads
=========
* Block-1.5 (stackblock)

  - Source code `block-1.5.0.gz <http://www.sunqm.net/pyscf/files/src/block-1.5.0.gz>`_.

  - Source code `block-1.5.0-serial.gz <http://www.sunqm.net/pyscf/files/src/block-1.5.0-serial.gz>`_.
    This serial version is used by DMRG-NEVPT2 compressed perturber method.

  - Precompiled binary `block.spin_adapted-1.5.0.gz <http://www.sunqm.net/pyscf/files/bin/block.spin_adapted-1.5.0.gz>`_
    (+ OpenMPI + Boost-1.55 + MKL-11) for Linux x86_64.

  - Precompiled binary `block.spin_adapted-1.5.0-serial.gz <http://www.sunqm.net/pyscf/files/bin/block.spin_adapted-1.5.0-serial.gz>`_
    (+ Boost-1.55 + MKL-11) for Linux x86_64.  This serial version is used by
    DMRG-NEVPT2 compressed perturber method.

* Block-1.1.1

  - Precompiled binary `block.spin_adapted-1.1.1.gz <http://www.sunqm.net/pyscf/files/bin/block.spin_adapted-1.1.1.gz>`_
    (+ OpenMPI + Boost-1.55 + MKL-11) for Linux x86_64

  - Precompiled binary `block.spin_adapted-1.1.1-serial.gz <http://www.sunqm.net/pyscf/files/bin/block.spin_adapted-1.1.1-serial.gz>`_
    (+ Boost-1.55 + MKL-11) for Linux x86_64.  This binary is used by
    DMRG-NEVPT2 compressed perturber method.


Use Block with PySCF package
============================

Block program supports two executing modes: running standalone through command
line or as a plugin of other quantum chemistry package.  The Python-based
quantum chemistry program package `PySCF <http://www.pyscf.org>`_ provides a
simple solution to run Block program.  It is the recommended way to use
Block program in most scenario.  Please see the userguide :ref:`dmrg_pyscf`.

.. Calling `BLOCK` as an external function
.. =======================================
.. 
.. The makefile distributed with `Block` code can be used to generate a library file called
.. libqcdmrg.a. 
.. To call `Block` as a subroutine from a C++ program, the library file has to be
.. linked to the program. 
.. A DMRG calculation can be performed using the function call ``calldmrg(inputf, outputf)``,
.. where ``inputf`` and ``outputf`` are C-style character arrays specifying the `Block` input and output fies respectively.

License and how to cite
=======================

`Block` is distributed under the GNU GPL license which is reproduced in the file LICENSE.
In addition, `Block` contains a full copy of the Newmat C++ matrix library by Robert Davies.

We would appreciate if you cite the following papers in publications resulting from the
use of `Block`:

* G. K.-L. Chan and M. Head-Gordon, J. Chem. Phys. 116, 4462 (2002),
* G. K.-L. Chan, J. Chem. Phys. 120, 3172 (2004),
* D. Ghosh, J. Hachmann, T. Yanai, and G. K.-L. Chan, J. Chem. Phys., 128, 144117 (2008),
* S. Sharma and G. K-.L. Chan, J. Chem. Phys. 136, 124121 (2012),
* R. Olivares-Amaya, W. Hu, N. Nakatani, S. Sharma, J. Yang and G. K.-L. Chan, J. Chem. Phys. 142, 034102 (2015).

In addition, a useful list of DMRG references relevant to quantum chemistry can be found
in the article above by Sharma and Chan. 

