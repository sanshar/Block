.. _benchmark:


Benchmark
*********

========= ==============================
Platform
========= ==============================
CPU       4 Intel E5-2670 @ 2.6 GB
Memory    128 GB DDR3
OS        Custom Redhat 6.6
BLAS      MKL 11.0
Compiler  GCC 4.8.2
========= ==============================


Computation cost
================

Block-1.1

================================= ============= ======================== =============
   Problem size                    Program             CPU                  Memory     
================================= ============= ======================== =============
 439 AOs, CAS(16e,20o), M = 1000   DMRG-CASCI    ~ 1 h (16 core)          < 1 GB/core 
 \                                 DMRG-CASSCF   ~ 1.5 h/iter (16 core)   < 1 GB/core               
 \                                 NEVPT2        48 h (8 core)            ~12 GB/core
 \                                 MPS-NEVPT2    5.5 h (16 core)          < 4 GB/core
 439 AOs, CAS(22e,27o), M = 1000   DMRG-CASCI    6 h (16 core)            < 2 GB/core 
 \                                 DMRG-CASSCF   9 h/iter (16 core)       < 2 GB/core 
 \                                 MPS-NEVPT2    29 h (16 core)           ~10 GB/core
 760 AOs, CAS(30e,36o), M = 1500   DMRG-CASCI    24 h (16 core)           < 2 GB/core
================================= ============= ======================== =============


Block-1.5 (stackblock)

================================= ============= ======================== =============
   Problem size                    Program             CPU                  Memory     
================================= ============= ======================== =============
 439 AOs, CAS(16e,20o), M = 1000   DMRG-CASCI    ~ .5 h (16 core)         < 1 GB/core
 \                                 DMRG-CASSCF   ~ 1  h/iter (16 core)    < 1 GB/core              
 \                                 MPS-NEVPT2    3   h (16 core)          < 4 GB/core
 439 AOs, CAS(22e,27o), M = 1000   DMRG-CASCI    3 h (16 core)            < 2 GB/core
 \                                 DMRG-CASSCF   5 h/iter (16 core)       < 2 GB/core
 \                                 MPS-NEVPT2    15 h (16 core)           ~10 GB/core
 760 AOs, CAS(30e,36o), M = 1500   DMRG-CASCI    12 h (16 core)           < 2 GB/core
================================= ============= ======================== =============

