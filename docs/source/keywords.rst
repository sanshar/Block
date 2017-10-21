.. _keywords_list:

Keywords
********

Block program offers many keywords to control the calculation.  In the Block
input file, each keyword (and its value) should take one line.  There are two
types of format to input the keywords.  One is the key-value pattern::

    keyword value1 value2 ...

The other is the `schedule` section::

    schedule
    0 50  1.0e-6 1e-06
    5 50  1.0e-7 0e-06
    end

Hamiltonian Types
=================

Default Hamiltonian type is the chemistry Hamiltonian, in which the
2-electron interaction term is a dense matrix.  Block also supports model
Hamiltonians:

heisenberg
  No associated value.

hubbard
  No associated value.

Algorithm Types
===============

If no algorithm is specified, the default is twodot\_to\_onedot algorithm.

onedot
  No associated value.

twodot
  No associated value.

twodot\_to\_onedot
  followed by a single number to indicate the iteration when to switch from the
  two-dot to the one-dot algorithm.

Warm-up Types
=============

If ``warmup`` is not specified, the default is ``local\_0site``.

warmup
  followed by one of the following strings.

  | local\_0site
  | local\_2site
  | local\_3site
  | local\_4site
  | wilson


Solver Types
============

If not specified, default is ``davidson`` solver.

davidson
  No associated value

lanczos 
  No associated value

Orbital Reorder Types
=====================

If not specified, ``fielder`` ordering is used by default.

fiedler
  No associated value

gaopt
  followed by the filename

reorder
  followed by the filename which saves the indices of the reordered orbitals.

noreorder
  No associated value

Calculation Types
=================

dmrg
  No associated value.  This is the default calculation type.

calcoverlap
  No associated value.  If set, Compute only the overlap of two DMRG wave function.

calchamiltonian
  No associated value.

fourpdm (Not support in Block-1.5)
  No associated value.  If set, besides the regular dmrg calculation, the
  program also compute the 4-particle density matrix.

fullrestart 
  No associated value

nevpt2\_npdm
  No associated value

onepdm
  No associated value, If set, the program also computes and stores the
  1-particle density matrix in :file:`spatial_onepdm.*.*.txt`.

restart\_fourpdm
  No associated value

restart\_nevpt2\_npdm
  No associated value

restart\_onepdm
  No associated value

restart\_threepdm
  No associated value

restart\_tran\_onepdm
  No associated value

restart\_tran\_twopdm
  No associated value

restart\_twopdm
  No associated value

threepdm
  No associated value, If set, the program also computes and stores the
  3-particle density matrix in :file:`spatial_threepdm.*.*.txt`.

transition\_onepdm
  No associated value, If set, the program also computes and stores the
  transition density matrix in :file:`spatial_onepdm.*.*.txt`.

transition\_twopdm
  No associated value, If set, the program also computes and stores the
  transition 2-particle density matrix in :file:`spatial_twopdm.*.*.txt`.

twopdm
  No associated value, If set, the program also computes and stores the
  2-particle density matrix in :file:`spatial_twopdm.*.*.txt`.
  

Schedule
========

DMRG sweep schedule.  It supports two types of input format.  One is::

    schedule default

to generate the sweep schedule by Block program.  This is the default value of
``schedule`` keyword.  The other is to explicitly input the schedule, eg::

    schedule
    0 50  1.0e-6 1e-06
    5 50  1.0e-7 0e-06
    end

The first column is the sweep-iteration, the second column is the corresponding
bond dimension, the third column is the convergence tolerance for Davidson
diagonalization, the fourth column is the noise.


Expert Keywords
===============

backward
  No associated value.  If set, the program starts with backward sweep.

hf\_occ
  The initial HF wave function occupancies, in spin orbital.  It can be one of
  the following values.  The recommended one is ``integral``.

  integral
    Generate occupancy pattern based on the 1-electron integrals.
    
  canonical
    Based on the input orbital ordering, take the first N orbitals as occupied
    orbitals

  manual
    followed by the HF occupancy in the same line.

irrep
  followed by one or two numbers. If one number is given, it determines the wave
  function symmetry.  If two numbers are given, it means that calculations of
  transition density matrix between wavefunctions with different irreps.  The
  irrep value follows the Molpro symmetry convention, see also
  https://github.com/sunqm/pyscf/blob/master/future/dmrgscf/dmrg_sym.py

lastM
  followed by an integer.  Default is 500

maxiter
  followed by an integer represents the max sweep iterations. Default is 10.

maxM
  followed by an integer, the maximum bond dimension.

nelec
  followed by an integer, for the total number of electrons.

.. new\_npdm\_code

nonspinadapted 
  No associated value.

nroots
  followed by one integer, for the number of states to solve. Default is 1.

.. occ nocc
..   followed by one integer, for the number of occupied orbitals.

orbitals
  followed by the orbital file

outputlevel **0** `|| 1 || 2 || 3`
  followed by an integer, range from 0 (less output) to 3 (very noise).  Default is 0.

.. pdm\_unsorted

scratch
  followed by the scratch directory to store the intermediates, the resultant
  wave function, and density matrices.

screen\_tol
  followed by a float number.  Default is 0.0

spin `2S`
  followed by an integer represents the difference of alpha and beta electron
  numbers. 

startM
  followed by an integer.  Default is 250.

statespecific 
  No associated value.  This option implies that a previous dmrg calculation has
  already been performed.  This calculation will take the previous wavefunctions
  and refine them.

sweep\_tol
  followed by a float number, for the convergence tolerance.

sym
  followed by the string for the symbol of point group symmetry, can be d2h and
  its subgroup.

weights
  followed by a list of float numbers for weight of each state.  You could chose
  to omit the keyword weights in which case the weights will be distributed
  uniformly bet ween the different roots.

num_thrds (new in Block-1.5)
  followed by an integer for the number of OpenMP threads to use.

memory (new in Block-1.5)
  followed by an integer and a letter as the unit (k, m, g).
  The default is ``2 g`` (2 GB memory).
