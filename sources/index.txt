.. BLOCK documentation master file, created by
   sphinx-quickstart on Wed Mar 18 11:45:47 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to BLOCK's documentation!
=================================
.. figure:: images/block_logo.jpg
   :align: left
   :scale: 15%

`BLOCK` implements the density matrix renormalization group (DMRG) algorithm for quantum chemistry.
The DMRG is a variational wavefunction method. Compared to other quantum chemical methods,
it efficiently describes strong, multi-reference correlation in a large number of active orbitals (occupancies far from 0 or 2).
The method is also provably optimal for correlation with a one-dimensional topology, that is,
where orbitals are arranged with a chain- or ring-like connectivity. 
However, with the possible exception of small molecules, for correlation that is dynamic in character, 
the DMRG may be less computationally efficient than other methods 
such as coupled cluster theory or multireference configuration interaction.
We recommend the use of the DMRG in problems requiring active spaces too large for
standard complete active space (CAS) techniques. Thus, if you are interested in:

* a CAS-like treatment of low-lying eigenstates in real problems with more than 50 active orbitals,
* or, one-dimensional orbital topologies with up to 100 active orbitals,
* and, standard chemical accuracy (1 kcal/mol in energy differences),

then the DMRG may be the right method for you.

Contents
--------

.. toctree::
   :maxdepth: 2

   overview.rst
   build.rst
   with-pyscf.rst
   examples.rst
   keywords.rst
   benchmark.rst
   CHANGELOG.rst

