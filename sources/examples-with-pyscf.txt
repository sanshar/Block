DMRG for Electronic Structure Calculations
******************************************

Using PySCF
===========

`PySCF <http://chemists.princeton.edu/chan/software/pyscf/>`_ package provided
an interface to call ``BLOCK`` code for DMRG-CASSCF and DMRG-NEVPT2 calculations.
See the :ref:`installation <pyscf-itrf>` section to setup the PySCF/BLOCK interface.
In the section, we will demonstrate how to use ``BLOCK`` and ``PySCF`` packages
to study static and dynamic correlations with DMRG-CASSCF/DMRG-CASCI and
DMRG-MRPT solvers for large active space problems.

DMRG-CASSCF
-----------

We start from a simple example::

        $ cat example1.py
        from pyscf import gto, scf, dmrgscf
        mf = gto.M(atom="C 0 0 0; C 0 0 1", basis="ccpvdz").apply(scf.RHF).run()
        mc = dmrgscf.dmrgci.DMRGSCF(mf, 6, 6)
        mc.run()

Executing this script in command line::

        $ python example1.py

will start BLOCK program with 4 processors, assuming that you have the
configuration ``dmrgscf.settings.MPIPREFIX = "mpirun -n 4"``.
The number of parallel processors can be dynamically adjusted using
``sys.argv``, eg::

        $ cat example2.py
        import sys
        from pyscf import gto, scf, dmrgscf
        dmrgscf.dmrgci.settings.MPIPREFIX = "mpirun -n %s" % sys.argv[1]
        mf = gto.M(atom="C 0 0 0; C 0 0 1", basis="ccpvdz").apply(scf.RHF).run()
        mc = dmrgscf.dmrgci.DMRGSCF(mf, 6, 6)
        mc.run()

        $ python example2.py 4

In the above examples, ``gto``, ``scf`` are the standard modules provided by
PySCF package.  For the use of PySCF package, we refer the reader to the
`PySCF documentation <http://www.pyscf.org>`_.  ``dmrgscf`` module is the code
where we put Block interface.  It is designed to control all Block input
parameters, access the results from Block, including but not limiting to
regular DMRG calculation, N-particle density matrices (up to 4-PDM) and
transition density matrices (up to 2-PDM), DMRG-NEVPT2 calculations.

The standard way to start a DMRG-CASSCF calculation needs to modify the
``fcisolver`` attribute of CASSCF or CASCI object::

        from pyscf import gto, scf, mcscf, dmrgscf
        mol = gto.M(atom="C 0 0 0; C 0 0 1", basis="ccpvdz")
        mf = scf.RHF(mol)
        mf.run()
        norb = 6
        nelec = 6
        mc = mcscf.CASSCF(mf, norb, nelec)
        dmrgsolver = dmrgscf.dmrgci.DMRGCI(mol)
        dmrgsolver.maxM = 50
        dmrgsolver.maxIter = 10
        mc.fcisolver = dmrgscfsolver
        mc.run()

        mc = mcscf.CASCI(mf, norb, nelec)
        mc.fcisolver = dmrgscf.dmrgci.DMRGCI(mol)
        mc.run()

``dmrgsolver = dmrgscf.dmrgci.DMRGCI(mol)`` created an object ``dmrgsolver`` to
hold Block input parameters and runtime environments.  By default,
``maxM=1000`` is applied.   One can control the DMRG calculation by changing
the settings of ``dmrgsolver`` object, eg to set the sweep schedule::

        dmrgsolver.scheduleSweeps = [0, 4, 8, 12, 16, 20, 24, 28, 30, 34]
        dmrgsolver.scheduleMaxMs  = [200, 400, 800, 1200, 2000, 4000, 3000, 2000, 1000, 500]
        dmrgsolver.scheduleTols   = [0.0001, 0.0001, 0.0001, 0.0001, 1e-5, 1e-6, 1e-7, 1e-7, 1e-7, 1e-7]
        dmrgsolver.scheduleNoises = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0, 0.0, 0.0, 0.0]
        dmrgsolver.twodot_to_onedot = 38

For more details of the default settings and control parameters, we refer to
the `PySCF source code <https://github.com/sunqm/pyscf/blob/master/future/dmrgscf/dmrgci.py>`_
and the corresponding Block code :ref:`keywords_list` list.

To make the embedded DMRG solver work more efficiently in CASSCF iteration, one
need carefully tune the DMRG runtime parameters.  It means more input arguments
and options to be specified for ``dmrgsolver`` object.  To simplify the input,
we provided a shortcut function ``DMRGSCF`` in the ``dmrgscf`` module, as shown
by the first example ``example1.py``.  In the ``DMRGSCF`` function, we created
a CASSCF object, and assigned ``dmrgsolver`` to its ``fcisolver``, then hooked
a function which is used to dynamically adjust sweep scheduler to CASSCF object.
Now, the DMRG-CASSCF calculation can be executed in one line of input::

        mc = dmrgscf.dmrgci.DMRGSCF(mf, norb, nelec).run()

The DMRG-CASSCF results such as orbital coefficients, natural occupancy etc.
will be held in the ``mc`` object.  The DMRG wave function will be stored on
disk, more precisely, in the directory specified by ``mc.fcisolver.scratchDirectory``.
Apparently, we can modify its value to change the place to save the DMRG
wave function.  The default directory is read from the dmrgscf configuration
parameter ``dmrgscf.settings.BLOCKSCRATCHDIR``.

.. note:: Be sure the ``mc.fcisolver.scratchDirectory`` is properly assigned.
  Since all DMRGCI object by default uses the same ``BLOCKSCRATCHDIR`` settings,
  it's easy to cause name conflicts on the scratch directory, especially when
  two DMRG-CASSCF calculations are executed on the same node.

.. note:: Usually, the DMRG wave function is very large.  Be sure that the
  disk which ``BLOCKSCRATCHDIR`` pointed to has enough space.

Due to the complexity of multi-configuration model, it's common that we need
interrupt the CASSCF calculation and restart the calculation with modified
parameters.  To restart the CASSCF calculation, we need the information such as
orbital coefficients and active space CI wave function from last simulation.
Although the orbital coefficients can be save/load through
`PySCF chkfile module <https://github.com/sunqm/pyscf/blob/master/mcscf/chkfile.py>`_,  the CI wave
function are not saved by PySCF.  Unlike the regular Full CI based CASSCF
calculation in which the Full CI wave function can be fast rebuilt by a fresh
running, the restart feature of DMRG-CASSCF calculation relies on the wave
function indicated by the ``mc.fcisolver.scratchDirectory`` attribute and the
``restart`` flag of DMRG solver::

        mc = dmrgscf.dmrgci.DMRGSCF(mol)
        mc.fcisolver.scratchDirectory = "/path/to/last/dmrg/scratch"
        mc.fcisolver.restart = True
        mc.run()

.. note:: A mismatched DMRG wave function (from wrong
  ``mc.fcisolver.scratchDirectory``) may cause DMRG-CASSCF crash.

Other common features like state-average DMRG-CASSCF or state-specific for
excited state can be easily called with the ``DMRGSCF`` wrapper function::

        from pyscf import gto, scf, mcscf, dmrgscf
        mol = gto.M(atom="C 0 0 0; C 0 0 1", basis="ccpvdz")
        mf = scf.RHF(mol)
        mf.run()
        mc = dmrgscf.dmrgci.DMRGSCF(mf, 6, 6)
        # half-half average over ground state and first excited state
        mc.state_average_([0.5, 0.5])
        mc.run()

        # Optimize the first excited state
        mc.state_specific_(state=1)
        mc.run()

More information of their usage can be found in PySCF examples
`10-state_average.py <https://github.com/sunqm/pyscf/blob/master/examples/dmrg/10-state_average.py>`_
and 
`11-excited_states.py <https://github.com/sunqm/pyscf/blob/master/examples/dmrg/11-excited_states.py>`_.


DMRG-NEVPT2
-----------

DMRG-NEVPT2 calculation is straightforward if the DMRG-CASCI or DMRG-CASSCF are
finished::

        from pyscf import gto, scf, dmrgscf, mrpt
        mol = gto.M(atom="C 0 0 0; C 0 0 1", basis="ccpvdz")
        mf = scf.RHF(mol).run()

        mc = dmrgscf.dmrgci.DMRGSCF(mf, 6, 6).run()
        mrpt.nevpt2.sc_nevpt(mc)

        mc = mcscf.CASCI(mf, 6, 6)
        mc.fcisolver = dmrgscf.dmrgci.DMRGCI(mol)
        mc.run()
        mrpt.nevpt2.sc_nevpt(mc)

However, the default DMRG-NEVPT2 calculation is extremely demanding on both CPU
and memory resources.  In Block code, there is an effective approximation
implemented based on compressed MPS perturber which can significantly
reduce the computation cost::

        from pyscf import gto, scf, dmrgscf, mrpt
        mol = gto.M(atom="C 0 0 0; C 0 0 1", basis="ccpvdz")
        mf = scf.RHF(mol).run()
        mc = dmrgscf.dmrgci.DMRGSCF(mf, 6, 6).run()

        mrpt.nevpt2.sc_nevpt(dmrgscf.compress_perturb(mc))

The efficient NEVPT2 needs be initialized with ``compress_perturb`` function,
in which the most demanding intermediates are precomputed and stored on disk.

.. note:: The efficient NEVPT2 algorithm is also very demanding, especially on
  the memory usage.  Please refer to the :ref:`benchmark` for approximate cost.

If the excitation energy is of interest, we can use DMRG-NEVPT2 to compute the
energy of excited state.  Note only the state-specific NEVPT2 calculation is
available in the current Block version::

        mc = mcscf.CASCI(mf, 6, 6)
        mc.fcisolver = dmrgscf.dmrgci.DMRGCI(mol)
        mc.fcisolver.nroots = 2
        mc.kernel()
        mps_nevpt_e1 = mrpt.nevpt2.sc_nevpt(dmrgscf.compress_perturb(mc, maxM=100, root=0))
        mps_nevpt_e2 = mrpt.nevpt2.sc_nevpt(dmrgscf.compress_perturb(mc, maxM=100, root=1))

In the above example, two NEVPT2 calculations are carried out separately for
two states which are indicated by the argument ``root=*``.

For DMRG-CASSCF and DMRG-NEVPT2 calculations, there are more examples available
in `PySCF source code <https://github.com/sunqm/pyscf/tree/master/examples/dmrg>`_.


Using Molpro
============

The examples of Block installation and DMRG-SCF calculation can be found in
`Molpro online manual <https://www.molpro.net/info/2015.1/doc/manual/node385.html>`_.


.. ORCA
.. ====
.. DMRG calculation within ORCA can be found in
.. https://sites.google.com/site/orcainputlibrary/cas-calculations/dmrg .


