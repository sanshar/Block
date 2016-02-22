.. _keywords_list:

Keywords
********

The keyword input syntax is simple::

 keyword value

The default keywords and values are **bolded**.

Hamiltonian Types
=================

* heisenberg
* hubbard
* **quantum\_chemistry**

Algorithm Types
===============

* onedot
* twodot
* **twodot\_to\_onedot**

Warm-up Types
=============

* warmup **local\_0site** `|| local\_2site || local\_3site || local\_4site || wilson`

Solver Types
============

* **davidson**
* lanczos 

Orbital Reorder Types
=====================

* **fiedler**
* gaopt `default`
* reorder  `reorder file`
* noreorder

Calculation Types
=================

* backward
* calchamiltonian
* calcoverlap
* **dmrg**
* fci
* fourpdm
* fullrestart 
* nevpt2\_npdm
* onepdm
* restart\_fourpdm
* restart\_nevpt2\_npdm
* restart\_onepdm
* restart\_threepdm
* restart\_tran\_onepdm
* restart\_tran\_twopdm
* restart\_twopdm
* threepdm
* transition\_onepdm
* transition\_twopdm
* twopdm


Expert Keywords
===============

* hf\_occ `integral || orbital || manual`
* irrep `isym`
* lastM **500** `||  lastM`
* maxiter **10** `|| max sweep iterations`
* maxM `maxM`
* nelec `nelec`
* new\_npdm\_code
* nonspinadapted 
* nroots **1** `|| nroots`
* occ nocc
* orbitals `orbital file`
* outputlevel **0** `|| 1 || 2 || 3`
* pdm\_unsorted
* **schedule default**
* schedule `sweep_iteration M davidson_tolerance noise` end
* scratch **current directory of input file** `|| scratch directory`
* screen\_tol **0.0** `|| ScreenTol`
* spin `2S`
* startM **250** `|| startM`
* statespecific 
* sweep\_tol **1.0e-5** `|| SweepTol`
* sym `point group`
* weights **1.0** `||` `W`\ :sub:`1`, `W`\ :sub:`2`, ..., `W`\ :sub:`nroots`
