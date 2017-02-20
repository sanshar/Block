<img src="https://raw.githubusercontent.com/sanshar/Block/master/README_Examples/block_logo.jpg" width="60px" height="60px" />

`BLOCK` implements the density matrix renormalization group (DMRG) algorithm for quantum chemistry.

This version Block-1.1.1 is not maintained anymore.  The DMRG code has been
rewritten for better memory and computational efficiency.  The project was moved
to bitbucket and released as Block-1.5.  The precompiled binary of Block-1.5 can
be [downloaded](https://sanshar.github.io/Block/build.html#precompiled-binary)

* [New features](#new-features)
* [How to cite Block](#how-to-cite-block)
* [Project web](http://chan.caltech.edu/software/block)
* [Documentation](https://sanshar.github.io/Block)
* [License GPLv3](../master/LICENSE.txt)


New features
------------

### Version 1.5

* Optimized memory usage and efficiency
* Supported (OpenMP) threads and shared memory

### [Version 1.1 (alpha)](../../releases/latest)

* DMRG-CASSCF: close integration with [PySCF](http://chemists.princeton.edu/chan/software/pyscf/)
  (up to about 40 orbitals and 3000 basis functions)
* DMRG-NEVPT2 (up to about 24 orbitals and 1200 basis functions with PySCF)
* DMRG-NEVPT2 with compressed MPS perturber (up to about 30 orbitals and 1200
  basis functions through PySCF)
* 1, 2, 3, 4-particle density matrices
* 1, 2-particle transition density matrices

How to cite Block
-----------------

`Block` is distributed under the GNU GPL license which is reproduced in the file LICENSE.
In addition, `Block` contains a full copy of the Newmat C++ matrix library by Robert Davies.

We would appreciate if you cite the following papers in publications resulting from the
use of `Block`:

* G. K.-L. Chan and M. Head-Gordon, J. Chem. Phys. 116, 4462 (2002),
* G. K.-L. Chan, J. Chem. Phys. 120, 3172 (2004),
* D. Ghosh, J. Hachmann, T. Yanai, and G. K.-L. Chan, J. Chem. Phys., 128, 144117 (2008),
* S. Sharma and G. K-.L. Chan, J. Chem. Phys. 136, 124121 (2012).

In addition, a useful list of DMRG references relevant to quantum chemistry can be found
in the article above by Sharma and Chan.


