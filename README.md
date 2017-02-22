<img src="https://raw.githubusercontent.com/sanshar/Block/master/README_Examples/block_logo.jpg" width="60px" height="60px" />

`BLOCK` implements the density matrix renormalization group (DMRG) algorithm for quantum chemistry.

This version Block-1.1.1 is not maintained anymore.  The DMRG code has been
rewritten for better memory and computational efficiency.  The project was moved
to bitbucket and released as Block-1.5.  Source code and manual of Block-1.5 can
be found [here](https://sanshar.github.io/Block).

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

* DMRG-CASSCF: close integration with PySCF (up to about 40 orbitals and 3000
  basis functions)
* DMRG-NEVPT2 (up to about 24 orbitals and 1200 basis functions with PySCF)
* DMRG-NEVPT2 with compressed MPS perturber (up to about 30 orbitals and 1200
  basis functions through PySCF)
* 1, 2, 3, 4-particle density matrices
* 1, 2-particle transition density matrices


FAQ
---

### Should I update to Block-1.5 (stackblock)?
We highly recommend you to move to Block-1.5 for many reasons:

* Bugfix.  We fixed many bugs in Block-1.5.  Calculations, especially NEVPT2,
  are more stable in Block-1.5.

* Performance.  Block-1.5 can be 1.5 to 5 times faster than Block-1.1.  The speed up
  is significant for small systems.

* Maintenance.  Block-1.5 is the long-term project that we'll invest energy to
  maintain.  We'll not spend time to fix problems that only exist in Block-1.1.

### Where can I get Block-1.5?
You can download the source code or the binary executable from online
[manual](https://sanshar.github.io/Block/build.html).  If you are using DMRG
for chemistry problems, we recommend to use Block-1.5 with its interface
in [PySCF](https://github.com/sunqm/pyscf) program.  The Block programs in
Molpro, ORCA and QChem are Block-1.1.  We are working to port Block-1.5
program to these quantum chemistry packages.

### Can I use the input of Block-1.1 in Block-1.5?
Yes.  Input files (FCIDUMP and dmrg.conf) of Block-1.1 can be used in Block-1.5
without any changes.  Block-1.5 provides new keywords num_thrds and memory to
control the number of threads and the total memory to use.
```
num_thrds 8
memory, 40, g
```
Note these keywords are not recognized in Block-1.1 and they will cause
Block-1.1 crashing.

### Is wave-function format compatible in block-1.1 and block-1.5?
No.  Block-1.1 and Block-1.5 have different data format for wave-function.
You're not able to use wfn files of Block-1.1 to restart calculation in
Block-1.5 and vice versa.

### Why does my Block-1.5 crash for multiple threads while it works all right with Block-1.1 and single-thread Block-1.5?
It is mostly due to the memory size you specified in the input.  In Block-1.5,
operators are held in memory within the limits you specified in the input (or
2GB/proc).  The memory is shared by all threads of the process.  If you specify
many threads (num_thrds keyword) in your calculation, you should increas
the memory size accordingly.


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
* S. Guo, M. A. Watson, W. Hu, Q. Sun, G. K.-L. Chan, J. Chem.  Theory Comput. 12, 1583 (2016)

In addition, a useful list of DMRG references relevant to quantum chemistry can be found
in the article above by Sharma and Chan.


