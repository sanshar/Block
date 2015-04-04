<img src="https://raw.githubusercontent.com/sanshar/Block/master/README_Examples/block_logo.jpg" width="60px" height="60px" />

`BLOCK` implements the density matrix renormalization group (DMRG) algorithm for quantum chemistry.

Build Block
-----------

### 1. Complile ```BLOCK```

```BLOCK``` requires BLAS, LAPACK and BOOST.
MPI library is needed for distributed-memory parallel compilation.
```BLOCK``` is compiled using the makefile supplied in the distribution. 
The following customizations need to be made to the makefile placed in the main directory ```./Block```. 

Choose compilers by specifying,

	CXX = g++  

For MPI-based parallel execution on distributed-memory machines,
  
        USE_MPI = yes
        MPICXX = mpicxx  

MPI library must be compiled using the same compiler as for compiling ```BLOCK```. 
Intel compiler such as ```icpc``` is also supported with approriate compiling flags chosen automatically.

To enable MKL library,

        USE_MKL = yes

And supply MKL and BOOST libraries by giving the locations,
    
	MKLLIB = /opt/intel/composer_xe_2013_sp1.0.080/mkl/lib/intel64/ 
	MKLFLAGS = /opt/intel/composer_xe_2013_sp1.0.080/mkl/include

	BOOSTLIB = /lib64/boost_1_55_0/lib/
	BOOSTINCLUDE = /lib64/boost_1_55_0/include/

When the makefile is configured, run in the directory ```./Block```

        $./make

The successful compilation generates the executable ```block.spin_adapted```, static and shared DMRG libraries ```libqcdmrg.a``` and ```libqcdmrg.so```.

### 2. Test ```BLOCK```

```BLOCK``` can be tested by executing the script in the directory ```./Block/dmrg_tests```,

        $cd dmrg_tests
        $./runtest

The tests require Python to be installed on the system.

### 3. Run ```BLOCK```

The standalone serial code can be executed running

        $block.spin_adapted input.dat > output.dat

```input.dat``` is the input file and the output of the program is piped into the output file ```output.dat```.

The MPI parallel mode can be called running

        $mpirun -np 4 block.spin_adapted input.dat > output.dat

How to cite `Block`
-------------------

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


Documentation
-------------

The online documentation is available at [https://sanshar.github.io/Block](https://sanshar.github.io/Block).

