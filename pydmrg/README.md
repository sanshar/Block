pydmrg
======

Python wrapper for Block code


Installation
------------

* Prerequisites
    - Cmake version 2.8 or higher
    - Python version 2.6 or higher
    - Numpy version 1.6 or higher
    - Cython version 0.20 or higher (optional, to generate _pydmrg.cpp)

* Make patches to Block
    - in spinblock.h, change 'class SpinBlock', make members public
    - in input.h, change 'class Input', make members public

* Compile Block code with flag  -fPIC  to get  libqcdmrg.so

    make libqcdmrg.so

* Compile and install pydmrg

    cd pydmrg/core
    mkdir build; cd build
    cmake ..
    make

* To make python be able to find pydmrg, add the path of Block code to
   PYTHONPATH, , e.g.

    echo 'export PYTHONPATH=/home/opt/Block:$PYTHONPATH' >> ~/.bashrc

* Use Intel MKL as BLAS library
    - add '-lmkl_avx' or '-lmkl_mc -lmkl_def' in Block makefile e.g.

        LAPACKBLAS = -L${MKLLIB} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_avx

      then compile  libqcdmrg.so  again

    - cmake with options -DBLA_VENDOR=Intel10_64lp

        cd build; cmake .. -DBLA_VENDOR=Intel10_64lp; make


Structure
---------
                        +-----------------------------------------------+
    Pydmrg layer        |   pydmrg.Wavefunction, pydmrg.SpinBlock ...   |
                        +-----------------------------------------------+
                                        A  |  A  |  A  |
                                        |  V  |  V  |  V
                        +-----------------------------------------------+
    Python C-API layer  | _dmrg.RawWavefunction, _dmrg.RawSpinBlock ... |
                        +-----------------------------------------------+
                                        A  |  A  |  A  |
                                        |  V  |  V  |  V
                        +-----------------------------------------------+
    Block code layer    |     Block.Wavefunction, Block.SpinBlock ...   |
                        +-----------------------------------------------+

* _dmrg is a lower interface layer which directly access Block code.
  It provides the most basic functions or class to represent the
  intrinsic data structure of Block.  The raw data are then wrapped in
  these xxx.py

* In _dmrg, a raw class use a pointer _this to save only one instance of
  the Block class.  So the memory can be managed by GC of python through
  raw class.

* sanity check or complicated stuff in xxx.py
