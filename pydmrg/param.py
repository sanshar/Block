#
# File: param.py
# Author: Qiming Sun <qimings@princeton.edu>
#

# guesswaveTypes
BASIC = 0
TRANSFORM = 1
TRANSPOSE = 2

#enum algorithmTypes
ONEDOT = 0
TWODOT = 1
TWODOT_TO_ONEDOT = 2

#enum hamTypes
QUANTUM_CHEMISTRY = 0
HUBBARD = 1

#enum solveTypes
LANCZOS = 0
DAVIDSON = 1

#enum noiseTypes
RANDOM = 0
EXCITEDSTATE = 1

#enum calcType
DMRG           = 0
ONEPDM         = 1
TWOPDM         = 2
RESTART_TWOPDM = 3
RESTART_ONEPDM = 4
TINYCALC       = 5
FCI            = 6

#enum orbitalFormat
MOLPROFORM = 0
DMRGFORM = 1

#enum reorderType
FIEDLER   = 0
GAOPT     = 1
MANUAL    = 2
NOREORDER = 3

#
NO_PARTICLE_SPIN_NUMBER_CONSTRAINT = 0
PARTICLE_SPIN_NUMBER_CONSTRAINT = 1
LOCAL_STORAGE = 0
DISTRIBUTED_STORAGE = 1

#
LOCAL_STORAGE = 0
DISTRIBUTED_STORAGE = 1
