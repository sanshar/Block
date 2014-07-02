#cython: boundscheck=False
#cython: wraparound=False
#distutils: language = c++

import os
#from libc.stdlib cimport malloc, free
from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.string cimport string
#from cpython cimport bool
from libc.string cimport memcpy
import numpy
cimport numpy
cimport cython


#cdef extern from 'Python.h':
#    char *PyString_AsString(object)

cdef extern from 'config.h':
    pass
cdef extern from 'MatrixBLAS.h':
    pass

cdef extern from 'newmat.h':
    cdef cppclass Matrix:
        int Nrows()
        int Ncols()
        double& element(int, int)

cdef extern from "BaseOperator.h" namespace 'SpinAdapted':
    cdef cppclass SparseMatrix:
        vector[int]& get_orbs()
        char& allowed(int i, int j)
        int nrows()
        int ncols()
        int get_sign()
        bint get_fermion()
        vector[SpinQuantum]& set_deltaQuantum()

cdef extern from 'wavefunction.h' namespace 'SpinAdapted':
    # Wavefunction class belongs to SparseMatrix, so the member vars of
    # SparseMatrix need to be tracked
    #   ObjectMatrix<Matrix> operatorMatrix;
    cdef cppclass Wavefunction(SparseMatrix):
        Wavefunction()
        Wavefunction(Wavefunction& wfn)
        bool& get_onedot()

cdef extern from 'IrrepSpace.h' namespace 'SpinAdapted':
    cdef cppclass IrrepSpace:
        IrrepSpace(int ir)

cdef extern from 'SpinSpace.h' namespace 'SpinAdapted':
    cdef cppclass SpinSpace:
        SpinSpace(int ir)
        int getirrep()

cdef extern from 'SpinQuantum.h' namespace 'SpinAdapted':
    cdef cppclass SpinQuantum:
        int particleNumber
        SpinSpace totalSpin
        # IrrepSpace orbitalSymmetry;
        SpinQuantum()
        SpinQuantum(int p, SpinSpace s, IrrepSpace orbS)


cdef extern from 'StateInfo.h' namespace 'SpinAdapted':
    # StateInfo class holds
    # some flags hasCollectedQuanta hasPreviousStateInfo
    # oldToNewState (maybe no use now)
    cdef cppclass StateInfo:
        bool hasAllocatedMemory
        bool initialised
        StateInfo()
        StateInfo(int n, SpinQuantum *q, const int *qS)
        StateInfo *leftStateInfo
        StateInfo *rightStateInfo
        int totalStates
        vector[int] quantaStates
        vector[SpinQuantum] quanta
        vector[int] newQuantaMap
        # allowedQuanta => get_StateInfo_allowedQuanta
        # quantaMap => get_StateInfo_quantaMap
        vector[int] leftUnMapQuanta
        vector[int] rightUnMapQuanta
        StateInfo unCollectedStateInfo
        void AllocateUnCollectedStateInfo()
        void quanta_distribution(vector[SpinQuantum]& qnumbers, vector[int]& distribution,
                                 bool complement)

cdef extern from 'spinblock.h' namespace 'SpinAdapted':
    # SpinBlock class holds
    # a list of ops
    # some flags complementary normal loopblock localstorage
    #            hasMemoryAllocated direct
    # name??
    cdef cppclass SpinBlock:
        SpinBlock()
        SpinBlock(int start, int finish, bool implicitTranspose, bool is_complement)
        SpinBlock(StateInfo& s)
        SpinBlock* leftBlock
        SpinBlock* rightBlock
        StateInfo braStateInfo
        StateInfo ketStateInfo
        vector[int] sites
        vector[int] complementary_sites
        vector[int]& get_sites()
        # complementary_sites = [all i not in sites], only op_component.C uses
        # it to search the sites which are connected to complementary by ops
        void printOperatorSummary()
        #TODO BaseOperator or SparseMatrix get_op_array() for ops
        void default_op_components(bool direct, SpinBlock& lBlock,
                                   SpinBlock& rBlock, bool haveNormops,
                                   bool haveCompops, bool implicitTranspose)
        void default_op_components(bool complementary, bool implicitTranspose)
        void build_iterators()
        #void build_operators(std::vector<Csf>& s, std::vector<<std::vector<Csf> >& ladders)
        void build_operators()
        void setstoragetype(int)
        #StateInfo& get_stateInfo() by x_SpinBlock_stateInfo
        void addAdditionalCompOps() #TODO: direct access ops
        void set_big_components() #TODO: direct access ops
        void transform_operators(vector[Matrix]& rotateMatrix)
        void BuildTensorProductBlock(vector[int]& new_sites)
        void set_loopblock(bool p_loopblock)
    cdef enum Storagetype:
        LOCAL_STORAGE, DISTRIBUTED_STORAGE
            
cdef extern from 'sweep_params.h' namespace 'SpinAdapted':
    cdef enum guessWaveTypes:
        BASIC, TRANSFORM, TRANSPOSE

cdef extern from *:
    # although TensorProduct is defined in namespace SpinAdapted, it can only
    # be correctly found by gcc in the global namespace
    # tests on clang is required
    void TensorProduct(StateInfo& a, StateInfo& b, StateInfo& c,
                       int constraint, StateInfo* compState)
cdef extern from 'solver.h' namespace 'SpinAdapted::Solver':
    void solve_wavefunction(vector[Wavefunction]& solution, vector[double]& energies,
                            SpinBlock& big, double tol, int guesswavetype,
                            bool& onedot, bool& dot_with_sys, bool& warmUp,
                            double additional_noise, int currentRoot,
                            vector[Wavefunction]& lowerStates)

cdef extern from 'guess_wavefunction.h' namespace 'SpinAdapted::GuessWave':
    void onedot_shufflesysdot(StateInfo& guessstateinfo, StateInfo& transposestateinfo,
                              Wavefunction& guesswf, Wavefunction& transposewf)

cdef extern from "input.h" namespace "SpinAdapted":
    cdef enum hamTypes: QUANTUM_CHEMISTRY, HUBBARD
    cdef enum solveTypes: LANCZOS, DAVIDSON
    cdef enum algorithmTypes: ONEDOT, TWODOT, TWODOT_TO_ONEDOT
    cdef enum noiseTypes: RANDOM, EXCITEDSTATE
    cdef enum calcType: DMRG, ONEPDM, TWOPDM, RESTART_TWOPDM, RESTART_ONEPDM, TINYCALC, FCI
    cdef enum orbitalFormat: MOLPROFORM, DMRGFORM

    cdef cppclass Input:
        #vector[int] m_thrds_per_node
        int m_norbs
        int m_alpha
        int m_beta
        int m_Sz
        #IrrepSpace m_total_symmetry_number
        SpinQuantum m_molecule_quantum
        bool m_spinAdapted
        bool m_Bogoliubov
        int m_total_spin
        int m_guess_permutations
#        bool m_stateSpecific
        vector[int] m_hf_occupancy
        string m_hf_occ_user
        vector[double] m_weights
        vector[int] m_sweep_iter_schedule
        vector[int] m_sweep_state_schedule
        vector[double] m_sweep_tol_schedule
        vector[double] m_sweep_noise_schedule
        #bool m_schedule_type_default
        #bool m_schedule_type_backward
        int m_lastM
        int m_startM
        int m_maxM
        int m_integral_disk_storage_thresh
        bool m_do_diis
        double m_diis_error
        int m_start_diis_iter
        int m_diis_keep_states
        double m_diis_error_tol
        calcType m_calc_type
        noiseTypes m_noise_type
        hamTypes m_ham_type
        int m_nroots
        solveTypes m_solve_type
        #bool m_do_deriv
        bool m_do_fci
        bool m_do_npdm_ops
        bool m_do_npdm_in_core
        bool m_new_npdm_code
        bool m_set_Sz
        int m_maxiter
        double m_oneindex_screen_tol
        double m_twoindex_screen_tol
        bool m_no_transform
        bool m_add_noninteracting_orbs
        int m_nquanta
        int m_sys_add
        int m_env_add
        int m_deflation_min_size
        int m_deflation_max_size
        algorithmTypes m_algorithm_type
        int m_twodot_to_onedot_iter
        #std::vector< std::map<SpinQuantum, int> > m_quantaToKeep
        string m_save_prefix
        string m_load_prefix
        bool m_direct
        vector[double] m_orbenergies
        int m_maxj
        int m_max_lanczos_dimension
        double m_sweep_tol
        bool m_restart
        bool m_fullrestart
        bool m_restart_warm
        bool m_reset_iterations
        bool m_implicitTranspose
        vector[int] m_spin_vector
        vector[int] m_spin_orbs_symmetry
        int m_num_spatial_orbs
        vector[int] m_spatial_to_spin
        vector[int] m_spin_to_spatial
        int m_outputlevel
        double m_core_energy
        orbitalFormat m_orbformat
        int m_reorderType
        string m_reorderfile
        vector[int] m_reorder
        string m_gaconffile

    Input dmrginp


cdef extern from 'itrf.h':
    int save_rotmat(char *filerotmat, vector[Matrix] *mat)
    int load_rotmat(char *filerotmat, vector[Matrix] *mat)
    int update_rotmat(vector[Matrix] *rotateMatrix,
                      Wavefunction *wfn, SpinBlock *sys, SpinBlock *big,
                      int keptstates, int keptqstates, double noise)
    int guess_rotmat(vector[Matrix] *rotateMatrix, SpinBlock *newSystem,
                     int keptstates)

    #void initialize_default_dmrginp(char *fcidump, string prefix, string sym)
    void init_dmrginp(char *conf)
    int get_last_site_id()
    void assign_deref_shared_ptr[T](T& dest, T& src)

    int save_wavefunction(char *filewave, Wavefunction *oldWave,
                          StateInfo *waveInfo)
    int load_wavefunction(char *filewave, Wavefunction *oldWave,
                          StateInfo *waveInfo)
    int x_SpinQuantum_irrep(SpinQuantum *sq)

    int save_spinblock(char *filespinblock, SpinBlock *b)
    int load_spinblock(char *filespinblock, SpinBlock *b)
    #StateInfo *x_SpinBlock_stateInfo(SpinBlock *b)
    #vector[int] *x_SpinBlock_complementary_sites(SpinBlock *b)
    void BuildSlaterBlock_with_stateinfo(SpinBlock& environ, StateInfo& si,
                                         vector[int]& envSites, bool haveNormops)
    #void set_SpinBlock_for_BuildSumBlock(SpinBlock *self, SpinBlock *lblock,
    #                                     SpinBlock *rblock, vector[int]& sites,
    #                                     StateInfo *si)
    void set_SpinBlock_twoInt(SpinBlock *self)

    int save_stateinfo(char *filesi, StateInfo *si)
    int load_stateinfo(char *filesi, StateInfo *si)
    vector[int] *x_StateInfo_quantaMap(StateInfo *s, int lquanta_id,
                                       int rquanta_id)
    char *x_StateInfo_allowedQuanta(StateInfo *s, int lquanta_id,
                                    int rquanta_id)
    int get_whole_StateInfo_allowedQuanta(StateInfo *s, char *tftab)
    void union_StateInfo_quanta(StateInfo *a, StateInfo *b)



# RawAclass does not allocate memory for Aclass._this.  Its pointer _this only
# refer to other objects.  It does not deallocate Aclass._this when the object
# is destroyed by python.
cdef class RawSpinQuantum:
    cdef SpinQuantum *_this
    property particleNumber:
        def __get__(self): return self._this.particleNumber
        #def __set__(self, x): self._this.particleNumber = x
    property totalSpin:
        def __get__(self): return self._this.totalSpin.getirrep()
        #def __set__(self, x): self._this.totalSpin = x
    def irrep(self): return x_SpinQuantum_irrep(self._this)
cdef class NewRawSpinQuantum(RawSpinQuantum):
    def __cinit__(self):
        self._this = new SpinQuantum()
    def __dealloc__(self):
        del self._this
    def init(self, nparticle, spin, irrep_id):
        del self._this
        cdef IrrepSpace *irrep = new IrrepSpace(irrep_id)
        self._this = new SpinQuantum(nparticle, SpinSpace(spin), irrep[0])
        del irrep


cdef class RawStateInfo:
    cdef StateInfo *_this
    property totalStates:
        def __get__(self): return self._this.totalStates
        #def __set__(self, x): self._this.totalStates = x
    property quantaStates:
        def __get__(self): return self._this.quantaStates
        #def __set__(self, x): self._this.quantaStates = x
    property newQuantaMap:
        def __get__(self): return self._this.newQuantaMap
    def get_quanta(self, i):
        cdef SpinQuantum *p = &self._this.quanta[i]
        rawq = RawSpinQuantum()
        rawq._this = p
        return rawq
    def get_quantaMap(self, lquanta_id, rquanta_id):
        cdef vector[int] *qmap = x_StateInfo_quantaMap(self._this, lquanta_id,
                                                       rquanta_id)
        return qmap[0]
    def get_allowedQuanta(self, lquanta_id, rquanta_id):
        cdef char *a = x_StateInfo_allowedQuanta(self._this, lquanta_id,
                                                 rquanta_id)
        return a[0]
    def get_whole_allowedQuanta(self):
        nrow = self._this.leftStateInfo.quanta.size()
        ncol = self._this.leftStateInfo.quanta.size()
        cdef numpy.ndarray tftab = numpy.zeros((nrow,ncol),dtype=numpy.bool8)
        get_whole_StateInfo_allowedQuanta(self._this, <char *>tftab.data)
        return tftab
    property leftUnMapQuanta:
        def __get__(self): return self._this.leftUnMapQuanta
    property rightUnMapQuanta:
        def __get__(self): return self._this.rightUnMapQuanta
    property leftStateInfo:
        def __get__(self):
            s = RawStateInfo()
            s._this = self._this.leftStateInfo
            return s
    property rightStateInfo:
        def __get__(self):
            s = RawStateInfo()
            s._this = self._this.rightStateInfo
            return s
    def save(self, filesi):
        save_stateinfo(filesi, self._this)
    def load(self, filesi):
        load_stateinfo(filesi, self._this)
cdef class NewRawStateInfo(RawStateInfo):
    def __cinit__(self):
        self._this = new StateInfo()
    def __dealloc__(self):
        del self._this
    def init_by_a_spinquantum(self, RawSpinQuantum sq):
        del self._this
        cdef int quantaStates = 1
        self._this = new StateInfo(1, sq._this, &quantaStates)
    def set_unCollectedStateInfo(self, RawStateInfo old):
        self._this.initialised = True
        self._this.AllocateUnCollectedStateInfo()
        self._this.leftStateInfo    = old._this.leftStateInfo
        self._this.rightStateInfo   = old._this.rightStateInfo
        self._this.leftUnMapQuanta  = old._this.leftUnMapQuanta
        self._this.rightUnMapQuanta = old._this.rightUnMapQuanta
        assign_deref_shared_ptr(self._this.unCollectedStateInfo, old._this[0])


cdef class RawSpinBlock:
    cdef SpinBlock *_this
    def get_braStateInfo(self):
        si = RawStateInfo()
        #si._this = x_SpinBlock_stateInfo(self._this)
        si._this = &self._this.braStateInfo
        return si
    def get_ketStateInfo(self):
        si = RawStateInfo()
        si._this = &self._this.ketStateInfo
        return si
    def get_stateInfo(self):
        return self.get_braStateInfo()
    #def get_sites(self): return self._this.get_sites()
    property sites:
        def __get__(self): return self._this.sites
    def printOperatorSummary(self):
        self._this.printOperatorSummary()
    def save(self, filespinblock):
        save_spinblock(filespinblock, self._this)
cdef class NewRawSpinBlock(RawSpinBlock):
    def __cinit__(self):
        self._this = new SpinBlock()
    def __dealloc__(self):
        del self._this
    def load(self, filespinblock):
        bra = RawStateInfo()
        bra._this = &self._this.braStateInfo
        ket = RawStateInfo()
        ket._this = &self._this.ketStateInfo
        load_spinblock(filespinblock, self._this)
        return bra, ket
    def init_by_dot_id(self, int start, int finish, implicitTranspose, is_complement=0):
        del self._this
        # FIXME: SpinBlock(start,finish) calls dmrginp
        self._this = new SpinBlock(start, finish, implicitTranspose, is_complement)
    def init_by_stateinfo(self, RawStateInfo si):
        del self._this
        self._this = new SpinBlock(si._this[0])
    def BuildTensorProductBlock(self, sites):
        self._this.BuildTensorProductBlock(sites)
    def default_op_components(self, direct,
                              RawSpinBlock lBlock, RawSpinBlock rBlock,
                              haveNormops, haveCompops,
                              implicitTranspose, storagetype):
        # storage type can be one of 0 = LOCAL_STORAGE, 1 = DISTRIBUTED_STORAGE
        self._this.default_op_components(direct,
                                         lBlock._this[0], rBlock._this[0],
                                         haveNormops, haveCompops,
                                         implicitTranspose)
        cdef Storagetype t = storagetype
        self._this.setstoragetype(t)
    def default_op_components_compl(self, complementary, implicitTranspose):
        self._this.default_op_components(complementary, implicitTranspose)
    def set_complementary_sites(self, vector[int] c_sites):
        #cdef vector[int] *csites = x_SpinBlock_complementary_sites(self._this)
        #for i in c_sites:
        #    csites[0].push_back(i)
        #self._this.complementary_sites.clear()
        #for i in c_sites:
        #    self._this.complementary_sites.push_back(i)
        self._this.complementary_sites = c_sites
    def set_twoInt(self):
        set_SpinBlock_twoInt(self._this)
    def build_ops(self): # TODO add csf for overloaded build_operators
        self._this.build_iterators()
        self._this.build_operators()
    def addAdditionalCompOps(self):
        self._this.addAdditionalCompOps()
    def set_big_components(self):
        self._this.set_big_components()
    def transform_operators(self, RawRotationMatrix rotatemat):
        self._this.transform_operators(rotatemat._this[0])
    def sync(self, RawSpinBlock lblock, RawSpinBlock rblock,
             vector[int] sites, RawStateInfo bra, RawStateInfo ket):
        #set_SpinBlock_for_BuildSumBlock(self._this, lblock._this, rblock._this,
        #                                sites, si._this)
        self._this.leftBlock = lblock._this
        self._this.rightBlock = rblock._this
        self._this.sites = sites
        self._this.braStateInfo = bra._this[0]
        self._this.ketStateInfo = ket._this[0]
    def set_loopblock(self, tf):
        self._this.set_loopblock(tf)


cdef class RawSparseMatrix:
    cdef SparseMatrix *_this
    def get_orbs(self): return self._this.get_orbs()
    def get_sign(self): return self._this.get_sign()
    def get_fermion(self): return self._this.get_fermion()
    def allowed(self, i, j): return <bint>self._this.allowed(i,j)
    def get_shape(self): return self._this.nrows(), self._this.ncols()

cdef class RawWavefunction:
    cdef Wavefunction *_this
    #cdef readonly NewRawStateInfo stateInfo
    #cdef public NewRawStateInfo stateInfo
    def get_deltaQuantum(self):
        #cdef SpinQuantum *p = &(self._this.set_deltaQuantum())
        #deltaQuantum = RawSpinQuantum()
        #deltaQuantum._this = p
        deltaQuantum = RawSpinQuantum()
        deltaQuantum._this = &(self._this.set_deltaQuantum()[0])
        return deltaQuantum
    def get_onedot(self): return self._this.get_onedot()
    def get_orbs(self): return self._this.get_orbs()
    def get_sign(self): return self._this.get_sign()
    def get_fermion(self): return self._this.get_fermion()
    def allowed(self, i, j): return <bint>self._this.allowed(i,j)
    def get_shape(self): return self._this.nrows(), self._this.ncols()
    def save(self, wfnfile, RawStateInfo stateInfo):
        save_wavefunction(wfnfile, self._this, stateInfo._this)
cdef class NewRawWavefunction(RawWavefunction):
    def __cinit__(self):
        self._this = new Wavefunction()
    def __dealloc__(self):
        del self._this
        #TODO: call stateInfo.Free to release leftStateInfo, rightStateInfo, ...
    def load(self, wfnfile):
        #self.stateInfo = NewRawStateInfo()
        stateInfo = NewRawStateInfo()
        stateInfo._this.hasAllocatedMemory = True
        left = NewRawStateInfo()
        leftleft = NewRawStateInfo()
        leftright = NewRawStateInfo()
        right = NewRawStateInfo()
        rightleft = NewRawStateInfo()
        rightright = NewRawStateInfo()
        left._this.leftStateInfo = leftleft._this
        left._this.rightStateInfo = leftright._this
        stateInfo._this.leftStateInfo = left._this
        right._this.leftStateInfo = rightleft._this
        right._this.rightStateInfo = rightright._this
        stateInfo._this.rightStateInfo = right._this
        load_wavefunction(wfnfile, self._this, stateInfo._this)
        return stateInfo, left, leftleft, leftright, \
                right, rightleft, rightright


cdef class RawMatrix:
    cdef Matrix *_this
    def get_shape(self):
        return self._this.Nrows(), self._this.Ncols()

cdef class RawRotationMatrix:
    cdef vector[Matrix] *_this
    def get_matrix_by_quanta_id(self, quanta_id):
        #mat = RawMatrix()
        #mat._this = &self._this.at(qid) # bug: vague return type?
        #mat.update_allprop()
        cdef Matrix *mati = &(self._this.at(quanta_id))
        cdef int nrow = mati.Nrows()
        cdef int ncol = mati.Ncols()
        cdef numpy.ndarray mat = numpy.empty((nrow,ncol))
        if nrow*ncol > 0:
            memcpy(<double *>mat.data, &mati.element(0,0), nrow*ncol*sizeof(int))
        return mat
    def get_size(self): return self._this.size()
    def save(self, filerotmat):
        save_rotmat(filerotmat, self._this)
cdef class NewRawRotationMatrix(RawRotationMatrix):
    def __cinit__(self):
        self._this = new vector[Matrix]()
    def __dealloc__(self):
        del self._this
    def load(self, filerotmat):
        load_rotmat(filerotmat, self._this)



#################################################
#
#################################################

def PyTensorProduct(RawStateInfo a, RawStateInfo b, int constraint):
    c = NewRawStateInfo()
    # constraint = 0 for NO_PARTICLE_SPIN_NUMBER_CONSTRAINT
    # constraint = 1 for PARTICLE_SPIN_NUMBER_CONSTRAINT
    # I didn't find any call with compState other than NULL
    TensorProduct(a._this[0], b._this[0], c._this[0], constraint, NULL)
    return c

def Pyupdate_rotmat(RawWavefunction wfn, RawSpinBlock sys, RawSpinBlock big,
                    keep_states, keep_qstates, noise):
    # rmat is resized in update_rotmat => makeRotateMatrix => assign_matrix_by_dm
    rmat = NewRawRotationMatrix()

    # TODO: add noise
    update_rotmat(rmat._this, wfn._this, sys._this, big._this,
                  keep_states, keep_qstates, noise)
    return rmat

def Pyunion_StateInfo_quanta(RawStateInfo dest, RawStateInfo source):
    union_StateInfo_quanta(dest._this, source._this)
    return dest

def PyBuildSlaterBlock_with_stateinfo(RawSpinBlock environ, RawStateInfo si,
                                      envSites, haveNormops):
    BuildSlaterBlock_with_stateinfo(environ._this[0], si._this[0], envSites,
                                    haveNormops)

def Pysolve_wavefunction(RawSpinBlock big, nroots, dot_with_sys, warmUp,
                         onedot, tol, guesstype, additional_noise):
    cdef vector[Wavefunction] solution
    solution.resize(nroots)
    cdef vector[double] energies
    energies.resize(nroots)
    cdef guessWaveTypes gt = guesstype
    cdef int currentRoot = -1
    cdef vector[Wavefunction] lowerStates
    solve_wavefunction(solution, energies, big._this[0], tol, gt,
                       onedot, dot_with_sys, warmUp, additional_noise,
                       currentRoot, lowerStates)
    wfn = NewRawWavefunction()
    wfn._this[0] = solution[0]
    return wfn, energies[0]

def Pyonedot_shufflesysdot(RawStateInfo sguess, RawStateInfo stranspose,
                           RawWavefunction wfguess):
    wftranspose = NewRawWavefunction()
    wftranspose._this[0] = wfguess._this[0]
    onedot_shufflesysdot(sguess._this[0], stranspose._this[0],
                         wfguess._this[0], wftranspose._this[0])
    return wftranspose

def Pyguess_rotmat(RawSpinBlock newsys, keep_states):
    rotmat = NewRawRotationMatrix()
    guess_rotmat(rotmat._this, newsys._this, keep_states)
    return rotmat

#def Pyinitialize_defaults(fcidump, prefix, sym):
#    initialize_default_dmrginp(fcidump, prefix, sym)
def Pyinitialize_defaults(inp_conf):
    init_dmrginp(inp_conf)

def Pyget_last_site_id():
    return get_last_site_id()



def Pysync2dmrginp(dmrgenv):
    #dmrginp.m_norbs =
    dmrginp.m_alpha = (dmrgenv.nelec + dmrgenv.spin) / 2
    dmrginp.m_beta = (dmrgenv.nelec - dmrgenv.spin) / 2
    #dmrginp.m_Sz =
    #dmrginp.m_total_spin = dmrgenv.spin
    #dmrginp.m_guess_permutations =
    dmrginp.m_hf_occupancy = dmrgenv.hf_occupancy
    dmrginp.m_hf_occ_user = dmrgenv.hf_occ_user
    dmrginp.m_spinAdapted = dmrgenv.spinAdapted
    dmrginp.m_Bogoliubov = dmrgenv.Bogoliubov
    dmrginp.m_weights = dmrgenv.weights
    dmrginp.m_sweep_iter_schedule = dmrgenv.sweep_iter_schedule
    dmrginp.m_sweep_state_schedule = dmrgenv.sweep_state_schedule
    dmrginp.m_sweep_tol_schedule = dmrgenv.davidson_tol_schedule
    dmrginp.m_sweep_noise_schedule = dmrgenv.noise_schedule
    #dmrginp.m_schedule_type_default =
    #dmrginp.m_schedule_type_backward =
    dmrginp.m_lastM = dmrgenv.lastM
    dmrginp.m_startM = dmrgenv.startM
    dmrginp.m_maxM = dmrgenv.maxM
#    dmrginp.m_integral_disk_storage_thresh = dmrgenv.integral_disk_storage_thresh
#    dmrginp.m_do_diis = dmrgenv.do_diis
#    dmrginp.m_diis_error = dmrgenv.diis_error
#    dmrginp.m_start_diis_iter = dmrgenv.start_diis_iter
#    dmrginp.m_diis_keep_states = dmrgenv.diis_keep_states
#    dmrginp.m_diis_error_tol = dmrgenv.diis_error_tol
    dmrginp.m_calc_type = dmrgenv.calc_type
    dmrginp.m_noise_type = dmrgenv.noise_type
    dmrginp.m_ham_type = dmrgenv.ham_type
    dmrginp.m_nroots = dmrgenv.nroots
    dmrginp.m_solve_type = dmrgenv.solve_type
    dmrginp.m_do_fci = dmrgenv.do_fci
    dmrginp.m_do_npdm_ops     = dmrgenv.do_npdm_ops
    dmrginp.m_do_npdm_in_core = dmrgenv.do_npdm_in_core
    dmrginp.m_new_npdm_code   = dmrgenv.new_npdm_code
#    dmrginp.m_set_Sz = dmrgenv.set_Sz
    dmrginp.m_maxiter = dmrgenv.maxiter
    dmrginp.m_oneindex_screen_tol = dmrgenv.oneindex_screen_tol
    dmrginp.m_twoindex_screen_tol = dmrgenv.twoindex_screen_tol
#    dmrginp.m_no_transform = dmrgenv.no_transform
    dmrginp.m_add_noninteracting_orbs = dmrgenv.add_noninteracting_orbs
    dmrginp.m_nquanta = dmrgenv.nquanta
    dmrginp.m_sys_add = dmrgenv.sys_add
    #dmrginp.m_env_add =
    dmrginp.m_deflation_min_size = dmrgenv.deflation_min_size
    dmrginp.m_deflation_max_size = dmrgenv.deflation_max_size
    dmrginp.m_algorithm_type = dmrgenv.algorithm_type
    dmrginp.m_twodot_to_onedot_iter = dmrgenv.onedot_start_cycle
    dmrginp.m_save_prefix = dmrgenv.scratch_prefix
    dmrginp.m_load_prefix = dmrgenv.cratch_prefix
    dmrginp.m_direct = dmrgenv.direct
#    dmrginp.m_orbenergies = dmrgenv.orbenergies
#    dmrginp.m_maxj = dmrgenv.m_maxj
    dmrginp.m_max_lanczos_dimension = dmrgenv.max_lanczos_dimension
    dmrginp.m_sweep_tol = dmrgenv.sweep_tol
    #dmrginp.m_restart =
    #dmrginp.m_fullrestart =
    #dmrginp.m_restart_warm =
    #dmrginp.m_reset_iterations =
#    dmrginp.m_spin_vector = dmrgenv.spin_vector
    dmrginp.m_implicitTranspose = dmrgenv.implicitTranspose
    dmrginp.m_spin_orbs_symmetry = dmrgenv.spin_orbs_symmetry
    dmrginp.m_num_spatial_orbs = dmrgenv.tot_sites
    dmrginp.m_spatial_to_spin = dmrgenv.spatial_to_spin
    dmrginp.m_spin_to_spatial = dmrgenv.spin_to_spatial
    dmrginp.m_outputlevel = dmrgenv.outputlevel
    dmrginp.m_core_energy = dmrgenv.core_energy
    dmrginp.m_orbformat   = dmrgenv.orbformat
    dmrginp.m_reorderType = dmrgenv.reorderType
    dmrginp.m_reorderfile = dmrgenv.reorderfile
    dmrginp.m_gaconffile = dmrgenv.gaconffile

def Pysync_from_dmrginp(dmrgenv):
    dmrgenv.nelec = dmrginp.m_alpha + dmrginp.m_beta
    dmrgenv.spin  = dmrginp.m_alpha - dmrginp.m_beta
    dmrgenv.spinAdapted = dmrginp.m_spinAdapted
    dmrgenv.hf_occupancy = dmrginp.m_hf_occupancy
    dmrgenv.hf_occ_user = dmrginp.m_hf_occ_user
    dmrgenv.weights = dmrginp.m_weights
    dmrgenv.Bogoliubov = dmrginp.m_Bogoliubov
    dmrgenv.sweep_iter_schedule = dmrginp.m_sweep_iter_schedule
    dmrgenv.sweep_state_schedule = dmrginp.m_sweep_state_schedule
    dmrgenv.davidson_tol_schedule = dmrginp.m_sweep_tol_schedule
    dmrgenv.noise_schedule = dmrginp.m_sweep_noise_schedule
    dmrgenv.maxM = dmrginp.m_maxM
    dmrgenv.lastM = dmrginp.m_lastM
    dmrgenv.startM = dmrginp.m_startM
#    dmrgenv.integral_disk_storage_thresh = dmrginp.m_integral_disk_storage_thresh
#    dmrgenv.do_diis = dmrginp.m_do_diis
#    dmrgenv.diis_error = dmrginp.m_diis_error
#    dmrgenv.start_diis_iter = dmrginp.m_start_diis_iter
#    dmrgenv.diis_keep_states = dmrginp.m_diis_keep_states
#    dmrgenv.diis_error_tol = dmrginp.m_diis_error_tol
    dmrgenv.calc_type = dmrginp.m_calc_type
    dmrgenv.noise_type = dmrginp.m_noise_type
    dmrgenv.ham_type = dmrginp.m_ham_type
    dmrgenv.nroots = dmrginp.m_nroots
    dmrgenv.solve_type = dmrginp.m_solve_type
    dmrgenv.do_fci = dmrginp.m_do_fci
    dmrgenv.do_npdm_ops     = dmrginp.m_do_npdm_ops
    dmrgenv.do_npdm_in_core = dmrginp.m_do_npdm_in_core
    dmrgenv.new_npdm_code   = dmrginp.m_new_npdm_code
#    dmrgenv.set_Sz = dmrginp.m_set_Sz
    dmrgenv.maxiter = dmrginp.m_maxiter
    dmrgenv.oneindex_screen_tol = dmrginp.m_oneindex_screen_tol
    dmrgenv.twoindex_screen_tol = dmrginp.m_twoindex_screen_tol
#    dmrgenv.no_transform = dmrginp.m_no_transform
    dmrgenv.add_noninteracting_orbs = dmrginp.m_add_noninteracting_orbs
    dmrgenv.nquanta = dmrginp.m_nquanta
    dmrgenv.deflation_min_size = dmrginp.m_deflation_min_size
    dmrgenv.deflation_max_size = dmrginp.m_deflation_max_size
    dmrgenv.algorithm_type = dmrginp.m_algorithm_type
    dmrgenv.onedot_start_cycle = dmrginp.m_twodot_to_onedot_iter
    dmrgenv.scratch_prefix = dmrginp.m_save_prefix
    #dmrgenv.load_prefix = dmrginp.m_load_prefix
    dmrgenv.direct = dmrginp.m_direct
#    dmrgenv.orbenergies = dmrginp.m_orbenergies
#    dmrgenv.m_maxj = dmrginp.m_maxj
    dmrgenv.max_lanczos_dimension = dmrginp.m_max_lanczos_dimension
    dmrgenv.sweep_tol = dmrginp.m_sweep_tol
#    dmrgenv.spin_vector = dmrginp.m_spin_vector
    dmrgenv.spin_orbs_symmetry = dmrginp.m_spin_orbs_symmetry
    dmrgenv.tot_sites = dmrginp.m_num_spatial_orbs
    dmrgenv.spatial_to_spin = dmrginp.m_spatial_to_spin
    dmrgenv.spin_to_spatial = dmrginp.m_spin_to_spatial
    dmrgenv.outputlevel = dmrginp.m_outputlevel
    dmrgenv.core_energy = dmrginp.m_core_energy
    dmrgenv.orbformat = dmrginp.m_orbformat
    dmrgenv.reorderType = dmrginp.m_reorderType
    dmrgenv.reorder = dmrginp.m_reorder
    dmrgenv.reorderfile = dmrginp.m_reorderfile
    dmrgenv.gaconffile = dmrginp.m_gaconffile

