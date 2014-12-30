/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_INPUT_HEADER_H
#define SPIN_INPUT_HEADER_H
#include <vector>
#include <string>
#include <map>
#include <boost/serialization/serialization.hpp>
#include <boost/shared_array.hpp>
#include "IrrepSpace.h"
#include "SpinQuantum.h"
#include "timer.h"
#include "couplingCoeffs.h"
#include <boost/tr1/unordered_map.hpp>


namespace SpinAdapted{
class SpinBlock;
class OneElectronArray;
class TwoElectronArray;
class PairArray;
class CCCCArray;
class CCCDArray;

enum WarmUpTypes {WILSON, LOCAL0, LOCAL2, LOCAL3, LOCAL4};
enum hamTypes {QUANTUM_CHEMISTRY, HUBBARD, BCS, HEISENBERG};
enum solveTypes {LANCZOS, DAVIDSON, CONJUGATE_GRADIENT};
enum algorithmTypes {ONEDOT, TWODOT, TWODOT_TO_ONEDOT};
enum noiseTypes {RANDOM, EXCITEDSTATE};
enum calcType {DMRG, ONEPDM, TWOPDM, THREEPDM, FOURPDM, NEVPT2PDM, RESTART_TWOPDM,
               RESTART_ONEPDM, RESTART_THREEPDM, RESTART_FOURPDM, TINYCALC, FCI,
               EXCITEDDMRG, CALCOVERLAP, CALCHAMILTONIAN, COMPRESS, RESPONSE,
               TRANSITION_ONEPDM, TRANSITION_TWOPDM, RESTART_T_ONEPDM, RESTART_T_TWOPDM,
               NEVPT2,RESTART_NEVPT2};
enum orbitalFormat{MOLPROFORM, DMRGFORM};
enum reorderType{FIEDLER, GAOPT, MANUAL, NOREORDER};
enum keywords{ORBS, LASTM, STARTM, MAXM,  REORDER, HF_OCC, SCHEDULE, SYM, NELECS, SPIN, IRREP,
	      MAXJ, PREFIX, NROOTS, DOCD, DEFLATION_MAX_SIZE, MAXITER, BASENERGY,
	      SCREEN_TOL, ODOT, SWEEP_TOL, OUTPUTLEVEL, NONSPINADAPTED, BOGOLIUBOV, TWODOT_NOISE, WARMUP, NPDM_INTERMEDIATE, NPDM_NO_INTERMEDIATE, NPDM_MULTINODE, NPDM_NO_MULTINODE, NUMKEYWORDS};

class Input {

 private:
  std::vector<int> m_thrds_per_node;
  int m_norbs;
  int m_alpha;
  int m_beta;
  int m_Sz;
  bool m_spinAdapted;
  bool m_Bogoliubov;

  IrrepSpace m_total_symmetry_number;
  IrrepSpace m_bra_symmetry_number;// This is used when bra and ket have different spatial symmetry irrep;
                                // It is only used for transition density matrix calculations.
  bool m_transition_diff_spatial_irrep=false;
  SpinQuantum m_molecule_quantum;
  int m_total_spin;
  int m_guess_permutations;
  bool m_stateSpecific; //when targetting excited states we switch from state 
  //average to statespecific 
  int m_occupied_orbitals;

  vector<int> m_baseState;
  vector<int> m_projectorState;
  int m_targetState;

  std::vector<int> m_hf_occupancy;
  std::string m_hf_occ_user;
  std::vector<double> m_weights;

  std::vector<int> m_sweep_iter_schedule;
  std::vector<int> m_sweep_state_schedule;
  std::vector<int> m_sweep_qstate_schedule;
  std::vector<double> m_sweep_tol_schedule;
  std::vector<double> m_sweep_noise_schedule;
  std::vector<double> m_sweep_additional_noise_schedule;
  bool m_schedule_type_default;
  bool m_schedule_type_backward;
  int m_lastM;
  int m_startM;
  int m_maxM;
  int m_integral_disk_storage_thresh;
  int m_num_Integrals;

  bool m_do_diis;
  double m_diis_error;
  int m_start_diis_iter;
  int m_diis_keep_states;
  double m_diis_error_tol;

  calcType m_calc_type;
  noiseTypes m_noise_type;
  hamTypes m_ham_type;
  WarmUpTypes m_warmup;
  int m_nroots;
  solveTypes m_solve_type;
  bool m_do_deriv;
  bool m_do_fci;
  bool m_do_pdm;
  bool m_do_npdm_ops;
  bool m_do_npdm_in_core;
  bool m_new_npdm_code;
  bool m_store_spinpdm;
  bool m_spatpdm_disk_dump;
  bool m_pdm_unsorted;
  bool m_store_nonredundant_pdm;
  bool m_npdm_intermediate;
  bool m_npdm_multinode;
  bool m_set_Sz;
  int m_maxiter;
  double m_oneindex_screen_tol;
  double m_twoindex_screen_tol;
  bool m_no_transform;
  bool m_add_noninteracting_orbs;

  int m_nquanta;
  int m_sys_add;
  int m_env_add;

  int m_deflation_min_size;
  int m_deflation_max_size;

  algorithmTypes m_algorithm_type;
  int m_twodot_to_onedot_iter;
  std::vector< std::map<SpinQuantum, int> > m_quantaToKeep;
  std::string  m_save_prefix;
  std::string m_load_prefix;
  bool m_direct;
  std::vector<double> m_orbenergies;

  int m_maxj;
  ninejCoeffs m_ninej;
  int m_max_lanczos_dimension;
  
  int n_twodot_noise;
  double m_twodot_noise;
  double m_twodot_gamma;

  double m_sweep_tol;
  bool m_restart;
  bool m_backward;
  bool m_fullrestart;
  bool m_restart_warm;
  bool m_reset_iterations;
  bool m_implicitTranspose;


  std::vector<int> m_spin_vector;
  std::vector<int> m_spin_orbs_symmetry;
  int m_num_spatial_orbs;
  std::vector<int> m_spatial_to_spin;
  std::vector<int> m_spin_to_spatial;

  int m_outputlevel;
  orbitalFormat m_orbformat;

  int m_reorderType;
  string m_reorderfile;
  std::vector<int> m_reorder;//this can be manual, fiedler, gaopt or noreorder
  string m_gaconffile;
  
  bool m_calc_ri_4pdm;
  bool m_store_ripdm_readable;
  bool m_nevpt2;
  bool m_conventional_nevpt2;
  int m_kept_nevpt2_states;
  pair<bool,int> NevPrint;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & m_thrds_per_node & m_spinAdapted & m_Bogoliubov & m_stateSpecific & m_implicitTranspose & m_num_Integrals;
    ar & m_norbs & m_alpha & m_beta & m_solve_type & m_Sz & m_set_Sz & m_baseState& m_projectorState& m_targetState;
    ar & m_spin_vector & m_spin_orbs_symmetry & m_guess_permutations & m_nroots & m_weights & m_hf_occ_user & m_hf_occupancy;
    ar & m_sweep_iter_schedule & m_sweep_state_schedule & m_sweep_qstate_schedule & m_sweep_tol_schedule & m_sweep_noise_schedule &m_sweep_additional_noise_schedule & m_reorder;
    ar & m_molecule_quantum & m_total_symmetry_number & m_total_spin & m_orbenergies & m_add_noninteracting_orbs;
    ar & m_bra_symmetry_number;
    ar & m_save_prefix & m_load_prefix & m_direct & m_max_lanczos_dimension;
    ar & m_deflation_min_size & m_deflation_max_size & m_outputlevel & m_reorderfile;
    ar & m_algorithm_type & m_twodot_to_onedot_iter & m_orbformat ;
    ar & m_nquanta & m_sys_add & m_env_add & m_do_fci & m_no_transform ;
    ar & m_do_pdm & m_do_npdm_ops & m_do_npdm_in_core & m_new_npdm_code & m_occupied_orbitals;
    ar &  m_store_spinpdm &m_spatpdm_disk_dump & m_pdm_unsorted & m_npdm_intermediate& m_npdm_multinode;
    ar & m_maxj & m_ninej & m_maxiter & m_do_deriv & m_oneindex_screen_tol & m_twoindex_screen_tol & m_quantaToKeep & m_noise_type;
    ar & m_sweep_tol & m_restart & m_backward & m_fullrestart & m_restart_warm & m_reset_iterations & m_calc_type & m_ham_type & m_warmup;
    ar & m_do_diis & m_diis_error & m_start_diis_iter & m_diis_keep_states & m_diis_error_tol & m_num_spatial_orbs;
    ar & m_spatial_to_spin & m_spin_to_spatial & m_maxM & m_schedule_type_backward & m_schedule_type_default &m_integral_disk_storage_thresh;
    ar & n_twodot_noise & m_twodot_noise & m_twodot_gamma;
    ar & m_calc_ri_4pdm & m_store_ripdm_readable & m_nevpt2 & m_conventional_nevpt2 & m_kept_nevpt2_states & NevPrint;
  }


  void initialize_defaults();

  std::vector<int> hfOccGenerator_ ();

 public:
  //Input() : m_ninej(ninejCoeffs::getinstance()){}
  Input() {}
  Input (const std::string& config_name);
  // ROA
  void initCumulTimer()
  {
  ddscreen       = boost::shared_ptr<cumulTimer> (new cumulTimer());
  cdscreen       = boost::shared_ptr<cumulTimer> (new cumulTimer());
  dscreen       = boost::shared_ptr<cumulTimer> (new cumulTimer());
  buildcsfops       = boost::shared_ptr<cumulTimer> (new cumulTimer());
  sysdotmake       = boost::shared_ptr<cumulTimer> (new cumulTimer());
  initnewsystem       = boost::shared_ptr<cumulTimer> (new cumulTimer());
  guessgenT       = boost::shared_ptr<cumulTimer> (new cumulTimer());
  multiplierT     = boost::shared_ptr<cumulTimer> (new cumulTimer());
  operrotT        = boost::shared_ptr<cumulTimer> (new cumulTimer());
  davidsonT       = boost::shared_ptr<cumulTimer> (new cumulTimer()); 
  rotmatrixT      = boost::shared_ptr<cumulTimer> (new cumulTimer()); 
  blockdavid      = boost::shared_ptr<cumulTimer> (new cumulTimer()); 
  datatransfer    = boost::shared_ptr<cumulTimer> (new cumulTimer());
  hmultiply       = boost::shared_ptr<cumulTimer> (new cumulTimer()); 
  oneelecT        = boost::shared_ptr<cumulTimer> (new cumulTimer());
  twoelecT        = boost::shared_ptr<cumulTimer> (new cumulTimer());
  makeopsT        = boost::shared_ptr<cumulTimer> (new cumulTimer());
  collectqT       = boost::shared_ptr<cumulTimer> (new cumulTimer());
  opallocate      = boost::shared_ptr<cumulTimer> (new cumulTimer());
  opcatenate      = boost::shared_ptr<cumulTimer> (new cumulTimer());
  oprelease       = boost::shared_ptr<cumulTimer> (new cumulTimer());
  opequateT       = boost::shared_ptr<cumulTimer> (new cumulTimer());
  justmultiply    = boost::shared_ptr<cumulTimer> (new cumulTimer());
  spinrotation    = boost::shared_ptr<cumulTimer> (new cumulTimer());
  otherrotation   = boost::shared_ptr<cumulTimer> (new cumulTimer());
  solvewf         = boost::shared_ptr<cumulTimer> (new cumulTimer());
  postwfrearrange = boost::shared_ptr<cumulTimer> (new cumulTimer());
  couplingcoeff   = boost::shared_ptr<cumulTimer> (new cumulTimer());
  buildsumblock   = boost::shared_ptr<cumulTimer> (new cumulTimer());
  buildblockops   = boost::shared_ptr<cumulTimer> (new cumulTimer());
  addnoise        = boost::shared_ptr<cumulTimer> (new cumulTimer());
  s0time          = boost::shared_ptr<cumulTimer> (new cumulTimer());
  s1time          = boost::shared_ptr<cumulTimer> (new cumulTimer());
  s2time          = boost::shared_ptr<cumulTimer> (new cumulTimer());
  blockintegrals  = boost::shared_ptr<cumulTimer> (new cumulTimer());
  blocksites      = boost::shared_ptr<cumulTimer> (new cumulTimer());
  statetensorproduct = boost::shared_ptr<cumulTimer> (new cumulTimer());
  statecollectquanta = boost::shared_ptr<cumulTimer> (new cumulTimer());
  builditeratorsT = boost::shared_ptr<cumulTimer> (new cumulTimer());
  diskio = boost::shared_ptr<cumulTimer> (new cumulTimer());
  }
  void writeSummary();
#ifdef MOLPRO
  void writeSummaryForMolpro();
#endif
  void performSanityTest();
  void generateDefaultSchedule();
  void readorbitalsfile(string& dumpFile, OneElectronArray& v1, TwoElectronArray& v2, double& coreEnergy, int integralIndex);
  void readorbitalsfile(string& dumpFile, OneElectronArray& v1, TwoElectronArray& v2, double& coreEnergy, PairArray& vcc, CCCCArray& vcccc, CCCDArray& vcccd);  
  int getNumIntegrals() { return m_num_Integrals;}
  void readreorderfile(ifstream& dumpFile, std::vector<int>& reorder);
  std::vector<int> getgaorder(ifstream& gaconfFile, string& orbitalfile, std::vector<int>& fiedlerorder);
  std::vector<int> get_fiedler(string& dumpname);
  std::vector<int> get_fiedler_bcs(string& dumpname);  
  void usedkey_error(string& key, string& line);
  void makeInitialHFGuess();
  static void ReadMeaningfulLine(ifstream&, string&, int);

  boost::shared_ptr<cumulTimer> dscreen; 
  boost::shared_ptr<cumulTimer> cdscreen; 
  boost::shared_ptr<cumulTimer> ddscreen; 
  boost::shared_ptr<cumulTimer> buildcsfops; 
  boost::shared_ptr<cumulTimer> sysdotmake; 
  boost::shared_ptr<cumulTimer> initnewsystem; 
  boost::shared_ptr<cumulTimer> diskio; 
  boost::shared_ptr<cumulTimer> builditeratorsT; 
  boost::shared_ptr<cumulTimer> blockintegrals; 
  boost::shared_ptr<cumulTimer> blocksites; 
  boost::shared_ptr<cumulTimer> statetensorproduct; 
  boost::shared_ptr<cumulTimer> statecollectquanta; 
  boost::shared_ptr<cumulTimer> guessgenT; 
  boost::shared_ptr<cumulTimer> multiplierT; 
  boost::shared_ptr<cumulTimer> operrotT; 
  boost::shared_ptr<cumulTimer> davidsonT; 
  boost::shared_ptr<cumulTimer> rotmatrixT; 
  boost::shared_ptr<cumulTimer> blockdavid; 
  boost::shared_ptr<cumulTimer> datatransfer; 
  boost::shared_ptr<cumulTimer> hmultiply; 
  boost::shared_ptr<cumulTimer> oneelecT; 
  boost::shared_ptr<cumulTimer> twoelecT; 
  boost::shared_ptr<cumulTimer> makeopsT; 
  boost::shared_ptr<cumulTimer> collectqT; 
  boost::shared_ptr<cumulTimer> opallocate; 
  boost::shared_ptr<cumulTimer> opcatenate; 
  boost::shared_ptr<cumulTimer> oprelease; 
  boost::shared_ptr<cumulTimer> opequateT; 
  boost::shared_ptr<cumulTimer> justmultiply; 
  boost::shared_ptr<cumulTimer> spinrotation; 
  boost::shared_ptr<cumulTimer> otherrotation; 
  boost::shared_ptr<cumulTimer> solvewf; 
  boost::shared_ptr<cumulTimer> postwfrearrange; 
  boost::shared_ptr<cumulTimer> couplingcoeff; 
  boost::shared_ptr<cumulTimer> buildsumblock; 
  boost::shared_ptr<cumulTimer> buildblockops; 
  boost::shared_ptr<cumulTimer> addnoise; 
  boost::shared_ptr<cumulTimer> s0time; 
  boost::shared_ptr<cumulTimer> s1time; 
  boost::shared_ptr<cumulTimer> s2time; 

  const std::vector<int>& baseStates() const {return m_baseState;}
  const int& targetState() const {return m_targetState;}
  const std::vector<int>& projectorStates() const {return m_projectorState;}

  std::vector<int>& baseStates() {return m_baseState;}
  int& targetState() {return m_targetState;}
  std::vector<int>& projectorStates() {return m_projectorState;}

  const int& num_occupied_orbitals() const {return m_occupied_orbitals;}
  const bool& doimplicitTranspose() const {return m_implicitTranspose;}
  bool& setimplicitTranspose() {return m_implicitTranspose;}
  const bool& setStateSpecific() const {return m_stateSpecific;}
  bool& setStateSpecific() {return m_stateSpecific;}
  const orbitalFormat& orbformat() const {return m_orbformat;}
  const int& outputlevel() const {return m_outputlevel;}
  int& setOutputlevel()  {return m_outputlevel;}
  const vector<int>& spatial_to_spin() const {return m_spatial_to_spin;}
  int spatial_to_spin(int i) const {return m_spatial_to_spin.at(i);}
  const vector<int>& spin_to_spatial() const {return m_spin_to_spatial;}
  const double& diis_error_tol() const {return m_diis_error_tol;}
  const bool& do_diis() const {return m_do_diis;}
  const double& diis_error() const {return m_diis_error;}
  const int& start_diis_iter() const {return m_start_diis_iter;}
  const int& diis_keep_states() const {return m_diis_keep_states;}

  bool use_partial_two_integrals() const {return (m_norbs/2 >= m_integral_disk_storage_thresh);}
  bool& set_fullrestart() {return m_fullrestart;}
  const bool& get_fullrestart() const {return m_fullrestart;}
  const bool& get_backward() const {return m_backward;}
  const double& get_sweep_tol() const {return m_sweep_tol;}
  const int& get_twodot_method() const {return n_twodot_noise;} 
  const double& get_twodot_noise() const {return m_twodot_noise;} 
  double& set_twodot_noise()  {return m_twodot_noise;}
  const double& get_twodot_gamma() const {return m_twodot_gamma;}
  double& set_twodot_gamma()  {return m_twodot_gamma;}
  const bool& get_restart() const {return m_restart;}
  const bool& get_restart_warm() const {return m_restart_warm;}
  const bool& get_reset_iterations() const {return m_reset_iterations;}
  const ninejCoeffs& get_ninej() const {return m_ninej;}
  const hamTypes &hamiltonian() const {return m_ham_type;}
  const WarmUpTypes &warmup() const {return m_warmup;}
  const int &guess_permutations() const { return m_guess_permutations; }
  const int &max_lanczos_dimension() const {return m_max_lanczos_dimension;}
  std::vector<int> thrds_per_node() const { return m_thrds_per_node; }
  const calcType &calc_type() const { return m_calc_type; }
  const solveTypes &solve_method() const { return m_solve_type; }
  const noiseTypes &noise_type() const {return m_noise_type;}
  const bool &set_Sz() const {return m_set_Sz;}
  const algorithmTypes &algorithm_method() const { return m_algorithm_type; }
  algorithmTypes &set_algorithm_method() { return m_algorithm_type; }
  int twodot_to_onedot_iter() const { return m_twodot_to_onedot_iter; }
  std::vector< std::map<SpinQuantum, int> >& get_quantaToKeep() { return m_quantaToKeep;}
  const std::vector<int> &hf_occupancy() const { return m_hf_occupancy; }
  const std::vector<int> &spin_orbs_symmetry() const { return m_spin_orbs_symmetry; }
  std::vector<double> weights(int sweep_iter) const;// { return m_weights; }
  std::vector<double> weights() const { return m_weights; }
  const std::vector<int> &sweep_iter_schedule() const { return m_sweep_iter_schedule; }
  const std::vector<int> &sweep_state_schedule() const { return m_sweep_state_schedule; }
  const std::vector<int> &sweep_qstate_schedule() const { return m_sweep_qstate_schedule; }
  const std::vector<double> &sweep_tol_schedule() const { return m_sweep_tol_schedule; }
  const std::vector<double> &sweep_noise_schedule() const { return m_sweep_noise_schedule; }
  const std::vector<double> &sweep_additional_noise_schedule() const { return m_sweep_additional_noise_schedule; }
  std::vector<double> &set_sweep_noise_schedule() { return m_sweep_noise_schedule; }
  std::vector<double> &set_sweep_additional_noise_schedule() { return m_sweep_additional_noise_schedule; }
  int& Sz() {return m_Sz;}
  int nroots(int sweep_iter) const;
  int nroots() const {return m_nroots;}
  int real_particle_number() const { return (m_alpha + m_beta);}
  int total_particle_number() const { if(!m_add_noninteracting_orbs) return (m_alpha + m_beta); else return (2*m_alpha); }
  bool calc_ri_4pdm() const {return m_calc_ri_4pdm;}
  bool store_ripdm_readable() const {return m_store_ripdm_readable;}
  bool nevpt2() const {return m_nevpt2;}
  bool read_higherpdm() const {return m_conventional_nevpt2;}
  int kept_nevpt2_states() const {return m_kept_nevpt2_states;}
  bool Print() const {return NevPrint.first;}
  int PrintIndex() const{return NevPrint.second;}
  void SetPrint(bool p, int i=0){NevPrint.first = p;NevPrint.second=i;}
  const SpinSpace total_spin_number() const { if (!m_add_noninteracting_orbs) return SpinSpace(m_alpha - m_beta); else return SpinSpace(0); }
  int last_site() const 
  { 
    if(m_spinAdapted) 
      return m_num_spatial_orbs; 
    else 
      return 2*m_num_spatial_orbs; 
  }
  const bool &no_transform() const { return m_no_transform; }
  const int &deflation_min_size() const { return m_deflation_min_size; }
  const bool &direct() const { return m_direct; }
  const int &deflation_max_size() const { return m_deflation_max_size; }
  const IrrepSpace &total_symmetry_number() const { return m_total_symmetry_number; }
  const IrrepSpace &bra_symmetry_number() const { return m_bra_symmetry_number; }
  const SpinQuantum &molecule_quantum() const { return m_molecule_quantum; }
  const int &sys_add() const { return m_sys_add; }
  const bool &add_noninteracting_orbs() const {return m_add_noninteracting_orbs;}
  bool &add_noninteracting_orbs() {return m_add_noninteracting_orbs;}
  const int &nquanta() const { return m_nquanta; }
  const int &env_add() const { return m_env_add; }
  const bool &do_fci() const { return m_do_fci; }
  const int &max_iter() const { return m_maxiter; }
  const double &oneindex_screen_tol() const { return m_oneindex_screen_tol; }
  double &oneindex_screen_tol() { return m_oneindex_screen_tol; }
  const double &twoindex_screen_tol() const { return m_twoindex_screen_tol; }
  double &twoindex_screen_tol() { return m_twoindex_screen_tol; }
  const int &total_spin() const {return m_total_spin;}
  const std::vector<int> &spin_vector() const { return m_spin_vector; }
  const std::string &save_prefix() const { return m_save_prefix; }
  const std::string &load_prefix() const { return m_load_prefix; }
  SpinQuantum& set_molecule_quantum() {return m_molecule_quantum;}
  SpinQuantum effective_molecule_quantum() {
    if (!m_add_noninteracting_orbs) 
      return m_molecule_quantum;
    else 
      return SpinQuantum(m_molecule_quantum.particleNumber + m_molecule_quantum.totalSpin.getirrep(), SpinSpace(0), m_molecule_quantum.orbitalSymmetry);
    //return SpinQuantum(total_particle_number() + total_spin_number().getirrep(), SpinSpace(0), total_symmetry_number());
  }
  SpinQuantum bra_quantum() {
    if (!m_add_noninteracting_orbs) 
      return SpinQuantum(total_particle_number(), SpinSpace(m_alpha - m_beta), bra_symmetry_number());
    else 
      return SpinQuantum(total_particle_number() + total_spin_number().getirrep(), SpinSpace(0), bra_symmetry_number());
  }
  vector<SpinQuantum> effective_molecule_quantum_vec() {
    vector<SpinQuantum> q;
    if (!m_Bogoliubov) q.push_back(effective_molecule_quantum());
    else {
      SpinQuantum q_max = effective_molecule_quantum();
      for (int n = 0; n <= q_max.get_n(); n+=2) {
        q.push_back(SpinQuantum(n, q_max.get_s(), q_max.get_symm()));
      }
    }
    return q;
  }
  vector<SpinQuantum> bra_quantum_vec() {
    vector<SpinQuantum> q;
    if (!m_Bogoliubov) q.push_back(bra_quantum());
    else {
      SpinQuantum q_max = bra_quantum();
      for (int n = 0; n <= q_max.get_n(); n+=2) {
        q.push_back(SpinQuantum(n, q_max.get_s(), q_max.get_symm()));
      }
    }
    return q;
  }
  bool transition_diff_irrep(){
    return m_transition_diff_spatial_irrep;
  }
  std::vector<double>& get_orbenergies() {return m_orbenergies;}
  int getHFQuanta(const SpinBlock& b) const;
  const bool &do_pdm() const {return m_do_pdm;}
  bool &do_pdm() {return m_do_pdm;}
  const bool &do_npdm_ops() const {return m_do_npdm_ops;}
  bool &do_npdm_ops() {return m_do_npdm_ops;}
  const bool &do_npdm_in_core() const {return m_do_npdm_in_core;}
  bool &do_npdm_in_core() {return m_do_npdm_in_core;}
  const bool new_npdm_code() const{
    if( m_do_pdm) return m_new_npdm_code;
    else return false;
  }
  void set_new_npdm_code(){ m_new_npdm_code= true;}
  const bool &store_spinpdm() const {return m_store_spinpdm;}
  bool &store_spinpdm() {return m_store_spinpdm;}
  const bool &spatpdm_disk_dump() const {return m_spatpdm_disk_dump;}
  bool &spatpdm_disk_dump() {return m_spatpdm_disk_dump;}
  const bool &pdm_unsorted() const {return m_pdm_unsorted;}
  bool &pdm_unsorted(){return m_pdm_unsorted;}
  const bool &store_nonredundant_pdm() const { return m_store_nonredundant_pdm;}
  bool &store_nonredundant_pdm() { return m_store_nonredundant_pdm;}
  int slater_size() const {return m_norbs;}
  const std::vector<int> &reorder_vector() {return m_reorder;}
  bool spinAdapted() {return m_spinAdapted;}
  bool &npdm_intermediate() { return m_npdm_intermediate; }
  const bool &npdm_intermediate() const { return m_npdm_intermediate; }
  bool &npdm_multinode() { return m_npdm_multinode; }
  const bool &npdm_multinode() const { return m_npdm_multinode; }
};
}
#endif
