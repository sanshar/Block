#ifndef SPINBLOCK_HEADER
#define SPINBLOCK_HEADER
#include "csf.h" 
#include <para_array.h>
#include <list>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include "StateInfo.h"
#include "op_components.h"
#include "sweep_params.h"
#include "multiarray.h"
#include "boost/variant.hpp"
#include "perturb.h"


namespace SpinAdapted{
class Wavefunction;
class DensityMatrix;

enum Storagetype {LOCAL_STORAGE, DISTRIBUTED_STORAGE, DISTRIBUTED_STORAGE_FOR_ONEPDM};
boost::shared_ptr<Op_component_base> make_new_op(const opTypes &optype, const bool &is_core);

class SpinBlock
{
 private:
  friend class boost::serialization::access;
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & localstorage;
      ar & name;
      ar & complementary;
      ar & hasMemoryAllocated;
      ar & normal; 
      ar & direct;
      ar & loopblock;
      ar & sites;
      ar & complementary_sites ;
      ar & integralIndex;
      ar & nonactive_orbs;
      //ar & nevpt_perturbation;
      //FIX ME!! remove register_type stuff and add BOOST_CLASS_EXPORT to op_components.h (will take longer to compile)                     
      ar.register_type(static_cast<Op_component<Cre> *>(NULL));
      ar.register_type(static_cast<Op_component<Des> *>(NULL));

      ar.register_type(static_cast<Op_component<CreDes> *>(NULL));
      ar.register_type(static_cast<Op_component<DesCre> *>(NULL));
      ar.register_type(static_cast<Op_component<CreCre> *>(NULL));
      ar.register_type(static_cast<Op_component<DesDes> *>(NULL));

      ar.register_type(static_cast<Op_component<CreDesComp> *>(NULL));
      ar.register_type(static_cast<Op_component<DesCreComp> *>(NULL));
      ar.register_type(static_cast<Op_component<DesDesComp> *>(NULL));
      ar.register_type(static_cast<Op_component<CreCreComp> *>(NULL));
      ar.register_type(static_cast<Op_component<CreCreDesComp> *>(NULL));
      ar.register_type(static_cast<Op_component<CreDesDesComp> *>(NULL));

      ar.register_type(static_cast<Op_component<Ham> *>(NULL));
      ar.register_type(static_cast<Op_component<Overlap> *>(NULL));

      // 3PDM
      ar.register_type(static_cast<Op_component<RI3index> *>(NULL));
      ar.register_type(static_cast<Op_component<RI4index> *>(NULL));
      ar.register_type(static_cast<Op_component<CreCreDes> *>(NULL));
      ar.register_type(static_cast<Op_component<CreDesDes> *>(NULL));
      ar.register_type(static_cast<Op_component<CreDesCre> *>(NULL));
      ar.register_type(static_cast<Op_component<CreCreCre> *>(NULL));
      // 4PDM
      ar.register_type(static_cast<Op_component<DesCreDes> *>(NULL));
      ar.register_type(static_cast<Op_component<DesDesCre> *>(NULL));
      ar.register_type(static_cast<Op_component<DesCreCre> *>(NULL));
      ar.register_type(static_cast<Op_component<DesDesDes> *>(NULL));
      ar.register_type(static_cast<Op_component<CreCreDesDes> *>(NULL));
      ar.register_type(static_cast<Op_component<CreDesCreDes> *>(NULL));
      ar.register_type(static_cast<Op_component<CreDesDesCre> *>(NULL));
      ar.register_type(static_cast<Op_component<CreDesDesDes> *>(NULL));
      ar.register_type(static_cast<Op_component<CreCreCreDes> *>(NULL));
      ar.register_type(static_cast<Op_component<CreCreDesCre> *>(NULL));
      ar.register_type(static_cast<Op_component<CreDesCreCre> *>(NULL));
      ar.register_type(static_cast<Op_component<CreCreCreCre> *>(NULL));
      //mps_nevpt
      ar.register_type(static_cast<Op_component<CDD_sum> *>(NULL));
      ar.register_type(static_cast<Op_component<CDD_CreDesComp> *>(NULL));
      ar.register_type(static_cast<Op_component<CDD_DesDesComp> *>(NULL));
      ar.register_type(static_cast<Op_component<CCD_sum> *>(NULL));
      ar.register_type(static_cast<Op_component<CCD_CreDesComp> *>(NULL));
      ar.register_type(static_cast<Op_component<CCD_CreCreComp> *>(NULL));

      ar & ops;
    }

 private:

  std::map<opTypes, boost::shared_ptr<Op_component_base> > ops;
  bool complementary;
  bool normal;
  bool loopblock;
  bool localstorage;
  bool hasMemoryAllocated;
  bool direct;
  int name;
  int integralIndex;
  SpinBlock* leftBlock;
  SpinBlock* rightBlock;
  boost::shared_ptr<TwoElectronArray> twoInt;
  
  StateInfo braStateInfo;
  StateInfo ketStateInfo;
  std::vector<int> sites;
  std::vector<int> complementary_sites;
  // In nevpt2, the integral between nonactive space and active space should be considered. 
  // nonactive_orbs is the vector of orbitals that electrons are excited to or from. 
  std::vector<int> nonactive_orbs;
  //perturber nevpt_perturbation;
 public: 
  SpinBlock();
  SpinBlock (const StateInfo& s, int integralIndex);
  SpinBlock (const SpinBlock& b);
  SpinBlock (int start, int finish, int integralIndex, bool implicitTranspose, bool is_complement = false);
  SpinBlock(int start, int finish, const std::vector<int>& nonactive_orbs_, bool is_complement = false);
  void BuildTensorProductBlock (std::vector<int>& new_sites);
  
  static std::string  restore (bool forward, const vector<int>& sites, SpinBlock& b, int left, int right, char* name=0);//left and right are the bra and ket states and the name is the type of the MPO (currently only H)
  static void store (bool forward, const vector<int>& sites, SpinBlock& b, int left, int right, char* name=0);//left and right are the bra and ket states and the name is the type of the MPO (currently only H) 
  void Save (std::ofstream &ofs);
  void Load (std::ifstream &ifs);

  const boost::shared_ptr<TwoElectronArray> get_twoInt() const {return twoInt;}
  int get_integralIndex() const {return integralIndex;}
  int& set_integralIndex() {return integralIndex;}
  double memoryUsed();
  void addAdditionalCompOps();
  const StateInfo& get_stateInfo() const {return ketStateInfo;}
  const StateInfo& get_braStateInfo() const {return braStateInfo;}
  const StateInfo& get_ketStateInfo() const {return ketStateInfo;}
  StateInfo& set_braStateInfo() {return braStateInfo;}
  StateInfo& set_ketStateInfo() {return ketStateInfo;}
  static std::vector<int> make_complement(const std::vector<int>& sites);
  void setstoragetype(Storagetype st);
  void default_op_components(bool complementary_, bool implicitTranspose);
  void default_op_components(bool direct, SpinBlock& lBlock, SpinBlock& rBlock, bool haveNormops, bool haveCompops, bool implicitTranspose);
  void set_big_components();
  void perturb_op_components(bool direct, SpinBlock& lBlock, SpinBlock& rBlock, const perturber& pb);
  void printOperatorSummary();

  int size() const { return sites.size(); }
  const vector<int>& nonactive_orb() const { return nonactive_orbs;}
  vector<int>& nonactive_orb() { return nonactive_orbs;}
  const int& nonactive_orb(int i) const { return nonactive_orbs[i];}
  int& nonactive_orb(int i) { return nonactive_orbs[i];}
  //const perturber nevpt_pb() const { return nevpt_perturbation;}
  //perturber& nevpt_pb() { return nevpt_perturbation;}
  //void build_comp_remove_normal_ops();
  void remove_normal_ops();
  int get_name() const {return name;}
  std::vector<int>& set_sites() {return sites;}
  const std::vector<int>& get_sites() const {return sites;}
  const std::vector<int>& get_complementary_sites() const {return complementary_sites;}
  bool is_normal() const {return normal;}
  bool is_complementary() const {return complementary;}
  bool is_loopblock() const {return loopblock;}
  bool is_direct() const {return direct;}
  bool getlocalstorage() const {return localstorage;}
  bool has(opTypes optype) const
  {
    if(ops.find(optype) != ops.end())
      return true;
    else
      return false;
  }
  SpinBlock* get_leftBlock() const {return leftBlock;}
  SpinBlock* get_rightBlock() const {return rightBlock;}
  
  //void makeCCD_comp_ops();

  void initialise_op_array(opTypes optype, bool is_core);
  void erase(opTypes optype) {assert(has(optype)); ops.erase(optype);}
  boost::shared_ptr<Op_component_base>& set_op_array(opTypes optype){assert(has(optype));return ops.find(optype)->second;}
  Op_component_base& get_op_array(opTypes optype){assert(has(optype));return *(ops.find(optype)->second);}
  const Op_component_base& get_op_array(opTypes optype) const {assert(has(optype));return *(ops.find(optype)->second);}
  

  boost::shared_ptr<SparseMatrix> get_op_rep(const opTypes &optypes, const SpinQuantum& s, int i=-1, int j=-1, int k=-1,int l = -1) {
    assert(has(optypes));
    Op_component_base& opbase = *ops.find(optypes)->second;
    vector<SpinQuantum> temp(1, s);
    return opbase.get_op_rep(temp, i, j, k, l)->getworkingrepresentation(this);
  }
  const boost::shared_ptr<SparseMatrix> get_op_rep(const opTypes &optypes, const SpinQuantum& s, int i=-1, int j=-1, int k=-1, int l=-1) const {
    assert(has(optypes));
    Op_component_base& opbase = *ops.find(optypes)->second;
    vector<SpinQuantum> temp(1, s);    
    return opbase.get_op_rep(temp, i, j, k, l)->getworkingrepresentation(this);
  }
  boost::shared_ptr<SparseMatrix> get_op_rep(const opTypes &optypes, const std::vector<SpinQuantum>& s, int i=-1, int j=-1, int k=-1, int l=-1) {
    assert(has(optypes));
    Op_component_base& opbase = *ops.find(optypes)->second; 
    return opbase.get_op_rep(s, i, j, k, l)->getworkingrepresentation(this);
  }
  const boost::shared_ptr<SparseMatrix> get_op_rep(const opTypes &optypes, const vector<SpinQuantum>& s, int i=-1, int j=-1, int k=-1, int l= -1) const {
    assert(has(optypes));
    Op_component_base& opbase = *ops.find(optypes)->second;
    return opbase.get_op_rep(s, i, j, k, l)->getworkingrepresentation(this);
  }
  boost::shared_ptr<SparseMatrix> get_op_rep(const opTypes &optypes, const std::map< std::string, std::vector<SpinQuantum> >& s, int i=-1, int j=-1, int k=-1, int l= -1) {
    assert(has(optypes));
    Op_component_base& opbase = *ops.find(optypes)->second;
    return opbase.get_op_rep(s, i, j, k, l)->getworkingrepresentation(this);
  }
  const boost::shared_ptr<SparseMatrix> get_op_rep(const opTypes &optypes, const std::map< std::string, std::vector<SpinQuantum> >& s, int i=-1, int j=-1, int k=-1, int l= -1) const {
    assert(has(optypes));
    Op_component_base& opbase = *ops.find(optypes)->second;
    return opbase.get_op_rep(s, i, j, k, l)->getworkingrepresentation(this);
  }

  void operator= (const SpinBlock& b);
  void build_iterators();
  void build_operators(std::vector<Csf >& s, std::vector< std::vector<Csf> >& ladders);
  void build_operators();
  void build_and_renormalise_operators(const std::vector<Matrix>& rotateMatrix, const StateInfo *newStateInfo);
  void build_and_renormalise_operators(const std::vector<Matrix>& leftMat, const StateInfo *bra, const std::vector<Matrix>& rightMat, const StateInfo *ket);
  void renormalise_transform(const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo);
  void renormalise_transform(const std::vector<Matrix>& leftMat, const StateInfo *bra, const std::vector<Matrix>& rightMat, const StateInfo *ket);

  void BuildSumBlock(int condition, SpinBlock& b_1, SpinBlock& b_2, StateInfo* compState=0);
  void BuildSumBlock(int condition, SpinBlock& b_1, SpinBlock& b_2, const std::vector<SpinQuantum>& braquantum, const std::vector<SpinQuantum>& ketquantum);
  void BuildSumBlockSkeleton(int condition, SpinBlock& lBlock, SpinBlock& rBlock, StateInfo* compState=0);
  void BuildSumBlockSkeleton(int condition, SpinBlock& lBlock, SpinBlock& rBlock, const std::vector<SpinQuantum>& braquantum, const std::vector<SpinQuantum>& ketquantum);
  void BuildSlaterBlock (std::vector<int> sts, std::vector<SpinQuantum> qnumbers, std::vector<int> distribution, bool random, 
			 const bool haveNormops);
  void BuildSingleSlaterBlock(std::vector<int> sts);
  void set_loopblock(bool p_loopblock){loopblock = p_loopblock;}
  friend ostream& operator<< (ostream& os, const SpinBlock& b);
  void multiplyH(Wavefunction& c, Wavefunction* v, int num_threads) const;
  void multiplyH_Q(Wavefunction& c, Wavefunction* v, int num_threads, SpinQuantum &Q) const;
  void multiplyOverlap(Wavefunction& c, Wavefunction* v, int num_threads) const;
  void multiplyCDD_sum(Wavefunction& c, Wavefunction* v, int num_threads) const;
  void multiplyCCD_sum(Wavefunction& c, Wavefunction* v, int num_threads) const;
  void diagonalH(DiagonalMatrix& e) const;
  void clear();
  void sendcompOps(Op_component_base& opcomp, int I, int J, int optype, int compsite);
  void recvcompOps(Op_component_base& opcomp, int I, int J, int optype);
  void sendOneindexOps(Op_component_base& opcomp, int I, int optype, int otherproc);
  void recvOneindexOps(Op_component_base& opcomp, int I, int optype);

  void RenormaliseFrom (std::vector<double> &energies, std::vector<double> &spins, double &error, std::vector<Matrix>& rotateMatrix,
                        const int keptstates, const int keptqstates, const double tol, SpinBlock& big,
                        const guessWaveTypes &guesswavetype, const double noise, const double additional_noise, const bool &onedot, SpinBlock& system, 
			SpinBlock& sysDot, SpinBlock& environment, const bool& dot_with_sys, const bool& warmUp, int sweepiter, 
			int currenroot, std::vector<Wavefunction>& lowerStates, DensityMatrix* d=0);

  void transform_operators(std::vector<Matrix>& rotateMatrix);
  void transform_operators(std::vector<Matrix>& leftrotateMatrix, std::vector<Matrix>& rightrotateMatrix, bool clearRightBlock = true, bool clearLeftBlock = true);
};

 double makeRotateMatrix(DensityMatrix& tracedMatrix, vector<Matrix>& rotateMatrix, const int& keptstates, const int& keptqstates, std::vector<DiagonalMatrix> *eigs =0);
}
#endif
