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


namespace SpinAdapted{
class Wavefunction;
class DensityMatrix;

enum Storagetype {LOCAL_STORAGE, DISTRIBUTED_STORAGE};

class SpinBlock
{
 private:
  friend class boost::serialization::access;
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & localstorage & name & complementary & hasMemoryAllocated & normal & direct & loopblock
	& sites & complementary_sites & stateInfo;
      //FIX ME!! remove register_type stuff and add BOOST_CLASS_EXPORT to op_components.h (will take longer to compile)                     
      ar.register_type(static_cast<Op_component<Cre> *>(NULL));
      ar.register_type(static_cast<Op_component<CreDes> *>(NULL));
      ar.register_type(static_cast<Op_component<CreCre> *>(NULL));
      ar.register_type(static_cast<Op_component<CreDesComp> *>(NULL));
      ar.register_type(static_cast<Op_component<DesDesComp> *>(NULL));
      ar.register_type(static_cast<Op_component<CreCreDesComp> *>(NULL));
      ar.register_type(static_cast<Op_component<Ham> *>(NULL));
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
  SpinBlock* leftBlock;
  SpinBlock* rightBlock;
  
  StateInfo stateInfo;
  std::vector<int> sites;
  std::vector<int> complementary_sites;
 public: 
  SpinBlock();
  SpinBlock (const StateInfo& s);
  SpinBlock (const SpinBlock& b);
  SpinBlock (int start, int finish, bool is_complement = false);
  void BuildTensorProductBlock (std::vector<int>& new_sites);
  
  static std::string  restore (bool forward, const vector<int>& sites, SpinBlock& b);
  static void store (bool forward, const vector<int>& sites, SpinBlock& b);
  void Save (std::ofstream &ofs);
  void Load (std::ifstream &ifs);

  double memoryUsed();
  void addAdditionalCompOps();
  const StateInfo& get_stateInfo() const {return stateInfo;}
  std::vector<int> make_complement(const std::vector<int>& sites);
  void setstoragetype(Storagetype st);
  void default_op_components(bool complementary_);
  void default_op_components(bool direct, SpinBlock& lBlock, SpinBlock& rBlock, bool haveNormops, bool haveCompops);
  void set_big_components();
  void printOperatorSummary();

  int size() const { return sites.size(); }
  //void build_comp_remove_normal_ops();
  void remove_normal_ops();
  int get_name() const {return name;}
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
  
  Op_component_base& set_op_array(opTypes optype){assert(has(optype));return *(ops.find(optype)->second);}
  Op_component_base& get_op_array(opTypes optype){assert(has(optype));return *(ops.find(optype)->second);}
  const Op_component_base& get_op_array(opTypes optype) const {assert(has(optype));return *(ops.find(optype)->second);}
  
  boost::shared_ptr<SparseMatrix> get_op_rep(const opTypes &optypes, const SpinQuantum& s, int i=-1, int j=-1, int k=-1) {assert(has(optypes)); Op_component_base& opbase = *ops.find(optypes)->second; 
    return opbase.get_op_rep(s, i, j, k)->getworkingrepresentation(this);}
  const boost::shared_ptr<SparseMatrix> get_op_rep(const opTypes &optypes, const SpinQuantum& s, int i=-1, int j=-1, int k=-1) const {assert(has(optypes)); Op_component_base& opbase = *ops.find(optypes)->second; return opbase.get_op_rep(s, i, j, k)->getworkingrepresentation(this);}

  void operator= (const SpinBlock& b);
  void build_iterators();
  void build_operators(std::vector<Csf >& s, std::vector< std::vector<Csf> >& ladders);
  void build_operators();
  void BuildSumBlock(int condition, SpinBlock& b_1, SpinBlock& b_2, StateInfo* compState=0);
  void BuildSumBlockSkeleton(int condition, SpinBlock& lBlock, SpinBlock& rBlock, StateInfo* compState=0);
  void BuildSlaterBlock (std::vector<int> sts, std::vector<SpinQuantum> qnumbers, std::vector<int> distribution, bool random, 
			 const bool haveNormops);
  void set_loopblock(bool p_loopblock){loopblock = p_loopblock;}
  friend ostream& operator<< (ostream& os, const SpinBlock& b);
  void multiplyH(Wavefunction& c, Wavefunction* v, int num_threads) const;
  void diagonalH(DiagonalMatrix& e) const;
  void clear();
  void sendcompOps(Op_component_base& opcomp, int I, int J, int optype, int compsite);
  void recvcompOps(Op_component_base& opcomp, int I, int J, int optype);

  void RenormaliseFrom (std::vector<double> &energies, std::vector<double> &spins, double &error, std::vector<Matrix>& rotateMatrix,
                        const int keptstates, const int keptqstates, const double tol, SpinBlock& big,
                        const guessWaveTypes &guesswavetype, const double noise, const double additional_noise, const bool &onedot, SpinBlock& system, 
			SpinBlock& sysDot, SpinBlock& envDot, SpinBlock& environment, const bool& dot_with_sys, const bool& warmUp, int sweepiter);

  double makeRotateMatrix(DensityMatrix& tracedMatrix, vector<Matrix>& rotateMatrix, const int& keptstates, const int& keptqstates);
  void transform_operators(std::vector<Matrix>& rotateMatrix);
};
}
#endif
