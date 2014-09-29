#ifndef FOURPDM_OPERATORS_HEADER
#define FOURPDM_OPERATORS_HEADER

#include "BaseOperator.h"

namespace SpinAdapted{


//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// This is a skeleton class only used to set up infrastructure, and the operator contents are never stored or built during the standard sweep
class RI4index: public SpinAdapted::SparseMatrix
{
  public:
    RI4index() { abort(); }
    void build(const SpinBlock& b) { abort(); }
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block) { abort(); }
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b) { abort(); }
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreCreDesDes: public SpinAdapted::SparseMatrix
{
  public:
    CreCreDesDes() { orbs.resize(4); fermion = false; build_pattern = "(((CC)(D))(D))";} // default build_pattern
    void build(const SpinBlock& b);
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreDesCreDes: public SpinAdapted::SparseMatrix
{
  public:
    CreDesCreDes() { orbs.resize(4); fermion = false; build_pattern = "(((CD)(C))(D))";} // default build_pattern
    void build(const SpinBlock& b);
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreDesDesCre: public SpinAdapted::SparseMatrix
{
  public:
    CreDesDesCre() { orbs.resize(4); fermion = false; build_pattern = "(((CD)(D))(C))";} // default build_pattern
    void build(const SpinBlock& b);
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreDesDesDes: public SpinAdapted::SparseMatrix
{
  public:
    CreDesDesDes() { orbs.resize(4); fermion = false; build_pattern = "(((CD)(D))(D))";} // default build_pattern
    void build(const SpinBlock& b);
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreCreCreDes: public SpinAdapted::SparseMatrix
{
  public:
    CreCreCreDes() { orbs.resize(4); fermion = false; build_pattern = "(((CC)(C))(D))";} // default build_pattern
    void build(const SpinBlock& b);
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreCreDesCre: public SpinAdapted::SparseMatrix
{
  public:
    CreCreDesCre() { orbs.resize(4); fermion = false; build_pattern = "(((CC)(D))(C))";} // default build_pattern
    void build(const SpinBlock& b);
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreDesCreCre: public SpinAdapted::SparseMatrix
{
  public:
    CreDesCreCre() { orbs.resize(4); fermion = false; build_pattern = "(((CD)(C))(C))";} // default build_pattern
    void build(const SpinBlock& b);
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreCreCreCre: public SpinAdapted::SparseMatrix
{
  public:
    CreCreCreCre() { orbs.resize(4); fermion = false; build_pattern = "(((CC)(C))(C))";} // default build_pattern
    void build(const SpinBlock& b);
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

}

#endif
