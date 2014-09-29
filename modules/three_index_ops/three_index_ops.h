#ifndef THREE_INDEX_OPS_H
#define THREE_INDEX_OPS_H

#include "BaseOperator.h"

namespace SpinAdapted{

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// This is a skeleton class only used to set up infrastructure, and the operator contents are never stored or built during the standard sweep
class RI3index: public SpinAdapted::SparseMatrix
{
  public:
    RI3index() { abort(); } // should never be instantiated ??FIXME
    void build(const SpinBlock& b) { abort(); }
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block) { abort(); }
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b) { abort(); }
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// This can be avoided with transposes

class DesDesDes: public SpinAdapted::SparseMatrix
{
  public:
    DesDesDes() { orbs.resize(3); fermion = true; build_pattern = "((DD)(D))";} // default build_pattern
    void build(const SpinBlock& b);
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// 3PDM operators
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreCreDes: public SpinAdapted::SparseMatrix
{
  public:
    CreCreDes() { orbs.resize(3); fermion = true; build_pattern = "((CC)(D))";} // default build_pattern
    void build(const SpinBlock& b);
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreDesDes: public SpinAdapted::SparseMatrix
{
  public:
    CreDesDes() { orbs.resize(3); fermion = true; build_pattern = "((CD)(D))";} // default build_pattern
    void build(const SpinBlock& b);
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreDesCre: public SpinAdapted::SparseMatrix
{
  public:
    CreDesCre() { orbs.resize(3); fermion = true; build_pattern = "((CD)(C))";} // default build_pattern
    void build(const SpinBlock& b);
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CreCreCre: public SpinAdapted::SparseMatrix
{
  public:
    CreCreCre() { orbs.resize(3); fermion = true; build_pattern = "((CC)(C))";} // default build_pattern
    void build(const SpinBlock& b);
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// 4PDM operators
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class DesCreDes: public SpinAdapted::SparseMatrix
{
  public:
    DesCreDes() { orbs.resize(3); fermion = true; build_pattern = "((DC)(D))";} // default build_pattern
    void build(const SpinBlock& b);
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class DesDesCre: public SpinAdapted::SparseMatrix
{
  public:
    DesDesCre() { orbs.resize(3); fermion = true; build_pattern = "((DD)(C))";} // default build_pattern
    void build(const SpinBlock& b); 
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class DesCreCre: public SpinAdapted::SparseMatrix
{
  public:
    DesCreCre() { orbs.resize(3); fermion = true; build_pattern = "((D)(CC))";} // default build_pattern
    void build(const SpinBlock& b);
//    void build_from_disk(SpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    boost::shared_ptr<SparseMatrix> getworkingrepresentation(const SpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const SpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

}

#endif
