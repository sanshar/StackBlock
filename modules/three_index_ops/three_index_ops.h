#ifndef THREE_INDEX_OPS_H
#define THREE_INDEX_OPS_H

#include "StackBaseOperator.h"

namespace SpinAdapted{

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// This is a skeleton class only used to set up infrastructure, and the operator contents are never stored or built during the standard sweep
class RI3index: public SpinAdapted::StackSparseMatrix
{
  public:
    RI3index() { abort(); } // should never be instantiated ??FIXME
    void build(const StackSpinBlock& b) { abort(); }
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    //boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block) { abort(); }
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b) { abort(); }
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// This can be avoided with transposes

class StackDesDesDes: public SpinAdapted::StackSparseMatrix
{
  public:
    StackDesDesDes() { orbs.resize(3); fermion = true; build_pattern = "((DD)(D))";} // default build_pattern
    void build(const StackSpinBlock& b);
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    //boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// 3PDM operators
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class StackCreCreDes: public SpinAdapted::StackSparseMatrix
{
  public:
    StackCreCreDes() { orbs.resize(3); fermion = true; build_pattern = "((CC)(D))";} // default build_pattern
    void build(const StackSpinBlock& b);
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    //boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class StackCreDesDes: public SpinAdapted::StackSparseMatrix
{
  public:
    StackCreDesDes() { orbs.resize(3); fermion = true; build_pattern = "((CD)(D))";} // default build_pattern
    void build(const StackSpinBlock& b);
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    //boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class StackCreDesCre: public SpinAdapted::StackSparseMatrix
{
  public:
    StackCreDesCre() { orbs.resize(3); fermion = true; build_pattern = "((CD)(C))";} // default build_pattern
    void build(const StackSpinBlock& b);
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    //boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class StackCreCreCre: public SpinAdapted::StackSparseMatrix
{
  public:
    StackCreCreCre() { orbs.resize(3); fermion = true; build_pattern = "((CC)(C))";} // default build_pattern
    void build(const StackSpinBlock& b);
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    //boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// 4PDM operators
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class StackDesCreDes: public SpinAdapted::StackSparseMatrix
{
  public:
    StackDesCreDes() { orbs.resize(3); fermion = true; build_pattern = "((DC)(D))";} // default build_pattern
    void build(const StackSpinBlock& b);
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    //boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class StackDesDesCre: public SpinAdapted::StackSparseMatrix
{
  public:
    StackDesDesCre() { orbs.resize(3); fermion = true; build_pattern = "((DD)(C))";} // default build_pattern
    void build(const StackSpinBlock& b); 
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    //boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class StackDesCreCre: public SpinAdapted::StackSparseMatrix
{
  public:
    StackDesCreCre() { orbs.resize(3); fermion = true; build_pattern = "((D)(CC))";} // default build_pattern
    void build(const StackSpinBlock& b);
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs);
    //boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

}

#endif
