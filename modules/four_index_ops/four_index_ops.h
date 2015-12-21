#ifndef FOURPDM_OPERATORS_HEADER
#define FOURPDM_OPERATORS_HEADER

#include "StackBaseOperator.h"

namespace SpinAdapted{


//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// This is a skeleton class only used to set up infrastructure, and the operator contents are never stored or built during the standard sweep
class RI4index: public SpinAdapted::StackSparseMatrix
{
  public:
    RI4index() { abort(); }
    void build(const StackSpinBlock& b) { abort(); }
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block) { abort(); }
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b) { abort(); }
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class StackCreCreDesDes: public SpinAdapted::StackSparseMatrix
{
  public:
    StackCreCreDesDes() { orbs.resize(4); fermion = false; build_pattern = "(((CC)(D))(D))";} // default build_pattern
    void build(const StackSpinBlock& b);
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class StackCreDesCreDes: public SpinAdapted::StackSparseMatrix
{
  public:
    StackCreDesCreDes() { orbs.resize(4); fermion = false; build_pattern = "(((CD)(C))(D))";} // default build_pattern
    void build(const StackSpinBlock& b);
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class StackCreDesDesCre: public SpinAdapted::StackSparseMatrix
{
  public:
    StackCreDesDesCre() { orbs.resize(4); fermion = false; build_pattern = "(((CD)(D))(C))";} // default build_pattern
    void build(const StackSpinBlock& b);
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class StackCreDesDesDes: public SpinAdapted::StackSparseMatrix
{
  public:
    StackCreDesDesDes() { orbs.resize(4); fermion = false; build_pattern = "(((CD)(D))(D))";} // default build_pattern
    void build(const StackSpinBlock& b);
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class StackCreCreCreDes: public SpinAdapted::StackSparseMatrix
{
  public:
    StackCreCreCreDes() { orbs.resize(4); fermion = false; build_pattern = "(((CC)(C))(D))";} // default build_pattern
    void build(const StackSpinBlock& b);
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class StackCreCreDesCre: public SpinAdapted::StackSparseMatrix
{
  public:
    StackCreCreDesCre() { orbs.resize(4); fermion = false; build_pattern = "(((CC)(D))(C))";} // default build_pattern
    void build(const StackSpinBlock& b);
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class StackCreDesCreCre: public SpinAdapted::StackSparseMatrix
{
  public:
    StackCreDesCreCre() { orbs.resize(4); fermion = false; build_pattern = "(((CD)(C))(C))";} // default build_pattern
    void build(const StackSpinBlock& b);
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class StackCreCreCreCre: public SpinAdapted::StackSparseMatrix
{
  public:
    StackCreCreCreCre() { orbs.resize(4); fermion = false; build_pattern = "(((CC)(C))(C))";} // default build_pattern
    void build(const StackSpinBlock& b);
//    void build_from_disk(StackSpinBlock& b, std::ifstream& sysfs, std::ifstream& dotfs) { abort(); }
    boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block);
    double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

}

#endif
