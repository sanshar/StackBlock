#ifndef SPIN_STACKOPXOP_HEADER
#define SPIN_STACKOPXOP_HEADER
 
#include <newmat.h>

namespace SpinAdapted{
  class StackWavefunction;
  class StackSpinBlock;
  class StackSparseMatrix;

namespace stackopxop
{
  void cxcddcomp(const StackSpinBlock* otherblock, const std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, StackSparseMatrix* o);
  void cdxcdcomp(const StackSpinBlock* otherblock, const std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, StackSparseMatrix* o);
  void ddxcccomp(const StackSpinBlock* otherblock, const std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, StackSparseMatrix* o);

  void cxcddcomp_Element(const StackSpinBlock* otherblock, const std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, StackSparseMatrix* o, StackMatrix& m, int row, int col);
  void cdxcdcomp_Element(const StackSpinBlock* otherblock, const std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, StackSparseMatrix* o, StackMatrix& m, int row, int col);
  void ddxcccomp_Element(const StackSpinBlock* otherblock, const std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, StackSparseMatrix* o, StackMatrix& m, int row, int col);


  void cdxcdcomp_d(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, DiagonalMatrix* e);
  void ham_d(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, DiagonalMatrix* e, int proc);
  
 
  void cxcdcomp(const StackSpinBlock* otherblock, const std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, int I, StackSparseMatrix* o, double factor);
  void dxcccomp(const StackSpinBlock* otherBlock, const std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, int I, StackSparseMatrix* o, double scale);
  void cxcdcompElement(const StackSpinBlock* otherblock, const std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, int I, StackSparseMatrix* o, StackMatrix& m, int row, int col, double factor);
  void dxcccompElement(const StackSpinBlock* otherBlock, const std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, int I, StackSparseMatrix* o, StackMatrix& m, int row, int col, double scale);

  void dxcdcomp(const StackSpinBlock* otherBlock, const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, int I, StackSparseMatrix* o, double scale);
  void cxddcomp(const StackSpinBlock* otherBlock, const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, int K, StackSparseMatrix* o, double scale);


  void hamandoverlap(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q, double scale, int proc);
  void cxcddcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void cdxcdcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void ddxcccomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q); 

  void CreonLeft(boost::shared_ptr<StackSparseMatrix> op1, int luncollectedQPrime, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void OverlaponLeft(boost::shared_ptr<StackSparseMatrix> op1, int luncollectedQPrime, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void CreDesonLeft(boost::shared_ptr<StackSparseMatrix> op1, int luncollectedQPrime, const StackSpinBlock* cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void CreCreonLeft(boost::shared_ptr<StackSparseMatrix> op1, int luncollectedQPrime, const StackSpinBlock* cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);

  void CreDesonRight(boost::shared_ptr<StackSparseMatrix> op1, int luncollectedQPrime, const StackSpinBlock* cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void CreCreonRight(boost::shared_ptr<StackSparseMatrix> op1, int luncollectedQPrime, const StackSpinBlock* cblock, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void CreonRight(boost::shared_ptr<StackSparseMatrix> op1, int luncollectedQPrime, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void OverlaponRight(boost::shared_ptr<StackSparseMatrix> op1, int luncollectedQPrime, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);


  void cxcddcomp_3index(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void cxcddcomp_3indexElement(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, int index, const SpinQuantum& q);
  void cdxcdcomp_3index(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void cdxcdcomp_3indexElement(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, int index, const SpinQuantum& q);
  void ddxcccomp_3indexElement(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, int index, const SpinQuantum& q); 
  void ddxcccomp_3index(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q); 

  
  //these are only used when left and right states are different
  void dcxdccomp(const StackSpinBlock* otherblock, const std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, StackSparseMatrix* o);
  void dxccdcomp(const StackSpinBlock* CDDblock, const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackSparseMatrix* o);
  void dxccdcomp(const StackSpinBlock* otherblock, const std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void dcxdccomp(const StackSpinBlock* otherblock, const std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);

  //**********************************************************************************************************
  
  void cdd_cxddcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> & op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void cdd_dxcdcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> & op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void ccd_dxcccomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> & op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void ccd_cxcdcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> & op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  
  //*********************************************************************************
  
//  void cdd_cxddcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> & op1, const StackSpinBlock* b, StackSparseMatrix* o);
//  void cdd_dxcdcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> & op1, const StackSpinBlock* b, StackSparseMatrix* o);
//  void ccd_dxcccomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> & op1, const StackSpinBlock* b, StackSparseMatrix* o);
//  void ccd_cxcdcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> & op1, const StackSpinBlock* b, StackSparseMatrix* o);

  void cdd_cxddcomp(const StackSpinBlock* otherblock, const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackSparseMatrix* o);
  void cdd_dxcdcomp(const StackSpinBlock* otherblock, const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackSparseMatrix* o);
  void ccd_dxcccomp(const StackSpinBlock* otherblock, const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackSparseMatrix* o);
  void ccd_cxcdcomp(const StackSpinBlock* otherblock, const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackSparseMatrix* o);
//  void cdd_cxddcomp(const StackSpinBlock* otherblock, const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
//  void cdd_dxcdcomp(const StackSpinBlock* otherblock, const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
//  void ccd_dxcccomp(const StackSpinBlock* otherblock, const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
//  void ccd_cxcdcomp(const StackSpinBlock* otherblock, const std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
//  
//  //*********************************************************************************
//  
//
  void CDDandoverlap(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void CCDandoverlap(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
 
}
}
#endif 
