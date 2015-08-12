#ifndef SPIN_STACKOPXOP_HEADER
#define SPIN_STACKOPXOP_HEADER
 
#include <newmat.h>

namespace SpinAdapted{
  class StackWavefunction;
  class StackSpinBlock;
  class StackSparseMatrix;

namespace stackopxop
{
  void cxcddcomp(const StackSpinBlock* otherblock, std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, StackSparseMatrix* o);
  void cdxcdcomp(const StackSpinBlock* otherblock, std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, StackSparseMatrix* o);
  void ddxcccomp(const StackSpinBlock* otherblock, std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, StackSparseMatrix* o);


  void cdxcdcomp_d(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, DiagonalMatrix* e);
  void ham_d(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, DiagonalMatrix* e, int proc);
  
 
  void cxcdcomp(const StackSpinBlock* otherblock, std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, int I, StackSparseMatrix* o, double factor);
  void dxcccomp(const StackSpinBlock* otherBlock, std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, int I, StackSparseMatrix* o, double scale);

  void dxcdcomp(const StackSpinBlock* otherBlock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, int I, StackSparseMatrix* o, double scale);
  void cxddcomp(const StackSpinBlock* otherBlock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, int K, StackSparseMatrix* o, double scale);


  void hamandoverlap(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q, double scale, int proc);
  void cxcddcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void cdxcdcomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void ddxcccomp(const StackSpinBlock* otherblock, boost::shared_ptr<StackSparseMatrix> op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q); 

  
  //these are only used when left and right states are different
  void dcxdccomp(const StackSpinBlock* otherblock, std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, StackSparseMatrix* o);
  void dxccdcomp(const StackSpinBlock* CDDblock, std::vector<boost::shared_ptr<StackSparseMatrix> >& opvec1, const StackSpinBlock* b, StackSparseMatrix* o);
  void dxccdcomp(const StackSpinBlock* otherblock, std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
  void dcxdccomp(const StackSpinBlock* otherblock, std::vector< boost::shared_ptr<StackSparseMatrix> >& op1, const StackSpinBlock* b, StackWavefunction& c, StackWavefunction* v, const SpinQuantum& q);
 
}
}
#endif 
