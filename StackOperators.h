#ifndef SPIN_STACKOPERATORS_HEADER
#define SPIN_STACKOPERATORS_HEADER
#include "StackBaseOperator.h"
#include <boost/function.hpp>
#include <boost/functional.hpp>
#include <boost/bind.hpp>

typedef boost::function<void (std::vector<boost::shared_ptr<SpinAdapted::StackSparseMatrix> >)> FUNCTOR;
typedef boost::function<void (boost::shared_ptr<SpinAdapted::StackSparseMatrix>)> FUNCTOR2;
typedef boost::function<void (boost::shared_ptr<SpinAdapted::StackSparseMatrix>, int)> FUNCTOR3;

namespace SpinAdapted{
  class StackSpinBlock;

class StackCre: public SpinAdapted::StackSparseMatrix
{
 public:
  StackCre() { conj='n';fermion = true;}
  void build(const StackSpinBlock& block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
  virtual string opName() const {return "CRE";}
  void build(StackMatrix &m, int row, int col, const StackSpinBlock& block) ;
};

class StackDes: public SpinAdapted::StackSparseMatrix
{
 public:
  StackDes() { conj='n'; fermion = true;}
  void build(const StackSpinBlock& block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

class StackCreDes: public SpinAdapted::StackSparseMatrix
{
 public:
  StackCreDes() { conj='n'; fermion = false;}
  void build(const StackSpinBlock& block);
  void buildUsingCre(const StackSpinBlock* b);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
  void build(StackMatrix &m, int row, int col, const StackSpinBlock& block) ;
};

class StackDesCre: public SpinAdapted::StackSparseMatrix
{
 public:
  StackDesCre() { conj='n'; fermion = false;}
  void build(const StackSpinBlock& block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

class StackCreCre: public SpinAdapted::StackSparseMatrix
{
 public:
  StackCreCre() { conj='n'; fermion = false;}
  void build(const StackSpinBlock& block);
  void buildUsingCre(const StackSpinBlock* b);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
  void build(StackMatrix &m, int row, int col, const StackSpinBlock& block) ;
};

class StackDesDes: public SpinAdapted::StackSparseMatrix
{
 public:
  StackDesDes() { conj='n'; fermion = false;}
  void build(const StackSpinBlock& block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

class StackCreDesComp: public SpinAdapted::StackSparseMatrix
{
 public:
  StackCreDesComp() { conj='n'; fermion = false;}
  void build(const StackSpinBlock& block);
  void buildfromCreDes(StackSpinBlock& block);
  void buildUsingCre(const StackSpinBlock* b);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
  void build(StackMatrix &m, int row, int col, const StackSpinBlock& block);
};

class StackCreCreComp: public SpinAdapted::StackSparseMatrix
{
 public:
  StackCreCreComp() { conj='n'; fermion = false;}
  void build(const StackSpinBlock& block);
  void buildfromCreCre(StackSpinBlock& block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

class StackDesDesComp: public SpinAdapted::StackSparseMatrix
{
 public:
  StackDesDesComp() { conj='n'; fermion = false;}
  void build(const StackSpinBlock& block);
  void buildfromDesDes(StackSpinBlock& block);
  //void buildUsingCre(const StackSpinBlock* b);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
  void build(StackMatrix &m, int row, int col, const StackSpinBlock& block) ;
};

class StackDesCreComp: public SpinAdapted::StackSparseMatrix
{
 public:
  StackDesCreComp() { conj='n'; fermion = false;}
  void build(const StackSpinBlock& block);
  void buildfromDesCre(StackSpinBlock& block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

class StackCreDesDesComp: public SpinAdapted::StackSparseMatrix
{
 public:
  StackCreDesDesComp() { conj='n'; fermion = true;}
  void build(const StackSpinBlock& block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
};

class StackCreCreDesComp: public SpinAdapted::StackSparseMatrix
{
 public:
  StackCreCreDesComp() { conj='n'; fermion = true;}
  void build(const StackSpinBlock& block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
  void build(StackMatrix &m, int row, int col, const StackSpinBlock& block) ;
};


class StackHam: public SpinAdapted::StackSparseMatrix
{
 public:
  StackHam() { conj='n'; fermion = false;}
  void build(const StackSpinBlock& block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
  void build(StackMatrix &m, int row, int col, const StackSpinBlock& block) ;
};

class StackOverlap: public SpinAdapted::StackSparseMatrix
{
 public:
  StackOverlap() { conj='n'; fermion = false;}
  void build(const StackSpinBlock& block);
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b);
  virtual string opName() const {return "OVERLAP";}
  void build(StackMatrix &m, int row, int col, const StackSpinBlock& block) ;
};

}


#endif
