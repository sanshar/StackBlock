/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_STACKOP_COMPONENTS_H
#define SPIN_STACKOP_COMPONENTS_H

#include <boost/function.hpp>
#include <boost/functional.hpp>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <para_array.h>
#include <boost/shared_ptr.hpp>
#include <list>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>
#include "StackOperators.h"
#include "operatorloops.h"
#include "StackMatrix.h"
#include <string>
#include "three_index_ops.h"
#include "four_index_ops.h"
#include "para_array_3d.h"
#include "para_array_4d.h"
#include <tuple>

namespace SpinAdapted{
class StackSpinBlock;
 class StateInfo;
//===========================================================================================================================================================
//choose the type of array for different types of Operators

template <class T> struct ChooseArray {
  typedef para_array_1d<std::vector<boost::shared_ptr<StackSparseMatrix> > > ArrayType;
};
template <> struct ChooseArray<StackCre> {
  typedef para_array_1d<std::vector<boost::shared_ptr<StackCre> > > ArrayType; // Cre, CreDes, etc. are sparse matrices: <a|a_i^\dagger|b>
};
template <> struct ChooseArray<StackDes> {
  typedef para_array_1d<std::vector<boost::shared_ptr<StackDes> > > ArrayType;
};
template <> struct ChooseArray<StackCreDes> {
  typedef para_array_triang_2d<std::vector<boost::shared_ptr<StackCreDes> > > ArrayType;
};
template <> struct ChooseArray<StackDesCre> {
  typedef para_array_triang_2d<std::vector<boost::shared_ptr<StackDesCre> > > ArrayType;
};
template <> struct ChooseArray<StackCreCre> {
  typedef para_array_triang_2d<std::vector<boost::shared_ptr<StackCreCre> > > ArrayType;
};
template <> struct ChooseArray<StackDesDes> {
  typedef para_array_triang_2d<std::vector<boost::shared_ptr<StackDesDes> > > ArrayType;
};

template <> struct ChooseArray<StackCreDesComp> {
    typedef para_array_triang_2d<std::vector<boost::shared_ptr<StackCreDesComp> > > ArrayType;
};
template <> struct ChooseArray<StackDesCreComp> {
    typedef para_array_triang_2d<std::vector<boost::shared_ptr<StackDesCreComp> > > ArrayType;
};
template <> struct ChooseArray<StackDesDesComp> {
  typedef para_array_triang_2d<std::vector<boost::shared_ptr<StackDesDesComp> > > ArrayType;
};
template <> struct ChooseArray<StackCreCreComp> {
  typedef para_array_triang_2d<std::vector<boost::shared_ptr<StackCreCreComp> > > ArrayType;
};
template <> struct ChooseArray<StackCreCreDesComp> {
  typedef para_array_1d<std::vector<boost::shared_ptr<StackCreCreDesComp> > > ArrayType;
};
template <> struct ChooseArray<StackCreDesDesComp> {
  typedef para_array_1d<std::vector<boost::shared_ptr<StackCreDesDesComp> > > ArrayType;
};

template <> struct ChooseArray<StackHam> {
  typedef para_array_0d<std::vector<boost::shared_ptr<StackHam> > > ArrayType;
};
template <> struct ChooseArray<StackOverlap> {
  typedef para_array_0d<std::vector<boost::shared_ptr<StackOverlap> > > ArrayType;
};


//the NPDM operators
template <> struct ChooseArray<RI3index> {
  typedef para_array_3d<std::vector<boost::shared_ptr<RI3index> > > ArrayType;
};  
template <> struct ChooseArray<RI4index> {
  typedef para_array_4d<std::vector<boost::shared_ptr<RI4index> > > ArrayType;
};
// 3PDM
template <> struct ChooseArray<StackCreCreDes> {
  typedef para_array_3d<std::vector<boost::shared_ptr<StackCreCreDes> > > ArrayType;
};
template <> struct ChooseArray<StackCreDesDes> {
  typedef para_array_3d<std::vector<boost::shared_ptr<StackCreDesDes> > > ArrayType;
};
template <> struct ChooseArray<StackCreDesCre> {
  typedef para_array_3d<std::vector<boost::shared_ptr<StackCreDesCre> > > ArrayType;
};
template <> struct ChooseArray<StackCreCreCre> {
  typedef para_array_3d<std::vector<boost::shared_ptr<StackCreCreCre> > > ArrayType;
};
// 4PDM
template <> struct ChooseArray<StackDesCreDes> {
  typedef para_array_3d<std::vector<boost::shared_ptr<StackDesCreDes> > > ArrayType;
};
template <> struct ChooseArray<StackDesDesCre> {
  typedef para_array_3d<std::vector<boost::shared_ptr<StackDesDesCre> > > ArrayType;
};
template <> struct ChooseArray<StackDesCreCre> {
  typedef para_array_3d<std::vector<boost::shared_ptr<StackDesCreCre> > > ArrayType;
};
template <> struct ChooseArray<StackDesDesDes> {
  typedef para_array_3d<std::vector<boost::shared_ptr<StackDesDesDes> > > ArrayType;
};
template <> struct ChooseArray<StackCreCreDesDes> {
  typedef para_array_4d<std::vector<boost::shared_ptr<StackCreCreDesDes> > > ArrayType;
};
template <> struct ChooseArray<StackCreDesCreDes> {
  typedef para_array_4d<std::vector<boost::shared_ptr<StackCreDesCreDes> > > ArrayType;
};
template <> struct ChooseArray<StackCreDesDesCre> {
  typedef para_array_4d<std::vector<boost::shared_ptr<StackCreDesDesCre> > > ArrayType;
};
template <> struct ChooseArray<StackCreDesDesDes> {
  typedef para_array_4d<std::vector<boost::shared_ptr<StackCreDesDesDes> > > ArrayType;
};
template <> struct ChooseArray<StackCreCreCreDes> {
  typedef para_array_4d<std::vector<boost::shared_ptr<StackCreCreCreDes> > > ArrayType;
};
template <> struct ChooseArray<StackCreCreDesCre> {
  typedef para_array_4d<std::vector<boost::shared_ptr<StackCreCreDesCre> > > ArrayType;
};
template <> struct ChooseArray<StackCreDesCreCre> {
  typedef para_array_4d<std::vector<boost::shared_ptr<StackCreDesCreCre> > > ArrayType;
};
template <> struct ChooseArray<StackCreCreCreCre> {
  typedef para_array_4d<std::vector<boost::shared_ptr<StackCreCreCreCre> > > ArrayType;
};




//===========================================================================================================================================================

class StackOp_component_base
{

 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & m_core & m_deriv;
  }

 protected:
  bool m_core;
  bool m_deriv;

 public:
  virtual long build_iterators(StackSpinBlock& b, bool calcMemory)=0;
  virtual void build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_pair, std::map< std::tuple<int, int, int>, int>& tuple)=0;
  virtual void build_operators(StackSpinBlock& b)=0;
  virtual double* allocateOperators(const StateInfo& sl, const StateInfo& sr, double* pData) = 0;
  virtual long getRequiredMemory(const StateInfo& sl, const StateInfo& sr) = 0;
  virtual void build_csf_operators(std::vector< Csf >& dets, std::vector< std::vector<Csf> >& ladders, StackSpinBlock& b) = 0;
  virtual void build_and_renormalise_operators(StackSpinBlock&b, const opTypes &ot, const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo) =0;
  virtual void build_and_renormalise_operators(StackSpinBlock&b, const opTypes &ot, const std::vector<Matrix>& leftMat, const StateInfo *bra, const std::vector<Matrix>& rightMat, const StateInfo *ket) =0;

  //virtual string type_name() = 0;
  virtual int get_size() const =0;
  virtual int size() const=0;
  virtual void clear() =0;
  virtual std::vector<boost::shared_ptr<StackSparseMatrix> > get_local_element(int i) =0;
  virtual std::vector<boost::shared_ptr<StackSparseMatrix> > get_global_element(int i)=0;
  const bool &is_core() const {return m_core;}
  const bool &is_deriv() const {return m_deriv;}
  void set_core(bool is_core) {m_core = is_core;}
  virtual void add_local_indices(int i, int j=-1, int k=-1) {};
  virtual void remove_local_indices(int i, int j=-1, int k=-1, int l=-1) = 0;
  virtual bool is_local() const = 0;
  virtual bool& set_local() = 0; 
  virtual void set_length(int length) {};
  virtual std::vector< std::vector<int> > get_array() const =0;
  virtual std::vector< int > get_global_array() const =0;
  virtual const std::vector<boost::shared_ptr<StackSparseMatrix> > get_element(int i, int j=-1, int k=-1, int l=-1) const = 0;
  virtual std::vector<boost::shared_ptr<StackSparseMatrix> > get_element(int i, int j=-1, int k=-1, int l=-1) = 0;
  virtual bool has(int i, int j=-1, int k=-1, int l=-1) const = 0;
  virtual bool has_local_index(int i, int j=-1, int k=-1, int l=-1) const = 0;
  virtual boost::shared_ptr<StackSparseMatrix> get_op_rep(const std::vector<SpinQuantum>& s, int i=-1, int j=-1, int k=-1, int l=-1) = 0;
  virtual const boost::shared_ptr<StackSparseMatrix> get_op_rep(const std::vector<SpinQuantum>& s, int i=-1, int j=-1, int k=-1, int l=-1) const = 0;
  virtual boost::shared_ptr<StackSparseMatrix> get_op_rep(const std::map< std::string, std::vector<SpinQuantum> >& s, int i=-1, int j=-1, int k=-1, int l=-1) = 0;
  virtual const boost::shared_ptr<StackSparseMatrix> get_op_rep(const std::map< std::string, std::vector<SpinQuantum> >& s, int i=-1, int j=-1, int k=-1, int l=-1) const = 0;
  virtual std::string get_op_string() const = 0;
  virtual std::string get_filename() const = 0;
  virtual ~StackOp_component_base() {}  

};

//===========================================================================================================================================================

template <class Op> class StackOp_component : public StackOp_component_base
{

 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<StackOp_component_base>(*this);
    ar.register_type(static_cast<Op *>(NULL));
    ar & m_op & uniqueID;
  }

 protected:
  typedef typename ChooseArray<Op>::ArrayType paraarray;
  typedef Op OpType; 
  paraarray m_op;
  int uniqueID;
  // Use for unique filename for NPDM disk-based operator storage 
  static int nIDgenerator;

 public:
  StackOp_component() { m_deriv=false; uniqueID = nIDgenerator++; }
  StackOp_component(bool core) { m_core=core; m_deriv=false; uniqueID = nIDgenerator++; }
  std::string get_op_string() const;
  bool& set_local() {return m_op.set_local();}
  bool is_local() const {return m_op.is_local();}
  int get_size() const {return m_op.local_nnz();}
  int size() const  {return m_op.global_nnz();}
  bool has(int i, int j=-1, int k=-1, int l=-1) const {return m_op.has(i, j, k, l);}
  bool has_local_index(int i, int j=-1, int k=-1, int l=-1) const {return m_op.has_local_index(i, j, k, l);}
  virtual void add_local_indices(int i, int j=-1, int k=-1){};
  virtual void set_length(int length) {m_op.set_length(length);}
  void remove_local_indices(int i, int j=-1, int k=-1, int l=-1) {m_op.remove_local_index(i,j,k,l);}
  void clear(){m_op.clear();}

  long build_iterators(StackSpinBlock& b, bool calcMemory);
  void build_iterators(StackSpinBlock& b, std::vector<int>& screened_c_ix, std::vector<std::pair<int, int> >& screened_pair, std::map< std::tuple<int, int, int>, int>& tuple);
  void build_operators(StackSpinBlock& b) 
    { singlethread_build(*this, b); }

  // Note for NPDM higher-index operators, there are template specializations for these functions
  void build_and_renormalise_operators(StackSpinBlock&b, const opTypes &ot, const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo)
    {for_all_operators_singlethread(*this, bind(&StackSparseMatrix::build_and_renormalise_transform, _1, &b, boost::ref(rotateMatrix), stateinfo));}

  void build_and_renormalise_operators(StackSpinBlock&b, const opTypes &ot, const std::vector<Matrix>& leftMat, const StateInfo *bra, const std::vector<Matrix>& rightMat, const StateInfo *ket)
    {for_all_operators_singlethread(*this, bind(&StackSparseMatrix::build_and_renormalise_transform, _1, &b, boost::ref(leftMat), bra, boost::ref(rightMat), ket));}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Use for unique filename for NPDM disk-based operator storage
  std::string get_filename() const
  {
    std::string file;
    file = str( boost::format("%s%s%s%s%d%d%s") % dmrginp.load_prefix() % "/" % get_op_string() % "_p" % mpigetrank() % get_size() % ".tmp" );
    return file;
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  void build_csf_operators(std::vector< Csf >& c, vector< vector<Csf> >& ladders, StackSpinBlock& b) 
  {
    for_all_operators_singlethread( *this, bind(&StackSparseMatrix::buildUsingCsf, _1, boost::ref(b), boost::ref(ladders), boost::ref(c)) );
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  double* allocateOperators(const StateInfo& sl, const StateInfo& sr, double* pData) {
    double* data  = pData;
    for (int i=0; i<get_size(); i++) {
      int vecsize = get_local_element(i).size();
      for (int j=0; j< vecsize; j++) {
	get_local_element(i)[j]->getrowCompressedForm().resize(0);
	double* nextData = get_local_element(i)[j]->allocate(sl, sr, data);
	data = nextData;
      }
    }
    return data;
  }


  long getRequiredMemory(const StateInfo& sl, const StateInfo& sr) {
    long mem = 0;
    for (int i=0; i<get_size(); i++) {
      int vecsize = get_local_element(i).size();
      for (int j=0; j< vecsize; j++) {
	mem += SpinAdapted::getRequiredMemory(sl, sr, get_local_element(i)[j]->get_deltaQuantum());
      }
    }
    return mem;
  }



//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  std::vector<boost::shared_ptr<StackSparseMatrix> > get_local_element(int i) 
  {
    std::vector<boost::shared_ptr<StackSparseMatrix> > vec(m_op.get_local_element(i).size());
    for (int l=0; l<vec.size(); l++)
      vec[l] = m_op.get_local_element(i)[l]; 
    return vec;
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  std::vector<boost::shared_ptr<StackSparseMatrix> > get_global_element(int i)
  {
    std::vector<boost::shared_ptr<StackSparseMatrix> > vec(m_op.get_global_element(i).size());
    for (int l=0; l<vec.size(); l++)
      vec[l] = m_op.get_global_element(i)[l]; 
    return vec;
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  std::vector< std::vector<int> > get_array() const 
  {
    std::vector<int> orbs(2);
    std::vector< std::vector<int> > ret_val(m_op.local_nnz());
    for (int i=0; i<m_op.local_nnz(); i++)
      ret_val[i] = m_op.unmap_local_index(i);
    return ret_val;
  }

  std::vector< int > get_global_array() const 
  {
    return m_op.get_global_array();
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  const std::vector<boost::shared_ptr<StackSparseMatrix> >  get_element(int i, int j=-1, int k=-1, int l=-1) const 
  {
    std::vector<boost::shared_ptr<StackSparseMatrix> > vec(m_op(i,j,k,l).size());
    for (int p=0; p<vec.size(); p++)
      vec[p] = m_op(i,j,k,l)[p]; 
    return vec;
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  std::vector<boost::shared_ptr<StackSparseMatrix> >  get_element(int i, int j=-1, int k=-1, int l=-1)
  {
    std::vector<boost::shared_ptr<StackSparseMatrix> > vec(m_op(i,j,k,l).size());
    for (int p=0; p<vec.size(); p++)
      vec[p] = m_op(i,j,k,l)[p]; 
    return vec;
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  boost::shared_ptr<StackSparseMatrix> get_op_rep(const std::vector<SpinQuantum>& s, int i=-1, int j=-1, int k=-1, int l=-1)
  {
    Op* o = 0;
    std::vector<boost::shared_ptr<Op> >& vec = m_op(i,j,k,l);
    for (int p=0; p<vec.size(); p++) {
      if (s == vec[p]->get_deltaQuantum())
	    return m_op(i,j,k,l)[p];
    }
    return boost::shared_ptr<Op>(o);
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  const boost::shared_ptr<StackSparseMatrix> get_op_rep(const std::vector<SpinQuantum>& s, int i=-1, int j=-1, int k=-1, int l=-1) const
  {
    Op* o = 0;
    const std::vector<boost::shared_ptr<Op> >& vec = m_op(i,j,k,l);
    for (int p=0; p<vec.size(); p++)
      if (s == vec[p]->get_deltaQuantum())
	    return m_op(i,j,k,l)[p];
    return boost::shared_ptr<Op>(o);
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  boost::shared_ptr<StackSparseMatrix> get_op_rep(const std::map< std::string, std::vector<SpinQuantum> >& s, int i=-1, int j=-1, int k=-1, int l=-1)
  {
    Op* o = 0;
    std::vector<boost::shared_ptr<Op> >& vec = m_op(i,j,k,l);
    for (int p=0; p<vec.size(); p++) {
      if (s == vec[p]->get_quantum_ladder())
	    return m_op(i,j,k,l)[p];
    }
    return boost::shared_ptr<Op>(o);
  }

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

  const boost::shared_ptr<StackSparseMatrix> get_op_rep(const std::map< std::string, std::vector<SpinQuantum> >& s, int i=-1, int j=-1, int k=-1, int l=-1) const
  {
    Op* o = 0;
    const std::vector<boost::shared_ptr<Op> >& vec = m_op(i,j,k,l);
    for (int p=0; p<vec.size(); p++)
      if (s == vec[p]->get_quantum_ladder())
	    return m_op(i,j,k,l)[p];
    return boost::shared_ptr<Op>(o);
  }

};

//===========================================================================================================================================================

template<>
  void StackOp_component<StackCreCreCre>::build_and_renormalise_operators(StackSpinBlock&b, const opTypes &ot, const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo) ;
template<>
  void StackOp_component<StackCreCreDes>::build_and_renormalise_operators(StackSpinBlock&b, const opTypes &ot, const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo) ;
template<>
  void StackOp_component<StackCreDesDes>::build_and_renormalise_operators(StackSpinBlock&b, const opTypes &ot, const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo) ;
template<>
  void StackOp_component<StackCreDesCre>::build_and_renormalise_operators(StackSpinBlock&b, const opTypes &ot, const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo) ;
 
template <class Op> int StackOp_component<Op>::nIDgenerator = 1;

}

#endif
