/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_STACKWAVEFUNCTION_HEADER
#define SPIN_STACKWAVEFUNCTION_HEADER
#include "StackBaseOperator.h"

namespace SpinAdapted{
  class StackSpinBlock;

class StackWavefunction : public StackSparseMatrix
{
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<StackSparseMatrix>(*this);
    ar & onedot;
  } 
  bool onedot;
  //ObjectMatrix<StackMatrix> operatorMatrix;  // The StackMatrix does not own its data
public:
 StackWavefunction() : onedot(false), StackSparseMatrix(){}
 StackWavefunction(const StackWavefunction& wf) : StackSparseMatrix(wf), onedot(wf.onedot){}

  void copyData(const StackWavefunction& a);
  void operator=(const StackWavefunction& wf) {StackSparseMatrix::operator=(wf); onedot=wf.onedot;}
  //ObjectMatrix<StackMatrix>& get_operatorMatrix() {return operatorMatrix;}
  //void OperatorMatrixReference(ObjectMatrix<StackMatrix*>& m, const std::vector<int>& oldToNewStateI, const std::vector<int>& oldToNewStateJ);
  //const ObjectMatrix<StackMatrix>& get_operatorMatrix() const {return operatorMatrix;}
  void initialise(const vector<SpinQuantum>& dQ, const StateInfo& sl, const StateInfo& sr, const bool &onedot_, double* pData, long ptotalMemory);  
  void initialise(const vector<SpinQuantum>& dQ, const StateInfo& sl, const StateInfo& sr, const bool &onedot_);  
  void initialise(const StackWavefunction& w);

  virtual void deepCopy(const StackWavefunction& o) ;
  virtual void deepClearCopy(const StackWavefunction& o) ;
  const bool &get_onedot() const {return onedot;}
  void set_onedot(bool p_onedot) {onedot = p_onedot;}
  bool& set_onedot() {return onedot;}

  void build(const StackSpinBlock& block) {return;}
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b=0) {return 0.0;}

  void CollectFrom(const RowVector& C);
  void FlattenInto(Matrix& C);

  static bool exists(int state);
  static void CopyState(int from, int to);
  static void ChangeLastSite(int newLast, int oldLast, int state);
  void LoadWavefunctionInfo (StateInfo &waveInfo, const std::vector<int>& sites, const int wave_num, bool allocateData=false);
  void SaveWavefunctionInfo (const StateInfo &waveInfo, const std::vector<int>& sites, const int wave_num);
  double* allocateWfnOperatorMatrix();

  void UnCollectQuantaAlongRows(const StateInfo& sRow, const StateInfo& sCol);
  void CollectQuantaAlongRows(const StateInfo& sRow, const StateInfo& sCol); // FIXME what does this function do?
  void CollectQuantaAlongColumns(const StateInfo& sRow, const StateInfo& sCol);
  void UnCollectQuantaAlongColumns(const StateInfo& sRow, const StateInfo& sCol);
  StackWavefunction& operator+=(const StackWavefunction& other);
};
}
#endif
