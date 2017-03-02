#ifndef STACKSPINBLOCK_HEADER
#define STACKSPINBLOCK_HEADER
#include <list>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include "Stack_op_components.h"
#include "multiarray.h"
#include "boost/variant.hpp"
#include "sweep_params.h"
#include <newmat.h>
#include "perturb.h"


namespace SpinAdapted{
class StackWavefunction;
class StackDensityMatrix;
class Csf;
class StateInfo;

boost::shared_ptr<StackOp_component_base> make_new_stackop(const opTypes &optype, const bool &is_core);

class StackSpinBlock
{
 private:
  friend class boost::serialization::access;
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & localstorage;
      ar & name;
      ar & complementary;
      ar & normal; 
      ar & direct;
      ar & loopblock;
      ar & sites;
      ar & complementary_sites ;
      ar & integralIndex;
      //FIX ME!! remove register_type stuff and add BOOST_CLASS_EXPORT to op_components.h (will take longer to compile)                     
      ar.register_type(static_cast<StackOp_component<StackCre> *>(NULL));
      ar.register_type(static_cast<StackOp_component<StackDes> *>(NULL));

      ar.register_type(static_cast<StackOp_component<StackCreDes> *>(NULL));
      ar.register_type(static_cast<StackOp_component<StackDesCre> *>(NULL));
      ar.register_type(static_cast<StackOp_component<StackCreCre> *>(NULL));
      ar.register_type(static_cast<StackOp_component<StackDesDes> *>(NULL));

      ar.register_type(static_cast<StackOp_component<StackCreDesComp> *>(NULL));
      ar.register_type(static_cast<StackOp_component<StackDesCreComp> *>(NULL));
      ar.register_type(static_cast<StackOp_component<StackDesDesComp> *>(NULL));
      ar.register_type(static_cast<StackOp_component<StackCreCreComp> *>(NULL));
      ar.register_type(static_cast<StackOp_component<StackCreCreDesComp> *>(NULL));
      ar.register_type(static_cast<StackOp_component<StackCreDesDesComp> *>(NULL));

      ar.register_type(static_cast<StackOp_component<StackHam> *>(NULL));
      ar.register_type(static_cast<StackOp_component<StackOverlap> *>(NULL));
      ar.register_type(static_cast<StackOp_component<StackCDD_CreDesComp> *>(NULL));
      ar.register_type(static_cast<StackOp_component<StackCDD_DesDesComp> *>(NULL));
      ar.register_type(static_cast<StackOp_component<StackCCD_sum> *>(NULL));
      ar.register_type(static_cast<StackOp_component<StackCCD_CreDesComp> *>(NULL));
      ar.register_type(static_cast<StackOp_component<StackCCD_CreCreComp> *>(NULL));

      ar & ops;
    }

 private:

  std::map<opTypes, boost::shared_ptr<StackOp_component_base> > ops;
  bool complementary;
  bool normal;
  bool loopblock;
  bool localstorage;
  bool direct;
  int name;
  int integralIndex;
  StackSpinBlock* leftBlock;
  StackSpinBlock* rightBlock;
  boost::shared_ptr<TwoElectronArray> twoInt;
  
  StateInfo braStateInfo;
  StateInfo ketStateInfo;
  std::vector<int> sites;
  std::vector<int> complementary_sites;
  std::vector<int> nonactive_orbs;
  long totalMemory;
  double* data;
  long additionalMemory;
  double* additionaldata;
 public: 

  StackSpinBlock (const StateInfo& s, int integralIndex);

  //These are the constructors
  StackSpinBlock();

  //makes a shallow copy
  StackSpinBlock (const StackSpinBlock& b);

  //can only be called after the data has been initialized
  StackSpinBlock (int start, int finish, int integralIndex, bool implicitTranspose, bool is_complement = false);
  StackSpinBlock (int start, int finish, int integralIndex, std::vector<int> non_active_orbs, bool implicitTranspose, bool is_complement = false);



  //can only be called before the data is initialized
  void BuildTensorProductBlock (std::vector<int>& new_sites);
  void recreateStateInfo(int condition);
  void collectQuanta();

  //These functions need to be called before the required memory can be calculated
  //These functions decide what operators are to be built and stored in memory
  void default_op_components(bool complementary_, bool implicitTranspose);
  void default_op_components(bool direct, bool haveNormops, bool haveCompops, bool implicitTranspose);
  void nevpt_op_components(bool direct, StackSpinBlock& lBlock, StackSpinBlock& rBlock, const perturber& pb);
  void setstoragetype(Storagetype st);
  void set_big_components();
  void initialise_op_array(opTypes optype, bool is_core);


  //these will be tricky to implement
  void addAdditionalOps() ;
  void addAllCompOps() ;
  void addOneIndexNormOps() ;
  void addOneIndexOps() ;
  void messagePassTwoIndexOps() ;
  void formTwoIndexOps() ;
  void removeAdditionalOps() ;
  void sendcompOps(StackOp_component_base& opcomp, int I, int J, int optype, int compsite) ;
  void recvcompOps(StackOp_component_base& opcomp, int I, int J, int optype) ;

  //simple functions
  const boost::shared_ptr<TwoElectronArray> get_twoInt() const {return twoInt;}
  void set_twoInt(int integralIndex);
  int get_integralIndex() const {return integralIndex;}
  int& set_integralIndex() {return integralIndex;}
  const StateInfo& get_stateInfo() const {return ketStateInfo;}
  const StateInfo& get_braStateInfo() const {return braStateInfo;}
  const StateInfo& get_ketStateInfo() const {return ketStateInfo;}
  StateInfo& set_braStateInfo() {return braStateInfo;}
  StateInfo& set_ketStateInfo() {return ketStateInfo;}
  static std::vector<int> make_complement(const std::vector<int>& sites);
  void printOperatorSummary();
  int size() const { return sites.size(); }
  int get_name() const {return name;}
  std::vector<int>& set_sites() {return sites;}
  const std::vector<int>& get_sites() const {return sites;}
  const std::vector<int>& get_complementary_sites() const {return complementary_sites;}
  bool is_normal() const {return normal;}
  bool is_complementary() const {return complementary;}
  bool is_loopblock() const {return loopblock;}
  bool is_direct() const {return direct;}
  bool getlocalstorage() const {return localstorage;}
  StackSpinBlock* get_leftBlock() const {return leftBlock;}
  StackSpinBlock* get_rightBlock() const {return rightBlock;}
  const vector<int>& nonactive_orb() const { return nonactive_orbs;}
  vector<int>& nonactive_orb() { return nonactive_orbs;}
  const int& nonactive_orb(int i) const { return nonactive_orbs[i];}
  int& nonactive_orb(int i) { return nonactive_orbs[i];}

  //related to retriving information from oparray
  void CleanUpOperators();
  void erase(opTypes optype) {assert(has(optype)); ops.erase(optype);}
  boost::shared_ptr<StackOp_component_base>& set_op_array(opTypes optype){assert(has(optype));return ops.find(optype)->second;}
  StackOp_component_base& get_op_array(opTypes optype){assert(has(optype));return *(ops.find(optype)->second);}
  const StackOp_component_base& get_op_array(opTypes optype) const {assert(has(optype));return *(ops.find(optype)->second);}
  bool has(opTypes optype) const
  {
    if(ops.find(optype) != ops.end())
      return true;
    else
      return false;
  }  
  boost::shared_ptr<StackSparseMatrix> get_op_rep(const opTypes &optypes, const SpinQuantum& s, int i=-1, int j=-1, int k=-1,int l = -1) {
    assert(has(optypes));
    StackOp_component_base& opbase = *ops.find(optypes)->second;
    vector<SpinQuantum> temp(1, s);
    return opbase.get_op_rep(temp, i, j, k, l);
  }
  const boost::shared_ptr<StackSparseMatrix> get_op_rep(const opTypes &optypes, const SpinQuantum& s, int i=-1, int j=-1, int k=-1, int l=-1) const {
    assert(has(optypes));
    StackOp_component_base& opbase = *ops.find(optypes)->second;
    vector<SpinQuantum> temp(1, s);    
    return opbase.get_op_rep(temp, i, j, k, l);
  }
  boost::shared_ptr<StackSparseMatrix> get_op_rep(const opTypes &optypes, const std::vector<SpinQuantum>& s, int i=-1, int j=-1, int k=-1, int l=-1) {
    assert(has(optypes));
    StackOp_component_base& opbase = *ops.find(optypes)->second; 
    return opbase.get_op_rep(s, i, j, k, l);
  }
  const boost::shared_ptr<StackSparseMatrix> get_op_rep(const opTypes &optypes, const vector<SpinQuantum>& s, int i=-1, int j=-1, int k=-1, int l= -1) const {
    assert(has(optypes));
    StackOp_component_base& opbase = *ops.find(optypes)->second;
    return opbase.get_op_rep(s, i, j, k, l);
  }
  boost::shared_ptr<StackSparseMatrix> get_op_rep(const opTypes &optypes, const std::map< std::string, std::vector<SpinQuantum> >& s, int i=-1, int j=-1, int k=-1, int l= -1) {
    assert(has(optypes));
    StackOp_component_base& opbase = *ops.find(optypes)->second;
    return opbase.get_op_rep(s, i, j, k, l);
  }
  const boost::shared_ptr<StackSparseMatrix> get_op_rep(const opTypes &optypes, const std::map< std::string, std::vector<SpinQuantum> >& s, int i=-1, int j=-1, int k=-1, int l= -1) const {
    assert(has(optypes));
    StackOp_component_base& opbase = *ops.find(optypes)->second;
    return opbase.get_op_rep(s, i, j, k, l);
  }

  long memoryUsed() {return totalMemory;}
  long additionalMemoryUsed() {return additionalMemory;}
  void deallocate();
  double* getdata() {return data;}
  void moveToNewMemory(double* pData);
  void operator= (const StackSpinBlock& b);
  long build_iterators();
  void build_operators(std::vector<Csf >& s, std::vector< std::vector<Csf> >& ladders);
  void build_operators();
  void deallocate_coreops();
  void build_and_renormalise_operators(const std::vector<Matrix>& rotateMatrix, const StateInfo *newStateInfo);
  void build_and_renormalise_operators(const std::vector<Matrix>& leftMat, const StateInfo *bra, const std::vector<Matrix>& rightMat, const StateInfo *ket);
  void renormalise_transform(const std::vector<Matrix>& rotateMatrix, const StateInfo *stateinfo);
  void renormalise_transform(const std::vector<Matrix>& leftMat, const StateInfo *bra, const std::vector<Matrix>& rightMat, const StateInfo *ket);

  void BuildSumBlock(int condition, StackSpinBlock& b_1, StackSpinBlock& b_2, bool collectQuanta = true, StateInfo* compState=0);
  void BuildSumBlock(int condition, StackSpinBlock& lBlock, StackSpinBlock& rBlock, const std::vector<SpinQuantum>& braquantum, const std::vector<SpinQuantum>& ketquantum, bool collectQuanta= true);
  void BuildSumBlockSkeleton(int condition, StackSpinBlock& lBlock, StackSpinBlock& rBlock, bool collectQuanta = true, StateInfo* compState=0);
  void BuildSumBlockSkeleton(int condition, StackSpinBlock& lBlock, StackSpinBlock& rBlock, const std::vector<SpinQuantum>& braquantum, const std::vector<SpinQuantum>& ketquantum, bool collectQuanta= true);
  void BuildSlaterBlock (std::vector<int> sts, std::vector<SpinQuantum> qnumbers, std::vector<int> distribution, bool random, 
			 const bool haveNormops);
  void BuildSingleSlaterBlock(std::vector<int> sts);
  void set_loopblock(bool p_loopblock){loopblock = p_loopblock;}
  friend ostream& operator<< (ostream& os, const StackSpinBlock& b);
  void multiplyH(StackWavefunction& c, StackWavefunction* v, int num_threads) const;
  void multiplyH_2index(StackWavefunction& c, StackWavefunction* v, int num_threads) const;
  void multiplyOverlap(StackWavefunction& c, StackWavefunction* v, int num_threads) const;
  void multiplyCDD_sum(StackWavefunction& c, StackWavefunction* v, int num_threads) const;
  void multiplyCCD_sum(StackWavefunction& c, StackWavefunction* v, int num_threads) const;
  void diagonalH(DiagonalMatrix& e) const;
  void clear();

  void RenormaliseFrom (std::vector<double> &energies, std::vector<double> &spins, double &error, std::vector<Matrix>& rotateMatrix,
                        const int keptstates, const int keptqstates, const double tol, StackSpinBlock& big,
                        const guessWaveTypes &guesswavetype, const double noise, const double additional_noise, const bool &onedot, StackSpinBlock& system, 
			StackSpinBlock& sysDot, StackSpinBlock& environment, const bool& dot_with_sys, const bool& warmUp, int sweepiter, 
			int currenroot, std::vector<StackWavefunction>& lowerStates, StackDensityMatrix* d=0);

  void transform_operators(std::vector<Matrix>& rotateMatrix);
  void transform_operators(std::vector<Matrix>& leftrotateMatrix, std::vector<Matrix>& rightrotateMatrix, bool clearRightBlock = true, bool clearLeftBlock = true);

  static StackSpinBlock buildBigEdgeBlock(int start, int end, bool haveNorm, bool haveComp, int p_integralIndex, bool implicitTranspose);
    //static StackSpinBlock buildBigEdgeBlock (int start, int finish, int p_integralIndex, bool implicitTranspose);
  static void restore (bool forward, const vector<int>& sites, StackSpinBlock& b, int left, int right, char* name=0);//left and right are the bra and ket states and the name is the type of the MPO (currently only H)
  static void make_iterator(StackSpinBlock& b, opTypes op, int* data, int oneIndex, int numIndices) ;
  static void store (bool forward, const vector<int>& sites, StackSpinBlock& b, int left, int right, char* name=0);//left and right are the bra and ket states and the name is the type of the MPO (currently only H) 
  void Save (std::ofstream &ofs);
  void Load (std::ifstream &ifs);


};

 double makeRotateMatrix(StackDensityMatrix& tracedMatrix, vector<Matrix>& rotateMatrix, const int& keptstates, const int& keptqstates, std::vector<DiagonalMatrix> *eigs =0);
 void  initialiseSingleSiteBlocks(std::vector<StackSpinBlock>& stackblocks, int integralIndex);
}
#endif
