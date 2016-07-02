/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef SPIN_STACKBASEOPERATOR_HEADER
#define SPIN_STACKBASEOPERATOR_HEADER
#include "StackMatrix.h"
#include "ObjectMatrix.h"
#include "csf.h"
#include "StateInfo.h"
#include "boostutils.h"

namespace SpinAdapted{
  class StackSpinBlock;

template<class T> class Baseoperator  // The abstract class of an operator
{
 public:
  virtual bool get_fermion() const = 0;
  virtual int nrows() const = 0;
  virtual int ncols() const = 0;
  virtual bool get_initialised() const = 0;
  virtual const std::vector<int>& get_orbs() const = 0;
  virtual int get_orbs(int i) const = 0;
  virtual const char& allowed(int i, int j) const = 0;
  virtual char& allowed(int i, int j) = 0;
  virtual T& operator_element(int i, int j) = 0;
  virtual const T& operator_element(int i, int j) const = 0;
  virtual T& operator()(int i, int j) = 0;
  virtual const T& operator()(int i, int j) const = 0;
  virtual const char& allowed(int i, int j, char conj) const = 0;
  virtual char& allowed(int i, int j, char conj) = 0;
  virtual T& operator_element(int i, int j, char conj) = 0;
  virtual const T& operator_element(int i, int j, char conj) const = 0;
  virtual T& operator()(int i, int j, char conj) = 0;
  virtual const T& operator()(int i, int j, char conj) const = 0;
  virtual int get_deltaQuantum_size() const = 0;
  virtual std::vector<SpinQuantum> get_deltaQuantum() const = 0;
  virtual SpinQuantum get_deltaQuantum(int i) const = 0;  
  virtual char conjugacy() const = 0;
  virtual ~Baseoperator() {};
  virtual SpinSpace get_spin(int i=0) const = 0;
  virtual IrrepSpace get_symm(int i=0) const = 0;
  virtual double get_scaling(SpinQuantum leftq, SpinQuantum rightq) const = 0;
  virtual const std::vector<int>& getActiveRows(int i) const =0;
  virtual const std::vector<int>& getActiveCols(int i) const =0;
  virtual std::vector<int>& getActiveRows(int i)  =0;
  virtual std::vector<int>& getActiveCols(int i) =0;
  Baseoperator() {};
};


class StackSparseMatrix : public Baseoperator<StackMatrix>  // the sparse matrix representation of the operator
{
 private:
  friend class boost::serialization::access;
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & deltaQuantum \
         & quantum_ladder \
         & build_pattern \
         & fermion \
         & initialised \
         & built \
	 & built_on_disk \
	 & allowedQuantaMatrix \
         & Sign \
	& orbs \
	& rowCompressedForm \
	& colCompressedForm \
	& nonZeroBlocks \
	& mapToNonZeroBlocks \
	& conj		     \
	& filename \
	& totalMemory;
    }

 protected:
  long totalMemory; //the length of the data
  double* data;    //the data which this object does not own

  char conj;
  std::vector<int> orbs;
  bool fermion;
  ObjectMatrix<char> allowedQuantaMatrix;  // some whether it is allowed?
  bool initialised;
  bool built;
  bool built_on_disk;
  std::vector<SpinQuantum> deltaQuantum;    // allowed quantum
  int Sign;
  // With N-index ops (N>2), there are several ways to build them...
  std::string build_pattern;
  // ...and for each way we record the spin ladder components
  std::map< std::string, std::vector<SpinQuantum> > quantum_ladder;

  string filename; //if the operator is stored on disk what is the filename
  std::vector<std::vector<int> > rowCompressedForm;  //the ith vector corresponds to all the non zero blocks in the ith row
  std::vector<std::vector<int> > colCompressedForm;  //the ith vector corresponds to all the non zero blocks in the ith column
  std::vector< std::pair<std::pair<int, int>, StackMatrix> > nonZeroBlocks; //all the nonzero blocks, the first pair is the row and col index
  std::map< std::pair<int, int>, int> mapToNonZeroBlocks; //the pair of indices will give the element of the nonzeroBlocks with the same pair

 public:
 StackSparseMatrix() : totalMemory(0), data(0), fermion(false), orbs(2), initialised(false), built(false), built_on_disk(false), Sign(1), conj('n'){};
 StackSparseMatrix(const StackSparseMatrix& a) : 
  orbs(a.get_orbs()), deltaQuantum(a.get_deltaQuantum()), fermion(a.get_fermion()), quantum_ladder(a.quantum_ladder), build_pattern(a.build_pattern),
    initialised(a.get_initialised()), allowedQuantaMatrix(a.get_allowedQuantaMatrix()), 
    Sign(a.get_sign()), totalMemory(a.totalMemory), data(a.data), conj('n'), built(a.built),
    rowCompressedForm(a.rowCompressedForm), built_on_disk(a.built_on_disk),
    colCompressedForm(a.colCompressedForm), nonZeroBlocks(a.nonZeroBlocks), mapToNonZeroBlocks(a.mapToNonZeroBlocks), filename(a.filename) 
    {};

 StackSparseMatrix(double* pData, long pTotalMemory) : totalMemory(pTotalMemory), data(pData), fermion(false), orbs(2), initialised(false), built(false), built_on_disk(false), Sign(1), conj('n'){};
  void SaveThreadSafe() const;
  void LoadThreadSafe(bool allocate);
  virtual long memoryUsed() const {return totalMemory;}
  void allocate (const StateInfo& s);
  void allocate (const StateInfo& sl, const StateInfo& sr);
  void allocate (const StateInfo& s, double* pData);
  double* allocate(const StateInfo& rowSI, const StateInfo& colSI, double* pData);
  void allocateShell(const StateInfo& rowSI, const StateInfo& colSI);
  void deallocate() ;
  double* allocateOperatorMatrix();
  string get_filename() {return filename;}
  string& set_filename() {return filename;}
  virtual void build(StackMatrix &m, int row, int col, const StackSpinBlock& block) {};
  virtual void build(const StackSpinBlock& block) {};
  double* get_data() {return data;}
  const double* get_data() const {return data;}
  long& set_totalMemory() {return totalMemory;}
  void set_data(double* pData) {data = pData;}
  virtual void deepCopy(const StackSparseMatrix& o) ;
  virtual void deepClearCopy(const StackSparseMatrix& o) ;
  virtual string opName() const {return "None";}
  //I cannot simply allow resize because resizing should be accompanied with appropriate data allocation first
  //void resize(int n, int c) { operatorMatrix.ReSize(n, c); allowedQuantaMatrix.ReSize(n, c); }
  std::vector<std::pair<std::pair<int, int>, StackMatrix> >& get_nonZeroBlocks() {return nonZeroBlocks;} 
  const std::vector<std::pair<std::pair<int, int>, StackMatrix> >& get_nonZeroBlocks() const {return nonZeroBlocks;} 
  const std::vector<int>& getActiveRows(int i) const {return colCompressedForm[i];}
  const std::vector<int>& getActiveCols(int i) const {return rowCompressedForm[i];}
  std::vector<int>& getActiveRows(int i)  {return colCompressedForm[i];}
  std::vector<int>& getActiveCols(int i) {return rowCompressedForm[i];}
  std::vector<std::vector<int> >& getrowCompressedForm() {return rowCompressedForm;}
  std::vector<std::vector<int> >& getcolCompressedForm() {return colCompressedForm;}
  const StackMatrix& operator_element(int i, int j) const { 
    const int index = mapToNonZeroBlocks.at( make_pair(i,j) );
    if (conj == 'n') return nonZeroBlocks[index].second;
    else return nonZeroBlocks[mapToNonZeroBlocks.at(std::pair<int,int>(j,i))].second;
    //if (conj == 'n') return operatorMatrix(i, j); 
    //else return operatorMatrix(j, i);
  }
  const StackMatrix& operator_element(int i, int j, char conj) const { 
    if (conj == 'n') return nonZeroBlocks[mapToNonZeroBlocks.at(std::pair<int,int>(i,j))].second;
    else return nonZeroBlocks[mapToNonZeroBlocks.at(std::pair<int,int>(j,i))].second;
    //if (conj == 'n') return operatorMatrix(i, j); 
    //else return operatorMatrix(j, i);
  }
  const StackMatrix& operator()(int i, int j) const { 
    if (conj == 'n') return nonZeroBlocks[mapToNonZeroBlocks.at(std::pair<int,int>(i,j))].second;
    else return nonZeroBlocks[mapToNonZeroBlocks.at(std::pair<int,int>(j,i))].second;
    //if (conj == 'n') return operatorMatrix(i, j); 
    //else return operatorMatrix(j, i);
  }
  const StackMatrix& operator()(int i, int j, char conj) const { 
    if (conj == 'n') return nonZeroBlocks[mapToNonZeroBlocks.at(std::pair<int,int>(i,j))].second;
    else return nonZeroBlocks[mapToNonZeroBlocks.at(std::pair<int,int>(j,i))].second;
    //if (conj == 'n') return operatorMatrix(i, j); 
    //else return operatorMatrix(j, i);
  }
  StackMatrix& operator_element(int i, int j) { 
    if (conj == 'n') return nonZeroBlocks[mapToNonZeroBlocks.at(std::pair<int,int>(i,j))].second;
    else return nonZeroBlocks[mapToNonZeroBlocks.at(std::pair<int,int>(j,i))].second;
    //if (conj == 'n') return operatorMatrix(i, j); 
    //else return operatorMatrix(j, i);
  }
  StackMatrix& operator_element(int i, int j, char conj) { 
    if (conj == 'n') return nonZeroBlocks[mapToNonZeroBlocks.at(std::pair<int,int>(i,j))].second;
    else return nonZeroBlocks[mapToNonZeroBlocks.at(std::pair<int,int>(j,i))].second;
    //if (conj == 'n') return operatorMatrix(i, j); 
    //else return operatorMatrix(j, i);
  }
  StackMatrix& operator()(int i, int j) { 
    if (conj == 'n') return nonZeroBlocks[mapToNonZeroBlocks.at(std::pair<int,int>(i,j))].second;
    else return nonZeroBlocks[mapToNonZeroBlocks.at(std::pair<int,int>(j,i))].second;
    //if (conj == 'n') return operatorMatrix(i, j); 
    //else return operatorMatrix(j, i);
  }
  StackMatrix& operator()(int i, int j, char conj) { 
    if (conj == 'n') return nonZeroBlocks[mapToNonZeroBlocks.at(std::pair<int,int>(i,j))].second;
    else return nonZeroBlocks[mapToNonZeroBlocks.at(std::pair<int,int>(j,i))].second;
    //if (conj == 'n') return operatorMatrix(i, j); 
    //else return operatorMatrix(j, i);
  }

  friend ostream& operator<<(ostream& os, const StackSparseMatrix& a);
  void operator=(const StackSparseMatrix& m);
  void Randomise();
  void SymmetricRandomise();
  void Normalise(int* success);
  void Clear();
  void CleanUp();
  StackSparseMatrix& operator+=(const StackSparseMatrix& other);

  void Save(std::ofstream &ofs) const;
  void Load(std::ifstream &ifs, bool allocateData);


  virtual ~StackSparseMatrix(){};
  int nrows() const { 
    if (conj == 'n') return allowedQuantaMatrix.nrows();
    else return allowedQuantaMatrix.ncols(); 
  }
  int ncols() const { 
    if (conj == 'n') return allowedQuantaMatrix.ncols(); 
    else return allowedQuantaMatrix.nrows();
  }
  int nrows(char conj) const { return allowedQuantaMatrix.Nrows(conj); }
  char conjugacy() const { return conj; }
  int get_sign() const {return Sign;}
  int &set_sign() {return Sign;}
  int ncols(char conj) const { return allowedQuantaMatrix.Ncols(conj); }
  bool get_initialised() const { return initialised; }
  bool &set_initialised() { return initialised; }
  bool get_fermion() const { return fermion; }
  bool &set_fermion() { return fermion; }
  ObjectMatrix<char>& set_allowedQuantaMatrix() {return allowedQuantaMatrix;}
  const ObjectMatrix<char>& get_allowedQuantaMatrix() const {return allowedQuantaMatrix;}
  const char& allowed(int i, int j) const { 
    if (conj == 'n') return allowedQuantaMatrix(i, j); 
    else return allowedQuantaMatrix(j, i);
  }
  char& allowed(int i, int j) { 
    if (conj == 'n') return allowedQuantaMatrix(i, j); 
    else return allowedQuantaMatrix(j, i);
  }
  const char& allowed(int i, int j, char conj) const { 
    if (conj == 'n') return allowedQuantaMatrix(i, j); 
    else return allowedQuantaMatrix(j, i);
  }
  char& allowed(int i, int j, char conj) { 
    if (conj == 'n') return allowedQuantaMatrix(i, j); 
    else return allowedQuantaMatrix(j, i);
  }
  std::vector<SpinQuantum> &set_deltaQuantum() { return deltaQuantum; }
  SpinQuantum &set_deltaQuantum(int i) { return deltaQuantum[i]; }
  void resize_deltaQuantum(int i) { deltaQuantum.resize(i); }
  void set_deltaQuantum(int i, const SpinQuantum s) { deltaQuantum.assign(i, s); }
  int get_deltaQuantum_size() const { return deltaQuantum.size(); }
  SpinSpace get_spin(int i=0) const  { return deltaQuantum[i].get_s();}
  SpinQuantum get_deltaQuantum(int i) const 
  { 
    if (conj == 'n')
      return deltaQuantum[i]; 
    else
      return -deltaQuantum[i];
  }
  IrrepSpace get_symm(int i=0) const  { 
    if (conj == 'n')
      {
	return deltaQuantum[i].get_symm();
      }
    else
      {
	return -deltaQuantum[i].get_symm();
      }
  }
  std::vector<SpinQuantum> get_deltaQuantum() const 
  {
    if (conj == 'n')
      return deltaQuantum; 
    else {
      std::vector<SpinQuantum> q;
      for (int i = 0; i < get_deltaQuantum_size(); ++i) {
	q.push_back(-get_deltaQuantum(i));
      }
      return q;
    }
  }

 
  std::map< std::string, std::vector<SpinQuantum> >  get_quantum_ladder() const { return quantum_ladder; }
  std::map< std::string, std::vector<SpinQuantum> >& set_quantum_ladder() { return quantum_ladder; }
  std::string  get_build_pattern() const { return build_pattern; }
  std::string& set_build_pattern() { return build_pattern; }
  void set_conjugacy(char c) { conj = c;}
  void reset_conjugacy() { conj = 'n';}
  int get_orbs(int i) const 
  { 
    if(i >= orbs.size())
      return -1;
    else
      return orbs[i]; 
  }
  const std::vector<int>& get_orbs() const { return orbs; }
  std::vector<int>& set_orbs() { return orbs; }
  const bool& get_built() const { return built; }
  bool& set_built() { return built; }  
  const bool& get_built_on_disk() const { return built_on_disk; }
  bool& set_built_on_disk() { return built_on_disk; }
  virtual double get_scaling(SpinQuantum leftq, SpinQuantum rightq) const ;

  void build_and_renormalise_transform(StackSpinBlock *big, const std::vector<Matrix>& rotate_matrix, 
				       const StateInfo *newStateInfo);
  void build_and_renormalise_transform(StackSpinBlock *big, const std::vector<Matrix>& leftrotate_matrix, const StateInfo *leftstateinfo, 
				       const std::vector<Matrix>& rightrotate_matrix,  const StateInfo *newStateInfo);
  void buildUsingCsf(const StackSpinBlock& b, vector< vector<Csf> >& ladders, std::vector< Csf >& s) ;
  virtual void buildUsingCre(const StackSpinBlock* b) {};
  virtual double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b=0) {return 0;}
  double calcCompfactor(TensorOp& Top1, TensorOp& op2, CompType comp, const TwoElectronArray& v_2, int integralIndex);
  double calcCompfactor(int i, int j, int k, int l, int spin, CompType comp, const TwoElectronArray& v_2, int integralIndex);
  double calcCompfactor(TensorOp& Top1, TensorOp& op2, CompType comp, const CCCCArray& vcccc);
  double calcCompfactor(TensorOp& Top1, TensorOp& op2, CompType comp, const CCCDArray& vcccd);
  double calcCompfactor(TensorOp& Top1, TensorOp& op2, CompType comp, int op2index, const TwoElectronArray& v_2, int integralIndex);
  double calcCompfactor(TensorOp& Top1, TensorOp& op2, CompType comp, int op2index, const CCCCArray& vcccc);
  double calcCompfactor(TensorOp& Top1, TensorOp& op2, CompType comp, int op2index, const CCCDArray& vcccd);
  bool nonZeroTensorComponent(Csf& c1, SpinQuantum& opsym, Csf& ladder, int& nonzeroindex, double& cleb);
  std::vector<double> calcMatrixElements(Csf& c1, TensorOp& Top, Csf& c2, std::vector<bool>& backup1, std::vector<bool>& backup2);
  double calcMatrixElements(Csf& c1, TensorOp& Top, Csf& c2, std::vector<bool>& backup1, std::vector<bool>& backup2, int index);
};


class StackTransposeview : public StackSparseMatrix
{
private:
  boost::shared_ptr<StackSparseMatrix> opdata; 
public:
 StackTransposeview(const boost::shared_ptr<StackSparseMatrix>& opptr) : opdata(opptr) {}
  StackTransposeview(StackSparseMatrix& op) { opdata = boost::shared_ptr<StackSparseMatrix>(&op, boostutils::null_deleter());}
  int get_deltaQuantum_size() const { return opdata->get_deltaQuantum_size(); }  
  SpinQuantum get_deltaQuantum(int i) const {return -opdata->get_deltaQuantum(i);}
  std::vector<SpinQuantum> get_deltaQuantum() const {
    std::vector<SpinQuantum> q;
    for (int i = 0; i < opdata->get_deltaQuantum_size(); ++i) {
      q.push_back(-opdata->get_deltaQuantum(i));
    }
    return q;
  }
  virtual long memoryUsed() const {return opdata->memoryUsed();}
  const std::vector<int>& getActiveRows(int i) const {return opdata->getActiveCols(i);}
  const std::vector<int>& getActiveCols(int i) const {return opdata->getActiveRows(i);}
  std::vector<int>& getActiveRows(int i)  {return opdata->getActiveCols(i);}
  std::vector<int>& getActiveCols(int i) {return opdata->getActiveRows(i);}
  bool get_fermion() const { return opdata->get_fermion(); }
  bool get_initialised() const { return opdata->get_initialised(); }
  int nrows() const { return opdata->ncols(); }
  int ncols() const { return opdata->nrows(); }
  virtual string opName() const {return opdata->opName();}

  virtual void build(StackMatrix &m, int row, int col, const StackSpinBlock& block) {opdata->build(m, col, row, block);}
  const char &allowed(int i, int j) const { return opdata->allowed(j, i); }
  char &allowed(int i, int j) { return opdata->allowed(j, i); }
  const StackMatrix& operator_element(int i, int j) const { return opdata->operator_element(j, i); }
  StackMatrix& operator_element(int i, int j) { return opdata->operator_element(j, i); }
  const StackMatrix& operator()(int i, int j) const { return opdata->operator()(j, i); }
  StackMatrix& operator()(int i, int j) { return opdata->operator()(j, i); }
  const StackMatrix& operator()(int i, int j, char conj) const { if (conj == 'n') return opdata->operator()(j, i); else return opdata->operator()(i,j);}
  StackMatrix& operator()(int i, int j, char conj) { if (conj == 'n') return opdata->operator()(j, i); else return opdata->operator()(i,j);}
 const char &allowed(int i, int j, char conj) const { if (conj == 'n') return opdata->allowed(j, i); else return opdata->allowed(i,j); }
  char &allowed(int i, int j, char conj) { if (conj =='n') return opdata->allowed(j, i); else return opdata->allowed(i,j);}
  const StackMatrix& operator_element(int i, int j, char conj) const { if (conj == 'n') return opdata->operator_element(j, i); else return opdata->operator_element(i,j);}
  StackMatrix& operator_element(int i, int j, char conj) { if (conj == 'n') return opdata->operator_element(j, i); else return opdata->operator_element(i,j);}
  //std::vector<std::pair<std::pair<int, int>, StackMatrix> >& get_nonZeroBlocks() {return nonZeroBlocks;} 
  //const std::vector<std::pair<std::pair<int, int>, StackMatrix> >& get_nonZeroBlocks() const {return nonZeroBlocks;} 

  SpinSpace get_spin(int i=0) const  { return -opdata->get_deltaQuantum(i).get_s();}
  IrrepSpace get_symm(int i=0) const  { return -opdata->get_deltaQuantum(i).get_symm();}
  int get_orbs(int i) const {return opdata->get_orbs(i);}
  const std::vector<int>& get_orbs() const { return opdata->get_orbs(); }
   char conjugacy() const { if (opdata->conjugacy() == 'n') return 't'; else return 'n';}
  //double get_scaling(SpinQuantum leftq, SpinQuantum rightq) const ;
  boost::shared_ptr<StackSparseMatrix> getworkingrepresentation(const StackSpinBlock* block) {return opdata;}
  void build(const StackSpinBlock& b){};
  double redMatrixElement(Csf c1, vector<Csf>& ladder, const StackSpinBlock* b){return 0.0;}
}; 

const StackTransposeview Transpose(StackSparseMatrix& op);


void assignloopblock(StackSpinBlock*& loopblock, StackSpinBlock*& otherblock, StackSpinBlock* leftBlock,
		     StackSpinBlock* rightBlock);
long getRequiredMemory(const StateInfo& sr, const StateInfo& sc, const std::vector<SpinQuantum>& q); 
long getRequiredMemoryForWavefunction(const StateInfo& sr, const StateInfo& sc, const std::vector<SpinQuantum>& q); 
long getRequiredMemory(const StackSpinBlock& b, const std::vector<SpinQuantum>& q); 
long getRequiredMemory(const StateInfo& s, const std::vector<SpinQuantum>& q);
void Normalise(StackSparseMatrix& a, int* success = 0);
void ScaleAdd(double d, const StackSparseMatrix& a, StackSparseMatrix& b);
double DotProduct(const StackSparseMatrix& lhs, const StackSparseMatrix& rhs);
double trace(const StackSparseMatrix& lhs);
void Scale(double d, StackSparseMatrix& a);
//void copy(const ObjectMatrix<StackMatrix>& a, ObjectMatrix<StackMatrix>& b);
//void copy(const ObjectMatrix<Matrix>& a, ObjectMatrix<StackMatrix>& b);
//void copy(const ObjectMatrix<StackMatrix>& a, ObjectMatrix<Matrix>& b);
void copy(const StackMatrix& a, StackMatrix& b);
void copy(const StackMatrix& a, Matrix& b);
void copy(const Matrix& a, StackMatrix& b);
void copy(const Matrix& a, Matrix& b);
double getStandAlonescaling(SpinQuantum opQ, SpinQuantum leftq, SpinQuantum rightq);

} ;


#endif
