/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "StackBaseOperator.h"
#include "MatrixBLAS.h"
#include "Stackspinblock.h"
#include "distribute.h"
#include "tensor_operator.h"
#include "blas_calls.h"
#include "pario.h"
#include "StackMatrix.h"
#include "csf.h"
#include "StateInfo.h"

namespace SpinAdapted{
  
void StackSparseMatrix::deepCopy(const StackSparseMatrix& o)
{
  *this = o;
  data = Stackmem[omprank].allocate(totalMemory);
  DCOPY(totalMemory, o.data, 1, data, 1);
  allocateOperatorMatrix();
}

void StackSparseMatrix::deepClearCopy(const StackSparseMatrix& o)
{
  *this = o;
  data = Stackmem[omprank].allocate(totalMemory);
  memset(data, 0, totalMemory*sizeof(double));
  allocateOperatorMatrix();
}

const StackTransposeview Transpose(StackSparseMatrix& op) { return StackTransposeview(op); };


void StackSparseMatrix::build_and_renormalise_transform(StackSpinBlock *big, const std::vector<Matrix>& rotate_matrix, 
							const StateInfo *newStateInfo) 
{
  //backup old data
  double* oldData = data;
  long oldtotalMemory = totalMemory;
  std::vector<std::vector<int> > oldrowCompressedForm = rowCompressedForm;
  std::vector<std::vector<int> > oldcolCompressedForm = colCompressedForm;
  std::vector< std::pair<std::pair<int, int>, StackMatrix> > oldnonZeroBlocks = nonZeroBlocks; 
  std::map< std::pair<int, int>, int> oldmapToNonZeroBlocks = mapToNonZeroBlocks; 
  //ObjectMatrix<StackMatrix> oldoperatorMatrix=operatorMatrix; 
  ObjectMatrix<char> oldallowedQuantaMatrix = allowedQuantaMatrix;

  //allocate new data and build the operator
  totalMemory = 0; data=0;
  colCompressedForm.clear();
  rowCompressedForm.clear();
  allocate(big->get_braStateInfo(), big->get_ketStateInfo());
  build(*big);

  //put the new operatorMatrix 
  StackSparseMatrix tmp; //tmp.operatorMatrix = operatorMatrix;
  tmp.rowCompressedForm = rowCompressedForm;
  tmp.colCompressedForm = colCompressedForm;
  tmp.nonZeroBlocks = nonZeroBlocks;
  tmp.mapToNonZeroBlocks = mapToNonZeroBlocks;
  tmp.data = data; tmp.totalMemory = totalMemory;
  tmp.allowedQuantaMatrix = allowedQuantaMatrix;
  tmp.initialised = true;


  //restore the data
  data = oldData;
  totalMemory = oldtotalMemory;
  rowCompressedForm = oldrowCompressedForm;
  colCompressedForm = oldcolCompressedForm;
  nonZeroBlocks = oldnonZeroBlocks;
  mapToNonZeroBlocks = oldmapToNonZeroBlocks;
  //operatorMatrix = oldoperatorMatrix;
  allowedQuantaMatrix = oldallowedQuantaMatrix;
  memset(data, 0, totalMemory*sizeof(double));

  const std::vector<int>& newQuantaMap = newStateInfo->newQuantaMap;


  int quanta_thrds = dmrginp.quanta_thrds();
#pragma omp parallel for schedule(dynamic) num_threads(quanta_thrds)
  for (int newQ = 0; newQ < newQuantaMap.size(); newQ++)
    for (int newQPrime = 0; newQPrime < newQuantaMap.size(); newQPrime++) {
      if (this->allowed(newQ, newQPrime)) {
	int Q = newQuantaMap[newQ], QPrime = newQuantaMap[newQPrime];
	MatrixRotate(rotate_matrix[Q], tmp(Q, QPrime), rotate_matrix[QPrime], this->operator()(newQ, newQPrime));
      }
    }
  tmp.deallocate();
}
  
void StackSparseMatrix::build_and_renormalise_transform(StackSpinBlock *big, const std::vector<Matrix>& leftrotate_matrix, const StateInfo *newleftStateInfo, 
							const std::vector<Matrix>& rightrotate_matrix,  const StateInfo *newrightStateInfo) 
{
  
  //backup old data
  double* oldData = data;
  long oldtotalMemory = totalMemory;
  std::vector<std::vector<int> > oldrowCompressedForm = rowCompressedForm;
  std::vector<std::vector<int> > oldcolCompressedForm = colCompressedForm;
  std::vector< std::pair<std::pair<int, int>, StackMatrix> > oldnonZeroBlocks = nonZeroBlocks; 
  std::map< std::pair<int, int>, int> oldmapToNonZeroBlocks = mapToNonZeroBlocks; 
  //ObjectMatrix<StackMatrix> oldoperatorMatrix=operatorMatrix; 
  ObjectMatrix<char> oldallowedQuantaMatrix = allowedQuantaMatrix;

  //allocate new data and build the operator
  totalMemory = 0; data=0;
  colCompressedForm.clear();
  rowCompressedForm.clear();
  allocate(big->get_braStateInfo(), big->get_ketStateInfo());
  initialised = true;
  build(*big);
  
  //put the new operatorMatrix 
  StackSparseMatrix tmp; //tmp.operatorMatrix = operatorMatrix;
  tmp.rowCompressedForm = rowCompressedForm;
  tmp.colCompressedForm = colCompressedForm;
  tmp.nonZeroBlocks = nonZeroBlocks;
  tmp.mapToNonZeroBlocks = mapToNonZeroBlocks;
  tmp.data = data; tmp.totalMemory = totalMemory;
  tmp.allowedQuantaMatrix = allowedQuantaMatrix;
  tmp.initialised = true;

  //restore the data
  data = oldData;
  totalMemory = oldtotalMemory;
  rowCompressedForm = oldrowCompressedForm;
  colCompressedForm = oldcolCompressedForm;
  nonZeroBlocks = oldnonZeroBlocks;
  mapToNonZeroBlocks = oldmapToNonZeroBlocks;
  //operatorMatrix = oldoperatorMatrix;
  allowedQuantaMatrix = oldallowedQuantaMatrix;
  memset(data, 0, totalMemory*sizeof(double));
  
  const std::vector<int>& lnewQuantaMap = newleftStateInfo->newQuantaMap;
  const std::vector<int>& rnewQuantaMap = newrightStateInfo->newQuantaMap;
  

  for (int newQ = 0; newQ < lnewQuantaMap.size(); newQ++)
    for (int newQPrime = 0; newQPrime < rnewQuantaMap.size(); newQPrime++) {
      if (this->allowed(newQ, newQPrime)) {
	int Q = lnewQuantaMap[newQ], QPrime = rnewQuantaMap[newQPrime];
	MatrixRotate(leftrotate_matrix[Q], tmp(Q, QPrime), rightrotate_matrix[QPrime], this->operator()(newQ, newQPrime));
      }
    }

  tmp.deallocate();
  tmp.CleanUp();
}



void StackSparseMatrix::buildUsingCsf(const StackSpinBlock& b, vector< vector<Csf> >& ladders, std::vector< Csf >& s) 
{
  StateInfo stateinfo = b.get_stateInfo();
  built = true;
  
  for (int index=0; index < nonZeroBlocks.size(); index++) {
    int i = nonZeroBlocks[index].first.first, j = nonZeroBlocks[index].first.second;
    for (int jq =stateinfo.unBlockedIndex[j]; jq < stateinfo.unBlockedIndex[j]+stateinfo.quantaStates[j]; jq++) 
      for (int iq =stateinfo.unBlockedIndex[i]; iq < stateinfo.unBlockedIndex[i]+stateinfo.quantaStates[i]; iq++) {
	      nonZeroBlocks[index].second.operator()(iq-stateinfo.unBlockedIndex[i]+1, jq-stateinfo.unBlockedIndex[j]+1) = redMatrixElement(s[iq], ladders[jq], &b);
      }
  }

}


void assignloopblock(StackSpinBlock*& loopblock, StackSpinBlock*& otherblock, StackSpinBlock* leftBlock,
			    StackSpinBlock* rightBlock)
{
  //if one of the blocks is a dot block then it should be the loopblock
  if (dmrginp.spinAdapted() ) {
    if (rightBlock->get_sites().size() == 1) {
      loopblock = rightBlock;
      otherblock = leftBlock;
      return;
    }
  }
  else {
    if (rightBlock->get_sites().size() == 2) {
      loopblock = rightBlock;
      otherblock = leftBlock;
      return;
    }
  }

  //if none are dot blocks then use the information provided
  if (leftBlock->is_loopblock()) {
    loopblock = leftBlock;
    otherblock = rightBlock;
  }
  else {
    loopblock = rightBlock; 
    otherblock = leftBlock;
  }
}


void StackSparseMatrix::operator=(const StackSparseMatrix& a) 
{
  filename = a.filename;
  deltaQuantum = a.get_deltaQuantum();
  quantum_ladder = a.quantum_ladder;
  build_pattern = a.build_pattern;
  fermion = a.get_fermion();
  initialised = a.get_initialised();
  built = a.built;
  built_on_disk = a.built_on_disk;
  allowedQuantaMatrix = a.get_allowedQuantaMatrix();
  //operatorMatrix = a.operatorMatrix;
  Sign = a.get_sign();
  orbs = a.get_orbs(); 
  rowCompressedForm = a.rowCompressedForm;
  colCompressedForm = a.colCompressedForm;
  nonZeroBlocks = a.nonZeroBlocks;
  mapToNonZeroBlocks = a.mapToNonZeroBlocks;
  conj = a.conj;
  totalMemory = a.totalMemory;
  data=a.data;
}

double StackSparseMatrix::get_scaling(SpinQuantum leftq, SpinQuantum rightq) const 
{
  if(!dmrginp.spinAdapted()) return 1.0;
  if (conjugacy() == 'n') {return 1.0;}

  int lspin = leftq.get_s().getirrep(), lirrep = leftq.get_symm().getirrep();
  int rspin = rightq.get_s().getirrep(), rirrep = rightq.get_symm().getirrep();
  int cspin = get_deltaQuantum(0).get_s().getirrep(), cirrep = get_deltaQuantum(0).get_symm().getirrep();

  int cirrepTranspose = (-get_deltaQuantum(0)).get_symm().getirrep();
  for (int lsz = -lspin; lsz<lspin+1; lsz+=2)
  for (int rsz = -rspin; rsz<rspin+1; rsz+=2)
  for (int ll = 0; ll<Symmetry::sizeofIrrep(lirrep); ll++)
  for (int rl = 0; rl<Symmetry::sizeofIrrep(rirrep); rl++)
  {
    double cleb = clebsch(lspin, lsz, cspin, -cspin, rspin, rsz);
    double clebspatial = Symmetry::spatial_cg(lirrep, cirrep, rirrep, ll, 0, rl);
    if (fabs(cleb) <= NUMERICAL_ZERO || fabs(clebspatial) <= NUMERICAL_ZERO)
      continue;
    else {
      ///CHANGE THE SPATIAL_CG cirrep,1 to cirrep,0 depending on how the transpose works out!!!
      double spinscale = pow(-1.0,cspin) * cleb/clebsch(rspin, rsz, cspin, cspin, lspin, lsz);
      double spatscale =  clebspatial/Symmetry::spatial_cg(rirrep, cirrepTranspose, lirrep, rl, Symmetry::sizeofIrrep(cirrep)-1, ll);  

      return spinscale*spatscale;
    }
  }
  pout << "Major trouble, inappropriate arguments to get_scaling!!!"<<endl;
  pout << leftq<<"  ";
  for (int i = 0; i < get_deltaQuantum_size(); ++i) {
    pout <<get_deltaQuantum(i)<<"  ";
  }
  pout <<rightq<<endl;
  exit(0);
  return 1.0;
}

double getStandAlonescaling(SpinQuantum opq, SpinQuantum leftq, SpinQuantum rightq) 
{
  if(!dmrginp.spinAdapted()) return 1.0;

  int lspin = leftq.get_s().getirrep(), lirrep = leftq.get_symm().getirrep();
  int rspin = rightq.get_s().getirrep(), rirrep = rightq.get_symm().getirrep();
  int cspin = opq.get_s().getirrep(), cirrep = opq.get_symm().getirrep();

  int cirrepTranspose = (-opq).get_symm().getirrep();
  for (int lsz = -lspin; lsz<lspin+1; lsz+=2)
  for (int rsz = -rspin; rsz<rspin+1; rsz+=2)
  for (int ll = 0; ll<Symmetry::sizeofIrrep(lirrep); ll++)
  for (int rl = 0; rl<Symmetry::sizeofIrrep(rirrep); rl++)
  {
    double cleb = clebsch(lspin, lsz, cspin, -cspin, rspin, rsz);
    double clebspatial = Symmetry::spatial_cg(lirrep, cirrep, rirrep, ll, 0, rl);
    if (fabs(cleb) <= NUMERICAL_ZERO || fabs(clebspatial) <= NUMERICAL_ZERO)
      continue;
    else {
      ///CHANGE THE SPATIAL_CG cirrep,1 to cirrep,0 depending on how the transpose works out!!!
      double spinscale = pow(-1.0,cspin) * cleb/clebsch(rspin, rsz, cspin, cspin, lspin, lsz);
      double spatscale =  clebspatial/Symmetry::spatial_cg(rirrep, cirrepTranspose, lirrep, rl, Symmetry::sizeofIrrep(cirrep)-1, ll);  

      return spinscale*spatscale;
    }
  }
  pout << "Major trouble, inappropriate arguments to get_scaling!!!"<<endl;
  pout << leftq<<"  "<<opq<<"  ";
  pout <<rightq<<endl;
  exit(0);
  return 1.0;
}


long getRequiredMemory(const StackSpinBlock& b, const std::vector<SpinQuantum>& q) {
  return getRequiredMemory(b.get_braStateInfo(), b.get_ketStateInfo(), q);
}
long getRequiredMemory(const StateInfo& s, const std::vector<SpinQuantum>& q) {return getRequiredMemory(s,s,q);} 

void StackSparseMatrix::deallocate() 
{ 
  //allowedQuantaMatrix.Clear();
  //operatorMatrix.Clear();
  Stackmem[omprank].deallocate(data, totalMemory);
  totalMemory = 0;
  data = 0;
}


void StackSparseMatrix::allocate(const StateInfo& sr, const StateInfo& sl) {
  if (totalMemory != 0) return; //allready built
  long requiredData = getRequiredMemory(sr, sl, get_deltaQuantum());
  double* data = Stackmem[omprank].allocate(requiredData);
  memset(data, 0, requiredData * sizeof(double));      
  allocate(sr, sl, data);
} 

void StackSparseMatrix::allocate(const StateInfo& s) {
  if (totalMemory != 0) return; //allready built
  long requiredData = getRequiredMemory(s, s, get_deltaQuantum());
  double* data = Stackmem[omprank].allocate(requiredData);
  memset(data, 0, requiredData * sizeof(double));      
  allocate(s, s, data);
} 

void StackSparseMatrix::allocate(const StateInfo& sr, double* pData) {  
  if (totalMemory != 0) return; //allready built
  allocate(sr, sr, pData);
}

double* StackSparseMatrix::allocate(const StateInfo& rowSI, const StateInfo& colSI, double* pData)
{
  if (totalMemory != 0) {
    perr << "Already have memory and allocating more memory"<<endl;
    perr << "Possible memory bug in StackSparseMatrix "<<endl;
    abort();
  }
  //the shell of the operator is already present
  //so just need to allocate the operatormatrix
  if (rowCompressedForm.size() != 0) {
    data = pData;
    return allocateOperatorMatrix();
  }
  rowCompressedForm.clear();rowCompressedForm.resize(rowSI.quanta.size(), vector<int>());
  colCompressedForm.clear();colCompressedForm.resize(colSI.quanta.size(), vector<int>());
  nonZeroBlocks.resize(0);
  mapToNonZeroBlocks.clear();

  data = pData;
  //operatorMatrix.resize(rowSI.quanta.size(), colSI.quanta.size());
  allowedQuantaMatrix.resize(rowSI.quanta.size(), colSI.quanta.size());

  long index = 0;
  for (int lQ = 0; lQ < rowSI.quanta.size (); ++lQ)
    for (int rQ = 0; rQ < colSI.quanta.size (); ++rQ)
    {
      bool allowedcoupling = false;
      for (int k = 0; k < deltaQuantum.size(); ++k) {
        if (rowSI.quanta[lQ].allow(deltaQuantum[k], colSI.quanta[rQ])) {
          allowedcoupling = true;
	        rowCompressedForm[lQ].push_back(rQ);
	        colCompressedForm[rQ].push_back(lQ);
          break;
        }
      }
      allowedQuantaMatrix (lQ,rQ) = allowedcoupling;
      if (allowedQuantaMatrix(lQ, rQ))
      {
	      StackMatrix m(&data[index], rowSI.quantaStates [lQ], colSI.quantaStates [rQ]);
	      index += rowSI.quantaStates [lQ]* colSI.quantaStates [rQ] + CACHEBUFFER;
	      nonZeroBlocks.push_back(std::pair< std::pair<int, int> , StackMatrix>( std::pair<int,int>(lQ,rQ), m));
	      mapToNonZeroBlocks.insert(std::pair< std::pair<int, int>, int>(std::pair<int, int>(lQ, rQ), nonZeroBlocks.size()-1));
      }
    }
  totalMemory = index;
  return &pData[index];
}

void StackSparseMatrix::allocateShell(const StateInfo& rowSI, const StateInfo& colSI)
{
  if (totalMemory != 0) {
    perr << "Already have memory and allocating more memory"<<endl;
    perr << "Possible memory bug in StackSparseMatrix "<<endl;
    abort();
  }
  //the shell of the operator is already present
  //so just need to allocate the operatormatrix
  if (rowCompressedForm.size() != 0) {
    return;
  }

  rowCompressedForm.clear();rowCompressedForm.resize(rowSI.quanta.size(), vector<int>());
  colCompressedForm.clear();colCompressedForm.resize(colSI.quanta.size(), vector<int>());
  nonZeroBlocks.resize(0);
  mapToNonZeroBlocks.clear();

  //operatorMatrix.resize(rowSI.quanta.size(), colSI.quanta.size());
  allowedQuantaMatrix.resize(rowSI.quanta.size(), colSI.quanta.size());

  long index = 0;
  for (int lQ = 0; lQ < rowSI.quanta.size (); ++lQ)
    for (int rQ = 0; rQ < colSI.quanta.size (); ++rQ)
    {
      bool allowedcoupling = false;
      for (int k = 0; k < deltaQuantum.size(); ++k) {
        if (rowSI.quanta[lQ].allow(deltaQuantum[k], colSI.quanta[rQ])) {
          allowedcoupling = true;
	  rowCompressedForm[lQ].push_back(rQ);
	  colCompressedForm[rQ].push_back(lQ);
          break;
        }
      }
      allowedQuantaMatrix (lQ,rQ) = allowedcoupling;
    }
  return;
}


  


double* StackSparseMatrix::allocateOperatorMatrix()
{
  long index = 0;
  for (int i=0; i<nonZeroBlocks.size(); i++) {
    int lQ = nonZeroBlocks[i].first.first, rQ = nonZeroBlocks[i].first.second;
    StackMatrix m(&data[index], operator()(lQ, rQ).Nrows(), operator()(lQ, rQ).Ncols());
    //operatorMatrix(lQ,rQ).allocate(&data[index], operatorMatrix(lQ, rQ).Nrows(), operatorMatrix(lQ, rQ).Ncols());
    index += operator()(lQ, rQ).Nrows()* operator()(lQ, rQ).Ncols()+ CACHEBUFFER;
    nonZeroBlocks[i].second= m;
    mapToNonZeroBlocks.insert(std::pair< std::pair<int, int>, int>(std::pair<int, int>(lQ, rQ), i));
  }
  totalMemory = index;
  return &data[index];
}



template<class T>  void write(std::ofstream &ofs, const T& t)
{
  ofs.write((char*)(&t), sizeof(t));
}
template<class T>  void read(std::ifstream &ifs, const T& t)
{
  ifs.read((char*)(&t), sizeof(t));
}
void SaveVecThreadSafe(std::ofstream &ofs, const std::vector<int>& v)
{
  write(ofs, v.size());
  for (int i=0; i<v.size(); i++)
    write(ofs, v[i]);

}

void LoadVecThreadSafe(std::ifstream &ifs, std::vector<int>& v)
{
  size_t n;
  read(ifs, n);
  v.resize(n);
  for (int i=0; i<v.size(); i++)
    read(ifs, v[i]);
}

void StackSparseMatrix::SaveThreadSafe() const
{
  std::ofstream ofs(filename.c_str(), ios::binary);
  //#pragma omp critical
  //{
  write(ofs, deltaQuantum.size());

  for (int i=0; i<deltaQuantum.size(); i++)
    deltaQuantum[i].SaveThreadSafe(ofs);

  size_t size = build_pattern.size();
  write(ofs, size);
  ofs.write(build_pattern.c_str(), build_pattern.size());
  write(ofs, fermion);
  write(ofs, initialised);
  write(ofs, built);
  write(ofs, built_on_disk);

  size_t nrows = allowedQuantaMatrix.nrows();
  size_t ncols = allowedQuantaMatrix.ncols();
  write(ofs, nrows);
  write(ofs, ncols);
  for (int i=0 ; i<allowedQuantaMatrix.nrows(); i++)
    for (int j=0; j<allowedQuantaMatrix.ncols(); j++)
      write(ofs, allowedQuantaMatrix(i,j));

  /*
  for (int i=0 ; i<allowedQuantaMatrix.nrows(); i++)
    for (int j=0; j<allowedQuantaMatrix.ncols(); j++)
      operatorMatrix(i,j).SaveThreadSafe(ofs);
  */
  write(ofs, Sign) ;
  size = orbs.size();
  write(ofs, size);
  for (int i=0; i< orbs.size(); i++)
    write(ofs, orbs[i]);

  write(ofs, quantum_ladder.size());
  for (std::map<std::string, std::vector<SpinQuantum> >::const_iterator it = quantum_ladder.begin(); it != quantum_ladder.end(); it++) {
    size = it->first.size();
    write(ofs, size);
    ofs.write(it->first.c_str(), it->first.size());

    size = it->second.size();
    write(ofs, size);
    for (int i=0; i<it->second.size(); i++)
      it->second[i].SaveThreadSafe(ofs);
  }

  size = rowCompressedForm.size();
  write(ofs, size);
  for (int i=0; i<rowCompressedForm.size(); i++)
    SaveVecThreadSafe(ofs, rowCompressedForm[i]);

  size = colCompressedForm.size();
  write(ofs, size);
  for (int i=0; i<colCompressedForm.size(); i++)
    SaveVecThreadSafe(ofs, colCompressedForm[i]);

  size = nonZeroBlocks.size();
  write(ofs, size);
  for (int i=0; i<nonZeroBlocks.size(); i++) {
    write(ofs, nonZeroBlocks[i].first.first);
    write(ofs, nonZeroBlocks[i].first.second);
    nonZeroBlocks[i].second.SaveThreadSafe(ofs);
  }

  size = mapToNonZeroBlocks.size();
  write(ofs, size);
  for (std::map<std::pair<int, int>, int >::const_iterator it = mapToNonZeroBlocks.begin(); it != mapToNonZeroBlocks.end(); it++) {
    write(ofs, it->first.first);
    write(ofs, it->first.second);
    write(ofs, it->second);
  }

  write(ofs, totalMemory);
  ofs.write((char*)(data), sizeof(*data)*totalMemory);
  //}
  ofs.close();
}


void StackSparseMatrix::LoadThreadSafe(bool allocate)
{
  std::ifstream ifs(filename.c_str(), ios::binary);
  //#pragma omp critical
  //{
  size_t size=0;
  read(ifs, size);
  deltaQuantum.resize(size);
  for (int i=0; i<deltaQuantum.size(); i++)
    deltaQuantum[i].LoadThreadSafe(ifs);

  read(ifs, size);
  char* temp = new char[size+1];
  ifs.read(temp, size);
  temp[size] = '\0';
  build_pattern = temp;
  delete [] temp;

  read(ifs, fermion);
  read(ifs, initialised);
  read(ifs, built);
  read(ifs, built_on_disk);

  size_t nrows, ncols;
  read(ifs, nrows);
  read(ifs, ncols);
  allowedQuantaMatrix.resize(nrows, ncols);
  //operatorMatrix.resize(nrows, ncols);
  for (int i=0 ; i<allowedQuantaMatrix.nrows(); i++)
    for (int j=0; j<allowedQuantaMatrix.ncols(); j++)
      read(ifs, allowedQuantaMatrix(i,j));
  /*
  for (int i=0 ; i<allowedQuantaMatrix.nrows(); i++)
    for (int j=0; j<allowedQuantaMatrix.ncols(); j++)
      operatorMatrix(i,j).LoadThreadSafe(ifs);
  */
  

  read(ifs, Sign) ;
  read(ifs, size);
  orbs.resize(size);
  for (int i=0; i< orbs.size(); i++)
    read(ifs, orbs[i]);

  read(ifs, size);
  for (int i=0; i<size; i++) {
    size_t size2;
    read(ifs, size2);
    char* temp = new char[size2+1];
    ifs.read(temp, size2);
    temp[size2] = '\0';
    string s = temp;
    delete [] temp;

    std::vector<SpinQuantum> qvec;

    read(ifs, size2);
    qvec.resize(size2);
    for (int j=0; j<size2; j++) 
      qvec[j].LoadThreadSafe(ifs);
    quantum_ladder[s] = qvec;
  }


  read(ifs, size);
  rowCompressedForm.resize(size);
  for (int i=0; i<rowCompressedForm.size(); i++)
    LoadVecThreadSafe(ifs, rowCompressedForm[i]);

  read(ifs, size);
  colCompressedForm.resize(size);
  for (int i=0; i<colCompressedForm.size(); i++)
    LoadVecThreadSafe(ifs, colCompressedForm[i]);

  read(ifs, size);
  nonZeroBlocks.resize(size);
  for (int i=0; i<nonZeroBlocks.size(); i++) {
    int first, second;
    read(ifs, first);
    read(ifs, second);
    StackMatrix m;
    m.LoadThreadSafe(ifs);
    nonZeroBlocks[i] = make_pair( make_pair(first, second), m);
  }

  read(ifs, size);
  for (int i=0; i<size; i++) {
    int first, second, third;
    read(ifs, first);
    read(ifs, second);
    read(ifs, third);
    mapToNonZeroBlocks[make_pair(first, second)] = third;
  }

  read(ifs, totalMemory);

  if (allocate) data = Stackmem[omprank].allocate(totalMemory);
  ifs.read((char*)(data), sizeof(*data)*totalMemory);

  ifs.close();
  //}
}

void StackSparseMatrix::Save(std::ofstream &ofs) const
{
    boost::archive::binary_oarchive save_op(ofs);
    save_op << *this;
    save_op << totalMemory;
    save_op << boost::serialization::make_array<double>(data, totalMemory);
}
 
void StackSparseMatrix::Load(std::ifstream &ifs, bool allocateData)
{
    boost::archive::binary_iarchive load_op(ifs);
    load_op >> *this;
    
    load_op >> totalMemory;
    if (allocateData) data = Stackmem[omprank].allocate(totalMemory);
    load_op >> boost::serialization::make_array<double>(data, totalMemory);
}

long getRequiredMemory(const StateInfo& sr, const StateInfo& sc, const std::vector<SpinQuantum>& q) {
  dmrginp.getreqMem->start();
  long memory = 0;

  //#pragma omp parallel for schedule(dynamic) reduction(+ : memory)
  for (int i=0; i < sr.quanta.size(); i++)
    for (int j=0; j<sc.quanta.size(); j++) {
      bool allowedcoupling = false;
      for (int k=0; k<q.size(); k++) {
	if (sr.quanta[i].allow(q[k], sc.quanta[j])) {
	  allowedcoupling = true;
	  break;
	}
      }
      if (allowedcoupling)
	memory += sr.quantaStates[i]*sc.quantaStates[j]+ CACHEBUFFER;
    }
  dmrginp.getreqMem->stop();

  return memory;
}


void StackSparseMatrix::CleanUp ()
{
  allowedQuantaMatrix.Clear();
  //operatorMatrix.Clear();
  nonZeroBlocks.resize(0);
  rowCompressedForm.resize(0);
  colCompressedForm.resize(0);
  mapToNonZeroBlocks.clear();
}

  //const StackTransposeview Transpose(StackSparseMatrix& op) { return StackTransposeview(op); };

ostream& operator<< (ostream& os, const StackSparseMatrix& a)
{
  os << " nonZeroBlocks "<<a.nonZeroBlocks.size()<<endl;
  for (int i=0; i<a.nonZeroBlocks.size(); i++)
    os << a.nonZeroBlocks[i].first.first<<"  "<<a.nonZeroBlocks[i].first.second<<endl;
  if (!a.initialised){
    os <<" not initialised"<<endl;
    //return os;
  };
  os<<"indices : ";
  for(int i=0; i<a.orbs.size(); i++)
    os<<a.orbs[i]<<"  ";
  os <<endl;
  for (int i = 0; i < a.get_deltaQuantum_size(); ++i) {
    os<<a.get_deltaQuantum(i)<<endl;
  }
  os<<a.totalMemory<<endl;
  for (int i = 0; i < a.nrows (); ++i)
	for (int j = 0; j < a.ncols (); ++j)
	{
	  if (a.allowed(i, j)) 
	    os << i << " " << j << endl << a.operator_element(i, j) << endl;
	}
  return os;
}


void StackSparseMatrix::Randomise ()
{
  for (int lQ = 0; lQ < nrows(); ++lQ)
    for (int rQ = 0; rQ < ncols(); ++rQ)
      if (allowed(lQ, rQ))
	SpinAdapted::Randomise(operator_element(lQ, rQ));
}

void StackSparseMatrix::SymmetricRandomise ()
{
  for (int lQ = 0; lQ < nrows(); ++lQ)
    for (int rQ = 0; rQ < ncols(); ++rQ)
      if (allowed(lQ, rQ)) 
	SpinAdapted::SymmetricRandomise(operator_element(lQ, rQ));
}

double trace(const StackSparseMatrix& lhs)
{
  assert(lhs.nrows() == lhs.ncols());
  double trace = 0.0;

  for(int lQ=0;lQ<lhs.nrows();++lQ)
    if(lhs.allowed(lQ,lQ))
      for(int i=0;i< lhs(lQ,lQ).Nrows();++i)
	trace += (lhs)(lQ,lQ)(i+1,i+1);

  return trace;
}

double DotProduct(const StackSparseMatrix& lhs, const StackSparseMatrix& rhs)
{
  double result = 0.;
  for (int lQ = 0; lQ < lhs.nrows(); ++lQ)
    for (int rQ = 0; rQ < lhs.ncols (); ++rQ)
      if (lhs.allowed(lQ, rQ) && rhs.allowed(lQ, rQ))
	    result += MatrixDotProduct(lhs.operator_element(lQ, rQ), rhs.operator_element(lQ, rQ));

  return result;
}

void Scale(double d, StackSparseMatrix& a)
{
  for (int lQ = 0; lQ < a.nrows(); ++lQ)
    for (int rQ = 0; rQ < a.ncols(); ++rQ)
      if (a.allowed(lQ, rQ))
        MatrixScale(d, a.operator_element(lQ, rQ));
}

void ScaleAdd(double d, const StackSparseMatrix& a, StackSparseMatrix& b)
{
  for (int lQ = 0; lQ < a.nrows(); ++lQ)
    for (int rQ = 0; rQ < a.ncols(); ++rQ)
      if (a.allowed(lQ, rQ))
      {
	    if (!b.allowed(lQ, rQ))
	      pout <<"Not a valid addition"<<endl;
        assert(b.allowed(lQ, rQ));
        MatrixScaleAdd(d, a.operator_element(lQ, rQ), b.operator_element(lQ, rQ));
      }
}

void Normalise(StackSparseMatrix& a, int* success)
{
  a.Normalise(success);
}

void StackSparseMatrix::Normalise (int* success)
{
  double normalisation = DotProduct(*this, *this);

  //if the norm is really small then dont normalize??
  if(normalisation > NUMERICAL_ZERO)
    Scale(1./sqrt(normalisation), *this);
  else {
    pout << "\t\t\t Warning :: Didn't Normalise, because norm is " << normalisation<<endl;
    *success = 1; //not successful in normlaising                                                                                  
  }
}

void StackSparseMatrix::Clear ()
{
  built = false;
  memset(data, 0, totalMemory*sizeof(double));
}
void copy(const ObjectMatrix<Matrix>& a, ObjectMatrix<Matrix>& b)
{
  b.resize(a.Nrows(), a.Ncols());
  for (int i = 0; i < a.Nrows(); ++i)
    for (int j = 0; j < a.Ncols(); ++j)
      copy(a(i, j), b(i, j));
}

void copy(const Matrix& a, Matrix& b)
{
  if ((b.Nrows() != a.Nrows()) || (b.Ncols() != a.Ncols()))
    b.ReSize(a.Nrows(), a.Ncols());

#ifdef BLAS
  DCOPY((FORTINT) a.Storage(), a.Store(), (FORTINT) 1, b.Store(), (FORTINT) 1);
#else
  b = a;
#endif
}

//Make sure that b is already initialized with data
void copy(const ObjectMatrix<StackMatrix>& a, ObjectMatrix<StackMatrix>& b)
{
  assert(a.Nrows() == b.Nrows() && a.Ncols() == b.Ncols());
  for (int i = 0; i < a.Nrows(); ++i)
    for (int j = 0; j < a.Ncols(); ++j)
      copy(a(i, j), b(i, j));
}

//Make sure that b is already initialized with data
void copy(const ObjectMatrix<Matrix>& a, ObjectMatrix<StackMatrix>& b)
{
  assert(a.Nrows() == b.Nrows() && a.Ncols() == b.Ncols());
  for (int i = 0; i < a.Nrows(); ++i)
    for (int j = 0; j < a.Ncols(); ++j)
      copy(a(i, j), b(i, j));
}

//Make sure that b is already initialized with data
void copy(const ObjectMatrix<StackMatrix>& a, ObjectMatrix<Matrix>& b)
{
  assert(a.Nrows() == b.Nrows() && a.Ncols() == b.Ncols());
  for (int i = 0; i < a.Nrows(); ++i)
    for (int j = 0; j < a.Ncols(); ++j)
      copy(a(i, j), b(i, j));
}

//Make sure that b is already initialized with data
void copy(const Matrix& a, StackMatrix& b)
{
  assert ((b.Nrows() == a.Nrows()) && (b.Ncols() == a.Ncols()));

#ifdef BLAS
  DCOPY((FORTINT) a.Storage(), a.Store(), (FORTINT) 1, b.Store(), (FORTINT) 1);
#else
  b = a;
#endif
}

//Make sure that b is already initialized with data
void copy(const StackMatrix& a, Matrix& b)
{
  if ((b.Nrows() != a.Nrows()) || (b.Ncols() != a.Ncols())) 
    b.ReSize(a.Nrows(), a.Ncols());

#ifdef BLAS
  DCOPY((FORTINT) a.Storage(), a.Store(), (FORTINT) 1, b.Store(), (FORTINT) 1);
#else
  b = a;
#endif
}



//Make sure that b is already initialized with data
void copy(const StackMatrix& a, StackMatrix& b)
{
  assert ((b.Nrows() == a.Nrows()) && (b.Ncols() == a.Ncols()));

#ifdef BLAS
  DCOPY((FORTINT) a.Storage(), a.Store(), (FORTINT) 1, b.Store(), (FORTINT) 1);
#else
  b = a;
#endif
}


StackSparseMatrix& StackSparseMatrix::operator+=(const StackSparseMatrix& other)
{
  for (int i = 0; i < nrows(); ++i)
    for (int j = 0; j < ncols(); ++j)
      if (allowed(i, j))
	{
	  assert(other.allowed(i, j));
	  MatrixScaleAdd(1., other.operator_element(i, j), operator_element(i, j));
	}
  return *this;
}


}
