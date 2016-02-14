/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_ROTATION_MAT_HEADER
#define SPIN_ROTATION_MAT_HEADER 
#include <vector>
#include "multiarray.h"

namespace SpinAdapted{
  class StackSparseMatrix;
  class StateInfo;
  void UnCollectQuantaAlongRows(std::vector<Matrix>& rotation, const StateInfo&  sRow);
  void CollectQuantaAlongRows(std::vector<Matrix>& rotation, const StateInfo&  sRow);
void SaveRotationMatrix (const std::vector<int>& sites, const std::vector<Matrix>& m1, int state =-1);
void LoadRotationMatrix (const std::vector<int>& sites, std::vector<Matrix>& m1, int state=-1);
void diagonalise_dm(StackSparseMatrix& tracedMatrix, std::vector<DiagonalMatrix>& eigenMatrix);
void svd_densitymat(StackSparseMatrix& tracedMatrix, std::vector<DiagonalMatrix>& eigenMatrix);
void svd_densitymat(StackSparseMatrix& tracedMatrix, std::vector<Matrix>& U, std::vector<DiagonalMatrix>& eigenMatrix, std::vector<Matrix>& V);
void sort_weights(std::vector<DiagonalMatrix>& eigenMatrix, vector<pair<int, int> >& inorderwts, vector<vector<int> >& weightsbyquanta);
double assign_matrix_by_dm(std::vector<Matrix>& rotatematrix, std::vector<DiagonalMatrix>& eigenmatrix, StackSparseMatrix& transformmatrix, vector<pair<int, int> >& inorderwts, vector<vector<int> >& wtsbyquanta, int totalstatesbydm, int totalstatesbyquanta, int left_block_size, int right_block_size);
int PickSingleVectorAtRandom(std::vector<DiagonalMatrix>& d);
 int PickSingleNumberAtRandom(std::vector<double>& d);
void makeRotationMatrix(std::vector<DiagonalMatrix>&d, std::vector<Matrix>& U, std::vector<Matrix>& V, std::vector<Matrix>& rotation);
int get_total_states(const int &this_size, const int &other_size);
bool can_connect(int n, int spin, int right_block_size);

void allocate(const StateInfo& row, const StateInfo& col, std::vector<Matrix>& rotations);


}
#endif
