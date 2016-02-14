/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#include <iostream>
#include <cmath>
#include <include/newmatutils.h>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include "MatrixBLAS.h"
#include "global.h"
#ifdef BLAS
#include "blas_calls.h"
#endif


double SpinAdapted::dotproduct(const ColumnVector& a, const ColumnVector& b)
{
  assert(a.Nrows() == b.Nrows());
#ifdef BLAS
  return DDOT(a.Storage(), a.Store(), 1, b.Store(), 1);
#else
  return a.t() * b;
#endif
}

double SpinAdapted::dotproduct(const RowVector& a, const RowVector& b)
{
  assert(a.Ncols() == b.Ncols());
#ifdef BLAS
  return DDOT(a.Storage(), a.Store(), 1, b.Store(), 1);
#else
  return a * b.t();
#endif
}

double SpinAdapted::rowdoubleproduct(Matrix& a, int rowa, Matrix& b, int rowb)
{
  assert(a.Ncols() == b.Ncols());
  double* aptr = a.Store() + a.Ncols() * rowa;
  double* bptr = b.Store() + b.Ncols() * rowb;
  return DDOT(a.Ncols(), aptr, 1, bptr, 1);
}




void SpinAdapted::xsolve_AxeqB(const Matrix& a, const ColumnVector& b, ColumnVector& x)
{
  FORTINT ar = a.Nrows();
  int bc = 1;
  int info=0;
  FORTINT* ipiv = new FORTINT[ar];
  double* bwork = new double[ar];
  for(int i = 0;i<ar;++i)
    bwork[i] = b.element(i);
  double* workmat = new double[ar*ar];
  for(int i = 0;i<ar;++i)
    for(int j = 0;j<ar;++j)
      workmat[i*ar+j] = a.element(j,i);

  GESV(ar, bc, workmat, ar, ipiv, bwork, ar, info);
  delete[] ipiv;
  delete[] workmat;

  for(int i = 0;i<ar;++i)
    x.element(i) = bwork[i];

  delete[] bwork;

  if(info != 0)
  {
     pout << "Xsolve failed with info error " << info << endl;
     abort();
  }
}

void SpinAdapted::svd(Matrix& M, DiagonalMatrix& d, Matrix& U, Matrix& V)
{
  int nrows = M.Nrows();
  int ncols = M.Ncols();

  assert(nrows >= ncols);

  int minmn = min(nrows, ncols);
  int maxmn = max(nrows, ncols);
  int eigenrows = min(minmn, minmn);
  d.ReSize(minmn);
  Matrix Ut;
  Ut.ReSize(nrows, nrows);
  V.ReSize(ncols, ncols);

  int lwork = maxmn * maxmn + 100;
  double* workspace = new double[lwork];

  // first transpose matrix
  Matrix Mt;
  Mt = M.t();
  int info = 0;
  DGESVD('A', 'A', nrows, ncols, Mt.Store(), nrows, d.Store(), 
	 Ut.Store(), nrows, V.Store(), ncols, workspace, lwork, info);

  U.ReSize(nrows, ncols);
  SpinAdapted::Clear(U);
  for (int i = 0; i < nrows; ++i)
    for (int j = 0; j < ncols; ++j)
      U(i+1,j+1) = Ut(j+1,i+1);
  delete[] workspace;
}

void SpinAdapted::svd(StackMatrix& M, DiagonalMatrix& d, Matrix& U, Matrix& V)
{
  int nrows = M.Nrows();
  int ncols = M.Ncols();

  assert(nrows >= ncols);

  int minmn = min(nrows, ncols);
  int maxmn = max(nrows, ncols);
  int eigenrows = min(minmn, minmn);
  d.ReSize(minmn);
  Matrix Ut;
  Ut.ReSize(nrows, nrows);
  V.ReSize(ncols, ncols);

  int lwork = maxmn * maxmn + 100;
  double* workspace = new double[lwork];

  // first transpose matrix
  Matrix Mt;
  //Mt = M.t();
  int info = 0;
  DGESVD('A', 'A', nrows, ncols, Mt.Store(), nrows, d.Store(), 
	 Ut.Store(), nrows, V.Store(), ncols, workspace, lwork, info);

  U.ReSize(nrows, ncols);
  SpinAdapted::Clear(U);
  for (int i = 0; i < nrows; ++i)
    for (int j = 0; j < ncols; ++j)
      U(i+1,j+1) = Ut(j+1,i+1);
  delete[] workspace;
}

void SpinAdapted::diagonalise_tridiagonal(std::vector<double>& diagonal, std::vector<double>& offdiagonal, int numelements, Matrix& vec)
{
  int nrows = numelements;
  int ncols = numelements;

  vec.ReSize(nrows, nrows);

  Matrix vec_transpose; vec_transpose = vec;
  vector<double> workarray(4*nrows-2,0);
  int info = 0;

  DSTEV('V', nrows, &(diagonal[0]), &(offdiagonal[0]), vec_transpose.Store(), nrows,  &(workarray[0]), info);

  if (info != 0)
    {
      pout << "failed to converge :: " <<info<< endl;
      abort();
    }

  for (int i = 0; i < nrows; ++i)
    for (int j = 0; j < ncols; ++j)
      vec(j+1,i+1) = vec_transpose(i+1,j+1);
}

void SpinAdapted::diagonalise(Matrix& sym, DiagonalMatrix& d, Matrix& vec)
{
  int nrows = sym.Nrows();
  int ncols = sym.Ncols();
  assert(nrows == ncols);
  d.ReSize(nrows);
  vec.ReSize(nrows, nrows);

  Matrix workmat;
  workmat = sym;
  vector<double> workquery(1);
  int info = 0;
  double* dptr = d.Store();

  int query = -1;
  DSYEV('V', 'L', nrows, workmat.Store(), nrows, dptr, &(workquery[0]), query, info); // do query to find best size
  
  int optlength = static_cast<int>(workquery[0]);
  vector<double> workspace(optlength);

  DSYEV('V', 'U', nrows, workmat.Store(), nrows, dptr, &(workspace[0]), optlength, info); // do query to find best size


  
  if (info > 0) 
    {
      pout << "failed to converge " << endl;
      abort(); 
    }
  
  for (int i = 0; i < nrows; ++i)
    for (int j = 0; j < ncols; ++j)
      vec(j+1,i+1) = workmat(i+1,j+1);
}



void SpinAdapted::diagonalise(StackMatrix& sym, DiagonalMatrix& d)
{
  int nrows = sym.Nrows();
  int ncols = sym.Ncols();
  assert(nrows == ncols);
  d.ReSize(nrows);

  std::vector<double> workmat(nrows*ncols, 0.0);
  DCOPY(nrows*ncols, sym.Store(), 1, &workmat[0], 1);
  vector<double> workquery(1);
  int info = 0;
  double* dptr = d.Store();

  int query = -1;
  DSYEV('V', 'L', nrows, &workmat[0], nrows, dptr, &(workquery[0]), query, info); // do query to find best size
  
  int optlength = static_cast<int>(workquery[0]);
  vector<double> workspace(optlength);

  DSYEV('V', 'U', nrows, &workmat[0], nrows, dptr, &(workspace[0]), optlength, info); // do query to find best size


  
  if (info > 0) 
    {
      pout << "failed to converge " << endl;
      abort(); 
    }
  
  for (int i = 0; i < nrows; ++i)
    for (int j = 0; j < ncols; ++j)
      sym(j+1,i+1) = workmat[i*ncols+j];
}




double SpinAdapted::CheckSum (Matrix& a)
{
  double val = 0.;
  for (int i = 0; i < a.Nrows (); ++i)
    for (int j = 0; j < a.Ncols (); ++j)
      val += a.element (i, j);
  return val;
}

void SpinAdapted::CatenateProduct (const ObjectMatrix<Matrix*>& a, Matrix& b, bool allocate)
{
  try
    {
      std::vector<int> indexRows (a.Nrows ());
      std::vector<int> indexCols (a.Ncols ());
      int rowLength = 0;
      int colLength = 0;
      for (int i = 0; i < indexRows.size (); ++i)
	{
	  indexRows [i] = (i > 0) ? a (i - 1,0)->Nrows () + indexRows [i - 1] : 1;
	  rowLength += a (i,0)->Nrows ();
	}
      for (int i = 0; i < indexCols.size (); ++i)
	{
	  indexCols [i] = (i > 0) ? a (0,i - 1)->Ncols () + indexCols [i - 1] : 1;
	  colLength += a (0,i)->Ncols ();
	}
      
      if (!allocate) 
	assert (b.Nrows () == rowLength && b.Ncols () == colLength); // precondition
      else
	b.ReSize (rowLength, colLength);

      for (int i = 0; i < a.Nrows (); ++i)
	for (int j = 0; j < a.Ncols (); ++j)
	  {
#ifdef BLAS
	    int bcols = b.Ncols();
	    double* bptr = b.Store() + bcols * (indexRows[i] - 1) + (indexCols[j] - 1);
	    Matrix* aij = a(i, j);
	    double* aptr = aij->Store();
	    int nrows = aij->Nrows();
	    int ncols = aij->Ncols();
	    for (int r = 0; r < nrows; ++r)
	      {
		DCOPY(ncols, aptr, 1, bptr, 1);
		aptr += ncols;
		bptr += bcols;
	      }
#else
	    b.SubMatrix (indexRows [i], indexRows [i] + a (i,j)->Nrows () - 1, indexCols [j], indexCols [j] + a (i,j)->Ncols () - 1) = *(a (i,j));
#endif
	  }
    }
  catch (Exception)
    {
      pout << Exception::what () << endl;
      abort ();
    }
}


void SpinAdapted::CatenateProduct (const ObjectMatrix<StackMatrix*>& a, StackMatrix& b)
{
  try
    {
      std::vector<int> indexRows (a.Nrows ());
      std::vector<int> indexCols (a.Ncols ());
      int rowLength = 0;
      int colLength = 0;
      for (int i = 0; i < indexRows.size (); ++i)
	{
	  indexRows [i] = (i > 0) ? a (i - 1,0)->Nrows () + indexRows [i - 1] : 1;
	  rowLength += a (i,0)->Nrows ();
	}
      for (int i = 0; i < indexCols.size (); ++i)
	{
	  indexCols [i] = (i > 0) ? a (0,i - 1)->Ncols () + indexCols [i - 1] : 1;
	  colLength += a (0,i)->Ncols ();
	}
      

      for (int i = 0; i < a.Nrows (); ++i)
	for (int j = 0; j < a.Ncols (); ++j)
	  {
#ifdef BLAS
	    int bcols = b.Ncols();
	    double* bptr = b.Store() + bcols * (indexRows[i] - 1) + (indexCols[j] - 1);
	    StackMatrix* aij = a(i, j);
	    double* aptr = aij->Store();
	    int nrows = aij->Nrows();
	    int ncols = aij->Ncols();
	    for (int r = 0; r < nrows; ++r)
	      {
		DCOPY(ncols, aptr, 1, bptr, 1);
		aptr += ncols;
		bptr += bcols;
	      }
#else
	    b.SubMatrix (indexRows [i], indexRows [i] + a (i,j)->Nrows () - 1, indexCols [j], indexCols [j] + a (i,j)->Ncols () - 1) = *(a (i,j));
#endif
	  }
    }
  catch (Exception)
    {
      pout << Exception::what () << endl;
      abort ();
    }
}

  

void SpinAdapted::Save(const Matrix& a, std::ofstream &ofs)
{
  boost::archive::binary_oarchive save_mat(ofs);
  save_mat << a;
}

void SpinAdapted::Load(Matrix& a, std::ifstream &ifs)
{
  boost::archive::binary_iarchive load_mat(ifs);
  load_mat >> a;
}

void SpinAdapted::DebugPrint (vector<int>& v)
{
  for (int i = 0; i < v.size(); ++i)
    pout << v[i] << endl;
}
void SpinAdapted::DebugPrint (vector<double>& v)
{
  for (int i = 0; i < v.size(); ++i)
    pout << v[i] << endl;
}

