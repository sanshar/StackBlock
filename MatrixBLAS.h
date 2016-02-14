/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef SPIN_MATRIX_BLAS_HEADER
#define SPIN_MATRIX_BLAS_HEADER
#define WANT_STREAM
#include "StackMatrix.h"
#include <newmat.h>
#include <newmatio.h>
#include <cassert>
#include <cstdio>
#include <vector>
#include <map>
#include <fstream>
#include "blas_calls.h"
#include "ObjectMatrix.h"

using namespace std;


namespace SpinAdapted{

  template<class T> void MatrixScale(double d, T& a)
  {
#ifdef BLAS
    DSCAL(a.Storage(), d, a.Store(), 1);
#else
    a *= d;
#endif  
  }
  template<class T1, class T2> double MatrixDotProduct(const T1& a, const T2& b)
  {
    assert((a.Nrows() == b.Nrows()) && (a.Ncols() == b.Ncols()));
#ifdef BLAS
    return DDOT(a.Storage(), a.Store(), 1, b.Store(), 1);
#else
    abort();
#endif
  }
  template<class T> void MatrixNormalise(T& a)
  {
    double norm = MatrixDotProduct(a, a);
    MatrixScale(1./sqrt(norm), a);
  }
  
  template<class T> void Randomise (T& a)
  {
    Real* val = a.Store ();
    for (int i = 0; i < a.Storage (); ++i)
    {
      *val = double(rand ()) / RAND_MAX;
      ++val;
    }
  }
  
  template<class T> void SymmetricRandomise (T& a)
  {
    assert(a.Nrows() == a.Ncols());
    
    for (int i=0; i<a.Nrows(); i++)
      for (int j=0; j<i+1; j++) {
	a(i+1, j+1) = double(rand())/RAND_MAX;
	a(j+1, i+1) = a(i+1, j+1);
      }
    
  }

  template<class T1, class T2> void MatrixScaleAdd (double d, const T1& aref, T2& b){
#ifdef BLAS
    //ugly hack, sometimes aref is a TransposeMatrix which does not have the function Store()
    //so I cast it as a Matrix and then use the usual Matrix functions
    const T1& a = static_cast<const T1&>(aref); 
    assert (a.Nrows () == b.Nrows () && a.Ncols () == b.Ncols ());
    
    int n = a.Nrows () * a.Ncols ();
    assert (n == (b.Nrows () * b.Ncols ()));
    DAXPY (n, d, a.Store (), 1, b.Store (), 1);
#else
    b += d * a;
#endif
  }
  
  template<class T1, class T2, class T3> void MatrixTensorProduct (const T1& a_ref, char conjA, Real scaleA, const T2& b_ref, char conjB, Real scaleB, T3& c, int rowstride, int colstride){
#ifndef BLAS
    T1 A;
    T2 B;
#endif
    T1& a = const_cast<T1&>(a_ref); // for BLAS calls
    T2& b = const_cast<T2&>(b_ref);
    
    int arows = a.Nrows();
    int acols = a.Ncols();
    
    // some specialisations
#ifdef FAST_MTP
    //  if ((brows == 1) && (bcols == 1))
    {
      double b00 = *b.Store();
      if (conjA == 'n')
	{
	  double* cptr = c.Store()+ rowstride*c.Ncols() + colstride;
	  for (int i=0; i< a.Nrows();i++) 
	    DAXPY(a.Ncols(), scaleA * scaleB * b00, a.Store()+i*a.Ncols(), 1, cptr + i*c.Ncols(), 1);
	  return;
	}
      else 	
	{
	  double* aptr = a.Store();
	  double* cptr = c.Store() + rowstride*c.Ncols() + colstride;
	  for (int col = 0; col < acols; ++col)
	    {
	      DAXPY(arows, scaleA * scaleB * b00, aptr, acols, cptr, 1);
	      ++aptr;
	      cptr += c.Ncols();//arows;
	    }
	  
	  return;
	}	
    }
    //  else
    //    abort();
#else 
    try
      {
	if (conjA == 'n' && conjB == 'n')
	  {
#ifdef BLAS
	    int aRows = a.Nrows ();
	    int aCols = a.Ncols ();
	    int bRows = b.Nrows ();
	    int bCols = b.Ncols ();
	    
	    for (int i = 0; i < aRows; ++i)
	      for (int j = 0; j < aCols; ++j)
		{
		  Real scale = scaleA * scaleB * a (i+1,j+1);
		  for (int k = 0; k < bRows; ++k)
		    DAXPY (bCols, scale, &b (k+1,1), 1, &c (i * bRows + k+1 +rowstride,j * bCols+1+colstride), 1);
		}
	    return;
#else
	    A = a;
	    B = b;
#endif
	  }
	else if (conjA == 't' && conjB == 'n')
	  {
#ifdef BLAS
	    int aRows = a.Ncols ();
	    int aCols = a.Nrows ();
	    int bRows = b.Nrows ();
	    int bCols = b.Ncols ();
	    
	    for (int i = 0; i < aRows; ++i)
	      for (int j = 0; j < aCols; ++j)
		{
		  Real scale = scaleA * scaleB * a (j+1,i+1);
		  for (int k = 0; k < bRows; ++k)
		    DAXPY (bCols, scale, &b (k+1,1), 1, &c (i * bRows + k+1+rowstride,j * bCols+1+colstride), 1);
		}
	    return;
#else	  
	    A = a.t ();
	    B = b;
#endif
	  }
	else if (conjA == 'n' && conjB == 't')
	  {
#ifdef BLAS
	    int aRows = a.Nrows ();
	    int aCols = a.Ncols ();
	    int bRows = b.Ncols ();
	    int bCols = b.Nrows ();
	    
	    for (int i = 0; i < aRows; ++i)
	      for (int j = 0; j < aCols; ++j)
		{
		  Real scale = scaleA * scaleB * a (i+1,j+1);
		  for (int k = 0; k < bRows; ++k)
		    DAXPY (bCols, scale, &b (1,k+1), bRows, &c (i * bRows + k+1+rowstride,j * bCols+1+colstride), 1);
		}
	    return;
#else
	    A = a;
	    B = b.t ();
#endif
	  }
	else if (conjA == 't' && conjB == 't')
	  {
#ifdef BLAS
	    int aRows = a.Ncols ();
	    int aCols = a.Nrows ();
	    int bRows = b.Ncols ();
	    int bCols = b.Nrows ();
	    
	    for (int i = 0; i < aRows; ++i)
	      for (int j = 0; j < aCols; ++j)
		{
		  Real scale = scaleA * scaleB * a (j+1,i+1);
		  for (int k = 0; k < bRows; ++k)
		    DAXPY (bCols, scaleA * scaleB * a (j+1,i+1), &b (1,k+1), bRows, &c (i * bRows + k+1+rowstride,j * bCols+1+colstride), 1);
		}
	    return;
#else
	    A = a.t ();
	    B = b.t ();
#endif
	  }
	else
	  abort ();
#ifndef BLAS
	for (int i = 1; i <= A.Nrows (); ++i)
	  for (int j = 1; j <= A.Ncols (); ++j)
	    c.SubMatrix ((i - 1) * B.Nrows () + 1, i * B.Nrows (), (j - 1) * B.Ncols () + 1, j * B.Ncols ()) += (scaleA * scaleB) * A (i,j) * B; 
#endif
	
      }
    catch (Exception)
    {
      pout << Exception::what () << endl;
      abort ();
    }   
#endif
  }
  
  template<class T1, class T2> void MatrixMultiply (double d, const T1& a, T2& b)
    {
      //  b += d * a;
#ifdef BLAS 
      assert ((a.Nrows () == b.Nrows ()) && (a.Ncols () == b.Ncols ()));
      int n = a.Nrows () * a.Ncols ();
      DAXPY (n, d, a.Store (), 1, b.Store (), 1);
#else
      b += d * a;
#endif
    }
  
  template<class T1, class T2, class T3> void MatrixMultiply (const T1& a, char conjA, const T2& b, char conjB, T3& c, Real scale, double cfactor = 1.)
  {
    //dmrginp.justmultiply.start();
    //dmrginp.justmultiply -> start(); //ROA
    T1& a_ref = const_cast<T1&>(a); // for BLAS calls
    T2& b_ref = const_cast<T2&>(b);
    try
      {
	int aRows = a_ref.Nrows ();
	int aCols = a_ref.Ncols ();
	int bRows = b_ref.Nrows ();
	int bCols = b_ref.Ncols ();
	int cRows = c.Nrows ();
	int cCols = c.Ncols ();
	dmrginp.matmultNum++;
	if (conjA == 'n' && conjB == 'n')
	  {	  
	    dmrginp.matmultFlops[omprank] += aCols*cRows*cCols;
	    assert ((aCols == bRows) && (cRows == aRows) && (cCols == bCols));
#ifdef BLAS
	    dgemm_ (&conjA, &conjB, &bCols, &aRows, &bRows, &scale, b.Store (), &bCols, a.Store (), &aCols, &cfactor, c.Store (), &bCols);
#else
	    c += (scale * a) * b;
#endif
	  }
	else if (conjA == 'n' && conjB == 't')
	  {
	    dmrginp.matmultFlops[omprank] += aCols*cRows*cCols;
	    assert ((aCols == bCols) && (cRows == aRows) && (cCols == bRows));
#ifdef BLAS
	    dgemm_ (&conjB, &conjA, &bRows, &aRows, &bCols, &scale, b.Store (), &bCols, a.Store (), &aCols, &cfactor, c.Store (), &bRows);
#else
	    c += (scale * a) * b.t ();
#endif
	  } 
	else if (conjA == 't' && conjB == 'n')
	  {
	    dmrginp.matmultFlops[omprank] += aRows*cRows*cCols;
	    assert ((aRows == bRows) && (cRows == aCols) && (cCols == bCols));
#ifdef BLAS
	    dgemm_ (&conjB, &conjA, &bCols, &aCols, &bRows, &scale, b.Store (), &bCols, a.Store (), &aCols, &cfactor, c.Store (), &bCols);
#else
	    c += (scale * a.t ()) * b;
#endif
	  }
	else if (conjA == 't' && conjB == 't')
	  {
	    dmrginp.matmultFlops[omprank] += aRows*cRows*cCols;
	    assert ((aRows == bCols) && (cRows == aCols) && (cCols == bRows));
#ifdef BLAS
	    dgemm_ (&conjB, &conjA, &bRows, &aCols, &bCols, &scale, b.Store (), &bCols, a.Store (), &aCols, &cfactor, c.Store (), &bRows);
#else
	    c += (scale * a.t ()) * b.t ();
#endif
	  }
	else
	  abort ();
      }
    catch (Exception)
    {
      pout << Exception::what () << endl;
      abort ();
    }
    //dmrginp.justmultiply.stop();
    //dmrginp.justmultiply -> stop(); //ROA
  }
  
  template<class T> void Clear (T& a)
    {
#ifdef BLAS
      memset(a.Store(), 0, a.Storage() * sizeof(double));
#else
      a = 0.;
#endif
    }

  template<class T1, class T2> void MatrixRotate (T1 const& a, T2& b, T1 const& c, T2& d)
  {
    try
      {
	assert (d.Nrows () == a.Ncols () && d.Ncols () == c.Ncols ());
#ifdef BLAS
	T1 work (b.Nrows (), c.Ncols ());  
	Clear (work);
	MatrixMultiply (b, 'n', c, 'n', work, 1.);
	MatrixMultiply (a, 't', work, 'n', d, 1.);
#else
	d = a.t () * b * c;
#endif
      }
    catch (Exception)
    {
      pout << Exception::what () << endl;
      abort ();
    }
  }

template<class T> void MatrixDiagonalScale(double d, const T& a, double* b)
{
  //assert (a.Nrows () == a.Ncols () && a.Nrows () == b.Ncols ());
#ifdef BLAS
  int n = a.Nrows ();
  DAXPY (n, d, a.Store (), n+1, b, 1);
#else
  //b += d * a; Should add the non-blas analogue
#endif
}


  /*  
  template<> void MatrixRotate<Matrix> (const Matrix& a, const Matrix& b, const Matrix& c, Matrix& d)
  {
    try
      {
	assert (d.Nrows () == a.Ncols () && d.Ncols () == c.Ncols ());
#ifdef BLAS
	Matrix work (b.Nrows (), c.Ncols ());  
	Clear (work);
	MatrixMultiply (b, 'n', c, 'n', work, 1.);
	MatrixMultiply (a, 't', work, 'n', d, 1.);
#else
	d = a.t () * b * c;
#endif
      }
    catch (Exception)
    {
      pout << Exception::what () << endl;
      abort ();
    }
  }
  */
  
  void svd(Matrix& M, DiagonalMatrix& e, Matrix& U, Matrix& V);
  void svd(StackMatrix& M, DiagonalMatrix& e, Matrix& U, Matrix& V);
  void xsolve_AxeqB(const Matrix& a, const ColumnVector& b, ColumnVector& x);
  //void MatrixDiagonalScale(double d, const Matrix& a, double* b);
  
  
  double dotproduct(const ColumnVector& a, const ColumnVector& b);
  double dotproduct(const RowVector& a, const RowVector& b);
  double rowdoubleproduct(Matrix& a, int rowa, Matrix& b, int rowb);
  void diagonalise(Matrix& sym, DiagonalMatrix& d, Matrix& vec);
  void diagonalise(StackMatrix& sym, DiagonalMatrix& d);
  void diagonalise_tridiagonal(std::vector<double>& diagonal, std::vector<double>& offdiagonal, int numelements, Matrix& vec);
  
  
  void CatenateProduct (const ObjectMatrix<Matrix*>& a, Matrix& b, bool allocate = false);  
  void CatenateProduct (const ObjectMatrix<StackMatrix*>& a, StackMatrix& b);  
  char TransposeOf (char c);
  inline char TransposeOf (char c) { assert (c == 'n' || c == 't'); return (c == 'n') ? 't' : 'n'; }
  void MatrixRotate (const Matrix& a, const Matrix& b, const Matrix& c, Matrix& d);
  void Save (const Matrix& a, std::ofstream &ofs);
  void Load (Matrix& a, std::ifstream &ifs);
  
  template<class T> void print(ostream& os, const std::vector<T>& v)
    {
      for (int i = 0; i < v.size(); ++i)
	os << v[i] << " ";
      os << endl;
    }
  double CheckSum (Matrix& a);
  void DebugPrint (std::vector<int>& v);
  void DebugPrint (std::vector<double>& v);
  
  
  template<class T> void CyclicPermute (std::vector<T>& v, std::vector<T>& result, int n);
  
  template<class T> void CyclicPermute (std::vector<T>& v, std::vector<T>& result, int n)
    {
      assert (n >= 0);
      result.resize (v.size ());
      for (int i = 0; i < v.size (); ++i)
	result [i] = v [(i + n) % v.size ()];
    }
  
  template<class T> bool find(vector<T>& v, const T& value) { return std::find(v.begin(), v.end(), value) != v.end(); }
  
  template<class T> void get_sorted_indices(const std::vector<T>& data, std::vector<int>& indices)
    {
      multimap<T, int> sorted_data;
      for (int i = 0; i < data.size(); ++i)
	sorted_data.insert(make_pair<T, int>(data[i], i));
      typename multimap<T, int>::iterator it = sorted_data.begin();
      indices.reserve(data.size());
      while (it != sorted_data.end())
	{
	  indices.push_back(it->second);
	  ++it;
	}
      
    }
}
#endif
