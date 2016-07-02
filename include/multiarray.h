/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/


#ifndef multi_array_h
#define multi_array_h
#include <vector>
#include <string.h>
#include <iostream>
#include <stdio.h>
#include <cassert>
#include <iostream>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/serialization.hpp>

using namespace std;
template<class T> class multiarray
{
public:
  virtual T& operator() (const vector<int>& indices) = 0;
  virtual T operator() (const vector<int>& indices) const = 0;
  virtual string reflect_type() const = 0;
};

template<class T> class array_2d : public vector<T>, public multiarray<T>
{
 public:
  array_2d () : dim1_d (0), dim2_d (0), vector<T> () { }
  array_2d (const int d1, const int d2) : dim1_d (d1), dim2_d (d2), vector<T> (d1 * d2) { }
  void Clear()
  {
    for (int i=0; i<this->size(); ++i)
      (*this)[i]=0.;
  }
  T& operator() (const int i, const int j) 
    { 
      assert ((0 <= i) && (i < dim1_d)); assert ((0 <= j) && (j < dim2_d));
      return vector<T>::operator[] (i * dim2_d + j); 
    }
  T operator() (const int i, const int j) const
    { 
      assert ((0 <= i) && (i < dim1_d)); assert ((0 <= j) && (j < dim2_d));
      return vector<T>::operator[] (i * dim2_d + j); 
    }

  T& operator() (const vector<int>& indices)
  {
    assert(indices.size() == 2);
    return operator()(indices[0], indices[1]);
  }
  T operator() (const vector<int>& indices) const
  {
    assert(indices.size() == 2);
    return operator()(indices[0], indices[1]);
  }
  void resize (const int i, const int j) { vector<T>::resize (i * j); dim1_d = i; dim2_d = j; }
  int dim1() const { return dim1_d; }
  int dim2() const { return dim2_d; }
  string reflect_type() const { return string("array_2d"); }

 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & dim1_d & dim2_d;
    ar & boost::serialization::base_object<std::vector<T> >(*this);
  }
  int dim1_d;
  int dim2_d;
};

template<class T> class array_3d : public vector<T>, public multiarray<T>
{
 public:
  array_3d() : dim1_d (0), dim2_d (0), dim3_d (0), dim2_times_dim3_d (0), vector<T> () { }
  array_3d(const int d1, const int d2, const int d3) : dim1_d (d1), dim2_d (d2), dim3_d (d3), dim2_times_dim3_d (d2 * d3), vector<T> (d1 * d2 * d3) { }
  T& operator() (const int i, const int j, const int k) 
    { 
      assert((0 <= i) && (i < dim1_d)); assert((0 <= j) && (j < dim2_d)); assert((0 <= k) && (k < dim3_d));
      return vector<T>::operator[](i * dim2_times_dim3_d + j * dim3_d + k); 
    }
  void Clear()
  {
    for (int i=0; i<this->size(); ++i)
      (*this)[i]=0.;
  }
  T operator() (const int i, const int j, const int k) const
    { 
      assert((0 <= i) && (i < dim1_d)); assert((0 <= j) && (j < dim2_d)); assert((0 <= k) && (k < dim3_d));
      return vector<T>::operator[](i * dim2_times_dim3_d + j * dim3_d + k); 
    }
  T& operator() (const vector<int>& indices)
  {
    assert(indices.size() == 3);
    return operator()(indices[0], indices[1], indices[2]);
  }
  T operator() (const vector<int>& indices) const
  {
    assert(indices.size() == 3);
    return operator()(indices[0], indices[1], indices[2]);
  }

  void resize(const int i, const int j, const int k) { vector<T>::resize(i * j * k); dim1_d = i; dim2_d = j; dim3_d = k; dim2_times_dim3_d = dim2_d * dim3_d; }
  int dim1() const { return dim1_d; }
  int dim2() const { return dim2_d; }
  int dim3() const { return dim3_d; }
  string reflect_type() const { return string("array_3d"); }

 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & dim1_d & dim2_d & dim3_d & dim2_times_dim3_d;
    ar & boost::serialization::base_object<vector<T> >(*this);
  }
  int dim1_d;
  int dim2_d;
  int dim3_d;
  int dim2_times_dim3_d;
};

template<class T> class array_4d : public vector<T>, public multiarray<T>
{
 public:
  array_4d() : dim1_d(0), dim2_d(0), dim3_d(0), dim4_d(0), dim3_times_dim4_d(0), dim2_times_dim3_times_dim4_d (0), vector<T>() { }
  array_4d(const int d1, const int d2, const int d3, const int d4) : dim1_d (d1), dim2_d (d2), dim3_d (d3), dim4_d(d4), 
    dim3_times_dim4_d(d3 * d4), dim2_times_dim3_times_dim4_d(d2 * d3 * d4), vector<T> (d1 * d2 * d3 * d4) { }
  T& operator()(const int i, const int j, const int k, const int l) 
    { 
      assert((0 <= i) && (i < dim1_d)); assert((0 <= j) && (j < dim2_d)); assert((0 <= k) && (k < dim3_d)); assert((0 <= l) && (l < dim4_d));
      return vector<T>::operator[](i * dim2_times_dim3_times_dim4_d + j * dim3_times_dim4_d + k * dim4_d + l); 
    }
  T operator()(const int i, const int j, const int k, const int l) const
    { 
      assert((0 <= i) && (i < dim1_d)); assert((0 <= j) && (j < dim2_d)); assert((0 <= k) && (k < dim3_d)); assert((0 <= l) && (l < dim4_d));
      return vector<T>::operator[](i * dim2_times_dim3_times_dim4_d + j * dim3_times_dim4_d + k * dim4_d + l); 
    }
  array_4d<T>& operator+=(const array_4d<T>& C)
  {
    assert(dim1() == C.dim1() &&
           dim2() == C.dim2() &&
           dim3() == C.dim3() &&
           dim4() == C.dim4());
    for (int i=0; i<C.size(); ++i)
      (*this)[i]+=C[i];
    return *this;
  }
  void Clear()
  {
    for (int i=0; i<this->size(); ++i)
      (*this)[i]=0.;
  }
  T& operator() (const vector<int>& indices)
  {
    assert(indices.size() == 4);
    return operator()(indices[0], indices[1], indices[2], indices[3]);
  }
  T operator() (const vector<int>& indices) const
  {
    assert(indices.size() == 4);
    return operator()(indices[0], indices[1], indices[2], indices[3]);
  }

  void resize (const int i, const int j, const int k, const int l) 
  { 
    vector<T>::resize (i * j * k * l); dim1_d = i; dim2_d = j; dim3_d = k; dim4_d = l; 
    dim2_times_dim3_times_dim4_d = dim2_d * dim3_d * dim4_d; dim3_times_dim4_d = dim3_d * dim4_d;
  }
  string reflect_type() const { return string("array_4d"); }

  int dim1() const { return dim1_d; }
  int dim2() const { return dim2_d; }
  int dim3() const { return dim3_d; }
  int dim4() const { return dim4_d; } 
 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & dim1_d & dim2_d & dim3_d & dim4_d & dim2_times_dim3_times_dim4_d & dim3_times_dim4_d;
    ar & boost::serialization::base_object<std::vector<T> >(*this);
  }
  int dim1_d;
  int dim2_d;
  int dim3_d;
  int dim4_d;
  int dim2_times_dim3_times_dim4_d;
  int dim3_times_dim4_d;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

template<class T> class array_6d : public multiarray<T>
{
 public:

  array_6d() : dim1_d(0), dim2_d(0), dim3_d(0), dim4_d(0), dim5_d(0), dim6_d(0), 
    dim5xdim6_d(0),
    dim4xdim5xdim6_d(0),
    dim3xdim4xdim5xdim6_d(0),
    dim2xdim3xdim4xdim5xdim6_d(0),
    data(0) { }

 array_6d(const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, T* pdata) : 
    dim1_d (d1), dim2_d (d2), dim3_d (d3), dim4_d(d4), dim5_d(d5), dim6_d(d6), 
    dim5xdim6_d(d5*d6),
    dim4xdim5xdim6_d(d4*d5*d6),
    dim3xdim4xdim5xdim6_d(d3*d4*d5*d6),
      dim2xdim3xdim4xdim5xdim6_d(d2*d3*d4*d5*d6),
      data(pdata)
      {};
//std::cout << "array_6d:  " << d1 << " "
//<< d2 << " "
//<< d3 << " "
//<< d4 << " "
//<< d5 << " "
//<< d6 << std::endl;
// }

  T& operator()(const size_t i, const size_t j, const size_t k, const size_t l, const size_t m, const size_t n) 
    { 
      assert((0 <= i) && (i < dim1_d)); 
      assert((0 <= j) && (j < dim2_d)); 
      assert((0 <= k) && (k < dim3_d)); 
      assert((0 <= l) && (l < dim4_d));
      assert((0 <= m) && (m < dim5_d));
      assert((0 <= n) && (n < dim6_d));
      size_t p = i * dim2xdim3xdim4xdim5xdim6_d 
            + j * dim3xdim4xdim5xdim6_d 
            + k * dim4xdim5xdim6_d 
            + l * dim5xdim6_d 
            + m * dim6_d 
            + n; 
      return data[p];//vector<T>::operator[](p);
    }

  T& operator()(const size_t i, const size_t j, const size_t k, const size_t l, const size_t m, const size_t n) const
    { 
      assert((0 <= i) && (i < dim1_d)); 
      assert((0 <= j) && (j < dim2_d)); 
      assert((0 <= k) && (k < dim3_d)); 
      assert((0 <= l) && (l < dim4_d));
      assert((0 <= m) && (m < dim5_d));
      assert((0 <= n) && (n < dim6_d));
      size_t p = i * dim2xdim3xdim4xdim5xdim6_d 
            + j * dim3xdim4xdim5xdim6_d 
            + k * dim4xdim5xdim6_d 
            + l * dim5xdim6_d 
            + m * dim6_d 
            + n; 
      return data[p];//vector<T>::operator[](p);
    }

  array_6d<T>& operator+=(const array_6d<T>& C)
  {
    assert(dim1() == C.dim1() &&
           dim2() == C.dim2() &&
           dim3() == C.dim3() &&
           dim4() == C.dim4() &&
           dim5() == C.dim5() &&
           dim6() == C.dim6());

    for (size_t i=0; i<C.size(); ++i)
      (*this)[i]+=C[i];
    return *this;
  }

  void Clear()
  {
    size_t len = get_size()*sizeof(T);
    //if (mpigetrank() == 0) {
    //cout << len <<"  "<<sizeof(T)<<endl;
    //cout << "zeroing out "<<len<<" bytes "<<endl;
      //}
    memset(data, 0, len);
  }

  size_t get_size()
  {
    //if (mpigetrank() == 0) 
    //cout << dim1_d<<" "<<dim2_d<<"  "<<dim3_d<<"  "<<dim4_d<<"  "<<dim5_d<<"  "<<dim6_d<<endl;
    return dim1_d*dim2_d*dim3_d*dim4_d*dim5_d*dim6_d;
  }

  T& operator() (const vector<int>& indices)
  {
    assert(indices.size() == 6);
    return operator()(indices[0], indices[1], indices[2], indices[3], indices[4], indices[5]);
  }

  T operator() (const vector<int>& indices) const
  {
    assert(indices.size() == 6);
    return operator()(indices[0], indices[1], indices[2], indices[3], indices[4], indices[5]);
  }

  void resize (const size_t i, const size_t j, const size_t k, const size_t l, const size_t m, const size_t n, T* pdata) 
  { 
    if (data != 0) {
      std::cout << "cannot resize a non-zero sized array."<<std::endl;
      abort();
    }
    data = pdata;
    //data = (T*)(Stackmem[omprank].allocate(i*j*k*l*m*n* (sizeof(T)/sizeof(double) + sizeof(T)%sizeof(double))));}
    //data = (T*)(SpinAdapted::Stackmem[omprank].allocate(i*j*k*l*m*n* (sizeof(T)/sizeof(double) + sizeof(T)%sizeof(double))));
    dim1_d = i; 
    dim2_d = j; 
    dim3_d = k; 
    dim4_d = l; 
    dim5_d = m; 
    dim6_d = n; 
    dim5xdim6_d = m*n;
    dim4xdim5xdim6_d = l*m*n;
    dim3xdim4xdim5xdim6_d = k*l*m*n;
    dim2xdim3xdim4xdim5xdim6_d = j*k*l*m*n;
  }
  string reflect_type() const { return string("array_6d"); }

  size_t dim1() const { return dim1_d; }
  size_t dim2() const { return dim2_d; }
  size_t dim3() const { return dim3_d; }
  size_t dim4() const { return dim4_d; } 
  size_t dim5() const { return dim5_d; } 
  size_t dim6() const { return dim6_d; } 


  T* data;
 private:
  /*
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & dim1_d & dim2_d & dim3_d & dim4_d & dim5_d & dim6_d 
       & dim2xdim3xdim4xdim5xdim6_d
       & dim3xdim4xdim5xdim6_d
       & dim4xdim5xdim6_d
       & dim5xdim6_d;
    ar & boost::serialization::base_object<std::vector<T> >(*this);
  }
  */
  size_t dim1_d;
  size_t dim2_d;
  size_t dim3_d;
  size_t dim4_d;
  size_t dim5_d;
  size_t dim6_d;
  size_t dim2xdim3xdim4xdim5xdim6_d;
  size_t dim3xdim4xdim5xdim6_d;
  size_t dim4xdim5xdim6_d;
  size_t dim5xdim6_d;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

//Internally the 3RDM is stored as (i,j,k,l,m,n) = E^{i,j,k}_{n,m,l}
//we will store the array as a 3 dimention array with composite indices
// (i,n)  (j,m)  and (k,l). We will only store terms where these three composite indices
//are arranged in ascending order. All the 6 permutations have same values
template<class T> class array_6d_3rdm : public multiarray<T>
{
 public:

  array_6d_3rdm() : dim1_d(0), dim2_d(0), dim3_d(0), dim4_d(0), dim5_d(0), dim6_d(0), 
    dimA_d(0),
    dimB_d(0),
    dimC_d(0),
    data(0) { }

 array_6d_3rdm(const size_t d, T* pdata) : 
    dim1_d (d), dim2_d (d), dim3_d (d), dim4_d(d), dim5_d(d), dim6_d(d), 
    dimA_d(d*d),
    dimB_d(d*d),
    dimC_d(d*d),
      data(pdata)
      {};
//std::cout << "array_6d:  " << d1 << " "
//<< d2 << " "
//<< d3 << " "
//<< d4 << " "
//<< d5 << " "
//<< d6 << std::endl;
// }

  T& operator()(const size_t i, const size_t j, const size_t k, const size_t l, const size_t m, const size_t n) 
    { 
      size_t a = i*dim1_d+n, b = j*dim1_d+m, c = k*dim1_d+l;
      size_t A=0, B=0, C=0;
      A = max(a,max(b,c));
      if (A==a) {
	B = max(b,c); C = min(b,c);
      }
      else if (A==b) {
	B = max(a,c); C = min(a,c);
      }
      else {
	B = max(a,b); C = min(a,b);
      }

      size_t p = (A*A*A + 3*A*A + 2*A)/6  +  (B*B + B)/2 + C ;

      return data[p];//vector<T>::operator[](p);
    }

  T& operator()(const size_t i, const size_t j, const size_t k, const size_t l, const size_t m, const size_t n) const
    { 
      size_t a = i*dim1_d+n, b = j*dim1_d+m, c = k*dim1_d+l;
      size_t A=0, B=0, C=0;
      A = max(a,max(b,c));
      if (A==a) {
	B = max(b,c); C = min(b,c);
      }
      else if (A==b) {
	B = max(a,c); C = min(a,c);
      }
      else {
	B = max(a,b); C = min(a,b);
      }

      size_t p = (A*A*A + 3*A*A + 2*A)/6  +  (B*B + B)/2 + C ;

      return data[p];//vector<T>::operator[](p);
    }

  array_6d_3rdm<T>& operator+=(const array_6d_3rdm<T>& C)
  {
    for (size_t i=0; i<C.size(); ++i)
      (*this)[i]+=C[i];
    return *this;
  }

  void Clear()
  {
    size_t len = get_size()*sizeof(T);
    //if (mpigetrank() == 0) {
    //cout << len <<"  "<<sizeof(T)<<endl;
    //cout << "zeroing out "<<len<<" bytes "<<endl;
      //}
    memset(data, 0, len);
  }

  size_t get_size()
  {
    //if (mpigetrank() == 0) 
    //cout << dim1_d<<" "<<dim2_d<<"  "<<dim3_d<<"  "<<dim4_d<<"  "<<dim5_d<<"  "<<dim6_d<<endl;
    size_t A = dim1_d*dim1_d;
    return (A*A*A +  3*A*A + 2*A)/6 ;
  }

  T& operator() (const vector<int>& indices)
  {
    assert(indices.size() == 6);
    return operator()(indices[0], indices[1], indices[2], indices[3], indices[4], indices[5]);
  }

  T operator() (const vector<int>& indices) const
  {
    assert(indices.size() == 6);
    return operator()(indices[0], indices[1], indices[2], indices[3], indices[4], indices[5]);
  }

  void resize (const size_t i, T* pdata) 
  { 
    if (data != 0) {
      std::cout << "cannot resize a non-zero sized array."<<std::endl;
      abort();
    }
    data = pdata;
    //data = (T*)(Stackmem[omprank].allocate(i*j*k*l*m*n* (sizeof(T)/sizeof(double) + sizeof(T)%sizeof(double))));}
    //data = (T*)(SpinAdapted::Stackmem[omprank].allocate(i*j*k*l*m*n* (sizeof(T)/sizeof(double) + sizeof(T)%sizeof(double))));
    dim1_d = i; 
    dim2_d = i; 
    dim3_d = i; 
    dim4_d = i; 
    dim5_d = i; 
    dim6_d = i; 
    dimA_d = i*i;
    dimB_d = i*i;
    dimC_d = i*i;
  }
  string reflect_type() const { return string("array_6d_3rdm"); }

  size_t dim1() const { return dim1_d; }
  size_t dim2() const { return dim2_d; }
  size_t dim3() const { return dim3_d; }
  size_t dim4() const { return dim4_d; } 
  size_t dim5() const { return dim5_d; } 
  size_t dim6() const { return dim6_d; } 


  T* data;
 private:
  /*
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & dim1_d & dim2_d & dim3_d & dim4_d & dim5_d & dim6_d 
       & dim2xdim3xdim4xdim5xdim6_d
       & dim3xdim4xdim5xdim6_d
       & dim4xdim5xdim6_d
       & dim5xdim6_d;
    ar & boost::serialization::base_object<std::vector<T> >(*this);
  }
  */
  size_t dim1_d;
  size_t dim2_d;
  size_t dim3_d;
  size_t dim4_d;
  size_t dim5_d;
  size_t dim6_d;
  size_t dimA_d;
  size_t dimB_d;
  size_t dimC_d;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

template<class T> class array_8d : public vector<T>, public multiarray<T>
{
 public:

  array_8d() : 
    dim1_d (0), dim2_d (0), dim3_d (0), dim4_d(0), dim5_d(0), dim6_d(0), dim7_d(0), dim8_d(0), 
    dim7xdim8_d(0),
    dim6xdim7xdim8_d(0),
    dim5xdim6xdim7xdim8_d(0),
    dim4xdim5xdim6xdim7xdim8_d(0),
    dim3xdim4xdim5xdim6xdim7xdim8_d(0),
    dim2xdim3xdim4xdim5xdim6xdim7xdim8_d(0),
    vector<T> () { }

  array_8d(const int d1, const int d2, const int d3, const int d4, const int d5, const int d6, const int d7, const int d8) : 
    dim1_d (d1), dim2_d (d2), dim3_d (d3), dim4_d(d4), dim5_d(d5), dim6_d(d6), dim7_d(d7), dim8_d(d8), 
    dim7xdim8_d( d7*d8 ),
    dim6xdim7xdim8_d( d6*d7*d8 ),
    dim5xdim6xdim7xdim8_d( d5*d6*d7*d8 ),
    dim4xdim5xdim6xdim7xdim8_d( d4*d5*d6*d7*d8 ),
    dim3xdim4xdim5xdim6xdim7xdim8_d( d3*d4*d5*d6*d7*d8 ),
    dim2xdim3xdim4xdim5xdim6xdim7xdim8_d( d2*d3*d4*d5*d6*d7*d8 ),
    vector<T> ( d1*d2*d3*d4*d5*d6*d7*d8 ) { }

  T& operator()(const int i, const int j, const int k, const int l, const int m, const int n, const int p, const int q) 
    { 
      assert((0 <= i) && (i < dim1_d)); 
      assert((0 <= j) && (j < dim2_d)); 
      assert((0 <= k) && (k < dim3_d)); 
      assert((0 <= l) && (l < dim4_d));
      assert((0 <= m) && (m < dim5_d));
      assert((0 <= n) && (n < dim6_d));
      assert((0 <= p) && (p < dim7_d));
      assert((0 <= q) && (q < dim8_d));
      int u = i * dim2xdim3xdim4xdim5xdim6xdim7xdim8_d 
            + j * dim3xdim4xdim5xdim6xdim7xdim8_d 
            + k * dim4xdim5xdim6xdim7xdim8_d 
            + l * dim5xdim6xdim7xdim8_d 
            + m * dim6xdim7xdim8_d 
            + n * dim7xdim8_d 
            + p * dim8_d 
            + q;
      return vector<T>::operator[](u);
    }

  T operator()(const int i, const int j, const int k, const int l, const int m, const int n, const int p, const int q) const
    { 
      assert((0 <= i) && (i < dim1_d)); 
      assert((0 <= j) && (j < dim2_d)); 
      assert((0 <= k) && (k < dim3_d)); 
      assert((0 <= l) && (l < dim4_d));
      assert((0 <= m) && (m < dim5_d));
      assert((0 <= n) && (n < dim6_d));
      assert((0 <= p) && (p < dim7_d));
      assert((0 <= q) && (q < dim8_d));
      int u = i * dim2xdim3xdim4xdim5xdim6xdim7xdim8_d 
            + j * dim3xdim4xdim5xdim6xdim7xdim8_d 
            + k * dim4xdim5xdim6xdim7xdim8_d 
            + l * dim5xdim6xdim7xdim8_d 
            + m * dim6xdim7xdim8_d 
            + n * dim7xdim8_d 
            + p * dim8_d 
            + q;
      return vector<T>::operator[](u);
    }

  array_8d<T>& operator+=(const array_8d<T>& C)
  {
    assert(dim1() == C.dim1() &&
           dim2() == C.dim2() &&
           dim3() == C.dim3() &&
           dim4() == C.dim4() &&
           dim5() == C.dim5() &&
           dim6() == C.dim6() &&
           dim7() == C.dim7() &&
           dim8() == C.dim8());
    for (int i=0; i<C.size(); ++i)
      (*this)[i]+=C[i];
    return *this;
  }

  void Clear()
  {
    for (int i=0; i<this->size(); ++i)
      (*this)[i]=0.;
  }

  int get_size()
  {
    return vector<T>::size();
  }

  T& operator() (const vector<int>& indices)
  {
    assert(indices.size() == 8);
    return operator()(indices[0], indices[1], indices[2], indices[3], indices[4], indices[5], indices[6], indices[7]);
  }

  T operator() (const vector<int>& indices) const
  {
    assert(indices.size() == 8);
    return operator()(indices[0], indices[1], indices[2], indices[3], indices[4], indices[5], indices[6], indices[7]);
  }

  void resize (const int i, const int j, const int k, const int l, const int m, const int n, const int p, const int q) 
  { 
    vector<T>::resize (i * j * k * l * m * n * p * q); 
    dim1_d = i; 
    dim2_d = j; 
    dim3_d = k; 
    dim4_d = l; 
    dim5_d = m; 
    dim6_d = n; 
    dim7_d = p; 
    dim8_d = q; 
    dim7xdim8_d = p*q;
    dim6xdim7xdim8_d = n*p*q;
    dim5xdim6xdim7xdim8_d = m*n*p*q;
    dim4xdim5xdim6xdim7xdim8_d = l*m*n*p*q;
    dim3xdim4xdim5xdim6xdim7xdim8_d = k*l*m*n*p*q;
    dim2xdim3xdim4xdim5xdim6xdim7xdim8_d = j*k*l*m*n*p*q;
  }
  string reflect_type() const { return string("array_8d"); }

  int dim1() const { return dim1_d; }
  int dim2() const { return dim2_d; }
  int dim3() const { return dim3_d; }
  int dim4() const { return dim4_d; } 
  int dim5() const { return dim5_d; } 
  int dim6() const { return dim6_d; } 
  int dim7() const { return dim7_d; } 
  int dim8() const { return dim8_d; } 

 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & dim1_d & dim2_d & dim3_d & dim4_d & dim5_d & dim6_d & dim7_d & dim8_d 
       & dim2xdim3xdim4xdim5xdim6xdim7xdim8_d
       & dim3xdim4xdim5xdim6xdim7xdim8_d
       & dim4xdim5xdim6xdim7xdim8_d
       & dim5xdim6xdim7xdim8_d
       & dim6xdim7xdim8_d
       & dim7xdim8_d;
    ar & boost::serialization::base_object<std::vector<T> >(*this);
  }
  int dim1_d;
  int dim2_d;
  int dim3_d;
  int dim4_d;
  int dim5_d;
  int dim6_d;
  int dim7_d;
  int dim8_d;
  int dim2xdim3xdim4xdim5xdim6xdim7xdim8_d;
  int dim3xdim4xdim5xdim6xdim7xdim8_d;
  int dim4xdim5xdim6xdim7xdim8_d;
  int dim5xdim6xdim7xdim8_d;
  int dim6xdim7xdim8_d;
  int dim7xdim8_d;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

#endif
