#ifndef STACK_ALLOCATOR
#define STACK_ALLOCATOR

#include <iostream>
#include <vector>
#include <stdlib.h>  
#include <memory>

template<class T> class StackAllocator
{
 public:
  static std::size_t size;
  static T* data ;
  static std::size_t memused;

  typedef T* pointer;
  typedef const T* const_pointer;

  typedef T& reference;
  typedef const T& const_reference;

  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;

  typedef T value_type;


  StackAllocator(pointer data_ptr, size_type max_size)  {size =max_size; data=data_ptr;}
  
  StackAllocator(const StackAllocator &other) {}
  StackAllocator() {}

  template<typename U>
    struct rebind {typedef StackAllocator<U> other;};

  void operator=(const StackAllocator& other) {size=other.size; data=other.data;}
  void construct(pointer p, const_reference val) {*p = val;}
  void destroy(pointer p) {*p = 0;}
  pointer allocate(size_type n, const void* hint = 0) 
  {
    if (memused+n >=size) {
      std::cout << "exceeding allowed memory"<<std::endl;
      abort();
    }
    else {
      memused = memused+n;
      return &data[memused-n];
    }
  }
  void deallocate(void* ptr, size_type n) {
    if (n == 0) return;
    if (memused < n || ptr != &data[memused-n]) {
      std::cout << "deallocation not happening in reverse order"<<std::endl;
      abort();
    }
    else {
      memused = memused - n;
    }
  }
  size_type max_size() const {return size;}
  friend std::ostream& operator<<(std::ostream& os, const StackAllocator& c) {
    os<<c.size<<"  "<<c.data<<"  "<<StackAllocator::memused<<std::endl;
    return os;
  }
};


template<class T> 
std::size_t StackAllocator<T>::memused = 0;
template<class T>
T* StackAllocator<T>::data = 0;
template<class T>
std::size_t StackAllocator<T>::size = 0;

#endif
