#ifndef STACK_ALLOCATOR
#define STACK_ALLOCATOR

#include <iostream>
#include <vector>
#include <stdlib.h>  
#include <memory>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

void print_trace(int sig);

template<class T> class StackAllocator
{
 public:
  std::size_t size;
  T* data ;
  std::size_t memused;


 StackAllocator(T* data_ptr, std::size_t max_size): memused(0)  {size =max_size; data=data_ptr;}
  
 StackAllocator() : size(0), data(0), memused(0) {}
  void clear() {size = 0;data=0; memused=0;}
  T* allocate(std::size_t n, const void* hint = 0) 
  {
    if (memused+n >=size)
      {
	std::cout << "exceeding allowed memory"<<std::endl;
	print_trace(11);
	return 0;
      }
    else
      {
	memused = memused+n;
	return &data[memused-n];
      }
  }
  void deallocate(void* ptr, std::size_t n) {
    if (n == 0) return;
    if (memused < n || ptr != &data[memused-n]) {
      std::cout << "deallocation not happening in reverse order"<<std::endl;
      print_trace(11);
    }
    else {
      memused = memused - n;
    }
  }
  std::size_t max_size() const {return size;}
  friend std::ostream& operator<<(std::ostream& os, const StackAllocator& c) {
    os<<c.size<<"  "<<c.data<<"  "<<c.memused<<std::endl;
    return os;
  }
};

#endif
