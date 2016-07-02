/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#ifndef SPIN_DAVIDSON_HEADER
#define SPIN_DAVIDSON_HEADER

namespace SpinAdapted{
  class StackWavefunction;
  class StackSpinBlock;

struct Davidson_functor
{
  virtual void operator()(StackWavefunction& c, StackWavefunction& v) = 0;
  virtual const StackSpinBlock& get_block() = 0;
  virtual ~Davidson_functor() {};
};

class multiply_h : public Davidson_functor
{
 private:
  const StackSpinBlock& block;
 public:
  multiply_h(const StackSpinBlock& b, const bool &onedot_);
  void operator()(StackWavefunction& c, StackWavefunction& v);
  const StackSpinBlock& get_block() {return block;}
};

class multiply_h_2Index : public Davidson_functor
{
 private:
  const StackSpinBlock& block;
 public:
  multiply_h_2Index(const StackSpinBlock& b, const bool &onedot_);
  void operator()(StackWavefunction& c, StackWavefunction& v);
  const StackSpinBlock& get_block() {return block;}
};

}

#endif
