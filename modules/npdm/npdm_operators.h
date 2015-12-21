#ifndef NPDM_OPERATORS_H
#define NPDM_OPERATORS_H

#include "npdm_spin_ops.h"

namespace SpinAdapted{
namespace Npdm{

//FIXME constructors / destructors??

//===========================================================================================================================================================
//  4-INDEX compound Ops (built using RI approximation, exact on dot block)
//===========================================================================================================================================================

class Npdm_op_compound_CCDD : public NpdmSpinOps {
  public:
    Npdm_op_compound_CCDD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_CCDD>(new Npdm_op_compound_CCDD(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CCCD : public NpdmSpinOps {
  public:
    Npdm_op_compound_CCCD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_CCCD>(new Npdm_op_compound_CCCD(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CCDC : public NpdmSpinOps {
  public:
    Npdm_op_compound_CCDC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_CCDC>(new Npdm_op_compound_CCDC(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CDCC : public NpdmSpinOps {
  public:
    Npdm_op_compound_CDCC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_CDCC>(new Npdm_op_compound_CDCC(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CDCD : public NpdmSpinOps {
  public:
    Npdm_op_compound_CDCD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_CDCD>(new Npdm_op_compound_CDCD(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CDDC : public NpdmSpinOps {
  public:
    Npdm_op_compound_CDDC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_CDDC>(new Npdm_op_compound_CDDC(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CDDD : public NpdmSpinOps {
  public:
    Npdm_op_compound_CDDD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_CDDD>(new Npdm_op_compound_CDDD(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CCCC : public NpdmSpinOps {
  public:
    Npdm_op_compound_CCCC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_CCCC>(new Npdm_op_compound_CCCC(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_DCCD : public NpdmSpinOps {
  public:
    Npdm_op_compound_DCCD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_DCCD>(new Npdm_op_compound_DCCD(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_DCDC : public NpdmSpinOps {
  public:
    Npdm_op_compound_DCDC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_DCDC>(new Npdm_op_compound_DCDC(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_DDCC : public NpdmSpinOps {
  public:
    Npdm_op_compound_DDCC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_DDCC>(new Npdm_op_compound_DDCC(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_4INDEX).get_local_indices(); }
};

//===========================================================================================================================================================
//  4-INDEX Ops
//===========================================================================================================================================================

class Npdm_op_wrapper_CCDD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CCDD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_CRE_DES_DES).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_CCDD>(new Npdm_op_wrapper_CCDD(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_CRE_DES_DES).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CDCD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CDCD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_DES_CRE_DES).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_CDCD>(new Npdm_op_wrapper_CDCD(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_DES_CRE_DES).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CDDC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CDDC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_DES_DES_CRE).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_CDDC>(new Npdm_op_wrapper_CDDC(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_DES_DES_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CDDD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CDDD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_DES_DES_DES).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_CDDD>(new Npdm_op_wrapper_CDDD(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_DES_DES_DES).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CCCD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CCCD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_CRE_CRE_DES).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_CCCD>(new Npdm_op_wrapper_CCCD(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_CRE_CRE_DES).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CCDC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CCDC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_CRE_DES_CRE).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_CCDC>(new Npdm_op_wrapper_CCDC(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_CRE_DES_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CDCC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CDCC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_DES_CRE_CRE).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_CDCC>(new Npdm_op_wrapper_CDCC(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_DES_CRE_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CCCC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CCCC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_CRE_CRE_CRE).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_CCCC>(new Npdm_op_wrapper_CCCC(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_CRE_CRE_CRE).get_array(); }
};

//===========================================================================================================================================================
//  3-INDEX compound Ops (built using RI approximation, exact on dot block)
//===========================================================================================================================================================

class Npdm_op_compound_CCD : public NpdmSpinOps {
  public:
    Npdm_op_compound_CCD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_CCD>(new Npdm_op_compound_CCD(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_3INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CDD : public NpdmSpinOps {
  public:
    Npdm_op_compound_CDD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_CDD>(new Npdm_op_compound_CDD(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_3INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CDC : public NpdmSpinOps {
  public:
    Npdm_op_compound_CDC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_CDC>(new Npdm_op_compound_CDC(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_3INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_CCC : public NpdmSpinOps {
  public:
    Npdm_op_compound_CCC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_CCC>(new Npdm_op_compound_CCC(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_3INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_DCD : public NpdmSpinOps {
  public:
    Npdm_op_compound_DCD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_DCD>(new Npdm_op_compound_DCD(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_3INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_DDC : public NpdmSpinOps {
  public:
    Npdm_op_compound_DDC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_DDC>(new Npdm_op_compound_DDC(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_3INDEX).get_local_indices(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_compound_DCC : public NpdmSpinOps {
  public:
    Npdm_op_compound_DCC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_compound_DCC>(new Npdm_op_compound_DCC(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(RI_3INDEX).get_local_indices(); }
};

//===========================================================================================================================================================
//  3-INDEX Ops
//===========================================================================================================================================================

class Npdm_op_wrapper_CCC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CCC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_CRE_CRE).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_CCC>(new Npdm_op_wrapper_CCC(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_CRE_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CCD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CCD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_CCD>(new Npdm_op_wrapper_CCD(*this));}
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_CRE_DES).get_local_indices(); }
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_CRE_DES).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CDD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CDD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_DES_DES).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_CDD>(new Npdm_op_wrapper_CDD(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_DES_DES).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CDC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CDC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_DES_CRE).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_CDC>(new Npdm_op_wrapper_CDC(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_DES_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DCD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DCD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(DES_CRE_DES).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_DCD>(new Npdm_op_wrapper_DCD(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(DES_CRE_DES).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DDC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DDC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(DES_DES_CRE).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_DDC>(new Npdm_op_wrapper_DDC(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(DES_DES_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DCC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DCC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(DES_DES_CRE).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_DCC>(new Npdm_op_wrapper_DCC(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(DES_DES_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DDD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DDD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_CRE_CRE).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_DDD>(new Npdm_op_wrapper_DDD(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_CRE_CRE).get_array(); }
};

//===========================================================================================================================================================
//  2-INDEX Ops
//===========================================================================================================================================================

class Npdm_op_wrapper_CC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_CRE).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_CC>(new Npdm_op_wrapper_CC(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_CD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_CD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE_DES).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_CD>(new Npdm_op_wrapper_CD(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE_DES).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DC : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DC( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(DES_CRE).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_DC>(new Npdm_op_wrapper_DC(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(DES_CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_DD : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_DD( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(DES_DES).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_DD>(new Npdm_op_wrapper_DD(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(DES_DES).get_array(); }
};

//===========================================================================================================================================================
//  1-INDEX Ops
//===========================================================================================================================================================

class Npdm_op_wrapper_C : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_C( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE).get_local_indices(); }
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_C>(new Npdm_op_wrapper_C(*this));}
    std::vector< std::vector<int> > get_indices() { return spinBlock_->get_op_array(CRE).get_array(); }
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------

class Npdm_op_wrapper_D : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_D( StackSpinBlock * spinBlock );
    bool set_local_ops( int idx );
//    const std::vector< int >& get_1d_indices() { return spinBlock_->get_op_array(CRE).get_local_indices(); }
//    FIXME

    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_D>(new Npdm_op_wrapper_D(*this));}
    std::vector< std::vector<int> > get_indices() { 
      if(dmrginp.doimplicitTranspose())
	return spinBlock_->get_op_array(CRE).get_array(); 
      else
	return spinBlock_->get_op_array(DES).get_array(); 
    }
};

//===========================================================================================================================================================
//  NULL case (for empty creation-destruction patterns)
//===========================================================================================================================================================

class Npdm_op_wrapper_NULL : public NpdmSpinOps {
  public:
    Npdm_op_wrapper_NULL(StackSpinBlock * spinBlock);
    virtual boost::shared_ptr<NpdmSpinOps> getcopy() {return boost::shared_ptr<Npdm_op_wrapper_NULL>(new Npdm_op_wrapper_NULL(*this));}
    bool set_local_ops( int idx );
};

//===========================================================================================================================================================

}
}

#endif

