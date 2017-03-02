/*                                                                           
Developed by Sandeep Sharma and Garnet K.-L. Chan, 2012                      
Copyright (c) 2012, Garnet K.-L. Chan                                        
                                                                             
This program is integrated in Molpro with the permission of 
Sandeep Sharma and Garnet K.-L. Chan
*/

#include "IntegralMatrix.h"
#include "StackBaseOperator.h"
#include "csf.h"
#include "couplingCoeffs.h"
#include "operatorfunctions.h"
#include "operatorloops.h"
#include "distribute.h"
#include "tensor_operator.h"
#include "SpinQuantum.h"
#include "pario.h"
#include "enumerator.h"
//using namespace SpinAdapted::operatorfunctions;

bool SpinAdapted::StackSparseMatrix::nonZeroTensorComponent(Csf& c1, SpinQuantum& opsym, Csf& ladder, int& nonzeroindex, double& cleb)
{
  if (!dmrginp.spinAdapted()) {
    double clebsp = Symmetry::spatial_cg(ladder.sym_is().getirrep(), opsym.get_symm().getirrep(), c1.sym_is().getirrep(), ladder.row(), 0, c1.row());    
    if(c1.n_is() == ladder.n_is() + opsym.get_n() && c1.S.getirrep() == ladder.S.getirrep()+opsym.get_s().getirrep() &&
       fabs(clebsp) >=1.0e-14) {
      nonzeroindex = 0;
      cleb = 1.0;
      return true;
    }
    else
      return false;
  }
  nonzeroindex = 0;
  cleb = 0.0;
  bool found = false;
  int spin = opsym.get_s().getirrep();

  for (int Lz = 0; Lz< Symmetry::sizeofIrrep(opsym.get_symm().getirrep())&&!found; Lz++)
  for (int sz=spin; sz>-spin-1&&!found; sz-=2) 
  {
    cleb = Symmetry::spatial_cg(ladder.sym_is().getirrep(), opsym.get_symm().getirrep(), c1.sym_is().getirrep(), ladder.row(), Lz, c1.row()); 
    cleb *= cg(ladder.S.getirrep(), spin, c1.S.getirrep(), ladder.Sz, sz, c1.Sz);
    if (fabs(cleb) >= 1.0e-14) 
      found = true;
    else
      nonzeroindex++;
  }
  return found;
}

std::vector<double> SpinAdapted::StackSparseMatrix::calcMatrixElements(Csf& c1, TensorOp& Top, Csf& c2, vector<bool>& backupSlater1, vector<bool>& backupSlater2)
{
  vector< vector<double> >& szops = Top.Szops; 
  vector<int>& optypes = Top.optypes;
  std::vector<double> elements(szops.size(), 0.0);
  
  int slatersize = backupSlater1.size();
  for (int isz=0; isz<szops.size(); isz++) 
  for (map<Slater, double>::iterator it1 = c1.det_rep.begin(); it1!= c1.det_rep.end(); it1++) 
  for (map<Slater, double>::iterator it2 = c2.det_rep.begin(); it2!= c2.det_rep.end(); it2++) {
    vector<double>& sz = szops[isz];
    
    for (int j=0; j<sz.size(); j++) {

      if (fabs(sz[j]) < 1.0e-14)
	continue;
      vector<int>& Ind = Top.opindices[j]; //Indices of c and d operators
      
      Slater &s1 = const_cast<Slater&>((*it1).first); 
      Slater &s2 = const_cast<Slater&>((*it2).first);
      const double &d1 = (*it1).second, &d2 = (*it2).second;
      
      int sign1 = s1.getSign(), sign2 = s2.getSign();
      
      for (int i=0; i<slatersize; i++) {
	//backupSlater1[i] = s1.get_orbstring().get_occ_rep()[i];
	backupSlater2[i] = s2.get_orbstring().get_occ_rep()[i];
      }
      
      for (int k=Ind.size()-1; k>=0; k--) { //go from high to low
	if (optypes[k] == 1) s2.c(Ind[k]);
	else s2.d(Ind[k]);
	if (s2.isempty())
	  break;
      }
      elements[isz] += sz[j]*d1*d2*(s1.trace(s2));
      
      for (int i=0; i<slatersize; i++) {
	//s1.get_orbstring().get_occ_rep()[i] = backupSlater1[i];
	s2.get_orbstring().get_occ_rep()[i] = backupSlater2[i];
      }
      s1.setSign(sign1); s2.setSign(sign2);
      s1.setempty(false); s2.setempty(false);
    }
  }
  
  return elements;
}

double SpinAdapted::StackSparseMatrix::calcMatrixElements(Csf& c1, TensorOp& Top, Csf& c2, vector<bool>& backupSlater1, vector<bool>& backupSlater2, int isz)
{
  vector< vector<double> >& szops = Top.Szops; 
  vector<int>& optypes = Top.optypes;
  double element = 0.0;
  
  int slatersize = backupSlater1.size();

  for (map<Slater, double>::iterator it1 = c1.det_rep.begin(); it1!= c1.det_rep.end(); it1++) 
  for (map<Slater, double>::iterator it2 = c2.det_rep.begin(); it2!= c2.det_rep.end(); it2++){
    vector<double>& sz = szops[isz];
    
    for (int j=0; j<sz.size(); j++) {

      if (fabs(sz[j]) < 1.0e-14)
	continue;
      vector<int>& Ind = Top.opindices[j]; //Indices of c and d operators
      
      Slater &s1 = const_cast<Slater&>((*it1).first); 
      Slater &s2 = const_cast<Slater&>((*it2).first);
      const double &d1 = (*it1).second, &d2 = (*it2).second;
      
      int sign1 = s1.getSign(), sign2 = s2.getSign();
      
      for (int i=0; i<slatersize; i++)
	backupSlater2[i] = s2.get_orbstring().get_occ_rep()[i];
      
      for (int k=Ind.size()-1; k>=0; k--) { //go from high to low
	if (optypes[k] == 1) s2.c(Ind[k]);
	else s2.d(Ind[k]);
	if (s2.isempty())
	  break;
      }
      element += sz[j]*d1*d2*(s1.trace(s2));
      
      for (int i=0; i<slatersize; i++) 
	s2.get_orbstring().get_occ_rep()[i] = backupSlater2[i];
      
      s1.setSign(sign1); s2.setSign(sign2);
      s1.setempty(false); s2.setempty(false);
    }
  }
  
  return element;
}

// s*s -> 0
double SpinAdapted::StackSparseMatrix::calcCompfactor(int i, int j, int k, int l, int spin, CompType comp, const TwoElectronArray& v_2, int integralIndex)
{    
  if (!dmrginp.spinAdapted()) {
    if (comp == CD) 
      return 0.5*(-v_2(i, k, l, j) - v_2(k, i, j, l) 
		     + v_2(k, i, l, j) + v_2(i, k, j, l));
    else if (comp == DD)
      return 0.5*(v_2(i, j, l, k));
  }
  double cleb = clebsch(spin, spin, spin, -spin, 0, 0);
  double factor = 0.0;
  if (fabs(cleb) <= 1.0e-14)
    return 0.0;
  if (comp == CD && spin==0) {
    int Ind10 = dmrginp.spatial_to_spin()[i], Ind11 = dmrginp.spatial_to_spin()[j],  Ind20 = dmrginp.spatial_to_spin()[k],  Ind21 = dmrginp.spatial_to_spin()[l];
    factor += 0.5*(-v_2(Ind10, Ind20, Ind21, Ind11) - v_2(Ind20, Ind10, Ind11, Ind21) 
		   + v_2(Ind20, Ind10, Ind21, Ind11) + v_2(Ind10, Ind20, Ind11, Ind21))/cleb/2.;
    Ind10 = dmrginp.spatial_to_spin()[i], Ind11 = dmrginp.spatial_to_spin()[j],  Ind20 = dmrginp.spatial_to_spin()[k]+1,  Ind21 = dmrginp.spatial_to_spin()[l]+1;
    factor += 0.5*(-v_2(Ind10, Ind20, Ind21, Ind11) - v_2(Ind20, Ind10, Ind11, Ind21) 
		   + v_2(Ind20, Ind10, Ind21, Ind11) + v_2(Ind10, Ind20, Ind11, Ind21))/cleb/2.;
    Ind10 = dmrginp.spatial_to_spin()[i]+1, Ind11 = dmrginp.spatial_to_spin()[j]+1,  Ind20 = dmrginp.spatial_to_spin()[k],  Ind21 = dmrginp.spatial_to_spin()[l];
    factor += 0.5*(-v_2(Ind10, Ind20, Ind21, Ind11) - v_2(Ind20, Ind10, Ind11, Ind21) 
		   + v_2(Ind20, Ind10, Ind21, Ind11) + v_2(Ind10, Ind20, Ind11, Ind21))/cleb/2.;
    Ind10 = dmrginp.spatial_to_spin()[i]+1, Ind11 = dmrginp.spatial_to_spin()[j]+1,  Ind20 = dmrginp.spatial_to_spin()[k]+1,  Ind21 = dmrginp.spatial_to_spin()[l]+1;
    factor += 0.5*(-v_2(Ind10, Ind20, Ind21, Ind11) - v_2(Ind20, Ind10, Ind11, Ind21) 
		   + v_2(Ind20, Ind10, Ind21, Ind11) + v_2(Ind10, Ind20, Ind11, Ind21))/cleb/2.;
  }
  else if (comp == CD && spin==2) {
    int Ind10 = dmrginp.spatial_to_spin()[i], Ind11 = dmrginp.spatial_to_spin()[j]+1,  Ind20 = dmrginp.spatial_to_spin()[k]+1,  Ind21 = dmrginp.spatial_to_spin()[l];
    factor += -0.5*(-v_2(Ind10, Ind20, Ind21, Ind11) - v_2(Ind20, Ind10, Ind11, Ind21) 
		   + v_2(Ind20, Ind10, Ind21, Ind11) + v_2(Ind10, Ind20, Ind11, Ind21))/cleb;
  }
  else if (comp == DD && spin==0) {
    int Ind10 = dmrginp.spatial_to_spin()[i], Ind11 = dmrginp.spatial_to_spin()[j]+1,  Ind20 = dmrginp.spatial_to_spin()[k],  Ind21 = dmrginp.spatial_to_spin()[l]+1;
    factor += 0.5*(v_2(Ind10, Ind11, Ind21, Ind20))/cleb/2.;

    Ind10 = dmrginp.spatial_to_spin()[i], Ind11 = dmrginp.spatial_to_spin()[j]+1,  Ind20 = dmrginp.spatial_to_spin()[k]+1,  Ind21 = dmrginp.spatial_to_spin()[l];
    factor += -0.5*(v_2(Ind10, Ind11, Ind21, Ind20))/cleb/2.;

    Ind10 = dmrginp.spatial_to_spin()[i]+1, Ind11 = dmrginp.spatial_to_spin()[j],  Ind20 = dmrginp.spatial_to_spin()[k],  Ind21 = dmrginp.spatial_to_spin()[l]+1;
    factor += -0.5*(v_2(Ind10, Ind11, Ind21, Ind20))/cleb/2.;

    Ind10 = dmrginp.spatial_to_spin()[i]+1, Ind11 = dmrginp.spatial_to_spin()[j],  Ind20 = dmrginp.spatial_to_spin()[k]+1,  Ind21 = dmrginp.spatial_to_spin()[l];
    factor += 0.5*(v_2(Ind10, Ind11, Ind21, Ind20))/cleb/2.;
  }
  else if (comp == DD && spin==2) {
    int Ind10 = dmrginp.spatial_to_spin()[i], Ind11 = dmrginp.spatial_to_spin()[j],  Ind20 = dmrginp.spatial_to_spin()[k],  Ind21 = dmrginp.spatial_to_spin()[l];
    //cout <<Ind10<<"  "<<Ind11<<"  "<<Ind20<<"  "<<Ind21<<endl;
    factor = 0.5*(v_2(Ind10, Ind11, Ind21, Ind20))/cleb;
  }
  return factor;
}


double SpinAdapted::StackSparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, const TwoElectronArray& v_2, int integralIndex)
{
  double factor = 0.0;
  vector<double>& iSz1 = op1.Szops[0];
  bool found = false;
  for (int ilz2=0; ilz2 <op2.rows; ilz2++) 
  for (int sz2=-op2.Spin; sz2< (dmrginp.spinAdapted() ? op2.Spin+1 : -op2.Spin+1); sz2+=2) {
    if (found) break;
    
    int ilz1 = 0;

    int sz2index = dmrginp.spinAdapted() ? ilz2*(op2.Spin+1)+(-sz2+op2.Spin)/2 : 0;
    std::vector<double>&  iSz2 = op2.Szops[sz2index];
    
    double cleb = clebsch(op1.Spin, op1.Spin, op2.Spin, sz2, 0, 0);

    cleb *= Symmetry::spatial_cg(op1.irrep, op2.irrep, 0, ilz1, ilz2, 0);
    if (fabs(cleb) <= 1.0e-14)
      continue;
    else 
      found = true;
    int i1, i2;

    for (i1=0; i1<iSz1.size(); i1++)
      for (i2 =0; i2<iSz2.size(); i2++) {
	vector<int>& Ind1 = op1.opindices[i1], &Ind2 = op2.opindices[i2]; 
	if (comp == CD) {
	  factor += 0.5*(-v_2(Ind1[0], Ind2[0], Ind2[1], Ind1[1]) - v_2(Ind2[0], Ind1[0], Ind1[1], Ind2[1]) 
	    		 + v_2(Ind2[0], Ind1[0], Ind2[1], Ind1[1]) + v_2(Ind1[0], Ind2[0], Ind1[1], Ind2[1]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
	else if (comp == DD) {
	  factor += 0.5*v_2(Ind1[0], Ind1[1], Ind2[1], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	}
	else if (comp == CCD) {
	  factor += 0.5*(v_2(Ind1[0], Ind1[1], Ind2[0], Ind1[2]) - v_2(Ind1[1], Ind1[0], Ind2[0], Ind1[2]) )*iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
	else if (comp == CDD) {
	  factor += 0.5*(v_2(Ind2[0], Ind1[0], Ind1[2], Ind1[1]) - v_2(Ind1[0], Ind2[0], Ind1[2], Ind1[1]) )*iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
	else if (comp == C) {
	  if (op1.dn() + op2.dn() == 0) {
	    factor += 0.5*v_1[integralIndex](Ind1[0], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
	  } 
	  else { // D
	    factor += 0.5*(v_cc[integralIndex](Ind2[0], Ind1[0]) - v_cc[integralIndex](Ind1[0], Ind2[0]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
	  }
	}
      }
  }
  return factor;
}

double SpinAdapted::StackSparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, int op2index, const TwoElectronArray& v_2, int integralIndex)
{
  if(!dmrginp.spinAdapted())
    return calcCompfactor(op1, op2, comp, v_2, integralIndex);
  double factor = 0.0;
  vector<double>& iSz2 = op2.Szops[op2index];
  bool found = false;
  for (int ilz1=0; ilz1 <op1.rows; ilz1++)	
  for (int sz1=-op1.Spin; sz1< op1.Spin+1; sz1+=2) {
    if (found) break;
    
    int ilz2 = op2index/(op2.Spin+1);
    int sz2index = (op2index - ilz2*(op2.Spin+1)), sz2 = op2.Spin - 2*sz2index;
    std::vector<double>&  iSz1 = op1.Szops[ilz1*(op1.Spin+1)+(-sz1+op1.Spin)/2];
    
    //double cleb = cleb_(op1.Spin, sz1, op2.Spin, sz2, 0, 0);
    double cleb = clebsch(op1.Spin, sz1, op2.Spin, sz2, 0, 0);
    cleb *= Symmetry::spatial_cg(op1.irrep, op2.irrep, 0, ilz1, ilz2, 0);
    if (fabs(cleb) <= 1.0e-14)
      continue;
    else 
      found = true;
    int i1, i2;

    for (i1=0; i1<iSz1.size(); ++i1)
      for (i2 =0; i2<iSz2.size(); ++i2) {
	vector<int>& Ind1 = op1.opindices[i1], &Ind2 = op2.opindices[i2]; 
	if (comp == CD) {
	  factor += (-v_2(Ind1[0], Ind2[0], Ind2[1], Ind1[1]) + 
	    	     v_2(Ind2[0], Ind1[0], Ind2[1], Ind1[1]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
	else if (comp == DD) {
	  factor += 0.5*(v_2(Ind1[0], Ind1[1], Ind2[1], Ind2[0]))*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	}
	else if (comp == CCD) {
	  factor += v_2(Ind1[0], Ind1[1], Ind2[0], Ind1[2])*iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
	else if (comp == CDD) {
	  factor += 0.5*(v_2(Ind2[0], Ind1[0], Ind1[2], Ind1[1]) - v_2(Ind1[0], Ind2[0], Ind1[2], Ind1[1]) )*iSz1.at(i1)*iSz2.at(i2)/cleb;
	  //factor += 0.5*(v_2(Ind2[0], Ind1[0], Ind1[1], Ind1[2]) - v_2(Ind1[0], Ind2[0], Ind1[1], Ind1[2]) )*iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
	else if (comp == C) {
	  if (op1.dn() + op2.dn() == 0) {
	    factor += 0.5*v_1[integralIndex](Ind1[0], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
	  } 
	  else { // D
	    factor += 0.5*(v_cc[integralIndex](Ind2[0], Ind1[0]) - v_cc[integralIndex](Ind1[0], Ind2[0]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
	  }
	}
	else {
          abort();
        }
      }
  }
  return factor;
}

double SpinAdapted::StackSparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, const CCCDArray& vcccd) {
  double factor = 0.0;
  vector<double>& iSz1 = op1.Szops[0];
  bool found = false;
  for (int ilz2=0; ilz2 <op2.rows; ilz2++) 
    for (int sz2=-op2.Spin; sz2< (dmrginp.spinAdapted() ? op2.Spin+1 : -op2.Spin+1); sz2+=2) {
      if (found) break;
      
      int ilz1 = 0;
      int sz2index = dmrginp.spinAdapted() ? ilz2*(op2.Spin+1)+(-sz2+op2.Spin)/2 : 0;
      std::vector<double>&  iSz2 = op2.Szops[sz2index];
      
      double cleb = clebsch(op1.Spin, op1.Spin, op2.Spin, sz2, 0, 0);

      cleb *= Symmetry::spatial_cg(op1.irrep, op2.irrep, 0, ilz1, ilz2, 0);
      if (fabs(cleb) <= 1.0e-14)
        continue;
      else 
        found = true;
      int i1, i2;

      for (i1=0; i1<iSz1.size(); ++i1)
        for (i2 =0; i2<iSz2.size(); ++i2) {
          vector<int>& Ind1 = op1.opindices[i1], &Ind2 = op2.opindices[i2]; 
	  
          if (comp == CD) {
            if (op2.dn() == 2) {
              factor += 0.5*vcccd(Ind1[0],Ind2[0],Ind2[1],Ind1[1])*iSz1.at(i1)*iSz2.at(i2)/cleb;
            } else {
              factor += 0.5*vcccd(Ind2[1],Ind2[0],Ind1[1],Ind1[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
            }
          } else if (comp == DD) {
            if (op1.dn() == 2) {
              factor += 0.5*vcccd(Ind1[0],Ind1[1],Ind2[0],Ind2[1])*iSz1.at(i1)*iSz2.at(i2)/cleb;
            } else { // CC_comp
              factor += 0.5*vcccd(Ind1[1],Ind1[0],Ind2[1],Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
            }
          } else if (comp == CCD) { // two cases CDD and CCC
            if (op1.dn() == 3) { // CCC
              factor += (1./6) * vcccd(Ind1[0], Ind1[1], Ind1[2], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
            } else { // CDD
              factor += 0.5 * vcccd(Ind2[0], Ind1[2], Ind1[1], Ind1[0]) *iSz1.at(i1)*iSz2.at(i2)/cleb;
            }
          } else if (comp == CDD) {
            if (op1.dn() == -3) {
              factor += (1./6) * vcccd(Ind1[2], Ind1[1], Ind1[0], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
            } else {
              factor += 0.5 * vcccd(Ind2[0], Ind1[0], Ind1[1], Ind1[2])*iSz1.at(i1)*iSz2.at(i2)/cleb;
            }
          } else {
            abort();
          }
        }
    }
  return factor;
}

double SpinAdapted::StackSparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, int op2index, const CCCDArray& vcccd) {
  if(!dmrginp.spinAdapted())
    return calcCompfactor(op1, op2, comp, vcccd);
  pout << "Sorry, SpinAdapted BCS calculation not implemented" << endl;
  abort();
  return 0.;
}

double SpinAdapted::StackSparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, const CCCCArray& vcccc) {
  double factor = 0.0;
  vector<double>& iSz1 = op1.Szops[0];
  bool found = false;
  for (int ilz2=0; ilz2 <op2.rows; ilz2++) 
    for (int sz2=-op2.Spin; sz2< (dmrginp.spinAdapted() ? op2.Spin+1 : -op2.Spin+1); sz2+=2) {
      if (found) break;
      
      int ilz1 = 0;
      int sz2index = dmrginp.spinAdapted() ? ilz2*(op2.Spin+1)+(-sz2+op2.Spin)/2 : 0;
      std::vector<double>&  iSz2 = op2.Szops[sz2index];
      
      double cleb = clebsch(op1.Spin, op1.Spin, op2.Spin, sz2, 0, 0);

      cleb *= Symmetry::spatial_cg(op1.irrep, op2.irrep, 0, ilz1, ilz2, 0);
      if (fabs(cleb) <= 1.0e-14)
        continue;
      else 
        found = true;
      int i1, i2;

      for (i1=0; i1<iSz1.size(); ++i1)
        for (i2 =0; i2<iSz2.size(); ++i2) {
          vector<int>& Ind1 = op1.opindices[i1], Ind2 = op2.opindices[i2]; 
          if (comp == DD) {
            factor += 0.25 * vcccc(Ind1[0], Ind1[1], Ind2[0], Ind2[1]) *iSz1.at(i1)*iSz2.at(i2)/cleb;
          } else if (comp == CCD) {
            factor += (1./6) * vcccc(Ind2[0], Ind1[2], Ind1[1], Ind1[0]) * iSz1.at(i1)*iSz2.at(i2)/cleb;
          } else if (comp == CDD) {
            factor += (1./6) * vcccc(Ind2[0], Ind1[0], Ind1[1], Ind1[2]) * iSz1.at(i1)*iSz2.at(i2)/cleb;
          } else {
            abort();
          }
        }
    }
  return factor;
}

double SpinAdapted::StackSparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, const PerturbTwoElectronArray& v_2, int integralIndex)
{
  double factor = 0.0;
  vector<double>& iSz1 = op1.Szops[0];
  bool found = false;
  for (int ilz2=0; ilz2 <op2.rows; ilz2++) 
  for (int sz2=-op2.Spin; sz2< (dmrginp.spinAdapted() ? op2.Spin+1 : -op2.Spin+1); sz2+=2) {
    if (found) break;
    
    int ilz1 = 0;

    int sz2index = dmrginp.spinAdapted() ? ilz2*(op2.Spin+1)+(-sz2+op2.Spin)/2 : 0;
    std::vector<double>&  iSz2 = op2.Szops[sz2index];
    
    double cleb = clebsch(op1.Spin, op1.Spin, op2.Spin, sz2, 0, 0);

    cleb *= Symmetry::spatial_cg(op1.irrep, op2.irrep, 0, ilz1, ilz2, 0);
    if (fabs(cleb) <= 1.0e-14)
      continue;
    else 
      found = true;
    int i1, i2;

    for (i1=0; i1<iSz1.size(); i1++)
      for (i2 =0; i2<iSz2.size(); i2++) {
	vector<int>& Ind1 = op1.opindices[i1], Ind2 = op2.opindices[i2]; 
	if (comp == CDD_CD) {
	  factor += (v_2(Ind1[0], Ind2[0], Ind1[1], Ind2[1])-v_2(Ind1[0], Ind2[0], Ind2[1], Ind1[1]) )*iSz1.at(i1)*iSz2.at(i2)/cleb;
	  //factor +=  (v_2(Ind1[0], Ind2[0], Ind1[1], Ind2[1])- v_2(Ind1[0], Ind2[0], Ind2[1], Ind1[1]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
  else if (comp == CCD_CD){
	  factor += (v_2(Ind2[0], Ind1[0], Ind2[1], Ind1[1])-v_2(Ind1[0], Ind2[0], Ind2[1], Ind1[1]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
  }
	else if (comp == DD) {
	  //factor += 0.5*v_2(Ind1[0], Ind1[1], Ind2[1], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	  factor += v_2(Ind1[0], Ind1[1], Ind2[1], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	}
  else if (comp == CDD) {
	  //factor += v_2(Ind1[0], Ind2[0], Ind2[2], Ind2[1])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	  factor += (v_2(Ind1[0], Ind2[0], Ind2[2], Ind2[1]))*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
  }
  else if (comp == CCD) {
	  //factor += v_2(Ind1[0], Ind2[0], Ind2[2], Ind2[1])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	  factor += (v_2(Ind1[0], Ind1[1], Ind2[0], Ind1[2]))*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
  }
	else {
          assert(false);
        }
      }
  }
  return factor;
}


double SpinAdapted::StackSparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, int op2index, const PerturbTwoElectronArray& v_2, int integralIndex)
{
  if(!dmrginp.spinAdapted())
    return calcCompfactor(op1, op2, comp, v_2, integralIndex);
  double factor = 0.0;
  vector<double>& iSz2 = op2.Szops[op2index];
  bool found = false;
  for (int ilz1=0; ilz1 <op1.rows; ilz1++)	
  for (int sz1=-op1.Spin; sz1< op1.Spin+1; sz1+=2) {
    if (found) break;
    
    int ilz2 = op2index/(op2.Spin+1);
    int sz2index = (op2index - ilz2*(op2.Spin+1)), sz2 = op2.Spin - 2*sz2index;
    std::vector<double>&  iSz1 = op1.Szops[ilz1*(op1.Spin+1)+(-sz1+op1.Spin)/2];
    
    //double cleb = cleb_(op1.Spin, sz1, op2.Spin, sz2, 0, 0);
    double cleb = clebsch(op1.Spin, sz1, op2.Spin, sz2, 0, 0);
    cleb *= Symmetry::spatial_cg(op1.irrep, op2.irrep, 0, ilz1, ilz2, 0);
    if (fabs(cleb) <= 1.0e-14)
      continue;
    else 
      found = true;
    int i1, i2;

    for (i1=0; i1<iSz1.size(); i1++)
      for (i2 =0; i2<iSz2.size(); i2++) {
	vector<int>& Ind1 = op1.opindices[i1], Ind2 = op2.opindices[i2]; 
	if (comp == CDD_CD) {
	  //factor += 0.5*(-v_2(Ind1[0], Ind2[0], Ind2[1], Ind1[1])  
	  //  		  + v_2(Ind1[0], Ind2[0], Ind1[1], Ind2[1]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
	  factor +=  (v_2(Ind1[0], Ind2[0], Ind1[1], Ind2[1])-v_2(Ind1[0], Ind2[0], Ind2[1], Ind1[1])) *iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
  else if (comp == CCD_CD){
	  factor += (v_2(Ind2[0], Ind1[0], Ind2[1], Ind1[1])-v_2(Ind1[0], Ind2[0], Ind2[1], Ind1[1]))*iSz1.at(i1)*iSz2.at(i2)/cleb;
  }
	else if (comp == DD) {
	  //factor += 0.5*v_2(Ind1[0], Ind1[1], Ind2[1], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	  factor += v_2(Ind1[0], Ind1[1], Ind2[1], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	}
  else if (comp == CDD) {
	  factor += v_2(Ind1[0], Ind2[0], Ind2[2], Ind2[1])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
  }
  else if (comp == CCD) {
	  //factor += v_2(Ind1[0], Ind2[0], Ind2[2], Ind2[1])*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
	  factor += (v_2(Ind1[0], Ind1[1], Ind2[0], Ind1[2]))*iSz1.at(i1)*iSz2.at(i2)/cleb;	  
  }
	else {
          assert(false);
        }
      }
  }
  return factor;
}

double SpinAdapted::StackSparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, const OneElectronArray& v_1, int integralIndex)
{
  double factor = 0.0;
  vector<double>& iSz1 = op1.Szops[0];
  bool found = false;
  for (int ilz2=0; ilz2 <op2.rows; ilz2++) 
  for (int sz2=-op2.Spin; sz2< (dmrginp.spinAdapted() ? op2.Spin+1 : -op2.Spin+1); sz2+=2) {
    if (found) break;
    
    int ilz1 = 0;

    int sz2index = dmrginp.spinAdapted() ? ilz2*(op2.Spin+1)+(-sz2+op2.Spin)/2 : 0;
    std::vector<double>&  iSz2 = op2.Szops[sz2index];
    
    double cleb = clebsch(op1.Spin, op1.Spin, op2.Spin, sz2, 0, 0);

    cleb *= Symmetry::spatial_cg(op1.irrep, op2.irrep, 0, ilz1, ilz2, 0);
    if (fabs(cleb) <= 1.0e-14)
      continue;
    else 
      found = true;
    int i1, i2;

    for (i1=0; i1<iSz1.size(); i1++)
      for (i2 =0; i2<iSz2.size(); i2++) {
	vector<int>& Ind1 = op1.opindices[i1], Ind2 = op2.opindices[i2]; 
	if (comp == C) {
	  factor += v_1(Ind1[0], Ind2[0])*iSz1.at(i1)*iSz2.at(i2)/cleb;
	}
	else {
          assert(false);
        }
      }
  }
  return factor;
}



double SpinAdapted::StackSparseMatrix::calcCompfactor(TensorOp& op1, TensorOp& op2, CompType comp, int op2index, const CCCCArray& vcccc) {
  if(!dmrginp.spinAdapted())
    return calcCompfactor(op1, op2, comp, vcccc);
  pout << "Sorry, SpinAdapted BCS calculation not implemented" << endl;
  abort();
  return 0.;
}
