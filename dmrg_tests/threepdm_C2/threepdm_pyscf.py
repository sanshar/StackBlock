#!/usr/bin/env python
#
# Author: Sandeep Sharma <sanshar@gmail.com>
#         Qiming Sun <osirpt.sun@gmail.com>
#

import os, sys
import numpy as np
import pyscf.tools
import pyscf.lib.logger as logger
from pyscf.mrpt.nevpt2 import sc_nevpt
from pyscf.dmrgscf.dmrgci import *
from pyscf import fci


def convert2(size, two, one ):
  #twopdm is the twopdm in block definition
  #twopdm[i,j,k,l] = \sum_{spin} <a^+_ia^+_ja_ka^l>
  #two is the twopdm in spin-free form. 
  #two[i,j,k,l] = <E^i_jE^k_l>
  # a^+_ia^+_ja_ka^l =  E^i_lE^j_k -\delta_{j,l} E^i_k
  twopdm = np.zeros((size,size,size,size))
  for i in xrange(0,size):
    for j in xrange(0,size):
      for k in xrange(0,size):
        for l in xrange(0,size):
          twopdm[i,j,k,l] +=two[i,l,j,k] 
          if(j==l):
            twopdm[i,j,k,l] -=one[i,k] 
  return twopdm


def convert3(size, three, two, one, twopdm): 
  # threepdm[i,j,k,l,m,n] = E^{i,j,k}_{l,m,n}
  # E^{i,j,k}_{l,m,n} = E^{i,j}_{m,n}E^k_l -\delta_{k,m}E^{i,j}_{l,n}- \delta_{k,n}E^{i,j}_{m,l}
  # = E^i_nE^j_mE^k_l -\delta_{j,n}E^i_mE^k_l -\delta_{k,m}E^{i,j}_{l,n} -\delta_{k,n}E^{i,j}_{m,l}
  threepdm = np.zeros((size,size,size,size,size,size))
  for i in xrange(0,size):
    for j in xrange(0,size):
      for k in xrange(0,size):
        for l in xrange(0,size):
          for m in xrange(0,size):
            for n in xrange(0,size):
              threepdm[i,j,k,l,m,n] +=three[i,n,j,m,k,l] 
              if( j== n):
                threepdm[i,j,k,l,m,n] -=two[i,m,k,l] 
              if(k==m):
                threepdm[i,j,k,l,m,n] -=twopdm[i,j,l,n] 
              if(k==n):
                threepdm[i,j,k,l,m,n] -=twopdm[i,j,m,l] 
  return threepdm

def convert4(size, four, three, two, one, threepdm, twopdm):
  # fourpdm[i,j,k,l,m,n,p,q] = E^{i,j,k,l}_{m,n,p,q}
  # E^{i,j,k,l}_{m,n,p,q} = E^{i,j,k}_{n,p,q}E^l_m -\delta_{l,n}E^{i,j,k}_{m,p,q} -\delta_{l,p}E^{i,j,k}_{n,m,p}-\delta_{l,q}E^{i,j,k}_{n,p,m}
  # = E^{i,j,k}_{n,p,q}(E^i_qE^j_pE^k_n -\delta_{j,q}E^i_pE^k_n-\delta_{k,p}E^{i,j}_{n,q}-\delta_{k,q}E^{i,j}_{p,n} )E^l_m -\delta_{l,n}E^{i,j,k}_{m,p,q} -\delta_{l,p}E^{i,j,k}_{n,m,q}-\delta_{l,q}E^{i,j,k}_{n,p,m}
  # = E^i_qE^j_pE^k_nE^l_m -\delta_{j,q}E^i_pE^k_nE^l_m-\delta_{k,p}( E^i_qE^j_n- \delta{j,q}E^i_n)E^l_m-\delta_{k,q}(E^i_nE^j_p -\delta_{j,n}E^i_p)E^l_m -\delta_{l,n}E^{i,j,k}_{m,p,q} -\delta_{l,p}E^{i,j,k}_{n,m,q}-\delta_{l,q}E^{i,j,k}_{n,p,m}
  fourpdm = np.zeros((size,size,size,size,size,size,size,size))
  for i in xrange(0,size):
    for j in xrange(0,size):
      for k in xrange(0,size):
        for l in xrange(0,size):
          for m in xrange(0,size):
            for n in xrange(0,size):
              for p in xrange(0,size):
                for q in xrange(0,size):
                  fourpdm[i,j,k,l,m,n,p,q] += four[i,q,j,p,k,n,l,m]
                  if( j==q):
                    fourpdm[i,j,k,l,m,n,p,q] -= three[i,p,k,n,l,m] 
                  if( k==p):
                    fourpdm[i,j,k,l,m,n,p,q] -= three[i,q,j,n,l,m]
                    if( j==q):
                      fourpdm[i,j,k,l,m,n,p,q] +=two[i,n,l,m] 

                  if( k==q):
                    fourpdm[i,j,k,l,m,n,p,q] -=three[i,n,j,p,l,m] 
                    if (j==n):
                      fourpdm[i,j,k,l,m,n,p,q] += two[i,p,l,m]

                  if( l==n):
                    fourpdm[i,j,k,l,m,n,p,q] -=threepdm[i,j,k,m,p,q]
                  if( l==p):
                    fourpdm[i,j,k,l,m,n,p,q] -=threepdm[i,j,k,n,m,q]
                  if( l==q):
                    fourpdm[i,j,k,l,m,n,p,q] -=threepdm[i,j,k,n,p,m]
  return fourpdm


if __name__ == '__main__':
    from pyscf import gto
    from pyscf import scf
    from pyscf import mcscf
    import pickle
    from pyscf.dmrgscf import settings
    settings.MPIPREFIX =''
    settings.BLOCKEXE = '/tigress/shengg/stackblock/block.spin_adapted'

    b = 1.4
    mol = gto.Mole()
    mol.build(
        verbose = 7,
        output = 'out-casscf',
        atom = [['C', 0.,0., -0.5], ['C', 0.,0., 0.5]],
        basis = {'C': 'sto-3g'},
        #symmetry = 'd2h',
        spin = 0
    )
    m = scf.RHF(mol)
    m.scf()

    mc = mcscf.CASCI(m, 6, 4)

    mc.fcisolver.nroots=2
    mc.casci()
    size = 6
    dm1, dm2, dm3 = fci.rdm.make_dm123('FCI3pdm_kern_sf',mc.ci[0],mc.ci[0], mc.ncas, mc.nelecas)
    twopdm = convert2(6, dm2, dm1)
    threepdm = convert3(6, dm3, dm2, dm1, twopdm)
    f = open("rdm3_0_0", 'w')
    print >>f, size
    for i in xrange(0,size):
      for j in xrange(0,size):
        for k in xrange(0,size):
          for l in xrange(0,size):
            for m in xrange(0,size):
              for n in xrange(0,size):
                print >>f, i,j,k,l,m,n,threepdm[i,j,k,l,m,n]
    f.close()
    dm1, dm2, dm3 = fci.rdm.make_dm123('FCI3pdm_kern_sf',mc.ci[1],mc.ci[0], mc.ncas, mc.nelecas)
    twopdm = convert2(6, dm2, dm1)
    threepdm = convert3(6, dm3, dm2, dm1, twopdm)
    f = open("rdm3_1_0", 'w')
    print >>f, size
    for i in xrange(0,size):
      for j in xrange(0,size):
        for k in xrange(0,size):
          for l in xrange(0,size):
            for m in xrange(0,size):
              for n in xrange(0,size):
                print >>f, i,j,k,l,m,n,threepdm[i,j,k,l,m,n]
    f.close()
    dm1, dm2, dm3 = fci.rdm.make_dm123('FCI3pdm_kern_sf',mc.ci[1],mc.ci[1], mc.ncas, mc.nelecas)
    twopdm = convert2(6, dm2, dm1)
    threepdm = convert3(6, dm3, dm2, dm1, twopdm)
    f = open("rdm3_1_1", 'w')
    print >>f, size
    for i in xrange(0,size):
      for j in xrange(0,size):
        for k in xrange(0,size):
          for l in xrange(0,size):
            for m in xrange(0,size):
              for n in xrange(0,size):
                print >>f, i,j,k,l,m,n,threepdm[i,j,k,l,m,n]
    f.close()


    m = scf.RHF(mol)
    m.scf()
    mc2 = mcscf.CASCI(m, 6, 4)
    mc2.fcisolver = DMRGCI(mol,tol=1e-8)
    mc2.fcisolver.nroots = 2
    mc2.fcisolver.twopdm = False
    mc2.fcisolver.extraline.append("tran_threepdm")
    #mc2.fcisolver.extraline.append("tran_threepdm")
 #   mc2.fcisolver.extraline.append("npdm_no_intermediate")
    #mc2.fcisolver.extraline.append("npdm_no_intermediate")
    mc2.casci()
