#ifndef SPIN_ENUMERATOR_HEADER
#define SPIN_ENUMERATOR_HEADER

namespace SpinAdapted{
  enum opTypes{ HAM, CRE, CRE_CRE, DES_DESCOMP, CRE_DES, CRE_DESCOMP, CRE_CRE_DESCOMP, 
		DES, DES_DES, CRE_CRECOMP, DES_CRE, DES_CRECOMP, CRE_DES_DESCOMP, OVERLAP,

    //mps_sc_nevpt
    CDD_SUM, CDD_CRE_DESCOMP, CDD_DES_DESCOMP, CCD_SUM, CCD_CRE_DESCOMP, CCD_CRE_CRECOMP,

		// Extra (empty) operator classes for RI-approx NPDM                                                                                                                 
		RI_3INDEX, RI_4INDEX,

		// Extra 3PDM operators                                                                                                                                              
		CRE_CRE_CRE, CRE_CRE_DES, CRE_DES_DES, CRE_DES_CRE,
		// Extra 4PDM operators                                                                                                                                              
		DES_CRE_DES, DES_DES_CRE, DES_CRE_CRE, DES_DES_DES,
		CRE_CRE_DES_DES, CRE_DES_CRE_DES, CRE_DES_DES_CRE, CRE_DES_DES_DES,
		CRE_CRE_CRE_DES, CRE_CRE_DES_CRE, CRE_DES_CRE_CRE, CRE_CRE_CRE_CRE };

  enum CompType{CD, DD, CCD, C, CDD, CDD_CD, CCD_CD};
  enum sweepType {FULL, PARTIAL};
  enum guessWaveTypes {BASIC, TRANSFORM, TRANSPOSE};
  enum Storagetype {LOCAL_STORAGE, DISTRIBUTED_STORAGE, DISTRIBUTED_STORAGE_FOR_ONEPDM};
  enum WarmUpTypes {WILSON, LOCAL0, LOCAL2, LOCAL3, LOCAL4};
  enum hamTypes {QUANTUM_CHEMISTRY, HUBBARD, BCS, HEISENBERG};
  enum solveTypes {LANCZOS, DAVIDSON, CONJUGATE_GRADIENT};
  enum algorithmTypes {ONEDOT, TWODOT, TWODOT_TO_ONEDOT, PARTIAL_SWEEP};
  enum noiseTypes {RANDOM, EXCITEDSTATE};
  enum calcType {DMRG, ONEPDM, TWOPDM, THREEPDM, FOURPDM, NEVPT2PDM, RESTART_TWOPDM,
		 RESTART_ONEPDM, RESTART_THREEPDM, RESTART_FOURPDM, RESTART_NEVPT2PDM, TINYCALC, FCI,
		 EXCITEDDMRG, CALCOVERLAP, CALCHAMILTONIAN, COMPRESS, RESPONSE, RESPONSEBW,
		 TRANSITION_ONEPDM, TRANSITION_TWOPDM, TRANSITION_THREEPDM, RESTART_T_ONEPDM, RESTART_T_TWOPDM, RESTART_T_THREEPDM,
		 NEVPT2,RESTART_NEVPT2, RESPONSELCC, RESPONSEAAAV, RESPONSEAAAC, MPS_NEVPT, RESTART_MPS_NEVPT};
  enum orbitalFormat{MOLPROFORM, DMRGFORM};
  enum reorderType{FIEDLER, GAOPT, MANUAL, NOREORDER};
  enum keywords{ORBS, LASTM, STARTM, MAXM,  REORDER, HF_OCC, SCHEDULE, SYM, NELECS, SPIN, IRREP,
		MAXJ, PREFIX, NROOTS, DOCD, DEFLATION_MAX_SIZE, MAXITER, BASENERGY,
		SCREEN_TOL, ODOT, SWEEP_TOL, OUTPUTLEVEL, NONSPINADAPTED, BOGOLIUBOV, TWODOT_NOISE, WARMUP, 
		NPDM_INTERMEDIATE, NPDM_NO_INTERMEDIATE, NPDM_MULTINODE, NPDM_NO_MULTINODE, NUMKEYWORDS};
  
  enum {
    NO_PARTICLE_SPIN_NUMBER_CONSTRAINT,
    PARTICLE_SPIN_NUMBER_CONSTRAINT,
    SPIN_NUMBER_CONSTRAINT,
    PARTICLE_NUMBER_CONSTRAINT,
    HOLE_NUMBER_CONSTRAINT
  };
  
  enum {
    AnyQ,  
    LessThanQ, 
    EqualQ,
    EqualS,
    LessThanN,
    LessThanH,
    WITH_LIST /*< states are added together if the are allowed by the quantaList */
  };

  enum OnePerturbType{ Va_1=0, Vi_1, Vai_1, OnePertEnd};
  enum TwoPerturbType{Va=0,Vi,Vab,Vij,Vai,Vabi,Vaij,Vabij,TwoPertEnd};


};

#endif
