!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
! ******************************************************************************
!
! Purpose: Enforce bounds on variables.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
! 
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_EnforceBounds.F90,v 1.3 2016/02/08 20:03:59 fred Exp $
!
! Copyright: (c) 2011 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_EnforceBounds(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region 
  USE ModGrid, ONLY: t_grid
  USE ModMPI
  
  USE ModInterfaces, ONLY: RFLU_DecidePrint, &
                           MixtPerf_Eo_DGPUVW, &
                           MixtPerf_G_CpR, &
                           MixtPerf_P_DEoGVm2, &
                           MixtPerf_R_M, &
                           MixtPerf_T_DPR, &
                           RFLU_PrintLocInfo


  USE RFLU_ModJWL

#ifdef SPEC 
   USE ModInterfacesSpecies, ONLY: SPEC_GetSpeciesIndex
   USE ModSpecies,           ONLY: t_spec_input
   USE RFLU_ModConvertCv
#endif  
    
  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Parameters
! ==============================================================================   

  TYPE(t_region), POINTER :: pRegion 

! ==============================================================================
! Local variables
! ==============================================================================
 
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: icg,indCp,indMol,negValCntr
  REAL(RFREAL) :: e,Eo,gamma,ir,negPMin,p,pTol,rgas,r,t,u,v,Vm2,w,a,b,c,Yproducts
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,dv,gv,pCv
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

#ifdef SPEC
  LOGICAL :: scalarConvFlag
  INTEGER :: iCvSpecProducts
  TYPE(t_spec_input), POINTER :: pSpecInput
#endif
  !FRED - 9/10/15


! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_EnforceBounds.F90,v $'
  
  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_EnforceBounds',__FILE__)

! ==============================================================================
! Set pointers and variables
! ==============================================================================

  cv => pRegion%mixt%cv
  dv => pRegion%mixt%dv
  gv => pRegion%mixt%gv

  pCv => pRegion%mixt%cv 
  pGrid => pRegion%grid
  
  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

#ifdef SPEC
  pSpecInput => pRegion%specInput
#endif

  negValCntr = 0

  pTol    = 1.0E-12_RFREAL
  negPMin = 0.0_RFREAL
  
! ******************************************************************************
! Enforce bounds
! ******************************************************************************

  DO icg = 1,pGrid%nCells
    rgas  = MixtPerf_R_M(gv(GV_MIXT_MOL,icg*indMol))
    gamma = MixtPerf_G_CpR(gv(GV_MIXT_CP,icg*indCp),rgas)

    r  = cv(CV_MIXT_DENS,icg)
    ir = 1.0_RFREAL/r
    u  = ir*cv(CV_MIXT_XMOM,icg)
    v  = ir*cv(CV_MIXT_YMOM,icg)
    w  = ir*cv(CV_MIXT_ZMOM,icg)        
    Eo = ir*cv(CV_MIXT_ENER,icg)

    Vm2 = u*u + v*v + w*w
    e   = Eo - 0.5_RFREAL*Vm2

    IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
#ifdef SPEC
     IF (global%specUsed .EQV. .TRUE.) THEN
       IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
         scalarConvFlag = .FALSE.
       ELSE
         scalarConvFlag = .TRUE.

CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
    END IF

     iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')
     Yproducts = pRegion%spec%cv(iCvSpecProducts,icg)
     CALL RFLU_JWL_ComputePressureMixt(pRegion,icg,gamma,rgas,e,r,Yproducts,a,b,c,p,t)

    IF (scalarConvFlag .EQV. .TRUE.) THEN
     CALL RFLU_ScalarConvertCvPrim2Cons(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
    END IF
   END IF
#endif
    ELSE
    p = MixtPerf_P_DEoGVm2(r,Eo,gamma,Vm2)
    t = MixtPerf_T_DPR(r,p,rgas)
    END IF !FRED - added potential JWL Calc here to enforce bound then -9/10/15

    IF ( (p <= 0.0_RFREAL) .AND. (ABS(p) <= pTol) ) THEN
      negValCntr = negValCntr + 1 
      negPMin    = MIN(negPMin,p)        

      p = pTol
! Fix 1 ------------------------------                   
      cv(CV_MIXT_ENER,icg) = r*MixtPerf_Eo_DGPUVW(r,gamma,p,u,v,w)
! Fix 2 ------------------------------                   
!      r = e*(gamma-1.0_RFREAL)/p
!      cv(CV_MIXT_XMOM,icg) = r*ir*cv(CV_MIXT_XMOM,icg)
!      cv(CV_MIXT_YMOM,icg) = r*ir*cv(CV_MIXT_YMOM,icg)
!      cv(CV_MIXT_ZMOM,icg) = r*ir*cv(CV_MIXT_ZMOM,icg)
!      cv(CV_MIXT_ENER,icg) = r*ir*cv(CV_MIXT_ENER,icg)
!      cv(CV_MIXT_DENS,icg) = r
! ------------------------------------      
    END IF ! p
  END DO ! icg

! ******************************************************************************
! Write information
! ******************************************************************************
           
  IF ( (global%verbLevel > VERBOSE_LOW) .AND. &
       (global%myProcid == MASTERPROC) ) THEN 
    IF ( (RFLU_DecidePrint(global) .EQV. .TRUE.) .AND. (negValCntr > 0) ) THEN
      global%warnCounter = global%warnCounter + 1
     
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
            '*** WARNING *** Detected negative pressure!'
            
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal 
      IF ( global%flowType == FLOW_UNSTEADY ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                            global%currentTime
      END IF ! global%flowType
      
      WRITE(STDOUT,'(A,3X,A,1X,I1)') SOLVER_NAME,'Runge-Kutta stage:', &
                                     pRegion%irkStep      

      WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME, & 
            'Number of locations with negative species mass fractions:', &
            negValCntr
      WRITE(STDOUT,'(A,3X,A,1X,E23.16)') SOLVER_NAME, &
            'Largest negative pressure:',negPMin
      WRITE(STDOUT,'(A,3X,A,E23.16)') SOLVER_NAME, &
            'Negative pressure set to:',pTol
    END IF ! RFLU_DecidePrint
  END IF ! global%verbLevel           
           
! ******************************************************************************
! End
! ******************************************************************************  
      
  CALL DeregisterFunction(global)  
    
END SUBROUTINE RFLU_EnforceBounds


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLU_EnforceBounds.F90,v $
!   Revision 1.3  2016/02/08 20:03:59  fred
!   Fixing non-SPEC flag compiling error
!
!   Revision 1.2  2016/02/04 19:59:42  fred
!   Adding iterative JWL EOS capabilities for the cylindrical detonation case
!
!   Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
!   merged rocflu micro and macro
!
!   Revision 1.1.1.1  2014/07/15 14:31:37  brollin
!   New Stable version
!
!
! ******************************************************************************

