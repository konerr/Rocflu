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
!******************************************************************************
!
! Purpose: Allocate memory for Lagrangian particle solution.
!
! Description: None.
!
! Input: 
!   pRegion     Region pointer
!   pPlag       Pointer to particle data structure
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: PLAG_RFLU_AllocMemSol.F90,v 1.2 2015/12/18 23:22:29 rahul Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_AllocMemSol(pRegion,pPlag)

  USE ModDataTypes
  USE ModError
  USE ModParameters  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModPartLag, ONLY: t_plag
  USE ModMPI
   
  USE PLAG_ModParameters  
   
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion
  TYPE(t_plag), POINTER :: pPlag

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,icg,iCont,iPcl,iT,iVar,nAiv,nArv,nCont,nCv,nDv, &
             nPclsMax,nTv,nUnsteadyData
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_AllocMemSol.F90,v $ $Revision: 1.2 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_AllocMemSol',__FILE__)

! ******************************************************************************
! Set pointers and variables 
! ******************************************************************************

  pGrid => pRegion%grid

  nPclsMax = pRegion%plag%nPclsMax
  nCont    = pRegion%plagInput%nCont

  nAiv = pPlag%nAiv
  nArv = pPlag%nArv

  nCv  = pPlag%nCv
  nDv  = pPlag%nDv
  nTv  = pPlag%nTv

  nUnsteadyData = pRegion%plagInput%nUnsteadyData

! ******************************************************************************
! Allocate memory
! ******************************************************************************

! ==============================================================================
! State vector
! ==============================================================================

  ALLOCATE(pPlag%cv(nCv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%cv')
  END IF ! global%error 

! rahul- required for vFrac gradient reconstruction, for simple limiter
! precisely
  ALLOCATE(pPlag%plagcvinfo(1),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%plagcvInfo')
  END IF ! global%error 
! end -rahul

! ==============================================================================
! Dependent variables  
! ==============================================================================
    
  ALLOCATE(pPlag%dv(nDv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%dv')
  END IF ! global%error 

  ALLOCATE(pPlag%dvOld(nDv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%dvOld')
  END IF ! global%error 

! ==============================================================================
! Transport variables  
! ==============================================================================
    
  ALLOCATE(pPlag%tv(nTv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%tv')
  END IF ! global%error 

! ==============================================================================
! Additional variables  
! ==============================================================================

  ALLOCATE(pPlag%aiv(nAiv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%aiv')
  END IF ! global%error 

  ALLOCATE(pPlag%arv(nArv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%arv')
  END IF ! global%error

! ==============================================================================
! Lagrangian particles mass indices
! ==============================================================================

  ALLOCATE(pPlag%cvPlagMass(nCont),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global, ERR_ALLOCATE,__LINE__ ,'pPlag%cvPlagMass')
  END IF ! global%error

  DO iCont = 1,nCont
    pPlag%cvPlagMass(iCont) = CV_PLAG_LAST + iCont
  END DO ! iCont

! ==============================================================================
! Variables for unsteady force  
! ==============================================================================

  ALLOCATE(pPlag%timeBH(nUnsteadyData),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%timeBH')
  END IF ! global%error 

  ALLOCATE(pPlag%dImpulseMax(3,pRegion%grid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%dImpulseMax')
  END IF ! global%error 

  ALLOCATE(pPlag%forceTotal(3,pRegion%grid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%forceTotal')
  END IF ! global%error 

  ALLOCATE(pPlag%pgMixt(3,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%pgMixt')
  END IF ! global%error 

  ALLOCATE(pPlag%dudtMixt(3,nUnsteadyData,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%dudtMixt')
  END IF ! global%error 

  ALLOCATE(pPlag%dudtPlag(3,nUnsteadyData,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%dudtPlag')
  END IF ! global%error 

  ALLOCATE(pPlag%coeffIU(pRegion%grid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%coeffIU')
  END IF ! global%error 

! ==============================================================================
! Particle volume fraction and gradient variables  
! ==============================================================================

  ALLOCATE(pPlag%vFracE(1,pRegion%grid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%vFracE')
  END IF ! global%error 

  ALLOCATE(pPlag%vFracL(1,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%vFracL')
  END IF ! global%error 

  ALLOCATE(pPlag%gradVFracE(3,1,pRegion%grid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%gradVFracE')
  END IF ! global%error 

! rahul - pointer for gradient of gas volume fraction
  ALLOCATE(pPlag%gradVFracEg(3,1,pRegion%grid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%gradVFracEg')
  END IF ! global%error 
! end rahul

  ALLOCATE(pPlag%gradVFracL(3,1,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%gradVFracL')
  END IF ! global%error 

! ******************************************************************************
! Initialize memory
! ******************************************************************************

  DO iPcl = 1,nPclsMax
    DO iVar = 1,nCv
      pPlag%cv(iVar,iPcl) = 0.0_RFREAL
    END DO ! iVar
    
    DO iVar = 1,nDv
      pPlag%dv(iVar,iPcl)    = 0.0_RFREAL
      pPlag%dvOld(iVar,iPcl) = 0.0_RFREAL
    END DO ! iVar
   
    DO iVar = 1,nTv
      pPlag%tv(iVar,iPcl) = 0.0_RFREAL
    END DO ! iVar
    
    DO iVar = 1,nAiv
      pPlag%aiv(iVar,iPcl) = 0
    END DO ! iVar
    
    DO iVar = 1,nArv   
      pPlag%arv(iVar,iPcl) = 0
    END DO ! iVar

    DO iVar = 1,3
      pPlag%pgMixt(iVar,iPcl) = 0.0_RFREAL

      DO iT = 1,nUnsteadyData
        pPlag%dudtMixt(iVar,iT,iPcl) = 0.0_RFREAL
        pPlag%dudtPlag(iVar,iT,iPcl) = 0.0_RFREAL
      END DO ! iT
    END DO ! iVar

    DO iVar = 1,1
      pPlag%vFracL(iVar,iPcl) = 0.0_RFREAL
    END DO ! iVar

    DO iVar = 1,1
      pPlag%gradVFracL(1,iVar,iPcl) = 0.0_RFREAL
      pPlag%gradVFracL(2,iVar,iPcl) = 0.0_RFREAL
      pPlag%gradVFracL(3,iVar,iPcl) = 0.0_RFREAL
    END DO ! iVar
  END DO ! iPcl

  DO iT = 1,nUnsteadyData
    pPlag%timeBH(iT) = REAL(ABS(CRAZY_VALUE_INT),KIND=RFREAL)
  END DO ! iT
  pPlag%timeBH(1) = 0.0_RFREAL ! Index 1 is for current time

  DO icg = 1,pRegion%grid%nCellsTot
    pPlag%coeffIU(icg) = 0.0_RFREAL

    DO iVar = 1,1
      pPlag%vFracE(iVar,icg) = 0.0_RFREAL
    END DO ! iVar

    DO iVar = 1,1
      pPlag%gradVFracE(1,iVar,icg) = 0.0_RFREAL
      pPlag%gradVFracE(2,iVar,icg) = 0.0_RFREAL
      pPlag%gradVFracE(3,iVar,icg) = 0.0_RFREAL
    END DO ! iVar

    DO iVar = 1,1
      pPlag%gradVFracEg(1,iVar,icg) = 0.0_RFREAL
      pPlag%gradVFracEg(2,iVar,icg) = 0.0_RFREAL
      pPlag%gradVFracEg(3,iVar,icg) = 0.0_RFREAL
    END DO ! iVar

    DO iVar = 1,3
      pPlag%dImpulseMax(iVar,icg) = 0.0_RFREAL
      pPlag%forceTotal(iVar,icg)  = 0.0_RFREAL
    END DO ! iVar
  END DO ! icg

! rahul AUSM+up  
  DO iVar = 1,1
    pPlag%plagcvInfo(iVar) = 1
  END DO
! end- rahul AUSM+up
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_AllocMemSol

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_AllocMemSol.F90,v $
! Revision 1.2  2015/12/18 23:22:29  rahul
! Allocate memory to gradvFracEg and plagVarInfo.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:51  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/04/26 20:46:05  fnajjar
! Removed dvPlagVolu from allocation and deallocation calls
!
! Revision 1.1  2007/04/09 18:50:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.5  2007/03/06 23:15:32  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.4  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.3  2004/07/28 18:58:13  fnajjar
! Included new definition for nPclsTot from dynamic memory reallocation
!
! Revision 1.2  2004/03/08 23:02:28  fnajjar
! Added initialization section within DO-loop construct
!
! Revision 1.1  2004/02/26 21:00:36  haselbac
! Initial revision
!
!******************************************************************************

