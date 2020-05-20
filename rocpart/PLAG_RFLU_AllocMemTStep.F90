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
! Purpose: Allocate memory for Lagrangian particles related to time stepping.
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
! $Id: PLAG_RFLU_AllocMemTStep.F90,v 1.2 2015/12/18 23:28:36 rahul Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_AllocMemTStep(pRegion,pPlag)

  USE ModDataTypes
  USE ModError
  USE ModParameters  
  USE ModGlobal, ONLY: t_global 
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModPartLag, ONLY: t_plag  
  USE ModMPI
   
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
  INTEGER :: errorFlag,iPcl,iVar,nAiv,nArv,nCv,nPclsMax
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_AllocMemTStep.F90,v $ $Revision: 1.2 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_AllocMemTStep',__FILE__)

! ******************************************************************************
! Set pointers and variables 
! ******************************************************************************

  pGrid  => pRegion%grid

  nPclsMax = pRegion%plag%nPclsMax

  nAiv = pPlag%nAiv
  nArv = pPlag%nArv

  nCv  = pPlag%nCv

! ******************************************************************************
! Allocate memory
! ******************************************************************************

! ==============================================================================
! Old state vector
! ==============================================================================

  ALLOCATE(pPlag%cvOld(nCv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%cvOld')
  END IF ! global%error 

! Rahul - begin
  ALLOCATE(pPlag%vFracEOld(1,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%vFracEOld')
  END IF ! global%error 
! rahul - end

! ==============================================================================
! Additional variables  
! ==============================================================================

  ALLOCATE(pPlag%aivOld(nAiv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%aivOld')
  END IF ! global%error 

  ALLOCATE(pPlag%arvOld(nArv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%arvOld')
  END IF ! global%error

! ==============================================================================     
! Residuals 
! ==============================================================================
  
  ALLOCATE(pPlag%rhs(nCv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%rhs')
  END IF ! global%error  
  
  ALLOCATE(pPlag%rhsSum(nCv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%rhsSum')
  END IF ! global%error  

! ******************************************************************************
! Initialize memory
! ******************************************************************************

  DO iPcl = 1,nPclsMax
    DO iVar = 1,nCv
      pPlag%cvOld(iVar,iPcl)  = 0.0_RFREAL
      pPlag%rhs(iVar,iPcl)    = 0.0_RFREAL
      pPlag%rhsSum(iVar,iPcl) = 0.0_RFREAL
    END DO ! iVar

    DO iVar = 1,nAiv
      pPlag%aivOld(iVar,iPcl) = 0
    END DO ! iVar
    
    DO iVar = 1,nArv   
      pPlag%arvOld(iVar,iPcl) = 0
    END DO ! iVar
  END DO ! iPcl
 
  DO iVar = 1,pGrid%nCellsTot
    pPlag%vFracEOld(1,iVar) = 0.0_RFREAL
  END DO 

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_AllocMemTStep

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_AllocMemTStep.F90,v $
! Revision 1.2  2015/12/18 23:28:36  rahul
! Allocate memory to vFracEOld array.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:51  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:50:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.6  2007/03/06 23:15:32  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.5  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.4  2004/07/28 18:58:13  fnajjar
! Included new definition for nPclsTot from dynamic memory reallocation
!
! Revision 1.3  2004/07/16 20:09:39  fnajjar
! Bug fix to initialize cvOld instead of cv
!
! Revision 1.2  2004/03/08 23:02:28  fnajjar
! Added initialization section within DO-loop construct
!
! Revision 1.1  2004/02/26 21:00:39  haselbac
! Initial revision
!
!******************************************************************************

