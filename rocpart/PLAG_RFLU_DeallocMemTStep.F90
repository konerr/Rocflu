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
! Purpose: Deallocate memory for Lagrangian particles related to time stepping.
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
! $Id: PLAG_RFLU_DeallocMemTStep.F90,v 1.2 2015/12/18 23:29:54 rahul Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_DeallocMemTStep(pRegion,pPlag)

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
  INTEGER :: errorFlag
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_DeallocMemTStep.F90,v $ $Revision: 1.2 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_DeallocMemTStep',__FILE__)

! ******************************************************************************
! Set pointers and variables 
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

! ==============================================================================
! Old state vector
! ==============================================================================

  DEALLOCATE(pPlag%cvOld,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPlag%cvOld')
  END IF ! global%error 

! Rahul - begin
  DEALLOCATE(pPlag%vFracEOld,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPlag%vFracEOld')
  END IF ! global%error 
! Rahul - end

! ==============================================================================
! Additional variables  
! ==============================================================================

  DEALLOCATE(pPlag%aivOld,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPlag%aivOld')
  END IF ! global%error 

  DEALLOCATE(pPlag%arvOld,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPlag%arvOld')
  END IF ! global%error

! ==============================================================================     
! Residuals 
! ==============================================================================
  
  DEALLOCATE(pPlag%rhs,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPlag%rhs')
  END IF ! global%error  
  
  DEALLOCATE(pPlag%rhsSum,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPlag%rhsSum')
  END IF ! global%error  
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_DeallocMemTStep

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_DeallocMemTStep.F90,v $
! Revision 1.2  2015/12/18 23:29:54  rahul
! Deallocate memory of the vFracEOld array.
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
! Revision 1.2  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.1  2004/02/26 21:00:44  haselbac
! Initial revision
!
!******************************************************************************

