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
! Purpose: Check grid speeds.
!
! Description: Compute norms of difference of LHS and RHS of GCL.
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
! $Id: RFLU_CheckGridSpeeds.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_CheckGridSpeeds(pRegion)

  USE ModDataTypes
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  USE ModSortSearch

  USE RFLU_ModGrid

  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region), POINTER :: pRegion

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: c1,c2,errorFlag,ic,ifc,indGs,iPatch
  REAL(RFREAL) :: fs,nm,term
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: checkSum,checkSumSorted
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pXyz,pXyzOld
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid,pGridOld
  TYPE(t_patch), POINTER :: pPatch

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_CheckGridSpeeds.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_CheckGridSpeeds',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking grid speeds...'
  END IF ! global%myProcid
  
! *****************************************************************************
! Set pointers and variables
! *****************************************************************************
    
  pGrid    => pRegion%grid
  pGridOld => pRegion%gridOld

  indGs = pGrid%indGs

! *****************************************************************************
! Allocate memory
! *****************************************************************************

  ALLOCATE(checkSum(pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'checkSum')
  END IF ! global%error

  ALLOCATE(checkSumSorted(pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'checkSumSorted')
  END IF ! global%error

  DO ic = 1,pGrid%nCellsTot ! Explicit loop because of ASCI White problem
    checkSum(ic) = 0.0_RFREAL
  END DO ! ic

! *****************************************************************************
! Interior faces
! *****************************************************************************

  DO ifc = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifc)
    c2 = pGrid%f2c(2,ifc)
    
    nm = pGrid%fn(XYZMAG,ifc)
    fs = pGrid%gs(indGs*ifc)      
    
    checkSum(c1) = checkSum(c1) + fs*nm
    checkSum(c2) = checkSum(c2) - fs*nm    
  END DO ! ifc

! *****************************************************************************
! Boundary faces
! *****************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
    
    DO ifc = 1,pPatch%nBFaces
      c1 = pPatch%bf2c(ifc)
      
      nm = pPatch%fn(XYZMAG,ifc)
      fs = pPatch%gs(indGs*ifc)
      
      checkSum(c1) = checkSum(c1) + fs*nm
    END DO ! ifc  
  END DO ! iPatch

! *****************************************************************************
! Check values and compute norms
! *****************************************************************************

  DO ic = 1,pGrid%nCells
    checkSum(ic) = checkSum(ic)*global%dtMin - (pGrid%vol(ic)-pGridOld%vol(ic))
    checkSumSorted(ic) = checkSum(ic)
  END DO ! ic

  CALL QuickSortRFREAL(checkSumSorted(1:pGrid%nCells),pGrid%nCells)

  term = 0.0_RFREAL

  DO ic = 1,pGrid%nCells
    term = term + checkSumSorted(ic)*checkSumSorted(ic)
  END DO ! ic

  term = SQRT(term/REAL(pGrid%nCells,RFREAL))
  
  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Error in GCL:'
    WRITE(STDOUT,'(A,5X,A,1X,E15.8)') SOLVER_NAME,'L2-norm:',term
    WRITE(STDOUT,'(A,5X,A,1X,E15.8,1X,I6)') SOLVER_NAME,'L8-norm:', & 
                                            MAXVAL(checkSum(1:pGrid%nCells)), & 
                                            MAXLOC(checkSum(1:pGrid%nCells))
  END IF ! global%myProcid

! *****************************************************************************
! End
! *****************************************************************************

  DEALLOCATE(checkSum,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'checkSum')
  END IF ! global%error
  
  DEALLOCATE(checkSumSorted,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'checkSum')
  END IF ! global%error  

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking grid speeds done.'
  END IF ! global%myProcid

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CheckGridSpeeds

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CheckGridSpeeds.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:47  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:59  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:56  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:01  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2003/03/27 14:32:18  haselbac
! Fixed bug: Must not use MINLOC/MAXLOC on sorted list
!
! Revision 1.3  2003/03/15 18:25:54  haselbac
! Changed loop limit
!
! Revision 1.2  2003/01/28 14:27:13  haselbac
! Cosmetic changes only, use parameters in fn
!
! Revision 1.1  2002/11/08 21:55:20  haselbac
! Initial revision
!
! Revision 1.1  2002/10/27 19:20:08  haselbac
! Initial revision
!
!******************************************************************************

