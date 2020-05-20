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
! Purpose: Check that have coupled boundaries defined for coupled simulations.
!
! Description: None.
!
! Input: 
!   regions     Regions data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_CheckCouplingInput.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_CheckCouplingInput(regions)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModError
  USE ModParameters
  USE ModMPI

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:) :: regions

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN)   :: RCSIdentString
  INTEGER :: cntrGlob,cntrLoc,errorFlag,iPatch,iReg
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_CheckCouplingInput.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_CheckCouplingInput',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking coupling input...' 
  END IF ! global%verbLevel

  cntrGlob = 0
  cntrLoc  = 0

! ******************************************************************************
! Count number of interacting patches on this processor
! ******************************************************************************
  
  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
    DO iPatch = 1,regions(iReg)%grid%nPatches
      pPatch => regions(iReg)%patches(iPatch)

      IF ( pPatch%bcCoupled /= BC_NOT_COUPLED ) THEN 
        cntrLoc = cntrLoc + 1
      END IF ! pPatch%bcCoupled
    END DO ! iPatch  
  END DO ! iReg
  
! ******************************************************************************
! Reduce to total number of interacting patches, exit if none
! ******************************************************************************

  CALL MPI_Reduce(cntrLoc,cntrGlob,1,MPI_INTEGER,MPI_SUM,MASTERPROC, & 
                  global%mpiComm,errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
  END IF ! global%errorFlag
  
  IF ( global%myProcid == MASTERPROC .AND. cntrGlob == 0 ) THEN 
    CALL ErrorStop(global,ERR_BCCOUPLED_NONE,__LINE__)
  END IF ! global%myProcid

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking coupling input done.' 
  END IF ! global%verbLevel   

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CheckCouplingInput

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CheckCouplingInput.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:30  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:47:51  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:58:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2005/04/15 15:05:58  haselbac
! Converted to MPI, cosmetics
!
! Revision 1.1  2003/05/01 14:04:23  haselbac
! Initial revision
!
! ******************************************************************************

