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
! Purpose: Compute the L2-norm of the density changes.
!
! Description: None.
!
! Input: 
!   regions     Array of regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ResidualNorm.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ResidualNorm(regions)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI
  USE ModParameters

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,icg,iReg
  REAL(RFREAL) :: dr,drSum,drSumTot
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pCvOld
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ResidualNorm.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_ResidualNorm',__FILE__)

! ******************************************************************************
! Sum density changes over local regions 
! ******************************************************************************

  drSum = 0.0_RFREAL

  DO iReg = 1,global%nRegionsLocal
    pCv    => regions(iReg)%mixt%cv
    pCvOld => regions(iReg)%mixt%cvOld

    DO icg = 1,regions(iReg)%grid%nCells
      dr    = pCv(CV_MIXT_DENS,icg) - pCvOld(CV_MIXT_DENS,icg)
      drSum = drSum + dr*dr
    END DO ! icg
  END DO ! iReg

! ******************************************************************************
! Reduce across processes 
! ******************************************************************************

  IF ( global%nRegions > 1 ) THEN 
    CALL MPI_AllReduce(drSum,drSumTot,1,MPI_RFREAL,MPI_SUM,global%mpiComm, &
                       errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
    END IF ! global%error
  ELSE
    drSumTot = drSum
  END IF ! global%nRegions

! ******************************************************************************
! Finalize
! ******************************************************************************

  global%residual = SQRT(drSumTot)

  IF ( global%currentIter == 1 ) THEN
    global%resInit = global%residual
  END IF ! global

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ResidualNorm

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ResidualNorm.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:48  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:57  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:01  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.6  2006/04/07 15:19:22  haselbac
! Removed tabs
!
! Revision 1.5  2005/04/15 15:07:24  haselbac
! Converted to MPI
!
! Revision 1.4  2002/09/09 15:51:56  haselbac
! global now under region
!
! Revision 1.3  2002/07/25 14:30:40  haselbac
! Added FEM call to find proper residual norm for parallel runs
!
! Revision 1.2  2002/06/14 20:19:59  haselbac
! Deleted ModLocal, changed local%nRegions to global%nRegionsLocal
!
! Revision 1.1  2002/05/04 17:01:59  haselbac
! Initial revision
!
! ******************************************************************************

