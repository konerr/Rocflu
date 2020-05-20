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
! Purpose: Wrapper for initializing solution from scratch.
!
! Description: None.
!
! Input:
!   pRegion     	Pointer to region
!   pRegionSerial	Pointer to serial region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_InitFlowSerialWrapper.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_InitFlowSerialWrapper(pRegion,pRegionSerial)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  USE RFLU_ModCopyData, ONLY: RFLU_COPY_CellDataS2P_R2D

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion,pRegionSerial

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_InitFlowSerialWrapper.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_InitFlowSerialWrapper',__FILE__)

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing from serial solution...'
  END IF ! global%verbLevel

  IF ( global%verbLevel > VERBOSE_LOW ) THEN 
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                     pRegion%iRegionGlobal
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Mixture
! ******************************************************************************

  pRegion%mixt%cvState = CV_MIXT_STATE_CONS

  CALL RFLU_COPY_CellDataS2P_R2D(global,pGrid,pRegion%mixt%cv, & 
                                 pRegionSerial%mixt%cv)

  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
    CALL RFLU_COPY_CellDataS2P_R2D(global,pGrid,pRegion%mixt%cvOld, & 
                                   pRegionSerial%mixt%cvOld)
  END IF ! global%solverType

! ******************************************************************************
! Physical modules
! ******************************************************************************

#ifdef SPEC
! ==============================================================================
! Species
! ==============================================================================

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL RFLU_COPY_CellDataS2P_R2D(global,pGrid,pRegion%spec%cv, & 
                                   pRegionSerial%spec%cv)
  END IF ! global%specUsed
#endif

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Initializing from serial solution done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InitFlowSerialWrapper


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_InitFlowSerialWrapper.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2010/03/15 00:54:47  mparmar
! Added support for auxiliary data
!
! Revision 1.3  2008/12/06 08:43:55  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:08  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:55:41  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2005/08/19 02:36:02  haselbac
! Adapted to changes in RFLU_ModCopyData
!
! Revision 1.1  2005/04/15 15:08:21  haselbac
! Initial revision
!
! ******************************************************************************

