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
! Purpose: Calculate max. allowable local/global time step in the case
!          of inviscid flow.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_TimeStepImpulse.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_TimeStepImpulse(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModParameters

#ifdef PLAG
  USE ModPartLag, ONLY: t_plag
  USE PLAG_ModParameters
#endif

  USE ModTools, ONLY: MakeNonZero

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: icg,iPcl
  INTEGER :: dtMinImpulseLoc(1) 
  REAL(RFREAL) :: dtMinImpulse,dt(3)
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

#ifdef PLAG
  TYPE(t_plag), POINTER :: pPlag
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_TimeStepImpulse.F90,v $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TimeStepImpulse',__FILE__)

! ******************************************************************************
! Get dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

#ifdef PLAG
  pPlag => pRegion%plag
#endif

! ******************************************************************************
! Initialize time step, set dt=0 only in cells with particles
! ******************************************************************************

  DO icg = 1,pGrid%nCellsTot
    pRegion%dt(icg) = REAL(ABS(CRAZY_VALUE_INT),KIND=RFREAL)
  END DO ! icg

! ******************************************************************************
! Local time step
! ******************************************************************************

#ifdef PLAG
  DO icg = 1,pGrid%nCells
    IF ( pPlag%dImpulseMax(1,icg) > 0.0_RFREAL ) THEN
      dt(1) = pPlag%dImpulseMax(XCOORD,icg) &
            /ABS(MakeNonZero(pPlag%forceTotal(XCOORD,icg)))
      dt(2) = pPlag%dImpulseMax(YCOORD,icg) &
            /ABS(MakeNonZero(pPlag%forceTotal(YCOORD,icg)))
      dt(3) = pPlag%dImpulseMax(ZCOORD,icg) &
            /ABS(MakeNonZero(pPlag%forceTotal(ZCOORD,icg)))

      pRegion%dt(icg) = MINVAL(dt(1:pRegion%mixtInput%dimens))

!IF (1==2) THEN
!WRITE(*,'(A,1X,I4,2X,I7,2X,E12.6)') "iReg,icg,timeStepImpulse:", &
!                                    pRegion%iRegionGlobal,icg,pRegion%dt(icg)
!WRITE(*,'(3(E12.6,2X))') pPlag%dImpulseMax(XCOORD,icg), &
!                         pPlag%forceTotal(XCOORD,icg), dt(1)
!WRITE(*,'(3(E12.6,2X))') pPlag%dImpulseMax(YCOORD,icg), &
!                         pPlag%forceTotal(YCOORD,icg), dt(2)
!WRITE(*,'(3(E12.6,2X))') pPlag%dImpulseMax(ZCOORD,icg), &
!                         pPlag%forceTotal(ZCOORD,icg), dt(3)
!END IF
    END IF
  END DO ! icg
#endif

! ******************************************************************************
! For unsteady flows, determine minimum step
! ******************************************************************************

  IF ( global%flowType == FLOW_UNSTEADY ) THEN 
    dtMinImpulse    = MINVAL(pRegion%dt(1:pGrid%nCells))
    dtMinImpulseLoc = MINLOC(pRegion%dt(1:pGrid%nCells))

    pRegion%dtMinImpulse    = dtMinImpulse
    pRegion%dtMinImpulseLoc = dtMinImpulseLoc(1)

!IF (1==2) THEN
!WRITE(*,*) "--------------------------------------------"
!WRITE(*,'(A,1X,I4,2X,E12.6)') "iReg, Min(timeStepImpulse):", &
!                                    pRegion%iRegionGlobal,dtMinImpulse
!WRITE(*,*) "--------------------------------------------"
!END IF
  END IF ! global%flowType

! ******************************************************************************
! Finalize 
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TimeStepImpulse

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_TimeStepImpulse.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
! ******************************************************************************

