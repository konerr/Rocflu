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
! Purpose: Compute number of cells containing particles
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
! ******************************************************************************
!
! $Id: PLAG_RFLU_ComputeLagPclsPerReg.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_ComputeLagPclsPerReg(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModMPI
  USE ModParameters
  
  USE PLAG_ModParameters    
     
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
  INTEGER :: errorFlag,icl
  REAL(RFREAL) :: frac
  TYPE(t_global), POINTER :: global
  TYPE(t_grid),   POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: PLAG_RFLU_ComputeLagPclsPerReg.F90,v $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_ComputeLagPclsPerReg', &
                        __FILE__)
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                          'Computing number of particles per region ...'
  END IF ! global%verbLevel
  
! ******************************************************************************
! Set pointers and values
! ******************************************************************************

  pGrid => pRegion%grid
  pRegion%plag%nPcls = 0

  IF (pGrid%initPclPresent .EQV. .TRUE.) THEN
    frac = global%nPclsTarget - INT(global%nPclsTarget)
    frac = 1 - frac
    DO icl = 1,pGrid%nCells
      IF (pGrid%nPclsPerCell(icl) .LT. 0.0_RFREAL) THEN
        pGrid%nPclsPerCell(icl) = INT(0.0_RFREAL)
      ELSE
        IF (pGrid%nPclsPerCell(icl) .LE. frac) THEN
          pGrid%nPclsPerCell(icl) = INT(global%nPclsTarget)
        ELSE
          pGrid%nPclsPerCell(icl) = INT(global%nPclsTarget) + 1
        END IF
        pRegion%plag%nPcls = pRegion%plag%nPcls + pGrid%nPclsPerCell(icl)
      END IF
    END DO
  END IF

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A,1X,I8,1X,A,1X,I16)') SOLVER_NAME, &
          'Number of particles in region', pRegion%iRegionGlobal, ':', &
           pRegion%plag%nPcls
  END IF ! global%verbLevel

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Computing number of particles per region done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_ComputeLagPclsPerReg

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_ComputeLagPclsPerReg.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
! ******************************************************************************

