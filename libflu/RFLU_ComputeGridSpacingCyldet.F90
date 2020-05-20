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
! Purpose: Compute grid spacing.
!
! Description: None.
!
! Input:
!   pRegion             Region pointer
!
! Output: None.
!
! Notes: 
!   1. This routine creates computes the radial, circumferential and axial grid 
!      spacing , i.e., drad, dthe and dz for a given region
!   2. Applicable only to cyldet case - CCMT problem
!   3. Note each region should contain cells of the same size
!   4. The central section/ section containing charge is a Cartesian grid and 
!      therefore the grid spacing drad, dthe and dz are not applicable (however
!      they are computed for the sake of completeness). 
!
! ******************************************************************************
!
! $Id: RFLU_ComputeGridSpacingCyldet.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ComputeGridSpacingCyldet(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModBndPatch, ONLY: t_patch
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModMPI

  IMPLICIT NONE
  
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,icg,nRegSum
  REAL(RFREAL) :: rad,the,tol,z
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ComputeGridSpacingCyldet.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeGridSpacingCyldet',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Computing grid spacing...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set grid pointer
! ******************************************************************************

  pGrid => pRegion%grid
  
  tol = 1.0E-08_RFREAL
  
  nRegSum = pRegion%mixtInput%prepIntVal5 * pRegion%mixtInput%prepIntVal6 * &
            pRegion%mixtInput%prepIntVal7 + &
            pRegion%mixtInput%prepIntVal8 * pRegion%mixtInput%prepIntVal9 * &
            pRegion%mixtInput%prepIntVal10

  IF (pRegion%iRegionGlobal .LE. nRegSum) THEN

    icg = 1
    pGrid%radMin = DSQRT(pGrid%cofg(XCOORD,icg)**2.0_RFREAL + &
                        pGrid%cofg(YCOORD,icg)**2.0_RFREAL)
    pGrid%theMin = DATAN2(pGrid%cofg(YCOORD,icg),pGrid%cofg(XCOORD,icg))
    pGrid%zMin   = pGrid%cofg(ZCOORD,icg)
    pGrid%SgntheMin = INT(DSIGN(1.0_RFREAL,pGrid%theMin))
 
    ! Compute drad
    DO icg = 2,pGrid%nCells
      rad = DSQRT(pGrid%cofg(XCOORD,icg)**2.0_RFREAL + &
                 pGrid%cofg(YCOORD,icg)**2.0_RFREAL)
      the = DATAN2(pGrid%cofg(YCOORD,icg),pGrid%cofg(XCOORD,icg))
      z   = pGrid%cofg(ZCOORD,icg)
      IF (DABS(the-pGrid%theMin) .LE. tol .AND. &
          DABS(z-pGrid%zMin) .LE. tol) THEN
        pGrid%Imax = icg - 1
        pGrid%drad = rad -pGrid%radMin
        EXIT
      END IF
    END DO

    ! Compute dthe
    DO icg = 2,pGrid%nCells
      rad = DSQRT(pGrid%cofg(XCOORD,icg)**2.0_RFREAL + &
                 pGrid%cofg(YCOORD,icg)**2.0_RFREAL)
      the = DATAN2(pGrid%cofg(YCOORD,icg),pGrid%cofg(XCOORD,icg))
      z   = pGrid%cofg(ZCOORD,icg)
      IF (DABS(rad-pGrid%radMin) .LE. tol .AND. &
          DABS(z-pGrid%zMin) .LE. tol) THEN
         pGrid%dthe = DABS(DABS(the) -DABS(pGrid%theMin))
         EXIT
      END IF
    END DO

    ! Compute dz
    IF (pRegion%mixtInput%dimens == 2) THEN
      pGrid%dz = 2.0_RFREAL*pGrid%zMin
      pGrid%Jmax = 1
      pGrid%Kmax = pGrid%nCells/(pGrid%Imax*pGrid%Jmax)
    ELSE
      DO icg = 2,pGrid%nCells
        rad = DSQRT(pGrid%cofg(XCOORD,icg)**2.0_RFREAL + &
                   pGrid%cofg(YCOORD,icg)**2.0_RFREAL)
        the = DATAN2(pGrid%cofg(YCOORD,icg),pGrid%cofg(XCOORD,icg))
        z   = pGrid%cofg(ZCOORD,icg)
        IF (DABS(the-pGrid%theMin) .LE. tol .AND. &
            DABS(rad-pGrid%radMin) .LE. tol) THEN
          pGrid%dz = DABS(DABS(z) - DABS(pGrid%zMin))
          pGrid%Kmax = (icg-1)/pGrid%Imax
          pGrid%Jmax = pGrid%nCells/(pGrid%Imax*pGrid%Kmax)
          EXIT
        END IF
      END DO
    END IF
   
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,1X,A,1X,3(E23.16))') SOLVER_NAME, &
            'Grid-spacing: dr,dtheta,dz = ', &
             pGrid%drad,pGrid%dthe,pGrid%dz
    END IF ! global%verbLevel

  END IF

! ******************************************************************************
! End
! ******************************************************************************  

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Computing grid spacing done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputeGridSpacingCyldet

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeGridSpacingCyldet.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! ******************************************************************************

