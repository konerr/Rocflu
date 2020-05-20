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
! Purpose: Compute grid spacing for "Shktb" or "Hass" case. 
!
! Description: None.
!
! Input:
!   pRegion             Region pointer
!
! Output: None.
!
! Notes: 
!   1. This routine creates computes the dx, dy and dz for a given region
!   2. Applicable only to shktb or hass case - CCMT problem
!   3. Note the computed dx, dy, dz are for zones containing particles, i.e.
!      the grid may be non-uniform outside of the particle bed
!
! ******************************************************************************
!
! $Id: RFLU_ComputeGridSpacingShktb.F90,v 1.2 2016/02/03 20:40:29 rahul Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ComputeGridSpacingShktb(pRegion)

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
  INTEGER :: errorFlag,icg,icgref
  REAL(RFREAL) :: X,Y,Z,tol
  REAL(RFREAL) :: xmin,xmax,ymin,ymax,zmin,zmax
  LOGICAL :: pclsinreg
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ComputeGridSpacingShktb.F90,v $ $Revision: 1.2 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeGridSpacingShktb',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Computing grid spacing...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set grid pointer
! ******************************************************************************

  pGrid => pRegion%grid
  
  tol = 1.0E-08_RFREAL
 
  pclsinreg = .FALSE.

  xMin = pRegion%plagInput%iniRandXMin
  xMax = pRegion%plagInput%iniRandXMax
  yMin = pRegion%plagInput%iniRandYMin
  yMax = pRegion%plagInput%iniRandYMax
  zMin = pRegion%plagInput%iniRandZMin
  zMax = pRegion%plagInput%iniRandZMax

  pGrid%dx = HUGE(1.0_RFREAL)
  pGrid%dy = HUGE(1.0_RFREAL)
  pGrid%dz = HUGE(1.0_RFREAL)  

  IF (pRegion%mixtInput%dimens .EQ. 2) THEN
    DO icg = 1,pGrid%nCells
      X = pGrid%cofg(XCOORD,icg)
      Y = pGrid%cofg(YCOORD,icg)
      IF ((xMin-X .LE. tol) .AND. (X-xMax .LE. tol) .AND. &
         (yMin-Y .LE. tol) .AND. (Y-yMax .LE. tol)) THEN
         pclsinreg = .TRUE.
         icgref = icg
         EXIT
      END IF
    END DO
  ELSEIF (pRegion%mixtInput%dimens .EQ. 3) THEN
    DO icg = 1,pGrid%nCells
      X = pGrid%cofg(XCOORD,icg)
      Y = pGrid%cofg(YCOORD,icg)
      Z = pGrid%cofg(ZCOORD,icg)  
      IF ((xMin-X .LE. tol) .AND. (X-xMax .LE. tol) .AND. &
         (yMin-Y .LE. tol) .AND. (Y-yMax .LE. tol) .AND. &
         (zMin-Z .LE. tol) .AND. (Z-zMax .LE. tol)) THEN
         pclsinreg = .TRUE.
         icgref = icg
         EXIT
      END IF
    END DO
  END IF

  IF (pclsinreg .EQV. .TRUE.) THEN 
    IF (pRegion%mixtInput%dimens .EQ. 2) THEN
      ! Compute dx
      DO icg = icgref+1,pGrid%nCells
        IF (DABS(pGrid%cofg(YCOORD,icg) - &
           pGrid%cofg(YCOORD,icgref)) .LT. tol) THEN
          pGrid%dx = MIN(pGrid%dx,ABS(pGrid%cofg(XCOORD,icg) &
                        -pGrid%cofg(XCOORD,icgref)))
        END IF
      END DO
      ! Compute dy
      DO icg = icgref+1,pGrid%nCells
        IF (DABS(pGrid%cofg(XCOORD,icg) - &
           pGrid%cofg(XCOORD,icgref)) .LT. tol) THEN
          pGrid%dy = MIN(pGrid%dy,ABS(pGrid%cofg(YCOORD,icg) &
                        -pGrid%cofg(YCOORD,icgref)))
        END IF
      END DO
      ! Compute dz
      pGrid%dz = 0.0_RFREAL
    ELSEIF (pRegion%mixtInput%dimens .EQ. 3) THEN
      ! Compute dx
      DO icg = icgref+1,pGrid%nCells
        IF ((DABS(pGrid%cofg(YCOORD,icg) - &
          pGrid%cofg(YCOORD,icgref)) .LT. tol) .AND. &
          (DABS(pGrid%cofg(ZCOORD,icg) - pGrid%cofg(ZCOORD,icgref)) &
          .LT. tol)) THEN
          pGrid%dx = MIN(pGrid%dx,ABS(pGrid%cofg(XCOORD,icg) &
                    -pGrid%cofg(XCOORD,icgref)))
        END IF
      END DO
      ! Compute dy
      DO icg = icgref+1,pGrid%nCells
        IF ((DABS(pGrid%cofg(XCOORD,icg) - &
          pGrid%cofg(XCOORD,icgref)) .LT. tol) .AND. &
          (DABS(pGrid%cofg(ZCOORD,icg) - pGrid%cofg(ZCOORD,icgref)) &
          .LT. tol) ) THEN
          pGrid%dy = MIN(pGrid%dy,ABS(pGrid%cofg(YCOORD,icg) &
                      -pGrid%cofg(YCOORD,icgref)))
        END IF
      END DO
      ! Compute dz
      DO icg = icgref+1,pGrid%nCells
        IF ((DABS(pGrid%cofg(XCOORD,icg) - &
          pGrid%cofg(XCOORD,icgref)) .LT. tol) .AND. &
          (DABS(pGrid%cofg(YCOORD,icg) - pGrid%cofg(YCOORD,icgref)) &
          .LT. tol) ) THEN
          pGrid%dz = MIN(pGrid%dz,ABS(pGrid%cofg(ZCOORD,icg) &
                        -pGrid%cofg(ZCOORD,icgref)))
        END IF
      END DO
    END IF
  ELSE
    IF (pRegion%mixtInput%dimens .EQ. 2) THEN
      ! Compute dx
      DO icg = 2,pGrid%nCells
        IF ( (DABS(pGrid%cofg(YCOORD,icg) - pGrid%cofg(YCOORD,1)) &
          .LT. tol) ) THEN
          pGrid%dx = ABS(pGrid%cofg(XCOORD,icg)-pGrid%cofg(XCOORD,1))
          EXIT
        END IF
      END DO
      ! Compute dy
      DO icg = 2,pGrid%nCells
        IF ( (DABS(pGrid%cofg(XCOORD,icg) - pGrid%cofg(XCOORD,1)) &
          .LT. tol) ) THEN
          pGrid%dy = ABS(pGrid%cofg(YCOORD,icg)-pGrid%cofg(YCOORD,1))
          EXIT
        END IF
      END DO
      ! Compute dz
      pGrid%dz = 0.0_RFREAL
    ELSEIF (pRegion%mixtInput%dimens .EQ. 3) THEN
      ! Compute dx
      DO icg = 2,pGrid%nCells  
        IF ( (DABS(pGrid%cofg(YCOORD,icg) - pGrid%cofg(YCOORD,1)) &
          .LT. tol) .AND. (DABS(pGrid%cofg(ZCOORD,icg) -  &
          pGrid%cofg(ZCOORD,1)) .LT. tol) ) THEN
          pGrid%dx = ABS(pGrid%cofg(XCOORD,icg)-pGrid%cofg(XCOORD,1))
          EXIT 
        END IF
      END DO  
      ! Compute dy
      DO icg = 2,pGrid%nCells
        IF ( (DABS(pGrid%cofg(XCOORD,icg) - pGrid%cofg(XCOORD,1)) &
          .LT. tol) .AND. (DABS(pGrid%cofg(ZCOORD,icg) -  &
          pGrid%cofg(ZCOORD,1)) .LT. tol) ) THEN
          pGrid%dy = ABS(pGrid%cofg(YCOORD,icg)-pGrid%cofg(YCOORD,1))
          EXIT
        END IF
      END DO
      ! Compute dz
      DO icg = 2,pGrid%nCells
        IF ( (DABS(pGrid%cofg(XCOORD,icg) - pGrid%cofg(XCOORD,1)) &
          .LT. tol) .AND. (DABS(pGrid%cofg(YCOORD,icg) -  &
          pGrid%cofg(YCOORD,1)) .LT. tol) ) THEN
          pGrid%dz = ABS(pGrid%cofg(ZCOORD,icg)-pGrid%cofg(ZCOORD,1))
          EXIT
        END IF
      END DO
    END IF
  END IF

! Rahul - temporary fix for this case
  pGrid%dx = pRegion%mixtInput%prepRealVal21 
  pGrid%dy = pRegion%mixtInput%prepRealVal22 
  pGrid%dz = pRegion%mixtInput%prepRealVal23 
! Rahul - end

  IF ( global%myProcid == MASTERPROC .AND. & 
    global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A,1X,3(E23.16))') SOLVER_NAME, &
         'Grid-spacing: dx,dy,dz = ', &
          pGrid%dx,pGrid%dy,pGrid%dz
  END IF ! global%verbLevel


! ******************************************************************************
! End
! ******************************************************************************  

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Computing grid spacing done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputeGridSpacingShktb

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeGridSpacingShktb.F90,v $
! Revision 1.2  2016/02/03 20:40:29  rahul
! Grid spacing (dx, dy and dz) is now read from the .inp file. This is a
! temporary fix.
!
! Revision 1.1  2015/08/12 03:58:23  brollin
! Computes grid spacing when using uniform grid throughout or uniform grid in a particle curtain in the shock tube.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! ******************************************************************************

