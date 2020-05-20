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
! Purpose: Collection of routines for absorbing boundary condition.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModGFM.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2009 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModGFM

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModError
  USE ModMPI

  USE ModSortSearch

  USE RFLU_ModDifferentiationCells

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModGFM.F90,v $'

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: RFLU_GFM_CreateLevelSet, &
            RFLU_GFM_DestroyLevelSet, &
            RFLU_GFM_InitLevelSet, &
            RFLU_GFM_SetLevelSet, &
            RFLU_GFM_SetGhostFluid

! ==============================================================================
! Private functions
! ==============================================================================

  PRIVATE :: RFLU_GFM_NullifyLevelSet

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS







! ******************************************************************************
!
! Purpose: Allocate memory for LevelSet. 
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_GFM_CreateLevelSet(pRegion)

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

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_GFM_CreateLevelSet',__FILE__)

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  ALLOCATE(pRegion%mixt%levelSet(global%gfmNParticles,pRegion%grid%nCellsTot), &
                                 STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%mixt%levelSet')
  END IF ! global%error

  ALLOCATE(pRegion%mixt%gradCellLevelSet(XCOORD:ZCOORD,global%gfmNParticles, &
                                         pRegion%grid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%mixt%gradCellLevelSet')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_GFM_CreateLevelSet









! ******************************************************************************
!
! Purpose: Destroy sigma variables.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_GFM_DestroyLevelSet(pRegion)

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

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_GFM_DestroyLevelSet',__FILE__)

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DEALLOCATE(pRegion%mixt%levelSet,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%levelSet')
  END IF ! global%error

  DEALLOCATE(pRegion%mixt%gradCellLevelSet,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%gradCellLevelSet')
  END IF ! global%error

  CALL RFLU_GFM_NullifyLevelSet(pRegion)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_GFM_DestroyLevelSet









! ******************************************************************************
!
! Purpose: Initialize LevelSet. 
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_GFM_InitLevelSet(pRegion)

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

  INTEGER :: icg,iParticle
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_GFM_InitLevelSet',__FILE__)

! ******************************************************************************
! Initialize sigma
! ******************************************************************************

  DO icg = 1,pRegion%grid%nCellsTot
    DO iParticle = 1,global%gfmNParticles
      pRegion%mixt%levelSet(iParticle,icg) = 0.0_RFREAL

      pRegion%mixt%gradCellLevelSet(XCOORD,iParticle,icg) = 0.0_RFREAL
      pRegion%mixt%gradCellLevelSet(YCOORD,iParticle,icg) = 0.0_RFREAL
      pRegion%mixt%gradCellLevelSet(ZCOORD,iParticle,icg) = 0.0_RFREAL
    END DO ! iParticle
  END DO ! icg

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_GFM_InitLevelSet








! ******************************************************************************
!
! Purpose: Nullify sigma variable.
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

SUBROUTINE RFLU_GFM_NullifyLevelSet(pRegion)

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

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_GFM_NullifyLevelSet',__FILE__)

! ******************************************************************************
! Nullify memory
! ******************************************************************************

  NULLIFY(pRegion%mixt%levelSet)
  NULLIFY(pRegion%mixt%gradCellLevelSet)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_GFM_NullifyLevelSet







! ******************************************************************************
!
! Purpose: Set LevelSet. 
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_GFM_SetLevelSet(pRegion)

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

  INTEGER :: errorFlag,icg,icgGuess,iParticle,Nx,Ny,Nz
  INTEGER, DIMENSION(:), ALLOCATABLE :: varInfoParticles
  REAL(RFREAL) :: distance,iMagnitude,magnitude,x1,x2,xc,y1,y2,yc
  REAL(RFREAL) :: dx,dy,dz,xCenter,yCenter,zCenter,radius,xMin,xMax,yMin,yMax, &
                  zMin,zMax
  REAL(RFREAL) :: vals(4)
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_GFM_SetLevelSet',__FILE__)

! ******************************************************************************
! Set varInfo
! TODO: Remove later
! ******************************************************************************

  ALLOCATE(varInfoParticles(global%gfmNParticles),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varInfoParticles')
  END IF ! global%error

  DO iParticle = 1,global%gfmNParticles
    varInfoParticles(iParticle) = V_PARTICLE_VAR1 + iParticle - 1
  END DO ! iParticle

! ******************************************************************************
! Set LevelSet Function depending on particle shape
! ******************************************************************************

  DO iParticle=1,global%gfmNParticles

! ==============================================================================
!   TODO: Check for particle shape
!   For now, assume there is only one rectangular particle of known
!   size and location
! ==============================================================================

! HARDCODE: Manoj: rectangle
    x1 = 0.45_RFREAL 
    x2 = 0.60_RFREAL

    y1 = 0.45_RFREAL
    y2 = 0.55_RFREAL
! --------------------------------------------

! HARDCODE: Manoj: circle
    xCenter = 0.5
    yCenter = 0.5
    radius  = 0.075
! --------------------------------------------

    DO icg = 1,pRegion%grid%nCellsTot
      xc = pRegion%grid%cofg(XCOORD,icg)
      yc = pRegion%grid%cofg(YCOORD,icg)

IF (1==2) THEN     
      IF ( (xc-x1)*(xc-x2)<0 .AND. (yc-y1)*(yc-y2)<0 ) THEN
        vals(1) = ABS(xc-x1)
        vals(2) = ABS(xc-x2)
        vals(3) = ABS(yc-y1)
        vals(4) = ABS(yc-y2)

        distance = -MINVAL(vals)
      ELSEIF ( (yc-y1)*(yc-y2)<0 ) THEN
        vals(1) = ABS(xc-x1)
        vals(2) = ABS(xc-x2)
        vals(3) = 2.0_RFREAL*vals(1) ! Dummy 
        vals(4) = 2.0_RFREAL*vals(2) ! Dummy

        distance = MINVAL(vals)
      ELSEIF ( (xc-x1)*(xc-x2)<0 ) THEN
        vals(1) = ABS(yc-y1)
        vals(2) = ABS(yc-y2)
        vals(3) = 2.0_RFREAL*vals(1) ! Dummy 
        vals(4) = 2.0_RFREAL*vals(2) ! Dummy

        distance = MINVAL(vals)
      ELSE
        vals(1) = SQRT((xc-x1)*(xc-x1)+(yc-y1)*(yc-y1))
        vals(2) = SQRT((xc-x1)*(xc-x1)+(yc-y2)*(yc-y2))
        vals(3) = SQRT((xc-x2)*(xc-x2)+(yc-y1)*(yc-y1))
        vals(4) = SQRT((xc-x2)*(xc-x2)+(yc-y2)*(yc-y2))

        distance = MINVAL(vals)
      END IF
END IF

IF (1==1) THEN     
      distance = SQRT((xc-xCenter)*(xc-xCenter)+(yc-yCenter)*(yc-yCenter)) &
                 - radius
END IF

      pRegion%mixt%levelSet(iParticle,icg) = distance
! TEMPORARY: Manoj: Testing LevelSet      
      pRegion%mixt%cv(iParticle,icg) = distance
! END TEMPORARY
    END DO ! icg
  END DO ! iParticle

! ==============================================================================
! Compute gradient of LevelSet function and normalize it
! TODO: Figure out what varInfoMixt does?
! ==============================================================================

  CALL RFLU_ComputeGradCellsWrapper(pRegion,1,global%gfmNParticles, &
                                    1,global%gfmNParticles, &
                                    varInfoParticles,pRegion%mixt%levelSet, &
                                    pRegion%mixt%gradCellLevelSet)

  DO iParticle=1,global%gfmNParticles
    DO icg = 1,pRegion%grid%nCellsTot
      magnitude = SQRT(pRegion%mixt%gradCellLevelSet(XCOORD,iParticle,icg) & 
                      *pRegion%mixt%gradCellLevelSet(XCOORD,iParticle,icg) &
                      +pRegion%mixt%gradCellLevelSet(YCOORD,iParticle,icg) &
                      *pRegion%mixt%gradCellLevelSet(YCOORD,iParticle,icg) &
                      +pRegion%mixt%gradCellLevelSet(ZCOORD,iParticle,icg) &
                      *pRegion%mixt%gradCellLevelSet(ZCOORD,iParticle,icg))
      iMagnitude = 1.0_RFREAL/magnitude

      pRegion%mixt%gradCellLevelSet(XCOORD,iParticle,icg) = iMagnitude &
                            *pRegion%mixt%gradCellLevelSet(XCOORD,iParticle,icg)
      pRegion%mixt%gradCellLevelSet(YCOORD,iParticle,icg) = iMagnitude &
                            *pRegion%mixt%gradCellLevelSet(YCOORD,iParticle,icg)
      pRegion%mixt%gradCellLevelSet(ZCOORD,iParticle,icg) = iMagnitude &
                            *pRegion%mixt%gradCellLevelSet(ZCOORD,iParticle,icg)
! TEMPORARY: Manoj: Testing LevelSet      
!      pRegion%mixt%cv(2,icg) = pRegion%mixt%gradCellLevelSet(XCOORD,iParticle,icg)
!      pRegion%mixt%cv(3,icg) = pRegion%mixt%gradCellLevelSet(YCOORD,iParticle,icg)
!      pRegion%mixt%cv(4,icg) = pRegion%mixt%gradCellLevelSet(ZCOORD,iParticle,icg)
! END TEMPORARY
    END DO ! icg
  END DO ! iParticle

! ******************************************************************************
! Deallocate varInfo
! ******************************************************************************

  DEALLOCATE(varInfoParticles,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varInfoParticles')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_GFM_SetLevelSet










! ******************************************************************************
!
! Purpose: Set GhostFluid. 
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_GFM_SetGhostFluid(pRegion)

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

  INTEGER :: errorFlag,icg,icgGuess,icgTarget,iLoc,iParticle,Nx,Ny,Nz
  REAL(RFREAL) :: distance,iMagnitude,magnitude,normalX,normalY,normalZ,x1,x2, &
                  xc,xTarget,y1,y2,yc,yTarget,z1,z2,zc,zTarget
  REAL(RFREAL) :: dx,dy,dz,xMin,xMax,yMin,yMax,zMin,zMax
  REAL(RFREAL) :: vals(4)
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_GFM_SetGhostFluid',__FILE__)

! ******************************************************************************
! Set LevelSet Function depending on particle shape
! ******************************************************************************

  DO iParticle=1,global%gfmNParticles

! ==============================================================================
!   TODO: Check for particle shape
!   For now, assume there is only one rectangular particle of known
!   size and location
! ==============================================================================

    Nx = global%gfmGridNx 
    Ny = global%gfmGridNy 
    NZ = global%gfmGridNz 

    xMin = global%gfmGridX1
    xMax = global%gfmGridX2
    yMin = global%gfmGridY1
    yMax = global%gfmGridY2
    zMin = global%gfmGridZ1
    zMax = global%gfmGridZ2

    dx = global%gfmGridDx
    dy = global%gfmGridDy
    dz = global%gfmGridDz

    DO icg = 1,pRegion%grid%nCellsTot
      xc = pRegion%grid%cofg(XCOORD,icg)
      yc = pRegion%grid%cofg(YCOORD,icg)

      IF ( pRegion%mixt%levelSet(iParticle,icg) < 0 ) THEN
        distance = -pRegion%mixt%levelSet(iParticle,icg)
        normalX  = pRegion%mixt%gradCellLevelSet(XCOORD,iParticle,icg)      
        normalY  = pRegion%mixt%gradCellLevelSet(YCOORD,iParticle,icg)      
        normalZ  = pRegion%mixt%gradCellLevelSet(ZCOORD,iParticle,icg)      

        xTarget = xc + 2.0_RFREAL*distance*normalX
        yTarget = yc + 2.0_RFREAL*distance*normalY
        zTarget = zc + 2.0_RFREAL*distance*normalZ

        icgGuess  = (Nx-1)*FLOOR(yc/dy)      + FLOOR(xc/dx)      + 1
        icgTarget = (Nx-1)*FLOOR(yTarget/dy) + FLOOR(xTarget/dx) + 1

        IF ( global%nRegionsLocal > 1 ) THEN
! ------- Find index of serial cell
          CALL BinarySearchInteger(pRegion%grid%sc2pc(1:1, &
                                   1:pRegion%grid%nCellsTot), &
                                   pRegion%grid%nCellsTot,icgTarget,iLoc)
                                   
! ------- Reset target cell to local index of serial cell
          IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN
            icgTarget = pRegion%grid%sc2pc(2,iLoc)
          ELSE
            CALL ErrorStop(global,ERR_CELL_NOT_FOUND,__LINE__)
          END IF ! iLoc
        END IF ! global%nRegionsLocal > 1

! TEMPORARY: Manoj: Testing LevelSet      
WRITE(*,'(2(E12.6,1X),3(I6,1X),E12.6)') xc,yc,icg,icgGuess,icgTarget,pRegion%mixt%levelSet(iParticle,icg)
! END TEMPORARY

! TEMPORARY: Manoj: GFM: For now just reverse vectors, to mimic no-slip
!                        This way normal velocity is reversed and tangential
!                        velocity anyway does not affect flux
        pRegion%mixt%cv(CV_MIXT_DENS,icg) = &
                                        pRegion%mixt%cv(CV_MIXT_DENS,icgTarget)
        pRegion%mixt%cv(CV_MIXT_XMOM,icg) = &
                                       -pRegion%mixt%cv(CV_MIXT_XMOM,icgTarget)
        pRegion%mixt%cv(CV_MIXT_YMOM,icg) = &
                                       -pRegion%mixt%cv(CV_MIXT_YMOM,icgTarget)
        pRegion%mixt%cv(CV_MIXT_ZMOM,icg) = &
                                       -pRegion%mixt%cv(CV_MIXT_ZMOM,icgTarget)
        pRegion%mixt%cv(CV_MIXT_ENER,icg) = &
                                        pRegion%mixt%cv(CV_MIXT_ENER,icgTarget)

        pRegion%mixt%cvOld(CV_MIXT_DENS,icg) = &
                                        pRegion%mixt%cvOld(CV_MIXT_DENS,icgTarget)
        pRegion%mixt%cvOld(CV_MIXT_XMOM,icg) = &
                                       -pRegion%mixt%cvOld(CV_MIXT_XMOM,icgTarget)
        pRegion%mixt%cvOld(CV_MIXT_YMOM,icg) = &
                                       -pRegion%mixt%cvOld(CV_MIXT_YMOM,icgTarget)
        pRegion%mixt%cvOld(CV_MIXT_ZMOM,icg) = &
                                       -pRegion%mixt%cvOld(CV_MIXT_ZMOM,icgTarget)
        pRegion%mixt%cvOld(CV_MIXT_ENER,icg) = &
                                        pRegion%mixt%cvOld(CV_MIXT_ENER,icgTarget)

!        pRegion%mixt%cv(CV_MIXT_DENS,icg) = &
!                                        pRegion%mixt%levelSet(iParticle,icg)
!        pRegion%mixt%cv(CV_MIXT_XMOM,icg) = &
!                            pRegion%mixt%gradCellLevelSet(XCOORD,iParticle,icg)
!        pRegion%mixt%cv(CV_MIXT_YMOM,icg) = &
!                            pRegion%mixt%gradCellLevelSet(YCOORD,iParticle,icg)
!        pRegion%mixt%cv(CV_MIXT_ZMOM,icg) = &
!                            pRegion%mixt%gradCellLevelSet(ZCOORD,iParticle,icg)
      END IF ! pRegion%mixt%levelSet(iParticle,icg)
    END DO ! icg
! TEMPORARY: Manoj: Testing LevelSet      
!WRITE(*,*) "Stopping here..."
!STOP
! END TEMPORARY
  END DO ! iParticle

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_GFM_SetGhostFluid











END MODULE RFLU_ModGFM

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModGFM.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
!
! ******************************************************************************
