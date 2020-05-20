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
!*******************************************************************************
!
! Purpose: Suite of routines to deallocate memory.
!
! Description: None.
!
! Notes: None.
!
!*******************************************************************************
!
! $Id: RFLU_ModDeallocateMemory.F90,v 1.2 2015/12/18 23:00:07 rahul Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!*******************************************************************************

MODULE RFLU_ModDeallocateMemory

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModMPI

  USE ModInterfaces, ONLY: RFLU_DecideNeedBGradFace

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_DeallocateMemoryAuxVars, &
            RFLU_DeallocateMemoryGSpeeds, &
            RFLU_DeallocateMemorySol, &
            RFLU_DeallocateMemorySolCv, &
            RFLU_DeallocateMemorySolDv, &
            RFLU_DeallocateMemorySolGv, &
            RFLU_DeallocateMemorySolTv, &
            RFLU_DeallocateMemoryTStep, & 
            RFLU_DeallocateMemoryTStep_C, & 
            RFLU_DeallocateMemoryTStep_I

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

  CHARACTER(CHRLEN) :: RCSIdentString = &
    '$RCSfile: RFLU_ModDeallocateMemory.F90,v $ $Revision: 1.2 $'

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS





! ******************************************************************************
!
! Purpose: Deallocate memory for auxiliary vars required in non-dissipative
!          solver by Hou-Mahesh.
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

SUBROUTINE RFLU_DeallocateMemoryAuxVars(pRegion)

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

  CALL RegisterFunction(global,'RFLU_DeallocateMemoryAuxVars',__FILE__)

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DEALLOCATE(pRegion%mixt%cvOld,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%cvOld')
  END IF ! global%error

  DEALLOCATE(pRegion%mixt%cvOld2,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%cvOld2')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemoryAuxVars






! ******************************************************************************
!
! Purpose: Deallocate memory for grid speeds.
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

SUBROUTINE RFLU_DeallocateMemoryGSpeeds(pRegion)

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

  INTEGER :: errorFlag,iPatch
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemoryGSpeeds',__FILE__)

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid      => pRegion%grid
  pMixtInput => pRegion%mixtInput

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

! ==============================================================================
! Grid motion active
! ==============================================================================

  IF ( pMixtInput%moveGrid .EQV. .TRUE. ) THEN

! ------------------------------------------------------------------------------
!   Interior faces
! ------------------------------------------------------------------------------

    DEALLOCATE(pGrid%gs,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%gs')
    END IF ! global%error

! ------------------------------------------------------------------------------
!   Patch faces
! ------------------------------------------------------------------------------

    IF ( pGrid%nPatches > 0 ) THEN
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

        DEALLOCATE(pPatch%gs,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%gs')
        END IF ! global%error
      END DO ! iPatch
    END IF ! pGrid%nPatches

! ==============================================================================
! Grid motion not active
! ==============================================================================

  ELSE

! ------------------------------------------------------------------------------
!   Interior faces
! ------------------------------------------------------------------------------

    DEALLOCATE(pGrid%gs,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%gs')
    END IF ! global%error

! ------------------------------------------------------------------------------
!   Patch faces
! ------------------------------------------------------------------------------

    IF ( pGrid%nPatches > 0 ) THEN
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

        DEALLOCATE(pPatch%gs,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%gs')
        END IF ! global%error
      END DO ! iPatch
    END IF ! pGrid%nPatches
  END IF ! pMixtInput%moveGrid

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemoryGSpeeds






! ******************************************************************************
!
! Purpose: Deallocate memory for mixture solution.
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

SUBROUTINE RFLU_DeallocateMemorySol(pRegion)

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
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemorySol',__FILE__)

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  CALL RFLU_DeallocateMemorySolCv(pRegion)
  CALL RFLU_DeallocateMemorySolDv(pRegion)
  CALL RFLU_DeallocateMemorySolGv(pRegion)
  CALL RFLU_DeallocateMemorySolTv(pRegion)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemorySol






! ******************************************************************************
!
! Purpose: Deallocate memory for conserved variables of mixture.
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

SUBROUTINE RFLU_DeallocateMemorySolCv(pRegion)

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
  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemorySolCv',__FILE__)

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DEALLOCATE(pRegion%mixt%cv,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%cv')
  END IF ! global%error

  DEALLOCATE(pRegion%mixt%cvInfo,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%cvInfo')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemorySolCv






! ******************************************************************************
!
! Purpose: Deallocate memory for dependent variables of mixture.
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

SUBROUTINE RFLU_DeallocateMemorySolDv(pRegion)

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
  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemorySolDv',__FILE__)

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  IF ( pRegion%mixtInput%nDv /= 0 ) THEN
    DEALLOCATE(pRegion%mixt%dv,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%dv')
    END IF ! global%error
  END IF ! solverType

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemorySolDv






! ******************************************************************************
!
! Purpose: Deallocate memory for gas variables of mixture.
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

SUBROUTINE RFLU_DeallocateMemorySolGv(pRegion)

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

  CALL RegisterFunction(global,'RFLU_DeallocateMemorySolGv',__FILE__)

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DEALLOCATE(pRegion%mixt%gv,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%gv')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemorySolGv





! ******************************************************************************
!
! Purpose: Deallocate memory for transport variables of mixture.
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

SUBROUTINE RFLU_DeallocateMemorySolTv(pRegion)

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
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemorySolTv',__FILE__)

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pMixtInput => pRegion%mixtInput

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  IF ( pMixtInput%computeTv .EQV. .TRUE. ) THEN
    DEALLOCATE(pRegion%mixt%tv,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%tv')
    END IF ! global%error
  END IF ! pMixtInput%computeTv

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemorySolTv






! ******************************************************************************
!
! Purpose: Deallocate memory for mixture time-stepping.
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

SUBROUTINE RFLU_DeallocateMemoryTStep(pRegion)

  USE RFLU_ModOLES

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

  INTEGER :: errorFlag,iPatch,nBFaces,nBFacesTot
  TYPE(t_grid), POINTER :: pGrid,pGridOld,pGridOld2
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemoryTStep',__FILE__)

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid      => pRegion%grid
  pGridOld   => pRegion%gridOld
  pGridOld2  => pRegion%gridOld2
  pMixtInput => pRegion%mixtInput

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

! ==============================================================================
! Reference solutions for sponge layer 
! ==============================================================================

  IF ( global%abcFlag .EQV. .TRUE. ) THEN
    IF ( global%abcKind == 0 ) THEN
      DEALLOCATE(pRegion%mixt%cvRef,STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%cvRef')
      END IF ! global%error
    END IF ! global%abcKind
  END IF ! global%abcFlag 

! ==============================================================================
! Old solutions
! ==============================================================================

  IF ( global%solverType /= SOLV_IMPLICIT_HM ) THEN
    DEALLOCATE(pRegion%mixt%cvOld,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%cvOld')
    END IF ! global%error
  END IF ! solverType

  SELECT CASE ( global%solverType )
    CASE ( SOLV_IMPLICIT_NK )
      DEALLOCATE(pRegion%mixt%cvOld1,STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%cvOld1')
      END IF ! global%error

      DEALLOCATE(pRegion%mixt%cvOld2,STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%cvOld2')
      END IF ! global%error
    CASE ( SOLV_IMPLICIT_HM )
      CALL RFLU_DeallocateMemoryAuxVars(pRegion)
  END SELECT

! ==============================================================================
! Time step
! ==============================================================================

  IF ( global%solverType /= SOLV_IMPLICIT_HM ) THEN
    DEALLOCATE(pRegion%dt,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%dt')
    END IF ! global%error
  END IF ! solverType

! ==============================================================================
! Residuals
! ==============================================================================

  IF ( global%solverType /= SOLV_IMPLICIT_HM ) THEN
    DEALLOCATE(pRegion%mixt%rhs,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%rhs')
    END IF ! global%error

    DEALLOCATE(pRegion%mixt%diss,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%diss')
    END IF ! global%error

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      DEALLOCATE(pRegion%mixt%rhsSum,STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%rhsSum')
      END IF ! global%error
    END IF ! global%flowType
  END IF ! solverType

! ==============================================================================
! Gradients
! ==============================================================================

! ------------------------------------------------------------------------------
! Cell gradients
! ------------------------------------------------------------------------------

  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
    DEALLOCATE(pRegion%mixt%gradCell,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%gradCell')
    END IF ! global%error

    DEALLOCATE(pRegion%mixt%gradCellOld,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%gradCellOld')
    END IF ! global%error

    DEALLOCATE(pRegion%mixt%gradCellOld2,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%gradCellOld2')
    END IF ! global%error
  ELSE
    IF ( (pMixtInput%spaceDiscr == DISCR_UPW_ROE     ) .OR. &
         (pMixtInput%spaceDiscr == DISCR_UPW_HLLC    ) .OR. & 
         (pMixtInput%spaceDiscr == DISCR_UPW_AUSMPLUS) .OR. &
         (pMixtInput%spaceDiscr == DISCR_UPW_AUSMPLUSUP) ) THEN
      IF ( pMixtInput%spaceOrder > 1 ) THEN
        DEALLOCATE(pRegion%mixt%gradCell,STAT=errorFlag)
        global%error = errorFlag
        IF (global%error /= ERR_NONE) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%gradCell')
        END IF ! global%error
      END IF ! pMixtInput%spaceOrder
    ELSE IF ( pMixtInput%spaceDiscr == DISCR_OPT_LES ) THEN
      DEALLOCATE(pRegion%mixt%gradCell,STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%gradCell')
      END IF ! global%error
    END IF ! pMixtInput%spaceDiscr
#ifdef PLAG
! - Gradient of density and momentum need in inviscid unsteady force
    DEALLOCATE(pRegion%mixt%gradCellE,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%gradCellE')
    END IF ! global%error
#endif
  END IF ! solverType

! ------------------------------------------------------------------------------
! Face gradients
! ------------------------------------------------------------------------------

  IF ( global%solverType /= SOLV_IMPLICIT_HM ) THEN
    IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN
      DEALLOCATE(pRegion%mixt%gradFace,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%gradFace')
      END IF ! global%error
    END IF ! pMixtInput%flowModel
  END IF ! solverType

  IF ( pGrid%nPatches > 0 ) THEN
    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)
  
      IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
        DEALLOCATE(pPatch%mixt%gradFace,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%mixt%gradFace')
        END IF ! global%error
      END IF ! RFLU_DecideNeedBGradFace
    END DO ! iPatch
  END IF ! pGrid%nPatches

! ==============================================================================
! Grid motion. NOTE grid speeds are allocated separately because they are
! written into grid file, and hence they need to be allocated in pre- and
! postprocessors also.
! ==============================================================================

  IF ( pMixtInput%moveGrid .EQV. .TRUE. ) THEN

! ------------------------------------------------------------------------------
!   Residual
! ------------------------------------------------------------------------------

    DEALLOCATE(pGrid%rhs,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%rhs')
    END IF ! global%error

! ------------------------------------------------------------------------------
!   Displacement
! ------------------------------------------------------------------------------

    IF ( pMixtInput%moveGridType /= MOVEGRID_TYPE_XYZ ) THEN
      DEALLOCATE(pGrid%disp,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%disp')
      END IF ! global%error
    END IF ! pMixtInput%moveGridType

! ------------------------------------------------------------------------------
!   Old coordinates
! ------------------------------------------------------------------------------

    DEALLOCATE(pGridOld%xyz,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%gridOld%xyz')
    END IF ! global%error

    IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN 
      DEALLOCATE(pGridOld2%xyz,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%gridOld2%xyz')
      END IF ! global%error
    END IF ! global%solverType

! ------------------------------------------------------------------------------
!   Old volume
! ------------------------------------------------------------------------------

    DEALLOCATE(pGridOld%vol,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%gridOld%vol')
    END IF ! global%error

    IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN 
      DEALLOCATE(pGridOld2%vol,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%gridOld2%vol')
      END IF ! global%error
    END IF ! global%solverType

#ifndef GENX
! ------------------------------------------------------------------------------
!   Patch displacements. NOTE allocate here only if not running inside GENX,
!   because when running inside GENX allocate displacements also for virtual
!   vertices.
! ------------------------------------------------------------------------------

    IF ( pGrid%nPatches > 0 ) THEN
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

        DEALLOCATE(pPatch%dXyz,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%dXyz')
        END IF ! global%error

      END DO ! iPatch
    END IF ! pGrid%nPatches
#endif
  END IF ! pMixtInput%moveGrid

#ifdef STATS
! ==============================================================================
! Time averaged statistics
! ==============================================================================

  IF ( (global%flowType == FLOW_UNSTEADY) .AND. &
       (global%doStat == ACTIVE) ) THEN
    DEALLOCATE(pRegion%mixt%tav,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%tav')
    END IF ! global%error
  END IF ! global%flowType
#endif

! ==============================================================================
! Optimal LES
! ==============================================================================

  IF ( pMixtInput%spaceDiscr == DISCR_OPT_LES ) THEN
    CALL RFLU_DestroyStencilsWeightsOLES(pRegion)
  END IF ! pMixtInput%spaceDiscr

! ==============================================================================
! Substantial derivative
! ==============================================================================

  DEALLOCATE(pRegion%mixt%sd,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%sd')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemoryTStep






! ******************************************************************************
!
! Purpose: Deallocate memory for mixture time-stepping for compressible fluid.
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

SUBROUTINE RFLU_DeallocateMemoryTStep_C(pRegion)

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

  INTEGER :: errorFlag,iPatch
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemoryTStep_C',__FILE__)

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

! ==============================================================================
! Mass fluxes
! ==============================================================================

  DEALLOCATE(pRegion%mixt%mfMixt,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%mfMixt')
  END IF ! global%error

  IF ( pRegion%grid%nPatches > 0 ) THEN
    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

      DEALLOCATE(pPatch%mfMixt,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%mfMixt')
      END IF ! global%error
    END DO ! iPatch
  END IF ! pRegion%grid%nPatches

! ==============================================================================
! Pressure correction, needed in non-dissipative implicit scheme of Hou-Mahesh.
! ==============================================================================

  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
    DEALLOCATE(pRegion%mixt%delP,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%delP')
    END IF ! global%error
  END IF ! solverType

! ==============================================================================
! Face-normal speeds, needed for non-dissipative implicit scheme by Hou-Mahesh.
! ==============================================================================

  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
    DEALLOCATE(pRegion%mixt%vfMixt,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%vfMixt')
    END IF ! global%error

    DEALLOCATE(pRegion%mixt%vfMixtOld,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%vfMixtOld')
    END IF ! global%error
  END IF ! solverType

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemoryTStep_C






! ******************************************************************************
!
! Purpose: Deallocate memory for mixture time-stepping for incompressible fluid.
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

SUBROUTINE RFLU_DeallocateMemoryTStep_I(pRegion)

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

  INTEGER :: errorFlag,iPatch
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemoryTStep_I',__FILE__)

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

! ==============================================================================
! Face velocities
! ==============================================================================

  DEALLOCATE(pRegion%mixt%vfMixt,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%vfMixt')
  END IF ! global%error

  IF ( pRegion%grid%nPatches > 0 ) THEN
    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

      DEALLOCATE(pPatch%vfMixt,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%vfMixt')
      END IF ! global%error
    END DO ! iPatch
  END IF ! pRegion%grid%nPatches

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemoryTStep_I





! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLU_ModDeallocateMemory


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModDeallocateMemory.F90,v $
! Revision 1.2  2015/12/18 23:00:07  rahul
! Added AUSM+up case to the memory deallocation IF statement.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:36  brollin
! New Stable version
!
! Revision 1.5  2009/07/08 19:11:45  mparmar
! Added deallocation of cvRef
!
! Revision 1.4  2008/12/06 08:43:39  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:16:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/11/28 23:05:20  mparmar
! Added deallocation of SOLV_IMPLICIT_HM specific arrays
!
! Revision 1.1  2007/04/09 18:49:24  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:40  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.9  2006/11/01 15:50:00  haselbac
! Changed so implicit-solver arrays dealt with properly
!
! Revision 1.8  2006/08/19 15:39:03  mparmar
! Moved region%mixt%bGradFace to patch%mixt%gradFace
!
! Revision 1.7  2006/02/08 21:03:47  hdewey2
! Added old2 quantities
!
! Revision 1.6  2005/09/22 17:11:04  hdewey2
! Added deallocation of cvOld1 and cvOld2 for transient implicit solver.
!
! Revision 1.5  2005/07/14 21:43:26  haselbac
! Added AUSM flux function to memory deallocation IF statement
!
! Revision 1.4  2004/12/19 15:47:08  haselbac
! Added memory deallocation for incompressible solver
!
! Revision 1.3  2004/10/19 19:38:48  haselbac
! Updated for GEN3
!
! Revision 1.2  2004/07/30 22:47:36  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.1  2004/03/19 21:15:19  haselbac
! Initial revision
!
! ******************************************************************************

