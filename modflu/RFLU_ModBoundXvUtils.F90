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
! Purpose: Collection of routines for boundary data.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModBoundXvUtils.F90,v 1.4 2016/02/08 20:03:20 fred Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModBoundXvUtils

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModError
  USE ModMPI

  USE RFLU_ModNSCBC, ONLY: RFLU_NSCBC_InitSW, &
                           RFLU_NSCBC_InitNSWHeat, &
                           RFLU_NSCBC_InitNSWTemp, &
                           RFLU_NSCBC_InitIFTotAng, &
                           RFLU_NSCBC_InitIFVelTemp, &
                           RFLU_NSCBC_InitOF, &
                           RFLU_NSCBC_InitFF, &
                           RFLU_NSCBC_InitIJ, &
                           RFLU_NSCBC_DecideHaveNSCBC
                           
  USE ModInterfaces, ONLY: MixtPerf_C_GRT, &
                           MixtPerf_Eo_DGPUVW, &
                           MixtPerf_Eo_GRTUVW, &
                           MixtPerf_G_CpR, &
                           MixtPerf_P_DEoGVm2, &
                           MixtPerf_R_M, &
                           MixtPerf_T_DPR

   USE RFLU_ModJWL
   USE RFLU_ModConvertCv,    ONLY: RFLU_ScalarConvertCvCons2Prim, &
                                   RFLU_ScalarConvertCvPrim2Cons
#ifdef SPEC
   USE ModInterfacesSpecies, ONLY: SPEC_GetSpeciesIndex
   USE ModSpecies,           ONLY: t_spec_input
#endif

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModBoundXvUtils.F90,v $ $Revision: 1.4 $'

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: RFLU_BXV_ComputeVarsCv, &
            RFLU_BXV_CreateVarsCv, &
            RFLU_BXV_CreateVarsDv, &
            RFLU_BXV_CreateVarsTStep, &            
            RFLU_BXV_DestroyVarsCv, &
            RFLU_BXV_DestroyVarsDv, &
            RFLU_BXV_DestroyVarsTStep, &
            RFLU_BXV_InitVars, &
            RFLU_BXV_NullifyVarsCv, &
            RFLU_BXV_NullifyVarsDv, &
            RFLU_BXV_NullifyVarsTStep, &
            RFLU_BXV_ReadVarsASCII, &
            RFLU_BXV_ReadVarsBinary, &            
            RFLU_BXV_ReadVarsWrapper, &
            RFLU_BXV_SetDependentVars, &
            RFLU_BXV_WriteVarsASCII, &
            RFLU_BXV_WriteVarsBinary, &            
            RFLU_BXV_WriteVarsWrapper

! ==============================================================================
! Private functions
! ==============================================================================

  PRIVATE :: RFLU_BXV_CompEnergyPatch, &             
             RFLU_BXV_CompMomEnergyPatch, & 
             RFLU_BXV_SetDependentVarsPatch

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS







! ******************************************************************************
!
! Purpose: Compute conserved variables on patches.
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

SUBROUTINE RFLU_BXV_ComputeVarsCv(pRegion)

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

  INTEGER :: iPatch
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_ComputeVarsCv',__FILE__)

! ******************************************************************************
! Loop over patches and set dependent variables
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      SELECT CASE( pPatch%bcType )
        CASE (BC_SLIPWALL)
        CASE (BC_NOSLIPWALL_HFLUX)
        CASE (BC_NOSLIPWALL_TEMP)
        CASE (BC_INFLOW_TOTANG)
        CASE (BC_INFLOW_VELTEMP)
          IF ( pPatch%mixt%switches(BCSWI_INFLOW_TYPE) == BCOPT_SUBSONIC ) THEN
            IF ( pPatch%reflect == BC_REFLECTING ) THEN
              CALL RFLU_BXV_CompMomEnergyPatch(pRegion,pPatch)
            END IF ! pPatch%reflect
          END IF ! pPatch%mixt%switches
        CASE (BC_OUTFLOW)
          IF ( pPatch%mixt%switches(BCSWI_OUTFLOW_TYPE) == BCOPT_SUBSONIC ) THEN
            IF ( pPatch%reflect == BC_REFLECTING ) THEN
              CALL RFLU_BXV_CompEnergyPatch(pRegion,pPatch)
            END IF ! pPatch%reflect
          END IF ! pPatch%mixt%switches
        CASE (BC_FARFIELD)
        CASE (BC_INJECTION)
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%bcType
    END IF ! pPatch%bcKind
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_ComputeVarsCv








! ******************************************************************************
!
! Purpose: Compute energy on patches.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_BXV_CompEnergyPatch(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: ifl,icg,indCp,indMol
  REAL(RFREAL) :: Eo,r,u,v,w,p,ru,rv,rw,mw,cp,gc,g,Yproducts,a,ir,t
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv
  TYPE(t_global), POINTER :: global

#ifdef SPEC
  LOGICAL :: scalarConvFlag
  INTEGER :: iCvSpecProducts
  TYPE(t_spec_input), POINTER :: pSpecInput
#endif
!FRED - Added JWL EOS capability - 9/11/15


! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

#ifdef SPEC
  pSpecInput => pRegion%specInput
#endif

  CALL RegisterFunction(global,'RFLU_BXV_CompEnergyPatch',__FILE__)

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pCv => pPatch%mixt%cv
  pDv => pPatch%mixt%dv
  pGv => pRegion%mixt%gv ! NOTE gv taken from cells for now

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

! ******************************************************************************
! Compute total energy for patch data
! ******************************************************************************

  DO ifl = 1,pPatch%nBFaces
    icg  = pPatch%bf2c(ifl) 

    r  = pCv(CV_MIXT_DENS,ifl)
    ru = pCv(CV_MIXT_XMOM,ifl)
    rv = pCv(CV_MIXT_YMOM,ifl)
    rw = pCv(CV_MIXT_ZMOM,ifl)
    p  = pDv(DV_MIXT_PRES,ifl)

    u = ru/r
    v = rv/r
    w = rw/r

    mw = pGv(GV_MIXT_MOL,indMol*icg)
    cp = pGv(GV_MIXT_CP ,indCp *icg)
    gc = MixtPerf_R_M(mw)
    g  = MixtPerf_G_CpR(cp,gc)

    IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
#ifdef SPEC
     IF (global%specUsed .EQV. .TRUE.) THEN
         IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
             scalarConvFlag = .FALSE.
         ELSE
             scalarConvFlag = .TRUE.
    CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
         END IF

    ir = 1.0_RFREAL/r
    iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')
    Yproducts = pRegion%spec%cv(iCvSpecProducts,icg)

    CALL RFLU_JWL_ComputeEnergyMixt(pRegion,icg,g,gc,p,r,Yproducts,a,Eo,t)

    IF (scalarConvFlag .EQV. .TRUE.) THEN
  CALL RFLU_ScalarConvertcvPrim2Cons(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
        END IF
      END IF
#endif
    ELSE

    Eo = MixtPerf_Eo_DGPUVW(r,g,p,u,v,w)
    END IF

    pCv(CV_MIXT_ENER,ifl) = r*Eo
  END DO ! ifl


! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_CompEnergyPatch









! ******************************************************************************
!
! Purpose: Compute momentum and energy on patches.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_BXV_CompMomEnergyPatch(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: ifl,icg,indCp,indMol,distrib
  REAL(RFREAL) :: Eo,r,u,v,w,ru,rv,rw,mw,cp,gc,g,T,P,Yproducts,a,ir,td
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,vals
  TYPE(t_global), POINTER :: global

#ifdef SPEC
  LOGICAL :: scalarConvFlag
  INTEGER :: iCvSpecProducts
  TYPE(t_spec_input), POINTER :: pSpecInput
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_CompMomEnergyPatch',__FILE__)

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  distrib = pPatch%mixt%distrib

  pCv  => pPatch%mixt%cv
  vals => pPatch%mixt%vals
  pGv => pRegion%mixt%gv ! NOTE gv taken from cells for now

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

! ******************************************************************************
! Compute Momentum variables for patch data
! ******************************************************************************

  DO ifl = 1,pPatch%nBFaces
    icg  = pPatch%bf2c(ifl) 

    r  = pCv(CV_MIXT_DENS,ifl)

    u = vals(BCDAT_INFLOW_U,distrib*ifl) 
    v = vals(BCDAT_INFLOW_V,distrib*ifl)
    w = vals(BCDAT_INFLOW_W,distrib*ifl)

    T = vals(BCDAT_INFLOW_T,distrib*ifl)
    P = vals(BCDAT_INFLOW_P,distrib*ifl) 
    mw = pGv(GV_MIXT_MOL,indMol*icg)
    cp = pGv(GV_MIXT_CP ,indCp *icg)
    gc = MixtPerf_R_M(mw)
    g  = MixtPerf_G_CpR(cp,gc)

    IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
#ifdef SPEC
     IF (global%specUsed .EQV. .TRUE.) THEN
         IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
             scalarConvFlag = .FALSE.
         ELSE
             scalarConvFlag = .TRUE.
    CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
         END IF

    ir = 1.0_RFREAL/r
    iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')
    Yproducts = pRegion%spec%cv(iCvSpecProducts,icg)

    CALL RFLU_JWL_ComputeEnergyMixt(pRegion,icg,g,gc,p,r,Yproducts,a,Eo,td)

    IF (scalarConvFlag .EQV. .TRUE.) THEN
  CALL RFLU_ScalarConvertcvPrim2Cons(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
        END IF
      END IF
#endif
    ELSE
    Eo = MixtPerf_Eo_GRTUVW(g,gc,T,u,v,w)
    END IF

    pCv(CV_MIXT_XMOM,ifl) = r*u
    pCv(CV_MIXT_YMOM,ifl) = r*v
    pCv(CV_MIXT_ZMOM,ifl) = r*w
    pCv(CV_MIXT_ENER,ifl) = r*Eo
  END DO ! ifl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_CompMomEnergyPatch








! ******************************************************************************
!
! Purpose: Allocate memory for conserved variables on boundary. 
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

SUBROUTINE RFLU_BXV_CreateVarsCv(pRegion)

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
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_CreateVarsCv',__FILE__)

! ******************************************************************************
! Loop over patches and allocate memory
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF (pPatch%bcKind == BC_KIND_NSCBC ) THEN
      ALLOCATE(pPatch%mixt%cv(pRegion%mixtInput%nCv,pPatch%nBFaces), &
               STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%cv')
      END IF ! global%error
    END IF ! pPatch%bcKind
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_CreateVarsCv







! ******************************************************************************
!
! Purpose: Allocate memory for dependent variables on boundary.
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

SUBROUTINE RFLU_BXV_CreateVarsDv(pRegion)

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

  INTEGER :: errorFlag, iPatch
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_CreateVarsDv',__FILE__)

! ******************************************************************************
! Loop over patches and allocate memory
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      ALLOCATE(pPatch%mixt%dv(pRegion%mixtInput%nDv,pPatch%nBFaces), &
               STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%dv')
      END IF ! global%error
    END IF ! pPatch%bcKind
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_CreateVarsDv





! ******************************************************************************
!
! Purpose: Allocate memory for time-stepping variables on boundary. 
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

SUBROUTINE RFLU_BXV_CreateVarsTStep(pRegion)

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
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_CreateVarsTStep',__FILE__)

! ******************************************************************************
! Loop over patches and Allocate memory
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF (pPatch%bcKind == BC_KIND_NSCBC ) THEN
      ALLOCATE(pPatch%mixt%cvOld(pRegion%mixtInput%nCv,pPatch%nBFaces), &
               STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%cvOld')
      END IF ! global%error

      ALLOCATE(pPatch%mixt%rhs(pRegion%mixtInput%nCv,pPatch%nBFaces), &
               STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%rhs')
      END IF ! global%error

      ALLOCATE(pPatch%mixt%rhsSum(pRegion%mixtInput%nCv,pPatch%nBFaces), &
               STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%rhsSum')
      END IF ! global%error
    END IF ! pPatch%bcKind
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_CreateVarsTStep






! ******************************************************************************
!
! Purpose: Destroy conserved variables on patches.
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

SUBROUTINE RFLU_BXV_DestroyVarsCv(pRegion)

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

  INTEGER :: errorFlag, iPatch
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_DestroyVarsCv',__FILE__)

! ******************************************************************************
! Loop over patches and deallocate memory
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      DEALLOCATE(pPatch%mixt%cv,STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%mixt%cv')
      END IF ! global%error
    END IF ! pPatch%bcKind
  END DO ! iPatch

  CALL RFLU_BXV_NullifyVarsCv(pRegion)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_DestroyVarsCv






! ******************************************************************************
!
! Purpose: Destroy dependent variables on patches.
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

SUBROUTINE RFLU_BXV_DestroyVarsDv(pRegion)

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
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_DestroyVarsDv',__FILE__)

! ******************************************************************************
! Loop over patches and deallocate memory
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      DEALLOCATE(pPatch%mixt%dv,STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%mixt%dv')
      END IF ! global%error
    END IF ! pPatch%bcKind
  END DO ! iPatch

  CALL RFLU_BXV_NullifyVarsDv(pRegion)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_DestroyVarsDv





! ******************************************************************************
!
! Purpose: Destroy variables on patches related to time-stepping.
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

SUBROUTINE RFLU_BXV_DestroyVarsTStep(pRegion)

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

  INTEGER :: errorFlag, iPatch
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_DestroyVarsTStep',__FILE__)

! ******************************************************************************
! Loop over patches and deallocate memory
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      DEALLOCATE(pPatch%mixt%rhs,STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%mixt%rhs')
      END IF ! global%error

      DEALLOCATE(pPatch%mixt%rhsSum,STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%mixt%rhs')
      END IF ! global%error
    END IF ! pPatch%bcKind
  END DO ! iPatch

  CALL RFLU_BXV_NullifyVarsTStep(pRegion)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_DestroyVarsTStep








! ******************************************************************************
!
! Purpose: Initialize patch variables.
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

SUBROUTINE RFLU_BXV_InitVars(pRegion)

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

  INTEGER :: errorFlag, iPatch
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_InitVars',__FILE__)

! ******************************************************************************
! Loop over patches and initialize patch variable arrays 
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      SELECT CASE( pPatch%bcType )
        CASE (BC_SLIPWALL)
          CALL RFLU_NSCBC_InitSW(pRegion,pPatch)
        CASE (BC_NOSLIPWALL_HFLUX)
          CALL RFLU_NSCBC_InitNSWHeat(pRegion,pPatch)
        CASE (BC_NOSLIPWALL_TEMP)
          CALL RFLU_NSCBC_InitNSWTemp(pRegion,pPatch)
        CASE (BC_INFLOW_TOTANG)
          CALL RFLU_NSCBC_InitIFTotAng(pRegion,pPatch)
        CASE (BC_INFLOW_VELTEMP)
          CALL RFLU_NSCBC_InitIFVelTemp(pRegion,pPatch)
        CASE (BC_OUTFLOW)
          CALL RFLU_NSCBC_InitOF(pRegion,pPatch)
        CASE (BC_FARFIELD)
          CALL RFLU_NSCBC_InitFF(pRegion,pPatch)
        CASE (BC_INJECTION)
          CALL RFLU_NSCBC_InitIJ(pRegion,pPatch)
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%bcType
    END IF ! pPatch%bcKind
  END DO ! iPatch

! ******************************************************************************
! Set dependent variables for each patch 
! ******************************************************************************

  CALL RFLU_BXV_SetDependentVars(pRegion)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_InitVars






! ******************************************************************************
!
! Purpose: Nullify conserved variables on patches.
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

SUBROUTINE RFLU_BXV_NullifyVarsCv(pRegion)

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

  INTEGER :: errorFlag, iPatch
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_NullifyVarsCv',__FILE__)

! ******************************************************************************
! Loop over patches and nullify memory
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      NULLIFY(pPatch%mixt%cv)
    END IF ! pPatch%bcKind
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_NullifyVarsCv







! ******************************************************************************
!
! Purpose: Nullify dependent variables on patches.
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

SUBROUTINE RFLU_BXV_NullifyVarsDv(pRegion)

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

  INTEGER :: errorFlag, iPatch
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_NullifyVarsDv',__FILE__)

! ******************************************************************************
! Loop over patches and nullify memory
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      NULLIFY(pPatch%mixt%dv)
    END IF ! pPatch%bcKind
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_NullifyVarsDv







! ******************************************************************************
!
! Purpose: Nullify time-stepping variables on patches.
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

SUBROUTINE RFLU_BXV_NullifyVarsTStep(pRegion)

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

  INTEGER :: errorFlag, iPatch
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_NullifyVarsTStep',__FILE__)

! ******************************************************************************
! Loop over patches and nullify memory
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      NULLIFY(pPatch%mixt%cvOld)
      NULLIFY(pPatch%mixt%rhs)
      NULLIFY(pPatch%mixt%rhsSum)
    END IF ! pPatch%bcKind
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_NullifyVarsTStep








! ******************************************************************************
!
! Purpose: Read patch variables in ASCII format.
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

SUBROUTINE RFLU_BXV_ReadVarsASCII(pRegion)

  USE ModBuildFileNames, ONLY: BuildFileNameUnsteady, &
                               BuildFileNameSteady

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: errorString,iFileName,sectionString
  INTEGER :: dummyInteger,errorFlag,iData,iFile,ifl,iPatch,iPatchGlobal, &
             loopCounter,nBFaces,nData
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_ReadVarsASCII',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII boundary-condition' // &
                                         ' variable file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Open file
! ******************************************************************************

  iFile = IF_DISTR

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.bcva', &
                               pRegion%iRegionGlobal,global%currentTime, &
                               iFileName)
  ELSE
    CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.bcva', &
                             pRegion%iRegionGlobal,global%currentIter, &
                             iFileName)
  END IF ! global%flowType

  OPEN(iFile,FILE=iFileName,FORM='FORMATTED',STATUS='OLD',IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! Read header
! ******************************************************************************

  READ(iFile,'(A)') sectionString
  IF ( TRIM(sectionString) /= '# ROCFLU boundary-condition variable file' ) THEN
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
  END IF ! TRIM

! ******************************************************************************
! Read data
! ******************************************************************************

  loopCounter = 0

  DO ! set up infinite loop
    loopCounter = loopCounter + 1

    READ(iFile,'(A)') sectionString

    SELECT CASE ( TRIM(sectionString) )

! ==============================================================================
!     Patch data
! ==============================================================================

      CASE ( '# Patch data' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Patch variable...'
        END IF ! global%verbLevel

        READ(iFile,*) iPatch,iPatchGlobal,nBFaces,nData

! ------------------------------------------------------------------------------
!       Check that input correct
! ------------------------------------------------------------------------------

        IF ( iPatch > pGrid%nPatches ) THEN
          WRITE(errorString,'(A,1X,I3)') 'Patch index invalid:',iPatch
          CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, &
                         TRIM(errorString))
        END IF ! iPatch

        pPatch => pRegion%patches(iPatch)

! ------------------------------------------------------------------------------
!       Check if this patch kind is NSCBC
! ------------------------------------------------------------------------------

        IF ( pPatch%bcKind /= BC_KIND_NSCBC ) THEN
          WRITE(errorString,'(A,1X,I3)') 'Error in writing NSCBC variables for patch :', &
                                         iPatchGlobal
          CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, &
                         TRIM(errorString))
        END IF ! pPatch%bcKind

        IF ( nBFaces /= pPatch%nBFaces ) THEN
          WRITE(errorString,'(A,1X,I6)') 'Number of faces invalid:',nBFaces
          CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, &
                           TRIM(errorString))
        END IF ! nBFaces

        IF ( iPatchGlobal /= pPatch%iPatchGlobal ) THEN
          WRITE(errorString,'(A,1X,I3)') 'Global patch index invalid:', &
                                         iPatchGlobal
          CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, &
                         TRIM(errorString))
        END IF ! iPatchGlobal

        IF ( nData /= pRegion%mixtInput%nCv ) THEN
          WRITE(errorString,'(A,1X,I3)') &
              'Number of pieces of data invalid:',nData
          CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, &
                         TRIM(errorString))
        END IF ! iPatchGlobal

! ------------------------------------------------------------------------------
!       Read data
! ------------------------------------------------------------------------------

        pPatch%mixt%cvState = CV_MIXT_STATE_CONS

        DO ifl = 1,pPatch%nBFaces
          DO iData = 1,pRegion%mixtInput%nCv
            READ(iFile,'(1X,I6,1X,I2,1X,E23.16)') dummyInteger,dummyInteger, &
                                                  pPatch%mixt%cv(iData,ifl)
          END DO ! iData
        END DO ! ifl

! ==============================================================================
!     End marker
! ==============================================================================

      CASE ( '# End' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
        END IF ! global%verbLevel

        EXIT

! ==============================================================================
!     Invalid section string
! ==============================================================================

      CASE DEFAULT
        IF ( global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
        END IF ! verbosityLevel

        CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END SELECT ! TRIM

! ==============================================================================
!   Guard against infinite loop - might be unnecessary because of read errors?
! ==============================================================================

    IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN
      CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
    END IF ! loopCounter
  END DO ! <empty>

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(iFile,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII boundary-condition' // &
                                         ' variable file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_ReadVarsASCII






! ******************************************************************************
!
! Purpose: Read patch variables in binary format.
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

SUBROUTINE RFLU_BXV_ReadVarsBinary(pRegion)

  USE ModBuildFileNames, ONLY: BuildFileNameUnsteady, &
                               BuildFileNameSteady

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: errorString,iFileName,sectionString
  INTEGER :: dummyInteger,errorFlag,iData,iFile,ifl,iPatch,iPatchGlobal, &
             loopCounter,nBFaces,nData
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_ReadVarsBinary',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary boundary-condition' // &
                                         ' variable file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Open file
! ******************************************************************************

  iFile = IF_DISTR

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.bcv', &
                               pRegion%iRegionGlobal,global%currentTime, &
                               iFileName)
  ELSE
    CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.bcv', &
                             pRegion%iRegionGlobal,global%currentIter, &
                             iFileName)
  ENDIF ! global%flowType

!  OPEN(iFile,FILE=iFileName,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=errorFlag)
!  global%error = errorFlag
! BBR - begin
    IF( global%solutFormat .EQ. FORMAT_BINARY )THEN
    OPEN(iFile,FILE=iFileName,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_L )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
         ACCESS="SEQUENTIAL",CONVERT="LITTLE_ENDIAN",IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_B )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
         ACCESS="SEQUENTIAL",CONVERT="BIG_ENDIAN",IOSTAT=errorFlag)
    END IF
    global%error = errorFlag
! BBR - end 

  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! Read header
! ******************************************************************************

  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# ROCFLU boundary-condition variable file' ) THEN
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
  END IF ! TRIM

! ******************************************************************************
! Read data
! ******************************************************************************

  loopCounter = 0

  DO ! set up infinite loop
    loopCounter = loopCounter + 1

    READ(iFile) sectionString

    SELECT CASE ( TRIM(sectionString) )

! ==============================================================================
!     Patch data
! ==============================================================================

      CASE ( '# Patch data' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Patch variable...'
        END IF ! global%verbLevel

        READ(iFile) iPatch,iPatchGlobal,nBFaces,nData

! ------------------------------------------------------------------------------
!       Check that input correct
! ------------------------------------------------------------------------------

        IF ( iPatch > pGrid%nPatches ) THEN
          WRITE(errorString,'(A,1X,I3)') 'Patch index invalid:',iPatch
          CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, &
                         TRIM(errorString))
        END IF ! iPatch

        pPatch => pRegion%patches(iPatch)

! ------------------------------------------------------------------------------
!       Check if this patch kind is NSCBC
! ------------------------------------------------------------------------------

        IF ( pPatch%bcKind /= BC_KIND_NSCBC ) THEN
          WRITE(errorString,'(A,1X,I3)') 'Error in writing NSCBC variables for patch :', &
                                         iPatchGlobal
          CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, &
                         TRIM(errorString))
        END IF ! pPatch%bcKind

        IF ( nBFaces /= pPatch%nBFaces ) THEN
          WRITE(errorString,'(A,1X,I6)') 'Number of faces invalid:',nBFaces
          CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, &
                           TRIM(errorString))
        END IF ! nBFaces

        IF ( iPatchGlobal /= pPatch%iPatchGlobal ) THEN
          WRITE(errorString,'(A,1X,I3)') 'Global patch index invalid:', &
                                         iPatchGlobal
          CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, &
                         TRIM(errorString))
        END IF ! iPatchGlobal

        IF ( nData /= pRegion%mixtInput%nCv ) THEN
          WRITE(errorString,'(A,1X,I3)') &
              'Number of pieces of data invalid:',nData
          CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, &
                         TRIM(errorString))
        END IF ! iPatchGlobal

! ------------------------------------------------------------------------------
!       Read data
! ------------------------------------------------------------------------------

        pPatch%mixt%cvState = CV_MIXT_STATE_CONS

        DO ifl = 1,pPatch%nBFaces
          DO iData = 1,pRegion%mixtInput%nCv
            READ(iFile) dummyInteger,dummyInteger,pPatch%mixt%cv(iData,ifl)
          END DO ! iData
        END DO ! ifl

! ==============================================================================
!     End marker
! ==============================================================================

      CASE ( '# End' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
        END IF ! global%verbLevel

        EXIT

! ==============================================================================
!     Invalid section string
! ==============================================================================

      CASE DEFAULT
        IF ( global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
        END IF ! verbosityLevel

        CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END SELECT ! TRIM

! ==============================================================================
!   Guard against infinite loop - might be unnecessary because of read errors?
! ==============================================================================

    IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN
      CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
    END IF ! loopCounter
  END DO ! <empty>

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(iFile,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary boundary-condition' // &
                                         ' variable file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_ReadVarsBinary








! ******************************************************************************
!
! Purpose: Wrapper for reading of boundary condition variables
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

SUBROUTINE RFLU_BXV_ReadVarsWrapper(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Local variables
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_ReadVarsWrapper',__FILE__)

! ******************************************************************************
! Read bcv files
! ******************************************************************************

  IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
    IF ( global%solutFormat == FORMAT_ASCII ) THEN
      CALL RFLU_BXV_ReadVarsASCII(pRegion)
!    ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
! BBR - begin
       ELSE IF ( global%solutFormat == FORMAT_BINARY .OR. & 
                global%solutFormat == FORMAT_BINARY_L .OR. &
                global%solutFormat == FORMAT_BINARY_B ) THEN
! BBR - end 
      CALL RFLU_BXV_ReadVarsBinary(pRegion)
    ELSE
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! global%solutFormat
  END IF ! RFLU_NSCBC_DecideHaveNSCBC(pRegion)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_ReadVarsWrapper







! ******************************************************************************
!
! Purpose: Set dependent variables on patches.
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

SUBROUTINE RFLU_BXV_SetDependentVars(pRegion)

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

  INTEGER :: iPatch
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_SetDependentVars',__FILE__)

! ******************************************************************************
! Loop over patches and set dependent variables if NSCBC active.
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF (pPatch%bcKind == BC_KIND_NSCBC ) THEN
      CALL RFLU_BXV_SetDependentVarsPatch(pRegion,pPatch)
    END IF ! pPatch%bcKind
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_SetDependentVars






! ******************************************************************************
!
! Purpose: Set dependent variables on patch.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pPatch	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_BXV_SetDependentVarsPatch(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: ifl,icg,indCp,indMol
  REAL(RFREAL) :: cpGas,Eo,gGas,ir,mm,r,u,v,w,p,rGas,Vm2,a,T,Y,eJWL,ePerf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv
  TYPE(t_global), POINTER :: global

#ifdef SPEC
  LOGICAL :: scalarConvFlag
  INTEGER :: iCvSpecProducts
  TYPE(t_spec_input), POINTER :: pSpecInput
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_SetDependentVarsPatch',__FILE__)

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pCv => pPatch%mixt%cv
  pDv => pPatch%mixt%dv
  pGv => pRegion%mixt%gv ! NOTE gv taken from cells for now

#ifdef SPEC
  pSpecInput => pRegion%specInput
#endif

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

! ******************************************************************************
! Compute dependent variables
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%fluidModel )

! ==============================================================================
!   Incompressible fluid model
! ==============================================================================

    CASE ( FLUID_MODEL_INCOMP )

! ==============================================================================
!   Compressible fluid model. NOTE check state of solution vector.
! ==============================================================================

    CASE ( FLUID_MODEL_COMP )
      SELECT CASE ( pRegion%mixtInput%gasModel )

! ------------------------------------------------------------------------------
!       Thermally and calorically perfect gas (single species)
! ------------------------------------------------------------------------------

        CASE ( GAS_MODEL_TCPERF )
          IF ( pPatch%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN
            CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
          END IF ! pPatch%mixt%cvState

          DO ifl = 1,pPatch%nBFaces
            icg  = pPatch%bf2c(ifl) 

            mm    = pGv(GV_MIXT_MOL,icg*indMol)
            cpGas = pGv(GV_MIXT_CP,icg*indCp)
            rGas  = MixtPerf_R_M(mm)
            gGas  = MixtPerf_G_CpR(cpGas,rGas)

            r  = pCv(CV_MIXT_DENS,ifl)
            ir = 1.0_RFREAL/r

            u  = ir*pCv(CV_MIXT_XMOM,ifl)
            v  = ir*pCv(CV_MIXT_YMOM,ifl)
            w  = ir*pCv(CV_MIXT_ZMOM,ifl)
            Eo = ir*pCv(CV_MIXT_ENER,ifl)
            
            Vm2 = u*u + v*v + w*w

            p = MixtPerf_P_DEoGVm2(r,Eo,gGas,Vm2)

            pDv(DV_MIXT_PRES,ifl) = p 
            pDv(DV_MIXT_TEMP,ifl) = MixtPerf_T_DPR(r,pDv(DV_MIXT_PRES,ifl),rGas)
            pDv(DV_MIXT_SOUN,ifl) = MixtPerf_C_GRT(gGas,rGas, &
                                                   pDv(DV_MIXT_TEMP,ifl))
          END DO ! ifl

         CASE ( GAS_MODEL_MIXT_JWL )
          IF ( pPatch%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN
            CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
          END IF ! pPatch%mixt%cvState
#ifdef SPEC
        IF (global%specUsed .EQV. .TRUE.) THEN
         IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
             scalarConvFlag = .FALSE.
         ELSE
             scalarConvFlag = .TRUE.
    CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
         END IF

          DO ifl = 1,pPatch%nBFaces
            icg  = pPatch%bf2c(ifl)

            mm    = pGv(GV_MIXT_MOL,icg*indMol)
            cpGas = pGv(GV_MIXT_CP,icg*indCp)
            rGas  = MixtPerf_R_M(mm)
            gGas  = MixtPerf_G_CpR(cpGas,rGas)

            r  = pCv(CV_MIXT_DENS,ifl)
            ir = 1.0_RFREAL/r

            u  = ir*pCv(CV_MIXT_XMOM,ifl)
            v  = ir*pCv(CV_MIXT_YMOM,ifl)
            w  = ir*pCv(CV_MIXT_ZMOM,ifl)
            Eo = ir*pCv(CV_MIXT_ENER,ifl)

            Vm2 = u*u + v*v + w*w

            iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')
            Y = pRegion%spec%cv(iCvSpecProducts,icg)

            CALL RFLU_JWL_ComputePressureMixt(pRegion,icg,ggas,rgas,Eo,r,Y,a,eJWL,ePerf,p,T)

            pDv(DV_MIXT_PRES,ifl) = p
            pDv(DV_MIXT_TEMP,ifl) = T
            pDv(DV_MIXT_SOUN,ifl) = a
          END DO!FRED - Added JWL EOS capabilities - 9/11/15

        IF (scalarConvFlag .EQV. .TRUE.) THEN
  CALL RFLU_ScalarConvertcvPrim2Cons(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
        END IF
      END IF
#endif
! ------------------------------------------------------------------------------
!       Default
! ------------------------------------------------------------------------------

        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pRegion%mixtInput%gasModel

! ==============================================================================
!   Default
! ==============================================================================

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixtInput%fluidModel

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_SetDependentVarsPatch







! ******************************************************************************
!
! Purpose: Write solution on patches in ASCII format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pPatch	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_BXV_WriteVarsASCII(pRegion)

  USE ModBuildFileNames, ONLY: BuildFileNameUnsteady, &
                               BuildFileNameSteady

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName,sectionString
  INTEGER :: dummyInteger,errorFlag,iData,iFile,ifl,iPatch
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_WriteVarsASCII',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII boundary-condition' // &
                                         ' variable file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Open file
! ******************************************************************************

  iFile = IF_DISTR

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.bcva', &
                               pRegion%iRegionGlobal,global%currentTime, &
                               iFileName)
  ELSE
    CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.bcva', &
                             pRegion%iRegionGlobal,global%currentIter, &
                             iFileName)
  ENDIF ! global%flowType

  OPEN(iFile,FILE=iFileName,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! Write header
! ******************************************************************************

  sectionString = '# ROCFLU boundary-condition variable file'
  WRITE(iFile,'(A)') sectionString

! ******************************************************************************
! Loop over patches
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)

! ==============================================================================
!   If bcKind = NSCBC, write data
! ==============================================================================

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      WRITE(iFile,'(A)') '# Patch data'
      WRITE(iFile,'(2(1X,I3),1X,I6,1X,I2)') iPatch,pPatch%iPatchGlobal, &
                                            pPatch%nBFaces, &
                                            pRegion%mixtInput%nCv

      DO ifl = 1,pPatch%nBFaces
        DO iData = 1,pRegion%mixtInput%nCv
          WRITE(iFile,'(1X,I6,1X,I2,1X,E23.16)') ifl,iData, &
                                                 pPatch%mixt%cv(iData,ifl)
        END DO ! iData
      END DO ! ifl
    END IF ! pPatch%bcKind 
  END DO ! iPatch

! ******************************************************************************
! Write footer
! ******************************************************************************

  sectionString = '# End'
  WRITE(iFile,'(A)') sectionString

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(iFile,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII boundary-condition' // &
                                         ' variable file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_WriteVarsASCII







! ******************************************************************************
!
! Purpose: Write solution on patches in binary format.
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

SUBROUTINE RFLU_BXV_WriteVarsBinary(pRegion)

  USE ModBuildFileNames, ONLY: BuildFileNameUnsteady, &
                               BuildFileNameSteady

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName,sectionString
  INTEGER :: dummyInteger,errorFlag,iData,iFile,ifl,iPatch
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_WriteVarsBinary',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary boundary-condition' // &
                                         ' variable file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Open file
! ******************************************************************************

  iFile = IF_DISTR

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.bcv', &
                               pRegion%iRegionGlobal,global%currentTime, &
                               iFileName)
  ELSE
    CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.bcv', &
                             pRegion%iRegionGlobal,global%currentIter, &
                             iFileName)
  ENDIF ! global%flowType

!  OPEN(iFile,FILE=iFileName,FORM='UNFORMATTED',STATUS='UNKNOWN', &
!       IOSTAT=errorFlag)
! BBR - begin
    IF( global%solutFormat .EQ. FORMAT_BINARY )THEN
    OPEN(iFile,FILE=iFileName,FORM='UNFORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_L )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
         ACCESS="SEQUENTIAL",CONVERT="LITTLE_ENDIAN",IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_B )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
         ACCESS="SEQUENTIAL",CONVERT="BIG_ENDIAN",IOSTAT=errorFlag)
    END IF
! BBR - end

  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! Write header
! ******************************************************************************

  sectionString = '# ROCFLU boundary-condition variable file'
  WRITE(iFile) sectionString

! ******************************************************************************
! Loop over patches
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)

! ==============================================================================
!   If bcKind = NSCBC, write data
! ==============================================================================

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      sectionString = '# Patch data'
      WRITE(iFile) sectionString
      WRITE(iFile) iPatch,pPatch%iPatchGlobal,pPatch%nBFaces, &
                   pRegion%mixtInput%nCv

      DO ifl = 1,pPatch%nBFaces
        DO iData = 1,pRegion%mixtInput%nCv
          WRITE(iFile) ifl,iData,pPatch%mixt%cv(iData,ifl)
        END DO ! iData
      END DO ! ifl
    END IF ! pPatch%bcKind 
  END DO ! iPatch

! ******************************************************************************
! Write footer
! ******************************************************************************

  sectionString = '# End'
  WRITE(iFile) sectionString

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(iFile,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary boundary-condition' // &
                                         ' variable file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_WriteVarsBinary







! ******************************************************************************
!
! Purpose: Wrapper for writing of boundary condition variables
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

SUBROUTINE RFLU_BXV_WriteVarsWrapper(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Local variables
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_WriteVarsWrapper',__FILE__)

! ******************************************************************************
! Read bcv files
! ******************************************************************************

  IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
    IF ( global%solutFormat == FORMAT_ASCII ) THEN
      CALL RFLU_BXV_WriteVarsASCII(pRegion)
!    ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
! BBR - begin
       ELSE IF ( global%solutFormat == FORMAT_BINARY .OR. & 
                global%solutFormat == FORMAT_BINARY_L .OR. &
                global%solutFormat == FORMAT_BINARY_B ) THEN
! BBR - end 
      CALL RFLU_BXV_WriteVarsBinary(pRegion)
    ELSE
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! global%solutFormat
  END IF ! RFLU_NSCBC_DecideHaveNSCBC(pRegion)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_WriteVarsWrapper






END MODULE RFLU_ModBoundXvUtils

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModBoundXvUtils.F90,v $
! Revision 1.4  2016/02/08 20:03:20  fred
! Fixing non-SPEC flag compiling error
!
! Revision 1.3  2016/02/04 19:58:26  fred
! Adding iterative JWL EOS capabilities for the cylindrical detonation case
!
! Revision 1.2  2015/07/23 23:11:18  brollin
! 1) The pressure coefficient of the  collision model has been changed back to its original form
! 2) New options in the format of the grid and solutions have been added. Now the user can choose the endianness, and convert from one to the over in rfluconv.
! 3) The solutions are now stored in folders named by timestamp or iteration number.
! 4) The address enty in the hashtable has been changed to an integer(8) for cases when the grid becomes very large.
! 5) RFLU_WritePM can now compute PM2 on the fly for the Macroscale problem
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:23  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:39  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/10/20 21:17:43  mparmar
! Modified read/write of binary boundary variable files
!
! Revision 1.3  2006/08/21 16:11:19  haselbac
! Clean-up: comments, order of routines, missing headers, indentation, etc
!
! Revision 1.2  2006/08/19 19:44:08  haselbac
! Significant clean-up and cosmetic changes
!
! Revision 1.1  2006/08/19 15:37:40  mparmar
! Initial revision
!
! ******************************************************************************

