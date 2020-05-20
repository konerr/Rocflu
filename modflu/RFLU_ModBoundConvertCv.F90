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
! Purpose: Suite of routines to convert state vector.
!
! Description: None.
!
! Notes: 
!
! ******************************************************************************
!
! $Id: RFLU_ModBoundConvertCv.F90,v 1.3 2016/02/08 20:02:41 fred Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModBoundConvertCv

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region

  USE RFLU_ModConvertCv,   ONLY: RFLU_ScalarConvertCvCons2Prim, &
                                 RFLU_ScalarConvertCvPrim2Cons

    
  IMPLICIT NONE
   
  PRIVATE
  PUBLIC :: RFLU_BXV_ConvertCvCons2Prim, & 
            RFLU_BXV_ConvertCvPrim2Cons
   
  SAVE    
     
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
  
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS




! ******************************************************************************
!
! Purpose: Convert conserved state vector to primitive variables.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to data of current region
!   pPatch              Pointer to data of current patch
!   cvStateFuture       Future state of conserved variables
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_BXV_ConvertCvCons2Prim(pRegion,pPatch,cvStateFuture)

  USE ModInterfaces, ONLY: MixtPerf_R_M, &
                           MixtPerf_T_DPR, &
                           MixtPerf_G_CpR

  USE RFLU_ModJWL

#ifdef SPEC
   USE ModInterfacesSpecies, ONLY: SPEC_GetSpeciesIndex
   USE ModSpecies,           ONLY: t_spec_input
#endif

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: cvStateFuture
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: icg,ifl,indMol
  REAL(RFREAL) :: gc,ir,mw,p,r,Yproducts,e,a,b,c,gamma,rgas,cpAir,mwAir,u,v,w
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv
  TYPE(t_global), POINTER :: global

#ifdef SPEC
  LOGICAL :: scalarConvFlag
  INTEGER :: iCvSpecProducts,iCvSpecAir
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
  TYPE(t_spec_input), POINTER :: pSpecInput
#endif
  !FRED - Added JWL EOS Capability 9/15/15

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_ConvertCvCons2Prim',__FILE__)

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pCv   => pPatch%mixt%cv
  pDv   => pPatch%mixt%dv

  pGv   => pRegion%mixt%gv ! pGv taken from interior domain

  indMol = pRegion%mixtInput%indMol

#ifdef SPEC
  pSpecInput => pRegion%specInput
  pCvSpec => pRegion%spec%cv
#endif


! ******************************************************************************
! Actual conversion
! ******************************************************************************

  SELECT CASE ( pPatch%mixt%cvState )

! ==============================================================================
!   Convert from conservative to primitive form
! ==============================================================================

    CASE ( CV_MIXT_STATE_CONS )
      SELECT CASE ( cvStateFuture )

! ------------------------------------------------------------------------------
!       Convert to duvwp form
! ------------------------------------------------------------------------------

        CASE ( CV_MIXT_STATE_DUVWP )
          pPatch%mixt%cvState = CV_MIXT_STATE_DUVWP

          DO ifl = 1,pPatch%nBFaces
            ir = 1.0_RFREAL/pCv(CV_MIXT_DENS,ifl)

            pCv(CV_MIXT_XVEL,ifl) = ir*pCv(CV_MIXT_XMOM,ifl)
            pCv(CV_MIXT_YVEL,ifl) = ir*pCv(CV_MIXT_YMOM,ifl)
            pCv(CV_MIXT_ZVEL,ifl) = ir*pCv(CV_MIXT_ZMOM,ifl)
            pCv(CV_MIXT_PRES,ifl) = pDv(DV_MIXT_PRES,ifl)
          END DO ! ifl

! ------------------------------------------------------------------------------
!       Convert to duvwt form
! ------------------------------------------------------------------------------

        CASE (CV_MIXT_STATE_DUVWT)
          pPatch%mixt%cvState = CV_MIXT_STATE_DUVWT

          SELECT CASE ( pRegion%mixtInput%fluidModel )

! --------- Compressible fluid model -----------------------------------------

            CASE ( FLUID_MODEL_COMP )
              SELECT CASE ( pRegion%mixtInput%gasModel )

! ------------- TC perfect gas or mixture thereof, pseudo-gas

                CASE ( GAS_MODEL_TCPERF, &
                       GAS_MODEL_MIXT_TCPERF, &
                       GAS_MODEL_MIXT_PSEUDO )
                  DO ifl = 1,pPatch%nBFaces
                    icg = pPatch%bf2c(ifl)

                    r  = pCv(CV_MIXT_DENS,ifl)
                    p  = pDv(DV_MIXT_PRES,ifl)
                    ir = 1.0_RFREAL/r

                    pCv(CV_MIXT_XVEL,ifl) = ir*pCv(CV_MIXT_XMOM,ifl)
                    pCv(CV_MIXT_YVEL,ifl) = ir*pCv(CV_MIXT_YMOM,ifl)
                    pCv(CV_MIXT_ZVEL,ifl) = ir*pCv(CV_MIXT_ZMOM,ifl)

                    mw = pGv(GV_MIXT_MOL,indMol*icg)
                    gc = MixtPerf_R_M(mw)

                    pCv(CV_MIXT_TEMP,ifl) = MixtPerf_T_DPR(r,p,gc)
                  END DO ! ifl

! ------------- Gas-liquid mixture

                CASE ( GAS_MODEL_MIXT_GASLIQ )
                  DO ifl = 1,pPatch%nBFaces
                    icg = pPatch%bf2c(ifl)

                    ir = 1.0_RFREAL/pCv(CV_MIXT_DENS,ifl)

                    pCv(CV_MIXT_XVEL,ifl) = ir*pCv(CV_MIXT_XMOM,ifl)
                    pCv(CV_MIXT_YVEL,ifl) = ir*pCv(CV_MIXT_YMOM,ifl)
                    pCv(CV_MIXT_ZVEL,ifl) = ir*pCv(CV_MIXT_ZMOM,ifl)

                    pCv(CV_MIXT_TEMP,ifl) = pDv(DV_MIXT_TEMP,ifl)
                  END DO ! ifl

!------------------Mixture of Detonation product and perfect gases---------

                 CASE (GAS_MODEL_MIXT_JWL)
#ifdef SPEC
        IF (global%specUsed .EQV. .TRUE.) THEN
         IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
             scalarConvFlag = .FALSE.
         ELSE
             scalarConvFlag = .TRUE.
    CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
         END IF

                  DO ifl = 1,pPatch%nBFaces

                  icg = pPatch%bf2c(ifl)
                  r = pCv(CV_MIXT_DENS,ifl)
                  ir = 1.0_RFREAL/r

                  pCv(CV_MIXT_XVEL,ifl) = ir*pCv(CV_MIXT_XMOM,ifl)
                  pCv(CV_MIXT_YVEL,ifl) = ir*pCv(CV_MIXT_YMOM,ifl)
                  pCv(CV_MIXT_ZVEL,ifl) = ir*pCv(CV_MIXT_ZMOM,ifl)

                  u = pCv(CV_MIXT_XVEL,ifl)
                  v = pCv(CV_MIXT_YVEL,ifl)
                  w = pCv(CV_MIXT_ZVEL,ifl)

                  e = (ir*pCv(CV_MIXT_ENER,ifl))-0.5_RFREAL*(u*u+v*v+w*w)

           iCvSpecProducts = SPEC_GetSpeciesIndex(global,pspecInput,'PRODUCTS')
           iCvSpecAir = SPEC_GetSpeciesIndex(global,pspecInput,'AIR')

           Yproducts = pCvSpec(iCvSpecProducts,icg)
           mwAir = pSpecInput%specType(iCvSpecAir)%pMaterial%molw
           cpAir = pSpecInput%specType(iCvSpecAir)%pMaterial%spht
           rgas = MixtPerf_R_M(mwAir)
           gamma = MixtPerf_G_CpR(cpAir,rgas)

    CALL RFLU_JWL_ComputePressureMixt(pRegion,icg,gamma,rgas,e,r,Yproducts,a,b,c,p,pCv(CV_MIXT_TEMP,ifl))

                  END DO

            IF (scalarConvFlag .EQV. .TRUE.) THEN
  CALL RFLU_ScalarConvertcvPrim2Cons(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
            END IF
         END IF
#endif
                CASE DEFAULT
                  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
              END SELECT ! pRegion%mixtInput%gasModel

! --------- Default ----------------------------------------------------------

            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pRegion%mixtInput%fluidModel

! ------------------------------------------------------------------------------
!       Default
! ------------------------------------------------------------------------------

        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! cvStateFuture

! ==============================================================================
!   Error - invalid input
! ==============================================================================

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pPatch%mixt%cvStateFuture

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_ConvertCvCons2Prim






! ******************************************************************************
!
! Purpose: Convert primitive state vector to consverved variables.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to data of current region
!   pPatch              Pointer to data of current patch
!   cvStateFuture       Future state of conserved variables
!
! Output: None.
!
! Notes:
!   1. Strictly speaking, cvStateFuture is not needed (there is only one
!      state for conserved variables), but kept for consistency with
!      RFLU_ConvertCvCons2Prim.
!
! ******************************************************************************

SUBROUTINE RFLU_BXV_ConvertCvPrim2Cons(pRegion,pPatch,cvStateFuture)

  USE ModInterfaces, ONLY: MixtPerf_Cv_CpR, &
                           MixtGasLiq_Eo_CvmTVm2, &
                           MixtPerf_Eo_DGPUVW, &
                           MixtPerf_G_CpR, &
                           MixtPerf_R_M

  USE RFLU_ModJWL

#ifdef SPEC
   USE ModInterfacesSpecies, ONLY: SPEC_GetSpeciesIndex
   USE ModSpecies,           ONLY: t_spec_input
#endif

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: cvStateFuture
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: icg,ifl,indCp,indMol
  REAL(RFREAL) :: cp,Cvm,cvg,cvl,cvv,g,gc,mw,p,r,Rg,Rv,rYg,rYl,rYv,t,u,v, &
                  Vm2,w,a,b,c,Yproducts,gamma,rgas,ir,e
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv
  TYPE(t_global), POINTER :: global

#ifdef SPEC
  LOGICAL :: scalarConvFlag
  INTEGER :: iCvSpecProducts
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
  TYPE(t_spec_input), POINTER :: pSpecInput
#endif
  !FRED - Added JWL EOS capabilities  9/15/15

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_ConvertCvPrim2Cons',__FILE__)

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pCv   => pPatch%mixt%cv
  pDv   => pPatch%mixt%dv

  pGv   => pRegion%mixt%gv ! pGv is taken from interior cells

#ifdef SPEC
  pSpecInput => pRegion%specInput
#endif


  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

! ******************************************************************************
! Actual conversion
! ******************************************************************************

  IF ( pPatch%mixt%cvState == CV_MIXT_STATE_DUVWP .OR. &
       pPatch%mixt%cvState == CV_MIXT_STATE_DUVWT ) THEN

! ==============================================================================
!   Convert from primitive to conservative form
! ==============================================================================

    SELECT CASE ( cvStateFuture )
      CASE ( CV_MIXT_STATE_CONS )
        pPatch%mixt%cvState = CV_MIXT_STATE_CONS

        SELECT CASE ( pRegion%mixtInput%gasModel )

! ------------------------------------------------------------------------------
!         Thermally and calorically perfect gas or pseudo-gas
! ------------------------------------------------------------------------------

          CASE ( GAS_MODEL_TCPERF, &
                 GAS_MODEL_MIXT_TCPERF, &
                 GAS_MODEL_MIXT_PSEUDO )
            DO ifl = 1,pPatch%nBFaces
              icg = pPatch%bf2c(ifl)

              r = pCv(CV_MIXT_DENS,ifl)
              u = pCv(CV_MIXT_XVEL,ifl)
              v = pCv(CV_MIXT_YVEL,ifl)
              w = pCv(CV_MIXT_ZVEL,ifl)
              p = pDv(DV_MIXT_PRES,ifl)

              pCv(CV_MIXT_XMOM,ifl) = r*u
              pCv(CV_MIXT_YMOM,ifl) = r*v
              pCv(CV_MIXT_ZMOM,ifl) = r*w

              cp = pGv(GV_MIXT_CP,indCp*icg)
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)

              pCv(CV_MIXT_ENER,ifl) = r*MixtPerf_Eo_DGPUVW(r,g,p,u,v,w)
            END DO ! icg

!------------------Mixture of detonation products and perfect gas--------------

          CASE(GAS_MODEL_MIXT_JWL)
#ifdef SPEC
          IF ( global%specUsed .EQV. .TRUE. ) THEN

           IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
             scalarConvFlag = .FALSE.
         ELSE
             scalarConvFlag = .TRUE.
    CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
           END IF

          DO ifl = 1,pPatch%nBFaces

              icg = pPatch%bf2c(ifl)
              r = pCv(CV_MIXT_DENS,ifl)
              ir = 1.0_RFREAL/r

              u = pCv(CV_MIXT_XVEL,ifl)
              v = pCv(CV_MIXT_YVEL,ifl)
              w = pCv(CV_MIXT_ZVEL,ifl)
              p = pDv(DV_MIXT_PRES,ifl)

              pCv(CV_MIXT_XMOM,icg) = r*u
              pCv(CV_MIXT_YMOM,icg) = r*v
              pCv(CV_MIXT_ZMOM,icg) = r*w

           iCvSpecProducts = SPEC_GetSpeciesIndex(global,pspecInput,'PRODUCTS')

           Yproducts = pCvSpec(iCvSpecProducts,icg)

  e = Yproducts*pDv(DV_MIXT_EJWL,ifl) + (1.0_RFREAL-Yproducts)*pDv(DV_MIXT_EPERF,ifl)

              pCv(CV_MIXT_ENER,ifl) = r *(e + 0.5_RFREAL*(u*u+v*v+w*w))

            END DO

           IF (scalarConvFlag .EQV. .TRUE.) THEN
  CALL RFLU_ScalarConvertcvPrim2Cons(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
           END IF

             END IF ! global%specUsed
#else
              CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, &
                             'Can only be used with species module.')
#endif

! ------------------------------------------------------------------------------
!         Mixture of gas, liquid, and vapor
! ------------------------------------------------------------------------------

          CASE ( GAS_MODEL_MIXT_GASLIQ )
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, &
                           'Can only be used with species module.')

! ------------------------------------------------------------------------------
!         Other or invalid gas models
! ------------------------------------------------------------------------------

          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pRegion%mixtInput

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! cvStateFuture

! ==============================================================================
! Error - invalid input
! ==============================================================================

  ELSE
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END IF ! pPatch%mixt%cvState

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_ConvertCvPrim2Cons






  
! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLU_ModBoundConvertCv


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModBoundConvertCv.F90,v $
! Revision 1.3  2016/02/08 20:02:41  fred
! Fixing non-SPEC flag compiling error
!
! Revision 1.2  2016/02/04 19:58:42  fred
! Adding iterative JWL EOS capabilities for the cylindrical detonation case
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:36  brollin
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
! Revision 1.1  2006/08/19 15:37:43  mparmar
! Initial revision
!
! ******************************************************************************

