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
! Purpose: Set dependent variables.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!  icgBeg       Beginning cell index
!  icgEnd       Ending cell index
!
! Output: None.
!
! Notes: 
!   1. Restricted to perfect gas or JWL EOS.
!
! ******************************************************************************
!
! $Id: RFLU_SetDependentVars.F90,v 1.2 2016/02/04 19:59:23 fred Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_SetDependentVars(pRegion,icgBeg,icgEnd)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid

#ifdef PLAG
  USE PLAG_ModParameters
#endif  
#ifdef SPEC
  USE ModSpecies, ONLY: t_spec_input, &
                        t_spec_type
  USE ModInterfacesSpecies, ONLY: SPEC_GetSpeciesIndex                      
#endif

  USE RFLU_ModConvertCv, ONLY: RFLU_ScalarConvertCvCons2Prim, &
                               RFLU_ScalarConvertCvPrim2Cons

  USE ModInterfaces, ONLY: MixtGasLiq_C, &
                           MixtGasLiq_P, &
                           MixtLiq_C2_Bp, &
                           MixtLiq_D_DoBpPPoBtTTo, &
                           MixtPerf_C_GRT, &
                           MixtPerf_C2_GRT, &
                           MixtPerf_Cv_CpR, &                       
                           MixtPerf_D_PRT, &
                           MixtPerf_G_CpR, & 
                           MixtPerf_HM_T_DGMP, &
                           MixtPerf_P_DEoGVm2, & 
                           MixtPerf_R_M, &
                           MixtPerf_T_CvEoVm2, & 
                           MixtPerf_T_DPR
                                 
  USE RFLU_ModJWL

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(IN) :: icgBeg,icgEnd
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: icg,indCp,indMol,indVFracE
  REAL(RFREAL) :: a,cp,e,eJWL,Eo,ePerf,g,gc,ir,mw,nTol,p,r,refGamma,refM,rg,T, &
                  term,Vm2,YExplosive,YProducts  
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv
  TYPE(t_global), POINTER :: global

#ifdef SPEC
  LOGICAL :: scalarConvFlag
  INTEGER :: iCvSpecAir,iCvSpecExplosive,iCvSpecProducts,iSpec
  REAL(RFREAL) :: Bg2,Bl2,Bp,Bt,Bv2,Cg2,Cl2,cpAir,cpg,cpProducts,cpm,Cv2, &
                  cvg,cvl,Cvm,cvv,gAir,gcAir,gProducts,gcProducts,gcg,gcm, &
                  gcv,gm,immg,mmg,mmm,mwAir,mwProducts,phip,po,rl,ro,rv, &
                  rYg,rYi,rYl,rYv,to,vfg,vfl,vfv,Yg
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
  TYPE(t_spec_input), POINTER :: pSpecInput  
  TYPE(t_spec_type), POINTER :: pSpecType
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_SetDependentVars.F90,v $ $Revision: 1.2 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SetDependentVars',__FILE__)

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pCv => pRegion%mixt%cv
  pDv => pRegion%mixt%dv 
  pGv => pRegion%mixt%gv 
          
#ifdef SPEC  
  pCvSpec => pRegion%spec%cv
  pSpecInput => pRegion%specInput
#endif

  refGamma = global%refGamma
  refM     = global%refVelocity/ &
              (refGamma*global%refPressure/global%refDensity)**0.5_RFREAL 

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

  indVFracE = 0
#ifdef PLAG    
  indVFracE = pRegion%plagInput%indVFracE
#endif

  nTol = 1.0E-14_RFREAL

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
          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN
              IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_DUVWP ) THEN 
                CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
              END IF ! pRegion%mixt%cvState    

              DO icg = icgBeg,icgEnd
                r  = pCv(CV_MIXT_DENS,icg)
                p  = pCv(CV_MIXT_PRES,icg)

                T  = MixtPerf_HM_T_DGMP(r,refGamma,refM,p)

                pDv(DV_MIXT_TEMP,icg) = T
              END DO ! icg
            END IF ! pRegion%mixtInput%flowModel
          ELSE
            IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN 
              CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
            END IF ! pRegion%mixt%cvState    

            DO icg = icgBeg,icgEnd
              gc = MixtPerf_R_M(pGv(GV_MIXT_MOL,icg*indMol))
              g  = MixtPerf_G_CpR(pGv(GV_MIXT_CP,icg*indCp),gc)

              r  = pCv(CV_MIXT_DENS,icg)
              ir = 1.0_RFREAL/r
              Eo = ir*pCv(CV_MIXT_ENER,icg)

              Vm2 = ir*ir*(pCv(CV_MIXT_XMOM,icg)*pCv(CV_MIXT_XMOM,icg) + &
                           pCv(CV_MIXT_YMOM,icg)*pCv(CV_MIXT_YMOM,icg) + &
                           pCv(CV_MIXT_ZMOM,icg)*pCv(CV_MIXT_ZMOM,icg))

              ! Subbu - Debug
              IF (global%plagUsed .EQV. .TRUE.) THEN
                 r = r/(1.0_RFREAL - pRegion%plag%vFracE(1,indVFracE*icg))
              END IF
              ! Subbu - End Debug

              pDv(DV_MIXT_PRES,icg) = MixtPerf_P_DEoGVm2(r,Eo,g,Vm2)
              pDv(DV_MIXT_TEMP,icg) = MixtPerf_T_DPR(r,pDv(DV_MIXT_PRES,icg),gc)
              pDv(DV_MIXT_SOUN,icg) = MixtPerf_C_GRT(g,gc,pDv(DV_MIXT_TEMP,icg))    
            END DO ! icg
          END IF ! global%solverType

! ------------------------------------------------------------------------------
!       Mixture of thermally and calorically perfect gases
! ------------------------------------------------------------------------------     
      
        CASE ( GAS_MODEL_MIXT_TCPERF ) 
          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN
              IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_DUVWP ) THEN 
                CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
              END IF ! pRegion%mixt%cvState    

              DO icg = icgBeg,icgEnd
                r  = pCv(CV_MIXT_DENS,icg)
                p  = pCv(CV_MIXT_PRES,icg)

                T  = MixtPerf_HM_T_DGMP(r,refGamma,refM,p)

                pDv(DV_MIXT_TEMP,icg) = T
              END DO ! icg
            END IF ! pRegion%mixtInput%flowModel
          ELSE
            IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN 
              CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
            END IF ! pRegion%mixt%cvState    

            IF ( (global%pbaFlag .EQV. .FALSE.) .OR. &
                 (global%pbaBurnFlag .EQV. .FALSE.) ) THEN
! ----------- No program burn --------------------------------------------------            
              DO icg = icgBeg,icgEnd
                gc = MixtPerf_R_M(pGv(GV_MIXT_MOL,icg))
                g  = MixtPerf_G_CpR(pGv(GV_MIXT_CP,icg),gc)

                r  = pCv(CV_MIXT_DENS,icg)
                ir = 1.0_RFREAL/r
                Eo = ir*pCv(CV_MIXT_ENER,icg)

                Vm2 = ir*ir*(pCv(CV_MIXT_XMOM,icg)*pCv(CV_MIXT_XMOM,icg) + &
                             pCv(CV_MIXT_YMOM,icg)*pCv(CV_MIXT_YMOM,icg) + &
                             pCv(CV_MIXT_ZMOM,icg)*pCv(CV_MIXT_ZMOM,icg))

                r = r/(1.0_RFREAL - pRegion%plag%vFracE(1,indVFracE*icg))          
      
                pDv(DV_MIXT_PRES,icg) = MixtPerf_P_DEoGVm2(r,Eo,g,Vm2)
                pDv(DV_MIXT_TEMP,icg) = MixtPerf_T_DPR(r,pDv(DV_MIXT_PRES,icg),gc)
                pDv(DV_MIXT_SOUN,icg) = MixtPerf_C_GRT(g,gc,pDv(DV_MIXT_TEMP,icg))    
              END DO ! icg
            ELSE
#ifdef SPEC
! ----------- Program burn -----------------------------------------------------          
              IF ( global%specUsed .EQV. .FALSE. ) THEN
                CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) ! Defensive coding
              END IF ! global%specUsed

!             RFLU_SetDependentVars is called from rfluinit and hence spec%cv is
!             in conservative form. Need to convert to primitive form
              IF ( pRegion%spec%cvState == CV_MIXT_STATE_PRIM ) THEN 
                scalarConvFlag = .FALSE.
              ELSE
                scalarConvFlag = .TRUE.

                CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv, & 
                                                   pRegion%spec%cvState)
              END IF ! pRegion%spec%cvState    

              iCvSpecExplosive = SPEC_GetSpeciesIndex(global,pSpecInput, &
                                                      'EXPLOSIVE')
              iCvSpecProducts  = SPEC_GetSpeciesIndex(global,pSpecInput, &
                                                      'PRODUCTS')

              mwProducts = pSpecInput%specType(iCvSpecProducts)%pMaterial%molw
              cpProducts = pSpecInput%specType(iCvSpecProducts)%pMaterial%spht
              gcProducts = MixtPerf_R_M(mwProducts)
              gProducts  = MixtPerf_G_CpR(cpProducts,gcProducts)

              DO icg = icgBeg,icgEnd
                gc = MixtPerf_R_M(pGv(GV_MIXT_MOL,icg))
                g  = MixtPerf_G_CpR(pGv(GV_MIXT_CP,icg),gc)

                r  = pCv(CV_MIXT_DENS,icg)
                ir = 1.0_RFREAL/r

                Eo = ir*pCv(CV_MIXT_ENER,icg)

                Vm2 = ir*ir*(pCv(CV_MIXT_XMOM,icg)*pCv(CV_MIXT_XMOM,icg) + &
                             pCv(CV_MIXT_YMOM,icg)*pCv(CV_MIXT_YMOM,icg) + &
                             pCv(CV_MIXT_ZMOM,icg)*pCv(CV_MIXT_ZMOM,icg))

                r = r/(1.0_RFREAL - pRegion%plag%vFracE(1,indVFracE*icg))

                YExplosive = pCvSpec(iCvSpecExplosive,icg)
                YProducts  = pCvSpec(iCvSpecProducts,icg)

! DEBUG: Manoj-PBA1D, Notes: To match Jianghui's implementation
! 1. Set all cv and dv variables in pure explosive region
!    Changing it to set on dv variables
! 2. Use g mixture, Since Jianghui uses same material properties for products and explosive
!    Changing it in reaction zone to only use product g,gc
! 3. Use standard definition of speed of sound
!    Changing definition of sound speed to that given in paper
! END DEBUG                
                IF ( ABS(YExplosive) < nTol ) THEN
                  pDv(DV_MIXT_PRES,icg) = MixtPerf_P_DEoGVm2(r,Eo,g,Vm2)
                  pDv(DV_MIXT_TEMP,icg) = MixtPerf_T_DPR(r,pDv(DV_MIXT_PRES,icg),gc)
                  pDv(DV_MIXT_SOUN,icg) = MixtPerf_C_GRT(g,gc,pDv(DV_MIXT_TEMP,icg))    
                ELSEIF ( ABS(1.0_RFREAL-YExplosive) < nTol ) THEN
!                  pCv(CV_MIXT_DENS,icg) = pRegion%mixtInput%prepRealVal6
!                  pCv(CV_MIXT_XMOM,icg) = 0.0_RFREAL
!                  pCv(CV_MIXT_YMOM,icg) = 0.0_RFREAL
!                  pCv(CV_MIXT_ZMOM,icg) = 0.0_RFREAL
!                  pCv(CV_MIXT_ENER,icg) = ((pSpecInput%specType(iCvSpecExplosive)%pMaterial%detonVel &
!                                           + pRegion%mixtInput%prepRealVal7/pCv(CV_MIXT_DENS,icg)/ &
!                                            pSpecInput%specType(iCvSpecExplosive)%pMaterial%detonVel)**2/ &
!                                          (2.0_RFREAL*(g**2-1.0_RFREAL)))*pCv(CV_MIXT_DENS,icg)

                  pDv(DV_MIXT_PRES,icg) = pRegion%mixtInput%prepRealVal7
                  pDv(DV_MIXT_TEMP,icg) = MixtPerf_T_DPR(r,pDv(DV_MIXT_PRES,icg),gc)
                  pDv(DV_MIXT_SOUN,icg) = SQRT(pDv(DV_MIXT_PRES,icg)/pCv(CV_MIXT_DENS,icg))
!                  pDv(DV_MIXT_SOUN,icg) = MixtPerf_C_GRT(g,gc,pDv(DV_MIXT_TEMP,icg))    
                ELSE
                  pDv(DV_MIXT_PRES,icg) = YProducts*MixtPerf_P_DEoGVm2(r,Eo,gProducts,Vm2)
                  pDv(DV_MIXT_TEMP,icg) = MixtPerf_T_DPR(r,pDv(DV_MIXT_PRES,icg),gcProducts)
                  pDv(DV_MIXT_SOUN,icg) = SQRT((pDv(DV_MIXT_PRES,icg)/pCv(CV_MIXT_DENS,icg))* &
                                          (1.0_RFREAL+(gProducts-1.0_RFREAL)*YProducts))
                END IF ! YExplosive
              END DO ! icg

              IF ( scalarConvFlag .EQV. .TRUE. ) THEN 
                CALL RFLU_ScalarConvertCvPrim2Cons(pRegion,pRegion%spec%cv, & 
                                                   pRegion%spec%cvState)
              END IF ! scalarConvFlag
#else              
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) ! Defensive coding
#endif
            END IF ! pbaFlag 
          END IF ! global%solverType

! ------------------------------------------------------------------------------
!       Pseudogas
! ------------------------------------------------------------------------------     
      
        CASE ( GAS_MODEL_MIXT_PSEUDO ) 
#ifdef SPEC        
          IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN 
            CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
          END IF ! pRegion%mixt%cvState    

          DO icg = icgBeg,icgEnd
            mmm = pGv(GV_MIXT_MOL,icg)
            cpm = pGv(GV_MIXT_CP ,icg)            
            gcm = MixtPerf_R_M(mmm)
            gm  = MixtPerf_G_CpR(cpm,gcm)

            r  = pCv(CV_MIXT_DENS,icg)
            ir = 1.0_RFREAL/r
            Eo = ir*pCv(CV_MIXT_ENER,icg)

            cpg  = 0.0_RFREAL
            Yg   = 0.0_RFREAL
            phip = 0.0_RFREAL

            DO iSpec = 1,pRegion%specInput%nSpecies
              pSpecType => pRegion%specInput%specType(iSpec)

              rYi = pCvSpec(iSpec,icg)

              IF ( pSpecType%discreteFlag .EQV. .TRUE. ) THEN 
                phip = phip + rYi/pSpecType%pMaterial%dens
              ELSE 
                cpg  = cpg  + ir*rYi*pSpecType%pMaterial%spht 
                Yg   = Yg   + ir*rYi                                                              
              END IF ! pSpecType%discreteFlag
            END DO ! iSpec

            cpg = cpg/Yg
            
            Vm2 = ir*ir*(pCv(CV_MIXT_XMOM,icg)*pCv(CV_MIXT_XMOM,icg) + &
                         pCv(CV_MIXT_YMOM,icg)*pCv(CV_MIXT_YMOM,icg) + &
                         pCv(CV_MIXT_ZMOM,icg)*pCv(CV_MIXT_ZMOM,icg))

            term = MixtPerf_P_DEoGVm2(r,Eo,gm,Vm2)
            pDv(DV_MIXT_PRES,icg) = cpm/cpg*term/(1.0_RFREAL-phip)
            
            term = MixtPerf_T_DPR(r,pDv(DV_MIXT_PRES,icg),gcm)            
            pDv(DV_MIXT_TEMP,icg) = term*(1.0_RFREAL-phip)
            
            term = MixtPerf_C_GRT(gm,gcm,pDv(DV_MIXT_TEMP,icg))
            pDv(DV_MIXT_SOUN,icg) = term/(1.0_RFREAL-phip)
          END DO ! icg
#else
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) ! Defensive coding
#endif

! ------------------------------------------------------------------------------
!       Mixture of gas, liquid, and vapor
! ------------------------------------------------------------------------------     
      
        CASE ( GAS_MODEL_MIXT_GASLIQ ) 
#ifdef SPEC
          IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN 
            CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
          END IF ! pRegion%mixt%cvState  

          ro  = global%refDensityLiq
          po  = global%refPressLiq
          to  = global%refTempLiq
          Bp  = global%refBetaPLiq
          Bt  = global%refBetaTLiq
          cvl = global%refCvLiq

          gcg = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 gcg)

          gcv = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 gcv)

          DO icg = icgBeg,icgEnd
            r   = pCv(CV_MIXT_DENS,icg)
            ir  = 1.0_RFREAL/r
            Eo  = ir*pCv(CV_MIXT_ENER,icg)
            
            rYg = r*pCvSpec(1,icg)
            rYv = r*pCvSpec(2,icg) 
            rYl = r - rYg - rYv

            Vm2 = ir*ir*(pCv(CV_MIXT_XMOM,icg)*pCv(CV_MIXT_XMOM,icg) + &
                         pCv(CV_MIXT_YMOM,icg)*pCv(CV_MIXT_YMOM,icg) + &
                         pCv(CV_MIXT_ZMOM,icg)*pCv(CV_MIXT_ZMOM,icg))
            Cvm = (rYl*cvl + rYv*cvv + rYg*cvg)/r

            pDv(DV_MIXT_TEMP,icg) = MixtPerf_T_CvEoVm2(Cvm,Eo,Vm2) 

            Cl2 = MixtLiq_C2_Bp(Bp) 
            Cv2 = MixtPerf_C2_GRT(1.0_RFREAL,gcv,pDv(DV_MIXT_TEMP,icg))
            Cg2 = MixtPerf_C2_GRT(1.0_RFREAL,gcg,pDv(DV_MIXT_TEMP,icg))
            pDv(DV_MIXT_PRES,icg) = MixtGasLiq_P(rYl,rYv,rYg,Cl2,Cv2,Cg2,r, &
                                                 ro,po,to,Bp,Bt, &
                                                 pDv(DV_MIXT_TEMP,icg))

            rl = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,pDv(DV_MIXT_PRES,icg),po, &
                                        pDv(DV_MIXT_TEMP,icg),to)
            rv = MixtPerf_D_PRT(pDv(DV_MIXT_PRES,icg),gcv, &
                                pDv(DV_MIXT_TEMP,icg))
            rg = MixtPerf_D_PRT(pDv(DV_MIXT_PRES,icg),gcg, &
                                pDv(DV_MIXT_TEMP,icg))

            vfg = rYg/rg
            vfv = rYv/rv
            vfl = rYl/rl    

            Bl2 = -Bt/Bp
            Bv2 = rv*gcv
            Bg2 = rg*gcg

            pDv(DV_MIXT_SOUN,icg) = MixtGasLiq_C(Cvm,r,pDv(DV_MIXT_PRES,icg), &
                                                 rl,rv,rg,vfl,vfv,vfg,Cl2, &
                                                 Cv2,Cg2,Bl2,Bv2,Bg2)
          END DO ! icg
#else
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) ! Defensive coding
#endif

! ------------------------------------------------------------------------------
!       Mixture of thermally and calorically perfect gas and detonation products
! ------------------------------------------------------------------------------     
      
        CASE ( GAS_MODEL_MIXT_JWL ) 
#ifdef SPEC        
          IF ( global%specUsed .EQV. .TRUE. ) THEN
            IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN 
              CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
            END IF ! pRegion%mixt%cvState    

!           RFLU_SetDependentVars is called from rfluinit and hence spec%cv is
!           in conservative form. Need to convert to primitive form
            IF ( pRegion%spec%cvState == CV_MIXT_STATE_PRIM ) THEN 
              scalarConvFlag = .FALSE.
            ELSE
              scalarConvFlag = .TRUE.

              CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv, & 
                                                 pRegion%spec%cvState)
            END IF ! pRegion%spec%cvState    

            iCvSpecAir = SPEC_GetSpeciesIndex(global,pSpecInput,'AIR')
            mwAir = pSpecInput%specType(iCvSpecAir)%pMaterial%molw
            cpAir = pSpecInput%specType(iCvSpecAir)%pMaterial%spht
            gcAir = MixtPerf_R_M(mwAir)
            gAir  = MixtPerf_G_CpR(cpAir,gcAir)

            iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput, &
                                                   'PRODUCTS')

            DO icg = icgBeg,icgEnd
              gc = MixtPerf_R_M(pGv(GV_MIXT_MOL,icg))
              g  = MixtPerf_G_CpR(pGv(GV_MIXT_CP,icg),gc)

              r  = pCv(CV_MIXT_DENS,icg)
              ir = 1.0_RFREAL/r
              Eo = ir*pCv(CV_MIXT_ENER,icg)

              Vm2 = ir*ir*(pCv(CV_MIXT_XMOM,icg)*pCv(CV_MIXT_XMOM,icg) + &
                           pCv(CV_MIXT_YMOM,icg)*pCv(CV_MIXT_YMOM,icg) + &
                           pCv(CV_MIXT_ZMOM,icg)*pCv(CV_MIXT_ZMOM,icg))

              e  = Eo - 0.5_RFREAL*Vm2

              YProducts = pCvSpec(iCvSpecProducts,icg)

              !CALL RFLU_JWL_ComputePressureMixt(gAir,gcAir,e,r,YProducts,a, &
              !                                  eJWL,ePerf,p,T)
              CALL RFLU_JWL_ComputePressureMixt(pRegion,icg,gAir,gcAir,e,r,YProducts,a, &
                                                eJWL,ePerf,p,T)
              pDv(DV_MIXT_PRES,icg)  = p
              pDv(DV_MIXT_TEMP,icg)  = T
              pDv(DV_MIXT_SOUN,icg)  = a
              pDv(DV_MIXT_EJWL,icg)  = eJWL
              pDv(DV_MIXT_EPERF,icg) = ePerf
            END DO ! icg

            IF ( scalarConvFlag .EQV. .TRUE. ) THEN 
              CALL RFLU_ScalarConvertCvPrim2Cons(pRegion,pRegion%spec%cv, & 
                                                 pRegion%spec%cvState)
            END IF ! scalarConvFlag
          END IF ! global%specUsed
#else
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) ! Defensive coding
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

END SUBROUTINE RFLU_SetDependentVars

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_SetDependentVars.F90,v $
! Revision 1.2  2016/02/04 19:59:23  fred
! Adding iterative JWL EOS capabilities for the cylindrical detonation case
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:36  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:16:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/11/28 23:05:03  mparmar
! Added computation of non-dimensional temperature
!
! Revision 1.1  2007/04/09 18:48:53  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.9  2006/08/16 19:16:39  haselbac
! Bug fix: bad check-in, declarations missing
!
! Revision 1.8  2006/08/15 15:23:19  haselbac
! Bug fix: dv was computed incorrectly for pseudo-gas
!
! Revision 1.7  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.6  2006/03/26 20:21:40  haselbac
! Added support for GL model
!
! Revision 1.5  2005/11/14 16:55:17  haselbac
! Added support for pseudo-gas model
!
! Revision 1.4  2005/11/10 02:14:49  haselbac
! Added SELECT CASE on gasModel
!
! Revision 1.3  2005/04/15 15:06:20  haselbac
! Added range arguments, removed pGrid declaration
!
! Revision 1.2  2004/11/14 19:41:20  haselbac
! Added setting of dv for incompressible fluid model
!
! Revision 1.1  2004/10/19 19:23:53  haselbac
! Initial revision
!
! ******************************************************************************

