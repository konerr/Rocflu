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
! Purpose: Set transport variables.
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
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_SetTransportVars.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_SetTransportVars(pRegion,icgBeg,icgEnd)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid

  USE RFLU_ModHouMahesh, ONLY: RFLU_HM_ConvTvD2ND
                                 
#ifdef SPEC
  USE ModSpecies, ONLY: t_spec_type
#endif

! Subbu - Modules for JWL EOS
#ifdef SPEC
  USE ModInterfaces, ONLY: MixtPerf_R_M, &
                           MixtPerf_G_CpR
  USE ModSpecies, ONLY: t_spec_input
  USE ModInterfacesSpecies, ONLY: SPEC_GetSpeciesIndex
  USE SPEC_RFLU_ModPBAUtils, ONLY: SPEC_RFLU_PBA_GetSolution
#endif
! Subbu - End modules

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
  INTEGER :: icg,indCp
  REAL(RFREAL) :: absTemp,iPrLam,iPrLamLiq,kLiq,kGas,kVap,muLiq,muGas,muVap,&
                  refTemp,refVisc,refViscLiq,suthCoef,s1,s2,s3,s4,term,y,z
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pTv
#ifdef SPEC
  INTEGER :: iSpec,jSpec
  REAL(RFREAL) :: Cond,condi,condj,cp,cp2,molw,molw2,mui,muj,phiij,phiijCond, &
                  pr,pr2,suthCoef2,suthTemp,suthTemp2,refVisc2,termCond,Tmixt, &
                  Visc,xi,xj,xjphiij
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
  TYPE(t_spec_type), POINTER :: pSpecType,pSpecType2
#endif

  TYPE(t_global), POINTER :: global

! Subbu - Vars for JWL EOS
#ifdef SPEC
  INTEGER :: i,iCvSpecAir,iCvSpecProducts,j  
  REAL(RFREAL) :: cpMixt,gammaMixt,mlfAir,mlfCOTNT,mlfCO2TNT,mlfH2OTNT,mlfN2TNT, & 
                  mlfTNT,mlfO2TNT,mlfN2Air,mlfO2Air,mwAir,mwMixt,mwProd,RMixt,RSpec, &
                  termPhiTimesMlf,YProducts
  REAL(RFREAL), DIMENSION(5) :: shCp,epsk,k,mlf,mu,mw,omega,sigma
  REAL(RFREAL), DIMENSION(5,5) :: curfitCoef,epskAB,Phi,omegaAB,sigmaAB
  TYPE(t_spec_input), POINTER :: pSpecInput
#endif
! Subbu - End Vars

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_SetTransportVars.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SetTransportVars',__FILE__)

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************
 
  pCv => pRegion%mixt%cv
  pDv => pRegion%mixt%dv 
  pGv => pRegion%mixt%gv 
  pTv => pRegion%mixt%tv

#ifdef SPEC
  pCvSpec => pRegion%spec%cv
  pSpecInput => pRegion%specInput
#endif

  indCp = pRegion%mixtInput%indCp

  iPrLam   = 1.0_RFREAL/pRegion%mixtInput%prLam
  refVisc  = pRegion%mixtInput%refVisc
  refTemp  = pRegion%mixtInput%refViscTemp
  suthCoef = pRegion%mixtInput%suthCoef

! ******************************************************************************  
! Compute dependent variables
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%fluidModel )
  
! ==============================================================================  
!   Incompressible fluid model
! ==============================================================================  
   
    CASE ( FLUID_MODEL_INCOMP ) 

! ==============================================================================  
!   Compressible fluid model
! ==============================================================================  

    CASE ( FLUID_MODEL_COMP )             
      SELECT CASE ( pRegion%mixtInput%gasModel ) 
      
! ------------------------------------------------------------------------------
!       Thermally and calorically perfect or thermally perfect gas
! ------------------------------------------------------------------------------      
            
        CASE ( GAS_MODEL_TCPERF, & 
               GAS_MODEL_TPERF)
          SELECT CASE ( pRegion%mixtInput%viscModel ) 
      
! -------- Sutherland model ----------------------------------------------------
      
            CASE ( VISC_SUTHR )
              DO icg = icgBeg,icgEnd
                term = SQRT(pDv(DV_MIXT_TEMP,icg)/refTemp) & 
                      *(1.0_RFREAL + suthCoef/refTemp)/ & 
                       (1.0_RFREAL + suthCoef/pDv(DV_MIXT_TEMP,icg))

                pTv(TV_MIXT_MUEL,icg) = term*refVisc
                pTv(TV_MIXT_TCOL,icg) = pGv(GV_MIXT_CP  ,icg*indCp)* & 
                                        pTv(TV_MIXT_MUEL,icg      )*iPrLam
              END DO ! icg

! --------- Constant viscosity model -------------------------------------------
                
            CASE ( VISC_FIXED ) 
              DO icg = icgBeg,icgEnd
                pTv(TV_MIXT_MUEL,icg) = refVisc
                pTv(TV_MIXT_TCOL,icg) = pGv(GV_MIXT_CP  ,icg*indCp) & 
                                       *pTv(TV_MIXT_MUEL,icg      )*iPrLam
              END DO ! icg           
        
! --------- Antibes model ------------------------------------------------------
              
            CASE ( VISC_ANTIB ) 
              s1 = 110.0_RFREAL
              s2 = 120.0_RFREAL
              s3 = 230.0_RFREAL
              s4 = 1.0_RFREAL/refTemp*SQRT(s2/refTemp)*(refTemp+s1)/s3

              IF ( refTemp <= s2 ) THEN
                DO icg = icgBeg,icgEnd
                  absTemp = ABS(pDv(DV_MIXT_TEMP,icg))

                  IF ( absTemp < s2 ) THEN
                    term = absTemp/refTemp
                  ELSE
                    term = SQRT(absTemp/s2)*(absTemp/refTemp)*(s3/(absTemp+s1))
                  END IF ! absTemp

                  pTv(TV_MIXT_MUEL,icg) = refVisc*term
                  pTv(TV_MIXT_TCOL,icg) = pGv(GV_MIXT_CP  ,icg*indCp) & 
                                         *pTv(TV_MIXT_MUEL,icg      )*iPrLam
                END DO ! icg
              ELSE 
                DO icg = icgBeg,icgEnd
                  absTemp = ABS(pDv(DV_MIXT_TEMP,icg))

                  IF ( absTemp < s2 ) THEN
                    term = s4
                  ELSE
                    term = SQRT(absTemp/refTemp)*(absTemp/refTemp) & 
                          *(refTemp+s1)/(absTemp+s1)
                  END IF ! absTemp

                  pTv(TV_MIXT_MUEL,icg) = term*refVisc
                  pTv(TV_MIXT_TCOL,icg) = pGv(GV_MIXT_CP  ,icg*indCp) & 
                                         *pTv(TV_MIXT_MUEL,icg      )*iPrLam
                END DO ! icg      
              END IF ! refTemp         

! --------- Default ------------------------------------------------------------

            CASE DEFAULT 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pRegion%mixtInput%viscModel
      
! ------------------------------------------------------------------------------
!       Mixture of thermally and calorically perfect gases
! ------------------------------------------------------------------------------     

        CASE ( GAS_MODEL_MIXT_TCPERF )
#ifdef SPEC
          IF ( pRegion%spec%cvState /= CV_MIXT_STATE_PRIM ) THEN
            CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
          END IF ! pRegion%spec%cvState  

          SELECT CASE ( pRegion%mixtInput%viscModel )

! TEMPORARY: Manoj: 2012-05-10: Added temporarily to work with constant visc          
! --------- Constant viscosity model -------------------------------------------
                
            CASE ( VISC_FIXED ) 
              DO icg = icgBeg,icgEnd
                pTv(TV_MIXT_MUEL,icg) = refVisc
                pTv(TV_MIXT_TCOL,icg) = pGv(GV_MIXT_CP  ,icg*indCp) & 
                                       *pTv(TV_MIXT_MUEL,icg      )*iPrLam
              END DO ! icg    
! END TEMPORARY: Manoj: 2012-05-10: Added temporarily to work with constant visc          

! -------- Wilke Sutherland model ----------------------------------------------

            CASE ( VISC_WILKE_SUTH )
              DO icg = icgBeg,icgEnd
                Tmixt = pDv(DV_MIXT_TEMP,icg)
                Visc = 0.0_RFREAL
                Cond = 0.0_RFREAL

                DO iSpec = 1,pRegion%specInput%nSpecies
                  pSpecType => pRegion%specInput%specType(iSpec)

                  xi = pGv(GV_MIXT_MOL,icg)*pCvSpec(iSpec,icg) &
                       /pSpecType%pMaterial%molw

                  refVisc  = pSpecType%pMaterial%refVisc
                  suthTemp = pSpecType%pMaterial%suthTemp
                  suthCoef = pSpecType%pMaterial%suthCoef
                  molw     = pSpecType%pMaterial%molw
                  pr       = pSpecType%pMaterial%pr
                  cp       = pSpecType%pMaterial%spht

                  mui = refVisc*SQRT(Tmixt/suthTemp) &
                        *(1.0_RFREAL + suthCoef/suthTemp)/ &
                         (1.0_RFREAL + suthCoef/Tmixt)

                  condi = mui*cp/pr

                  term     = 0.0_RFREAL
                  termCond = 0.0_RFREAL

                  DO jSpec = 1,pRegion%specInput%nSpecies
                    pSpecType2 => pRegion%specInput%specType(jSpec)

                    xj = pGv(GV_MIXT_MOL,icg)*pCvSpec(jSpec,icg) &
                         /pSpecType2%pMaterial%molw

                    refVisc2  = pSpecType2%pMaterial%refVisc
                    suthTemp2 = pSpecType2%pMaterial%suthTemp
                    suthCoef2 = pSpecType2%pMaterial%suthCoef
                    molw2     = pSpecType2%pMaterial%molw
                    pr2       = pSpecType%pMaterial%pr
                    cp2       = pSpecType%pMaterial%spht

                    muj = refVisc2*SQRT(Tmixt/suthTemp2) &
                          *(1.0_RFREAL + suthCoef2/suthTemp2)/ &
                           (1.0_RFREAL + suthCoef2/Tmixt)

                    condj = muj*cp2/pr2

                    phiij = ((8.0_RFREAL*(1.0_RFREAL+molw/molw2) &
                             )**(-0.5_RFREAL)) &
                           *((1.0_RFREAL+((mui/muj)**0.50_RFREAL) &
                                        *((molw2/molw)**0.25_RFREAL) &
                             )**2.0_RFREAL)

! TEMPORARY: Using the formula which reduces to unity for single species
!                    phiijCond =
!                    1.065_RFREAL*((8.0_RFREAL*(1.0_RFREAL+molw/molw2) &
!                                              )**(-0.5_RFREAL)) &
!                              *((1.0_RFREAL+((condi/condj)**0.50_RFREAL) &
!                                           *((molw2/molw)**0.25_RFREAL) &
!                               )**2.0_RFREAL)
                    phiijCond = ((8.0_RFREAL*(1.0_RFREAL+molw/molw2) &
                                 )**(-0.5_RFREAL)) &
                              *((1.0_RFREAL+((condi/condj)**0.50_RFREAL) &
                                           *((molw2/molw)**0.25_RFREAL) &
                               )**2.0_RFREAL)
! END TEMPORARY

                    term     = term     + xj*phiij
                    termCond = termCond + xj*phiijCond
                  END DO ! jSpec

                  Visc = Visc + xi*mui/term
                  Cond = Cond + xi*condi/termCond
                END DO ! iSpec

                pTv(TV_MIXT_MUEL,icg) = Visc
                pTv(TV_MIXT_TCOL,icg) = Cond
              END DO ! icg

! --------- Default ------------------------------------------------------------

            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pRegion%mixtInput%viscModel
#else
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) ! Defensive coding
#endif     

! ------------------------------------------------------------------------------
!       Gas-liquid mixture fluid model
! ------------------------------------------------------------------------------

        CASE ( GAS_MODEL_MIXT_GASLIQ )
#ifdef SPEC        
          DO icg = icgBeg,icgEnd

! --------- Use Sutherland formula for gas -------------------------------------

            term = SQRT(pDv(DV_MIXT_TEMP,icg)/refTemp) & 
                   *(1.0_RFREAL + suthCoef/refTemp)/ & 
                    (1.0_RFREAL + suthCoef/pDv(DV_MIXT_TEMP,icg))

            muGas = term*refVisc
            kGas  = pRegion%specInput%specType(1)%pMaterial%spht*muGas*iPrLam
   
! --------- Use Sutherland formula for vapor -----------------------------------

            term = SQRT(pDv(DV_MIXT_TEMP,icg)/refTemp) & 
                   *(1.0_RFREAL + suthCoef/refTemp)/ & 
                    (1.0_RFREAL + suthCoef/pDv(DV_MIXT_TEMP,icg))

            muVap = term*refVisc   
            kVap  = pRegion%specInput%specType(2)%pMaterial%spht*muVap*iPrLam

! --------- Compute liquid viscosity and conductivity --------------------------

            refViscLiq = 1.0E-02_RFREAL ! NOTE hard-coded for now
            z = refTemp/pDv(DV_MIXT_TEMP,icg)
            y = -1.704_RFREAL -5.306_RFREAL*z + 7.003_RFREAL*z*z

            muLiq = refViscLiq*EXP(y) 
            kliq  = 0.585_RFREAL

! --------- Compute mixture viscosity and conductivity --------------------------

            pTv(TV_MIXT_MUEL,icg) = & 
              ((pCv(CV_MIXT_DENS,icg)- pCvSpec(1,icg) &
              - pCvSpec(2,icg))*muLiq + pCvSpec(1,icg)*muGas & 
              + pCvSpec(2,icg)*muVap)/pCv(CV_MIXT_DENS,icg)
            pTv(TV_MIXT_TCOL,icg) = & 
              ((pCv(CV_MIXT_DENS,icg)- pCvSpec(1,icg) &
              - pCvSpec(2,icg))*kLiq + pCvSpec(1,icg)*kGas & 
              + pCvSpec(2,icg)*kVap)/pCv(CV_MIXT_DENS,icg)
          END DO ! icg
#else
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
#endif          

! Subbu - Added JWL case
! ------------------------------------------------------------------------------
!       Mixture of detonation products and air
! ------------------------------------------------------------------------------
        CASE ( GAS_MODEL_MIXT_JWL )
#ifdef SPEC 
          DO icg = icgBeg,icgEnd
            Tmixt = pDv(DV_MIXT_TEMP,icg)

            iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput, &
                                                       'PRODUCTS')
            pSpecType => pSpecInput%specType(iCvSpecProducts)
            YProducts = pCvSpec(iCvSpecProducts,icg)/pCv(CV_MIXT_DENS,icg)
            mwProd  = pSpecType%pMaterial%molw

            iCvSpecAir = SPEC_GetSpeciesIndex(global,pSpecInput, &
                                                       'AIR')
            pSpecType => pSpecInput%specType(iCvSpecAir)
            mwAir  = pSpecType%pMaterial%molw
 
! ======================================================================================
! Compute mole fraction for each component
! ======================================================================================
            mlfTNT = YProducts/mwProd/(YProducts/mwProd + (1.0_RFREAL-YProducts)/mwAir)
            mlfAir = 1.0_RFREAL - mlfTNT
         
            mlfN2TNT  = pRegion%mixtInput%prepRealVal8
            mlfO2TNT  = pRegion%mixtInput%prepRealVal9
            mlfCOTNT  = pRegion%mixtInput%prepRealVal10
            mlfCO2TNT = pRegion%mixtInput%prepRealVal11
            mlfH2OTNT = pRegion%mixtInput%prepRealVal12

            mlfN2Air  = 0.79_RFREAL
            mlfO2Air  = 0.21_RFREAL

            mlf(1) = mlfTNT*mlfN2TNT  + mlfAir*mlfN2Air   ! N2
            mlf(2) = mlfTNT*mlfO2TNT  + mlfAir*mlfO2Air   ! O2
            mlf(3) = mlfTNT*mlfCOTNT                      ! C0
            mlf(4) = mlfTNT*mlfCO2TNT                     ! C02
            mlf(5) = mlfTNT*mlfH2OTNT                     ! H20

! =========================================================================================
! Define constants for each component and Curve fit coeffecients
! =========================================================================================
            mw(1:5) = (/28.013_RFREAL, 31.999_RFREAL, 28.010_RFREAL, 44.010_RFREAL, 18.015_RFREAL /)

            ! The following properties are given by BSL (2nd.Ed.) Appendix E, pp.864-866:
            ! Lennard-Jones radius [Angstrom] of [N2; O2; CO; CO2; H2O]
            sigma(1:5) = (/3.667_RFREAL, 3.433_RFREAL, 3.590_RFREAL, 3.996_RFREAL, 2.605_RFREAL /)
            ! Lennard-Jones temp. normalization (epsilon/kappa) [K] of [N2; O2; CO; CO2; H2O]
            epsk (1:5) = (/99.8_RFREAL, 113.0_RFREAL, 110.0_RFREAL, 190.0_RFREAL, 572.0_RFREAL /)
            omega(1:5) = 1.16145_RFREAL/(Tmixt/epsk(:))**0.14874_RFREAL   &
                       + 0.52487_RFREAL/exp(0.77320_RFREAL*Tmixt/epsk(:)) &
                       + 2.16178_RFREAL/exp(2.43787_RFREAL*Tmixt/epsk(:))

            ! pairwise-combined Lennard-Jones radius [Angstrom] and temp [K] (BSL 2nd Ed. p.527)
            DO i = 1,5
              DO j = 1,5
                sigmaAB(i,j) = 0.5_RFREAL*(sigma(i) + sigma(j))
                epskAB (i,j) = SQRT(epsk(i)*epsk(j))
                omegaAB(i,j) = 1.06036_RFREAL/(Tmixt/epskAB(i,j))**0.15610_RFREAL   &
                             + 0.19300_RFREAL/exp(0.47635_RFREAL*Tmixt/epskAB(i,j)) &
                             + 1.03587_RFREAL/exp(1.52996_RFREAL*Tmixt/epskAB(i,j)) &
                             + 1.76474_RFREAL/exp(3.89411_RFREAL*Tmixt/epskAB(i,j))
              END DO ! j
            END DO ! i

            ! Shomate equation fit coefficients (NIST)

            curfitCoef(1,:) = (/26.09200,  8.218801, -1.976141,  0.159274,  0.044434/)
            curfitCoef(2,:) = (/29.65900,  6.137261, -1.186521,  0.095780, -0.219663/)
            IF ( Tmixt <= 1300.0_RFREAL ) THEN
              curfitCoef(3,:) = (/25.56759,  6.096130,  4.054656, -2.671301,  0.131021/)
            ELSE
              curfitCoef(3,:) = (/35.15070,  1.300095, -0.205921,  0.013550, -3.282780/)
            END IF 

            IF ( Tmixt <= 1200.0_RFREAL ) THEN
              curfitCoef(4,:) = (/24.99735,  55.18696, -33.69137,  7.948387, -0.136638/)
            ELSE
              curfitCoef(4,:) = (/58.16639,  2.720074, -0.492289,  0.038844, -6.447293/)
            END IF

            IF ( Tmixt <= 1700.0_RFREAL ) THEN
              curfitCoef(5,:) = (/30.09200,  6.832514,  6.793435, -2.534480,  0.082139/)
            ELSE
              curfitCoef(5,:) = (/41.96426,  8.622053, -1.499780,  0.098119, -11.15764/)
            END IF 

! ==============================================================================
! Compute Gas properties 
! ==============================================================================
            mwMixt = 0.0_RFREAL
            DO i = 1,5
              mwMixt = mwMixt + mlf(i)*mw(i)
            END DO ! i 
           
            ! Heat capacity [J/mol*K] of a pure-species, low pressure (ideal) gas - Shomate equation fit:
            DO i = 1,5
              shCp(i) = curfitCoef(i,1)                                   &
                      + curfitCoef(i,2)* Tmixt/1000.0_RFREAL                 &
                      + curfitCoef(i,3)*(Tmixt/1000.0_RFREAL)**( 2.0_RFREAL) &
                      + curfitCoef(i,4)*(Tmixt/1000.0_RFREAL)**( 3.0_RFREAL) &
                      + curfitCoef(i,5)*(Tmixt/1000.0_RFREAL)**(-2.0_RFREAL)
            END DO  !i
            
            ! Heat capacity [J/mol*K] of a low pressure (ideal_mix = sum(x.*k./sum(Phi*x));l) gas mixture:
            cpMixt = 0.0_RFREAL
            DO i = 1,5
              cpMixt = cpMixt + mlf(i)*shCp(i)
            END DO ! i

! ************* Currently the gas variables cp & molwt are calculated here. Ideally
! this needs to be done in the subroutine RFLU_SetGasVars. ********************** !

            pGv(GV_MIXT_MOL,icg) = mwMixt
            pGv(GV_MIXT_CP ,icg) = cpMixt*1000.0_RFREAL/mwMixt

            RMixt = MixtPerf_R_M(mwMixt) ! Rmixt in J/kgK
            ! Convert cpMixt from J/molK to J/kgK
            gammaMixt = MixtPerf_G_CpR(cpMixt*1000.0_RFREAL/mwMixt,RMixt)

            ! Viscosity [10*g/cm*s = 10*Poise = kg/m*s] of a pure-species, 
            ! low-pressure gas-kinetic theory (BSL 2nd.Ed. p.26):
            DO i = 1,5
              mu(i) = 2.6693E-6_RFREAL*SQRT(mw(i)*Tmixt)/(sigma(i)**2.0_RFREAL*omega(i))
            END DO !i
            DO i = 1,5
              Rspec = MixtPerf_R_M(mw(i))
              ! Convert shCp from J/molK to J/kgK
              ! Thermal conductivity [W/m*K] of a prure-species, polyatomic,
              ! low-pressure gas - Eucken approximation (BSL 2nd Ed. p.276):
              k(i)  = (shCp(i)*1000.0_RFREAL/mw(i) + 1.25_RFREAL*Rspec)*mu(i)
            END DO !
            DO i = 1,5
              DO j = 1,5
                Phi(i,j) = 1.0_RFREAL/SQRT(8.0_RFREAL*(1.0_RFREAL+mw(i)/mw(j))) &
                         *(1.0_RFREAL + SQRT(mu(i)/mu(j))*(mw(j)/mw(i))**(0.25_RFREAL))**2.0_RFREAL ;
              END DO !j
            END DO ! i
 
            ! Viscosity [kg/m*s] of a low-pressure gas mixure - Wilke formula (BSL 2nd Ed.p.27)
            ! Thermal conductivity [W/m*K] of a low-pressure gas mixture - Mason & Saxena
            ! formula (BSL 2nd. Ed. p.276):

            Visc = 0.0_RFREAL
            Cond = 0.0_RFREAL
            DO i = 1,5
              termPhiTimesMlf = 0.0_RFREAL
              DO j = 1,5
                termPhiTimesMlf = termPhiTimesMlf + Phi(i,j)*mlf(j)
              END DO ! j
              Visc = Visc + mlf(i)*mu(i)/termPhiTimesMlf
              Cond = Cond + mlf(i)*k(i) /termPhiTimesMlf
            END DO ! i

            pTv(TV_MIXT_MUEL,icg) = Visc
            pTv(TV_MIXT_TCOL,icg) = Cond
              
          END DO

#else
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
#endif          
! Subbu - End JWL case
            
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

! ==============================================================================
! Non-dimensionalize transport variable for SOLV_IMPLICIT_HM
! ==============================================================================

  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
    IF ( pRegion%mixtInput%computeTv ) THEN
      CALL RFLU_HM_ConvTvD2ND(pRegion)
    END IF ! pRegion%mixtInput%computeTv
  END IF

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SetTransportVars

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_SetTransportVars.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.6  2009/07/09 20:44:02  mparmar
! Added support for GAS_MODEL_MIXT_TCPERF
!
! Revision 1.5  2009/05/02 23:45:24  mparmar
! Corrected usage of suthCoef and refViscTemp in Sutherland formula
!
! Revision 1.4  2008/12/06 08:43:36  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:16:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/11/28 23:05:04  mparmar
! Added non-dimensionalization of tv for SOLV_IMPLICIT_HM
!
! Revision 1.1  2007/04/09 18:48:53  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:50  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.3  2006/03/26 20:21:45  haselbac
! Added support for GL model
!
! Revision 1.2  2005/04/15 15:06:22  haselbac
! Added range arguments, removed pGrid declaration
!
! Revision 1.1  2004/11/14 20:02:49  haselbac
! Initial revision
!
! ******************************************************************************

