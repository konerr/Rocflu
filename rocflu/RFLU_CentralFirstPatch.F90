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
! Purpose: Compute central convective fluxes using first-order accurate
!   approximation through a patch by using an average of variables.
!
! Description: None.
!
! Input:
!   pRegion      Pointer to region
!   pPatch       Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_CentralFirstPatch.F90,v 1.5 2016/04/06 19:15:57 rahul Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_CentralFirstPatch(pRegion,pPatch)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters
  USE ModMPI

  USE RFLU_ModGridSpeedUtils, ONLY: RFLU_DescaleGridSpeed
  USE RFLU_ModRindStates, ONLY: RFLU_SetRindStateFarfieldPerf, &
                                RFLU_SetRindStateInjectPerf, & 
                                RFLU_SetRindStateSlipWallPerf, &
                                RFLU_SetRindStateSlipWallJWL, &
                                RFLU_SetRindStateSlipWallPerf1

  USE RFLU_ModRoeFlux, ONLY: RFLU_ROE_FluxFunction

#ifdef SPEC
   USE ModInterfacesSpecies, ONLY: SPEC_GetSpeciesIndex
   USE ModSpecies, ONLY: t_spec_type,&
                         t_spec_input
#endif

  USE ModInterfaces, ONLY: BcondFarfieldPerf, &
                           BcondInflowPerf, &
                           BcondOutflowPerf, &
                           MixtPerf_D_PRT, &                           
                           MixtPerf_G_CpR, &
                           MixtPerf_Eo_DGPUVW, &
                           MixtPerf_P_DEoGVm2, &
                           MixtPerf_P_DRT, &                           
                           MixtPerf_R_M, &
                           ReflectVector, &
                           RFLU_DecidePrint, &
                           RFLU_PrintLocInfo


  USE RFLU_ModConvertCv,   ONLY: RFLU_ScalarConvertCvCons2Prim, &
                                 RFLU_ScalarConvertCvPrim2Cons
  USE RFLU_ModJWL

#ifdef PLAG
  USE ModPartLag, ONLY: t_plag,t_plag_input
#endif

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

  LOGICAL :: corrFlag,decidePrintFlag
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER, PARAMETER :: MAX_INOUTFLOW_LOCS = 10
  INTEGER :: c1,bcOptFixed,bcOptType,distrib,gasModel,ifl,indCp,indGs,indMf, &
             indSd,indMol,nLocs
  INTEGER :: loc(MAX_INOUTFLOW_LOCS,MIN_VAL:MAX_VAL)
  REAL(RFREAL) :: aoa,aos,betah,betav,corr,cp,er,fs,fsu,g,h,Hl,iCmassRef, &
                  iCpRef,irl,irr,liftCoef,mach,mf,minj,mm,nm,nx,ny,nz,pa,pl, &
                  pf,pr,pRef,ptot,ql,qr,rgas,rhoa,rhoea,rhoua,rhova,rhowa,rl, &
                  rel,rer,rul,rur,rvl,rvr,rwl,rwr,rr,rRef,term,tf,tinj,tl,tr, &
                  ttot,uinj,ul,ur,vcont,vinj,vl,vm2,vr,vRef,winj,wl,wr,xc,yc,zc

  REAL(RFREAL) :: pg,ratio,rg,vFracGl
  REAL(RFREAL) :: flx(5)
  REAL(RFREAL), DIMENSION(:), POINTER :: mfMixt
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,dv,gv,rhs,vals,sd,pCv,pDv,pGv
  TYPE(t_global), POINTER :: global

#ifdef SPEC
  LOGICAL :: scalarConvFlag
  INTEGER :: iSpec,iCvSpecProducts,iCvSpecAir
  REAL(RFREAL) :: imm,Yproducts
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
  TYPE(t_spec_type), POINTER :: pSpecType
  TYPE(t_spec_input), POINTER :: pSpecInput
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_CentralFirstPatch.F90,v $ $Revision: 1.5 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_CentralFirstPatch',__FILE__)

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  indCp    = pRegion%mixtInput%indCp
  indMf    = pRegion%mixtInput%indMfMixt
  indSd    = pRegion%mixtInput%indSd
  indMol   = pRegion%mixtInput%indMol
  gasModel = pRegion%mixtInput%gasModel

  indGs = pRegion%grid%indGs

  distrib = pPatch%mixt%distrib

  cv  => pRegion%mixt%cv
  dv  => pRegion%mixt%dv
  gv  => pRegion%mixt%gv
  pCv  => pRegion%mixt%cv
  pDv  => pRegion%mixt%dv
  pGv  => pRegion%mixt%gv
  rhs => pRegion%mixt%rhs
  sd  => pRegion%mixt%sd

  mfMixt => pPatch%mfMixt

  nLocs = 0

  decidePrintFlag = RFLU_DecidePrint(global)

  pRef = global%refPressure
  rRef = global%refDensity
  vRef = global%refVelocity

  iCpRef    = 2.0_RFREAL/(rRef*vRef*vRef)
  iCmassRef = 1.0_RFREAL/(rRef*vRef)

#ifdef SPEC 
  pCvSpec => pRegion%spec%cv
  pSpecInput => pRegion%specInput
#endif

! ******************************************************************************
! Select boundary type
! ******************************************************************************

  SELECT CASE ( pPatch%bcType )

! ==============================================================================
!   Inflow (based on total quantities and flow angles)
! ==============================================================================

    CASE ( BC_INFLOW_TOTANG )
      vals       => pPatch%mixt%vals
      bcOptType  =  pPatch%mixt%switches(BCSWI_INFLOW_TYPE)
      bcOptFixed =  pPatch%mixt%switches(BCSWI_INFLOW_FIXED)

      DO ifl = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifl)

        nx = pPatch%fn(XCOORD,ifl)
        ny = pPatch%fn(YCOORD,ifl)
        nz = pPatch%fn(ZCOORD,ifl)
        nm = pPatch%fn(XYZMAG,ifl)

        fs = pPatch%gs(indGs*ifl)

        rl  = cv(CV_MIXT_DENS,c1)
        rul = cv(CV_MIXT_XMOM,c1)
        rvl = cv(CV_MIXT_YMOM,c1)
        rwl = cv(CV_MIXT_ZMOM,c1)
        rel = cv(CV_MIXT_ENER,c1)
        pl  = dv(DV_MIXT_PRES,c1)

        ptot  = vals(BCDAT_INFLOW_PTOT, distrib*ifl)
        ttot  = vals(BCDAT_INFLOW_TTOT, distrib*ifl)
        betah = vals(BCDAT_INFLOW_BETAH,distrib*ifl)
        betav = vals(BCDAT_INFLOW_BETAV,distrib*ifl)

        IF ( bcOptType /= BCOPT_SUBSONIC ) THEN
          mach = vals(BCDAT_INFLOW_MACH,distrib*ifl)
        ELSE
          mach = 0.0_RFREAL
        END IF ! bcOptType

        IF ( gasModel == GAS_MODEL_TCPERF ) THEN
          cp = gv(GV_MIXT_CP ,indCp *c1)
          mm = gv(GV_MIXT_MOL,indMol*c1)

          CALL BcondInflowPerf(bcOptType,bcOptFixed,ptot,ttot,betah,betav, &
                               mach,nx,ny,nz,cp,mm,rl,rul,rvl,rwl,rr,rur, &
                               rvr,rwr,rer,pr)
        ELSE 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! gasModel

        ql = (rul*nx + rvl*ny + rwl*nz)/rl - fs
        qr = (rur*nx + rvr*ny + rwr*nz)/rr - fs

        ul = rul/rl
        vl = rvl/rl
        wl = rwl/rl
        ur = rur/rr
        vr = rvr/rr
        wr = rwr/rr

        vFracGl = 1.0_RFREAL

        flx(1) = 0.5_RFREAL*(ql* rl                + qr* rr               )*nm
        flx(2) = 0.5_RFREAL*(ql* rul       + pl*nx + qr* rur       + pr*nx)*nm
        flx(3) = 0.5_RFREAL*(ql* rvl       + pl*ny + qr* rvr       + pr*ny)*nm
        flx(4) = 0.5_RFREAL*(ql* rwl       + pl*nz + qr* rwr       + pr*nz)*nm       
        flx(5) = 0.5_RFREAL*(ql*(rel + pl) + pl*fs + qr*(rer + pr) + pr*fs)*nm
        
        mfMixt(indMf*ifl) = flx(1)

        pPatch%cp(ifl)          = iCpRef*(0.5_RFREAL*(pl + pr) - pRef)
        pPatch%cmass(ifl)       = iCmassRef*flx(1)/nm
        pPatch%cmom(XCOORD,ifl) = iCpRef*0.5_RFREAL*(ul+ur)*flx(1)/nm
        pPatch%cmom(YCOORD,ifl) = iCpRef*0.5_RFREAL*(vl+vr)*flx(1)/nm
        pPatch%cmom(ZCOORD,ifl) = iCpRef*0.5_RFREAL*(wl+wr)*flx(1)/nm

        rhs(CV_MIXT_DENS,c1) = rhs(CV_MIXT_DENS,c1) + flx(1)
        rhs(CV_MIXT_XMOM,c1) = rhs(CV_MIXT_XMOM,c1) + flx(2)
        rhs(CV_MIXT_YMOM,c1) = rhs(CV_MIXT_YMOM,c1) + flx(3)
        rhs(CV_MIXT_ZMOM,c1) = rhs(CV_MIXT_ZMOM,c1) + flx(4)
        rhs(CV_MIXT_ENER,c1) = rhs(CV_MIXT_ENER,c1) + flx(5)

        sd(SD_XMOM,c1*indSd) = sd(SD_XMOM,c1*indSd) + 0.5_RFREAL*(rul/rl + rur/rr)*flx(1)
        sd(SD_YMOM,c1*indSd) = sd(SD_YMOM,c1*indSd) + 0.5_RFREAL*(rvl/rl + rvr/rr)*flx(1)
        sd(SD_ZMOM,c1*indSd) = sd(SD_ZMOM,c1*indSd) + 0.5_RFREAL*(rwl/rl + rwr/rr)*flx(1)

        IF ( pRegion%irkStep == 1 ) THEN
          global%massIn  = global%massIn  - MIN(flx(1),0.0_RFREAL)
          global%massOut = global%massOut + MAX(flx(1),0.0_RFREAL)
        END IF ! pRegion


        IF ( (global%checkLevel == CHECK_HIGH) .AND. &
             (global%verbLevel == VERBOSE_HIGH) .AND. &
             (global%myProcid == MASTERPROC) .AND. &
             (decidePrintFlag .EQV. .TRUE.) ) THEN
          IF ( flx(1) > 0.0_RFREAL ) THEN
            nLocs = nLocs + 1

            IF ( nLocs == 1 ) THEN
              global%warnCounter = global%warnCounter + 1

              WRITE(STDOUT,'(A,1X,A,1X,I9)') SOLVER_NAME, &
                    '*** WARNING *** Outflow detected at inflow boundary!'
              WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                               pRegion%iRegionGlobal
              IF ( global%flowType == FLOW_UNSTEADY ) THEN
                WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME, &
                                                    'Current time:', &
                                                    global%currentTime
              ELSE
                WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME, &
                                                 'Current iteration number:', &
                                                 global%currentIter
              END IF ! global%flowType

              WRITE(STDOUT,'(A,3X,A,1X,I1)') SOLVER_NAME, &
                                             'Runge-Kutta stage:', &
                                             pRegion%irkStep
            END IF ! nLocs

            IF ( nLocs <= MAX_INOUTFLOW_LOCS ) THEN
              loc(nLocs,MIN_VAL:MAX_VAL) = c1
            END IF ! nLocs
          END IF ! flx(1)
        END IF ! global%checkLevel
      END DO ! ifl

! ------------------------------------------------------------------------------
!     Write info on outflow at inflow boundary
! ------------------------------------------------------------------------------

      IF ( (global%checkLevel == CHECK_HIGH) .AND. &
           (global%verbLevel == VERBOSE_HIGH) .AND. &
           (global%myProcid == MASTERPROC) .AND. &
           (decidePrintFlag .EQV. .TRUE.) .AND. &
           (nLocs > 0) ) THEN
        IF ( nLocs > MAX_INOUTFLOW_LOCS ) THEN
           WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A,1X,I9,1X,A)') SOLVER_NAME, &
                 'Only wrote the first',MAX_INOUTFLOW_LOCS,'of',nLocs, &
                 'inflow faces with outflow.'
          CALL RFLU_PrintLocInfo(pRegion,loc,MAX_INOUTFLOW_LOCS, &
                                 LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
        ELSE
          CALL RFLU_PrintLocInfo(pRegion,loc(1:nLocs,MIN_VAL:MAX_VAL),nLocs, &
                                 LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
        END IF ! nLocs
      END IF ! global%verbLevel

! ==============================================================================
!   Inflow (based on velocities and temperature)
! ==============================================================================

    CASE ( BC_INFLOW_VELTEMP )
      vals       => pPatch%mixt%vals
      bcOptType  =  pPatch%mixt%switches(BCSWI_INFLOW_TYPE)

      DO ifl = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifl)

        nx = pPatch%fn(XCOORD,ifl)
        ny = pPatch%fn(YCOORD,ifl)
        nz = pPatch%fn(ZCOORD,ifl)
        nm = pPatch%fn(XYZMAG,ifl)

        fs = pPatch%gs(indGs*ifl)

        rl  = cv(CV_MIXT_DENS,c1)
        rul = cv(CV_MIXT_XMOM,c1)
        rvl = cv(CV_MIXT_YMOM,c1)
        rwl = cv(CV_MIXT_ZMOM,c1)
        rel = cv(CV_MIXT_ENER,c1)
        pl  = dv(DV_MIXT_PRES,c1)

        ur = vals(BCDAT_INFLOW_U,distrib*ifl)
        vr = vals(BCDAT_INFLOW_V,distrib*ifl)
        wr = vals(BCDAT_INFLOW_W,distrib*ifl)
        tr = vals(BCDAT_INFLOW_T,distrib*ifl)

        IF ( gasModel == GAS_MODEL_TCPERF ) THEN
          cp   = gv(GV_MIXT_CP ,indCp *c1)
          mm   = gv(GV_MIXT_MOL,indMol*c1)
          rgas = MixtPerf_R_M(mm)
          g    = MixtPerf_G_CpR(cp,rgas)

          IF ( bcOptType /= BCOPT_SUBSONIC ) THEN
            pr = vals(BCDAT_INFLOW_P,distrib*ifl)
            rr = MixtPerf_D_PRT(pr,rgas,tr)
          ELSE 
            rr = rl
            pr = MixtPerf_P_DRT(rr,rgas,tr)
          END IF ! bcOptType
          
          rer = rr*MixtPerf_Eo_DGPUVW(rr,g,pr,ur,vr,wr)
        ELSE IF ( gasModel == GAS_MODEL_MIXT_TCPERF ) THEN
#ifdef SPEC
          IF ( pRegion%spec%cvState /= CV_MIXT_STATE_CONS ) THEN
            CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
          END IF ! pRegion%spec%cvState  

          irl = 1.0_RFREAL/rl
          imm = 0.0_RFREAL
          cp  = 0.0_RFREAL

          DO iSpec = 1,pRegion%specInput%nSpecies
            pSpecType => pRegion%specInput%specType(iSpec)

            imm = imm + irl*pCvSpec(iSpec,c1)/pSpecType%pMaterial%molw
            cp  = cp  + irl*pCvSpec(iSpec,c1)*pSpecType%pMaterial%spht
          END DO ! iSpec

          mm   = 1.0_RFREAL/imm
          rgas = MixtPerf_R_M(mm)
          g    = MixtPerf_G_CpR(cp,rgas)

          IF ( bcOptType /= BCOPT_SUBSONIC ) THEN
            pr = vals(BCDAT_INFLOW_P,distrib*ifl)
            rr = MixtPerf_D_PRT(pr,rgas,tr)
          ELSE
            rr = rl
            pr = MixtPerf_P_DRT(rr,rgas,tr)
          END IF ! bcOptType

          rer = rr*MixtPerf_Eo_DGPUVW(rr,g,pr,ur,vr,wr)
! DEBUG: Manoj-PBA1D
IF (1==2) THEN
          rr  = rl
          ur  = rul/rl
          vr  = rvl/rl
          wr  = rwl/rl
          pr  = pl
!          tr  = tl
          rer = rel
!          WRITE(*,*) "Stopping here from inflow bc"
!          STOP
END IF ! 1==2
! END DEBUG
#else     
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) ! Defensive coding
#endif
        ELSE
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! gasModel

        rur = rr*ur
        rvr = rr*vr
        rwr = rr*wr

        ql = (rul*nx + rvl*ny + rwl*nz)/rl - fs
        qr = (rur*nx + rvr*ny + rwr*nz)/rr - fs

        ul = rul/rl
        vl = rvl/rl
        wl = rwl/rl

        vFracGl = 1.0_RFREAL

        flx(1) = 0.5_RFREAL*(ql* rl                + qr* rr               )*nm
        flx(2) = 0.5_RFREAL*(ql* rul       + pl*nx + qr* rur       + pr*nx)*nm
        flx(3) = 0.5_RFREAL*(ql* rvl       + pl*ny + qr* rvr       + pr*ny)*nm
        flx(4) = 0.5_RFREAL*(ql* rwl       + pl*nz + qr* rwr       + pr*nz)*nm       
        flx(5) = 0.5_RFREAL*(ql*(rel + pl) + pl*fs + qr*(rer + pr) + pr*fs)*nm

        mfMixt(indMf*ifl) = flx(1)

        pPatch%cp(ifl)          = iCpRef*(0.5_RFREAL*(pl + pr) - pRef)
        pPatch%cmass(ifl)       = iCmassRef*flx(1)/nm
        pPatch%cmom(XCOORD,ifl) = iCpRef*0.5_RFREAL*(ul+ur)*flx(1)/nm
        pPatch%cmom(YCOORD,ifl) = iCpRef*0.5_RFREAL*(vl+vr)*flx(1)/nm
        pPatch%cmom(ZCOORD,ifl) = iCpRef*0.5_RFREAL*(wl+wr)*flx(1)/nm

        rhs(CV_MIXT_DENS,c1) = rhs(CV_MIXT_DENS,c1) + flx(1)
        rhs(CV_MIXT_XMOM,c1) = rhs(CV_MIXT_XMOM,c1) + flx(2)
        rhs(CV_MIXT_YMOM,c1) = rhs(CV_MIXT_YMOM,c1) + flx(3)
        rhs(CV_MIXT_ZMOM,c1) = rhs(CV_MIXT_ZMOM,c1) + flx(4)
        rhs(CV_MIXT_ENER,c1) = rhs(CV_MIXT_ENER,c1) + flx(5)

        sd(SD_XMOM,c1*indSd) = sd(SD_XMOM,c1*indSd) + 0.5_RFREAL*(rul/rl + rur/rr)*flx(1)
        sd(SD_YMOM,c1*indSd) = sd(SD_YMOM,c1*indSd) + 0.5_RFREAL*(rvl/rl + rvr/rr)*flx(1)
        sd(SD_ZMOM,c1*indSd) = sd(SD_ZMOM,c1*indSd) + 0.5_RFREAL*(rwl/rl + rwr/rr)*flx(1)

        IF ( pRegion%irkStep == 1 ) THEN
          global%massIn  = global%massIn  - MIN(flx(1),0.0_RFREAL)
          global%massOut = global%massOut + MAX(flx(1),0.0_RFREAL)
        END IF ! pRegion


        IF ( (global%checkLevel == CHECK_HIGH) .AND. &
             (global%verbLevel == VERBOSE_HIGH) .AND. &
             (global%myProcid == MASTERPROC) .AND. &
             (decidePrintFlag .EQV. .TRUE.) ) THEN
          IF ( flx(1) > 0.0_RFREAL ) THEN
            nLocs = nLocs + 1

            IF ( nLocs == 1 ) THEN
              global%warnCounter = global%warnCounter + 1

              WRITE(STDOUT,'(A,1X,A,1X,I9)') SOLVER_NAME, &
                    '*** WARNING *** Outflow detected at inflow boundary!'
              WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                               pRegion%iRegionGlobal
              IF ( global%flowType == FLOW_UNSTEADY ) THEN
                WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME, &
                                                    'Current time:', &
                                                    global%currentTime
              ELSE
                WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME, &
                                                 'Current iteration number:', &
                                                 global%currentIter
              END IF ! global%flowType

              WRITE(STDOUT,'(A,3X,A,1X,I1)') SOLVER_NAME, &
                                             'Runge-Kutta stage:', &
                                             pRegion%irkStep
            END IF ! nLocs

            IF ( nLocs <= MAX_INOUTFLOW_LOCS ) THEN
              loc(nLocs,MIN_VAL:MAX_VAL) = c1
            END IF ! nLocs
          END IF ! flx(1)
        END IF ! global%checkLevel
      END DO ! ifl

! ------------------------------------------------------------------------------
!     Write info on outflow at inflow boundary
! ------------------------------------------------------------------------------

      IF ( (global%checkLevel == CHECK_HIGH) .AND. &
           (global%verbLevel == VERBOSE_HIGH) .AND. &
           (global%myProcid == MASTERPROC) .AND. &
           (decidePrintFlag .EQV. .TRUE.) .AND. &
           (nLocs > 0) ) THEN
        IF ( nLocs > MAX_INOUTFLOW_LOCS ) THEN
           WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A,1X,I9,1X,A)') SOLVER_NAME, &
                 'Only wrote the first',MAX_INOUTFLOW_LOCS,'of',nLocs, &
                 'inflow faces with outflow.'
          CALL RFLU_PrintLocInfo(pRegion,loc,MAX_INOUTFLOW_LOCS, &
                                 LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
        ELSE
          CALL RFLU_PrintLocInfo(pRegion,loc(1:nLocs,MIN_VAL:MAX_VAL),nLocs, &
                                 LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
        END IF ! nLocs
      END IF ! global%checkLevel

! ==============================================================================
!   Outflow
! ==============================================================================

    CASE ( BC_OUTFLOW )
      bcOptType = pPatch%mixt%switches(BCSWI_OUTFLOW_TYPE)

      IF ( bcOptType /= BCOPT_SUPERSONIC ) THEN
        vals => pPatch%mixt%vals
      END IF ! bcOptType

      DO ifl = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifl)

        nx = pPatch%fn(XCOORD,ifl)
        ny = pPatch%fn(YCOORD,ifl)
        nz = pPatch%fn(ZCOORD,ifl)
        nm = pPatch%fn(XYZMAG,ifl)

        fs = pPatch%gs(indGs*ifl)

        rl  = cv(CV_MIXT_DENS,c1)
        rul = cv(CV_MIXT_XMOM,c1)
        rvl = cv(CV_MIXT_YMOM,c1)
        rwl = cv(CV_MIXT_ZMOM,c1)
        rel = cv(CV_MIXT_ENER,c1)
        pl  = dv(DV_MIXT_PRES,c1)

        IF ( bcOptType /= BCOPT_SUPERSONIC ) THEN
          pr = vals(BCDAT_OUTFLOW_PRESS,distrib*ifl)
        ELSE
          pr = pl
        END IF ! bcOptType

        IF ( gasModel == GAS_MODEL_TCPERF ) THEN
          cp = gv(GV_MIXT_CP ,indCp *c1)
          mm = gv(GV_MIXT_MOL,indMol*c1)

          CALL BcondOutflowPerf(bcOptType,pr,nx,ny,nz,cp,mm,rl,rul,rvl,rwl, &
                                rel,pl,rr,rur,rvr,rwr,rer)
        ! Include JWL Gas model for outflow BC - Subbu
        !ELSE IF ( gasModel == GAS_MODEL_MIXT_TCPERF ) THEN
        
         !Fred - Added preliminary treatment for JWL EOS, need new scheme that
         !is
        !valid for JWL EOS as new treatment is the perfect gas one with the JWL
        !eqns, 9/22/15
        ELSE IF ( gasModel == GAS_MODEL_MIXT_TCPERF .OR. &
                  gasModel == GAS_MODEL_MIXT_JWL) THEN
#ifdef SPEC

         IF (global%specUsed .EQV. .TRUE.) THEN
         IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
             scalarConvFlag = .FALSE.
         ELSE
             scalarConvFlag = .TRUE.
    CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
         END IF
  
          cp = gv(GV_MIXT_CP ,indCp *c1)
          mm = gv(GV_MIXT_MOL,indMol*c1)

          iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')
          Yproducts =  pCvSpec(iCvSpecProducts,c1)

        ! IF (Yproducts .LE. 1.0E-06_RFREAL) THEN
          CALL BcondOutflowPerf(bcOptType,pr,nx,ny,nz,cp,mm,rl,rul,rvl,rwl, &
                                rel,pl,rr,rur,rvr,rwr,rer)

         !ELSE
         ! CALL BcondOutflowJWL(bcOptType,pr,nx,ny,nz,cp,mm,rl,rul,rvl,rwl, &
          !                      rel,pl,Yproducts,pRegion,rr,rur,rvr,rwr,rer)
         !ENDIF

!Fred - Add this in eventually...don't think a big enough mass fraction will
!reach the outflow boundary to make this worthwhile 
 
          IF (scalarConvFlag .EQV. .TRUE.) THEN
  CALL RFLU_ScalarConvertcvPrim2Cons(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
          END IF
      END IF

#else     
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) ! Defensive coding
#endif
        ELSE 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! gasModel

        ql = (rul*nx + rvl*ny + rwl*nz)/rl - fs
        qr = (rur*nx + rvr*ny + rwr*nz)/rr - fs

        ul = rul/rl
        vl = rvl/rl
        wl = rwl/rl
        ur = rur/rr
        vr = rvr/rr
        wr = rwr/rr

! DEBUG: Manoj-PBA1D
IF (1==2) THEN
          rr  = rl
!          ur  = ul
!          vr  = vl
!          wr  = wl
          pr  = pl
!          tr  = tl
          rur = rul
          rvr = rvl
          rwr = rwl
          rer = rel
          qr  = ql
!          WRITE(*,*) "Stopping here from outflow bc"
!          STOP
END IF ! 1==2
! END DEBUG

        flx(1) = 0.5_RFREAL*(ql* rl                + qr* rr               )*nm
        flx(2) = 0.5_RFREAL*(ql* rul       + pl*nx + qr* rur       + pr*nx)*nm
        flx(3) = 0.5_RFREAL*(ql* rvl       + pl*ny + qr* rvr       + pr*ny)*nm
        flx(4) = 0.5_RFREAL*(ql* rwl       + pl*nz + qr* rwr       + pr*nz)*nm       
        flx(5) = 0.5_RFREAL*(ql*(rel + pl) + pl*fs + qr*(rer + pr) + pr*fs)*nm
        
        mfMixt(indMf*ifl) = flx(1)

        pPatch%cp(ifl)          = iCpRef*(0.5_RFREAL*(pl + pr) - pRef)
        pPatch%cmass(ifl)       = iCmassRef*flx(1)/nm
        pPatch%cmom(XCOORD,ifl) = iCpRef*0.5_RFREAL*(ul+ur)*flx(1)/nm
        pPatch%cmom(YCOORD,ifl) = iCpRef*0.5_RFREAL*(vl+vr)*flx(1)/nm
        pPatch%cmom(ZCOORD,ifl) = iCpRef*0.5_RFREAL*(wl+wr)*flx(1)/nm

        rhs(CV_MIXT_DENS,c1) = rhs(CV_MIXT_DENS,c1) + flx(1)
        rhs(CV_MIXT_XMOM,c1) = rhs(CV_MIXT_XMOM,c1) + flx(2)
        rhs(CV_MIXT_YMOM,c1) = rhs(CV_MIXT_YMOM,c1) + flx(3)
        rhs(CV_MIXT_ZMOM,c1) = rhs(CV_MIXT_ZMOM,c1) + flx(4)
        rhs(CV_MIXT_ENER,c1) = rhs(CV_MIXT_ENER,c1) + flx(5)

        sd(SD_XMOM,c1*indSd) = sd(SD_XMOM,c1*indSd) + 0.5_RFREAL*(rul/rl + rur/rr)*flx(1)
        sd(SD_YMOM,c1*indSd) = sd(SD_YMOM,c1*indSd) + 0.5_RFREAL*(rvl/rl + rvr/rr)*flx(1)
        sd(SD_ZMOM,c1*indSd) = sd(SD_ZMOM,c1*indSd) + 0.5_RFREAL*(rwl/rl + rwr/rr)*flx(1)

        IF ( pRegion%irkStep == 1 ) THEN
          global%massIn  = global%massIn  - MIN(flx(1),0.0_RFREAL)
          global%massOut = global%massOut + MAX(flx(1),0.0_RFREAL)
        END IF ! pRegion

        IF ( (global%checkLevel == CHECK_HIGH) .AND. &
             (global%verbLevel == VERBOSE_HIGH) .AND. &
             (global%myProcid == MASTERPROC) .AND. &
             (decidePrintFlag .EQV. .TRUE.) ) THEN
          IF ( flx(1) < 0.0_RFREAL ) THEN
            nLocs = nLocs + 1

            IF ( nLocs == 1 ) THEN
              global%warnCounter = global%warnCounter + 1

              WRITE(STDOUT,'(A,1X,A,1X,I9)') SOLVER_NAME, &
                    '*** WARNING *** Inflow detected at outflow boundary!'
              WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                               pRegion%iRegionGlobal
              IF ( global%flowType == FLOW_UNSTEADY ) THEN
                WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME, &
                                                    'Current time:', &
                                                    global%currentTime
              ELSE
                WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME, &
                                                 'Current iteration number:', &
                                                 global%currentIter
              END IF ! global%flowType

              WRITE(STDOUT,'(A,3X,A,1X,I1)') SOLVER_NAME, &
                                             'Runge-Kutta stage:', &
                                             pRegion%irkStep
            END IF ! nLocs

            IF ( nLocs <= MAX_INOUTFLOW_LOCS ) THEN
              loc(nLocs,MIN_VAL:MAX_VAL) = c1
            END IF ! nLocs
          END IF ! flx(1)
        END IF ! global%checkLevel
      END DO ! ifl

! ------------------------------------------------------------------------------
!     Write info on inflow at outflow boundary
! ------------------------------------------------------------------------------

      IF ( (global%checkLevel == CHECK_HIGH) .AND. &
           (global%verbLevel == VERBOSE_HIGH) .AND. &
           (global%myProcid == MASTERPROC) .AND. &
           (decidePrintFlag .EQV. .TRUE.) .AND. &
           (nLocs > 0) ) THEN
        IF ( nLocs > MAX_INOUTFLOW_LOCS ) THEN
           WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A,1X,I9,1X,A)') SOLVER_NAME, &
                 'Only wrote the first',MAX_INOUTFLOW_LOCS,'of',nLocs, &
                 'outflow faces with inflow.'
          CALL RFLU_PrintLocInfo(pRegion,loc,MAX_INOUTFLOW_LOCS, &
                                 LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
        ELSE
          CALL RFLU_PrintLocInfo(pRegion,loc(1:nLocs,MIN_VAL:MAX_VAL),nLocs, &
                                 LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
        END IF ! nLocs
      END IF ! global%verbLevel

! ==============================================================================
!   Slip wall (weak imposition)
! ==============================================================================

    CASE ( BC_SLIPWALL )
      DO ifl = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifl)

        nx = pPatch%fn(XCOORD,ifl)
        ny = pPatch%fn(YCOORD,ifl)
        nz = pPatch%fn(ZCOORD,ifl)
        nm = pPatch%fn(XYZMAG,ifl)

        fs  = pPatch%gs(indGs*ifl)
        fsu = RFLU_DescaleGridSpeed(pRegion,fs)

        rl  = cv(CV_MIXT_DENS,c1)
        rul = cv(CV_MIXT_XMOM,c1)
        rvl = cv(CV_MIXT_YMOM,c1)
        rwl = cv(CV_MIXT_ZMOM,c1)
        pl  = dv(DV_MIXT_PRES,c1)

        SELECT CASE ( pRegion%mixtInput%gasModel )
          CASE ( GAS_MODEL_TCPERF, &
                 GAS_MODEL_MIXT_TCPERF, & 
                 GAS_MODEL_MIXT_PSEUDO )
            cp = gv(GV_MIXT_CP ,indCp *c1)
            mm = gv(GV_MIXT_MOL,indMol*c1)

        ! TEMPORARY:Subbu - RIND State causes the time-step to be reduced a
        ! great deal when secondary shock reaches internal wall in a cyl shock
        ! tube problem
            
            !CALL RFLU_SetRindStateSlipWallPerf(cp,mm,nx,ny,nz,rl,rul,rvl,rwl, &
            !                                   fsu,pl)

        ! END TEMPORARY
! TEMPORARY : using new implementation of slipwall
!             using roe flux solver rather than solving coupled eqn iteratively
!            h = 2.0_RFREAL* &
!                ABS((pRegion%grid%cofg(XCOORD,c1) - pPatch%fc(XCOORD,ifl))*nx &
!                  + (pRegion%grid%cofg(YCOORD,c1) - pPatch%fc(YCOORD,ifl))*ny &
!                  + (pRegion%grid%cofg(ZCOORD,c1) - pPatch%fc(ZCOORD,ifl))*nz)
!
!            rl  = cv(CV_MIXT_DENS,c1)
!            irl = 1.0_RFREAL/rl
!            ul  = cv(CV_MIXT_XMOM,c1)*irl
!            vl  = cv(CV_MIXT_YMOM,c1)*irl
!            wl  = cv(CV_MIXT_ZMOM,c1)*irl
!            pl  = dv(DV_MIXT_PRES,c1)
!            tl  = dv(DV_MIXT_TEMP,c1)
!
!            CALL RFLU_SetRindStateSlipWallPerf1(cp,mm,nx,ny,nz,h,rl,ul,vl,wl, &
!                                                fsu,pl,rg,pg)
!
!            cp = global%refCp
!            g  = global%refGamma
!
!            rr = rg
!            ratio = rl/rr
!            ratio = 1.0_RFREAL
!
!            ur = ul*ratio
!            vr = vl*ratio
!            wr = wl*ratio
!            CALL ReflectVector(nx,ny,nz,ur,vr,wr)
!            pr = pg
!            tr = tl*(rl*pr)/(rr*pl)
!
!            CALL RFLU_ROE_FluxFunction(nx,ny,nz,nm,fs,cp,g,rl,ul,vl,wl,pl,tl, &
!                                       rr,ur,vr,wr,pr,tr,flx)
! END TEMPORARY

          CASE ( GAS_MODEL_MIXT_JWL )
#ifdef SPEC
           IF (global%specUsed .EQV. .TRUE.) THEN
         IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
             scalarConvFlag = .FALSE.
         ELSE
             scalarConvFlag = .TRUE.
    CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
         END IF

            cp = gv(GV_MIXT_CP ,indCp *c1)
            mm = gv(GV_MIXT_MOL,indMol*c1)

           iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')
           Yproducts =  pCvSpec(iCvSpecProducts,c1)

         IF (Yproducts .LE. 1.0E-06_RFREAL) THEN
            CALL RFLU_SetRindStateSlipWallPerf(cp,mm,nx,ny,nz,rl,rul,rvl,rwl, &
                                               fsu,pl)
         ELSE
             CALL RFLU_SetRindStateSlipWallJWL(pRegion,c1,cp,mm,nx,ny,nz,rl,rul,rvl,rwl, &
                                               fsu,Yproducts,pl)
         ENDIF

          IF (scalarConvFlag .EQV. .TRUE.) THEN
  CALL RFLU_ScalarConvertcvPrim2Cons(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
          END IF
      END IF
#endif
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pRegion%mixtInput%gasModel

        flx(2) = pl*nx*nm 
        flx(3) = pl*ny*nm
        flx(4) = pl*nz*nm
        flx(5) = pl*fs*nm

        mfMixt(indMf*ifl) = 0.0_RFREAL

        pPatch%cp(ifl)          = iCpRef*(pl - pRef)

        pPatch%cmass(ifl)       = 0.0_RFREAL 
        pPatch%cmom(XCOORD,ifl) = 0.0_RFREAL
        pPatch%cmom(YCOORD,ifl) = 0.0_RFREAL
        pPatch%cmom(ZCOORD,ifl) = 0.0_RFREAL

        rhs(CV_MIXT_XMOM,c1) = rhs(CV_MIXT_XMOM,c1) + flx(2)
        rhs(CV_MIXT_YMOM,c1) = rhs(CV_MIXT_YMOM,c1) + flx(3)
        rhs(CV_MIXT_ZMOM,c1) = rhs(CV_MIXT_ZMOM,c1) + flx(4)
        rhs(CV_MIXT_ENER,c1) = rhs(CV_MIXT_ENER,c1) + flx(5)

      END DO ! ifl

! ==============================================================================
!   No-slip wall
! ==============================================================================

    CASE ( BC_NOSLIPWALL_HFLUX,BC_NOSLIPWALL_TEMP )
      DO ifl = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifl)

        nx = pPatch%fn(XCOORD,ifl)
        ny = pPatch%fn(YCOORD,ifl)
        nz = pPatch%fn(ZCOORD,ifl)
        nm = pPatch%fn(XYZMAG,ifl)

        fs = pPatch%gs(indGs*ifl)

        pl = dv(DV_MIXT_PRES,c1)

        flx(2) = pl*nx*nm 
        flx(3) = pl*ny*nm
        flx(4) = pl*nz*nm
        flx(5) = pl*fs*nm

        mfMixt(indMf*ifl) = 0.0_RFREAL

        pPatch%cp(ifl)          = iCpRef*(pl - pRef)
        pPatch%cmass(ifl)       = 0.0_RFREAL 
        pPatch%cmom(XCOORD,ifl) = 0.0_RFREAL
        pPatch%cmom(YCOORD,ifl) = 0.0_RFREAL
        pPatch%cmom(ZCOORD,ifl) = 0.0_RFREAL

        rhs(CV_MIXT_XMOM,c1) = rhs(CV_MIXT_XMOM,c1) + flx(2)
        rhs(CV_MIXT_YMOM,c1) = rhs(CV_MIXT_YMOM,c1) + flx(3)
        rhs(CV_MIXT_ZMOM,c1) = rhs(CV_MIXT_ZMOM,c1) + flx(4)
        rhs(CV_MIXT_ENER,c1) = rhs(CV_MIXT_ENER,c1) + flx(5)
      END DO ! ifl

! ==============================================================================
!   Farfield
! ==============================================================================

    CASE ( BC_FARFIELD )
      corr = pPatch%mixt%switches(BCSWI_FARF_CORR)
      
      IF ( corr == BCOPT_CORR_YES ) THEN 
        corrFlag = .TRUE. 
      ELSE 
        corrFlag = .FALSE.        
        liftCoef = 0.0_RFREAL
      END IF ! corr
      
      vals => pPatch%mixt%vals

      DO ifl = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifl)

        nx = pPatch%fn(XCOORD,ifl)
        ny = pPatch%fn(YCOORD,ifl)
        nz = pPatch%fn(ZCOORD,ifl)
        nm = pPatch%fn(XYZMAG,ifl)

        xc = pPatch%fc(XCOORD,ifl)
        yc = pPatch%fc(YCOORD,ifl)
        zc = pPatch%fc(ZCOORD,ifl)

        fs = pPatch%gs(indGs*ifl)

        rl  = cv(CV_MIXT_DENS,c1)
        rul = cv(CV_MIXT_XMOM,c1)
        rvl = cv(CV_MIXT_YMOM,c1)
        rwl = cv(CV_MIXT_ZMOM,c1)
        rel = cv(CV_MIXT_ENER,c1)
        pl  = dv(DV_MIXT_PRES,c1)

        mf  = vals(BCDAT_FARF_MACH  ,distrib*ifl)
        aoa = vals(BCDAT_FARF_ATTACK,distrib*ifl)
        aos = vals(BCDAT_FARF_SLIP  ,distrib*ifl)
        pf  = vals(BCDAT_FARF_PRESS ,distrib*ifl)
        tf  = vals(BCDAT_FARF_TEMP  ,distrib*ifl)

        IF ( gasModel == GAS_MODEL_TCPERF ) THEN
          cp = gv(GV_MIXT_CP ,indCp *c1)
          mm = gv(GV_MIXT_MOL,indMol*c1)
                                           
          CALL RFLU_SetRindStateFarfieldPerf(global,cp,mm,nx,ny,nz,mf,pf,tf, &
                                             aoa,aos,corrFlag,liftCoef,xc,yc, &
                                             zc,rl,rul,rvl,rwl,rel,rr,rur, & 
                                             rvr,rwr,rer,pr)      
        ELSE
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! gasModel

        ql = (rul*nx + rvl*ny + rwl*nz)/rl - fs
        qr = (rur*nx + rvr*ny + rwr*nz)/rr - fs

        ul = rul/rl
        vl = rvl/rl
        wl = rwl/rl
        ur = rur/rr
        vr = rvr/rr
        wr = rwr/rr

        flx(1) = 0.5_RFREAL*(ql* rl                + qr* rr               )*nm
        flx(2) = 0.5_RFREAL*(ql* rul       + pl*nx + qr* rur       + pr*nx)*nm
        flx(3) = 0.5_RFREAL*(ql* rvl       + pl*ny + qr* rvr       + pr*ny)*nm
        flx(4) = 0.5_RFREAL*(ql* rwl       + pl*nz + qr* rwr       + pr*nz)*nm       
        flx(5) = 0.5_RFREAL*(ql*(rel + pl) + pl*fs + qr*(rer + pr) + pr*fs)*nm
        
        mfMixt(indMf*ifl) = flx(1)

        pPatch%cp(ifl)          = iCpRef*(0.5_RFREAL*(pl + pr) - pRef)
        pPatch%cmass(ifl)       = iCmassRef*flx(1)/nm
        pPatch%cmom(XCOORD,ifl) = iCpRef*0.5_RFREAL*(ul+ur)*flx(1)/nm
        pPatch%cmom(YCOORD,ifl) = iCpRef*0.5_RFREAL*(vl+vr)*flx(1)/nm
        pPatch%cmom(ZCOORD,ifl) = iCpRef*0.5_RFREAL*(wl+wr)*flx(1)/nm

        rhs(CV_MIXT_DENS,c1) = rhs(CV_MIXT_DENS,c1) + flx(1)
        rhs(CV_MIXT_XMOM,c1) = rhs(CV_MIXT_XMOM,c1) + flx(2)
        rhs(CV_MIXT_YMOM,c1) = rhs(CV_MIXT_YMOM,c1) + flx(3)
        rhs(CV_MIXT_ZMOM,c1) = rhs(CV_MIXT_ZMOM,c1) + flx(4)
        rhs(CV_MIXT_ENER,c1) = rhs(CV_MIXT_ENER,c1) + flx(5)

        sd(SD_XMOM,c1*indSd) = sd(SD_XMOM,c1*indSd) + 0.5_RFREAL*(rul/rl + rur/rr)*flx(1)
        sd(SD_YMOM,c1*indSd) = sd(SD_YMOM,c1*indSd) + 0.5_RFREAL*(rvl/rl + rvr/rr)*flx(1)
        sd(SD_ZMOM,c1*indSd) = sd(SD_ZMOM,c1*indSd) + 0.5_RFREAL*(rwl/rl + rwr/rr)*flx(1)

        IF ( pRegion%irkStep == 1 ) THEN
          global%massIn  = global%massIn  - MIN(flx(1),0.0_RFREAL)
          global%massOut = global%massOut + MAX(flx(1),0.0_RFREAL)
        END IF ! pRegion
      END DO ! ifl

! ==============================================================================
!   Injection
! ==============================================================================

    CASE ( BC_INJECTION )
      vals => pPatch%mixt%vals

      DO ifl = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifl)

        nx = pPatch%fn(XCOORD,ifl)
        ny = pPatch%fn(YCOORD,ifl)
        nz = pPatch%fn(ZCOORD,ifl)
        nm = pPatch%fn(XYZMAG,ifl)

        fs  = pPatch%gs(indGs*ifl)
        fsu = RFLU_DescaleGridSpeed(pRegion,fs)

        IF ( gasModel == GAS_MODEL_TCPERF ) THEN
          cp = gv(GV_MIXT_CP ,indCp *c1)
          mm = gv(GV_MIXT_MOL,indMol*c1)

          minj = vals(BCDAT_INJECT_MFRATE,distrib*ifl)
                
          IF ( minj > 0.0_RFREAL ) THEN ! Surface burning
            tinj = vals(BCDAT_INJECT_TEMP,distrib*ifl)

            pl = dv(DV_MIXT_PRES,c1)

            CALL RFLU_SetRindStateInjectPerf(cp,mm,nx,ny,nz,minj,tinj,pl, &
                                             fsu,rl,ul,vl,wl,Hl)
                                             
            flx(1) =  -minj            *nm
            flx(2) = (-minj*ul + pl*nx)*nm
            flx(3) = (-minj*vl + pl*ny)*nm
            flx(4) = (-minj*wl + pl*nz)*nm
            flx(5) = (-minj*Hl + pl*fs)*nm                                   
          ELSE ! Surface NOT burning
            rl  = cv(CV_MIXT_DENS,c1)     
            rul = cv(CV_MIXT_XMOM,c1)
            rvl = cv(CV_MIXT_YMOM,c1)
            rwl = cv(CV_MIXT_ZMOM,c1) 
            pl  = dv(DV_MIXT_PRES,c1)
                
            CALL RFLU_SetRindStateSlipWallPerf(cp,mm,nx,ny,nz,rl,rul,rvl,rwl, &
                                               fsu,pl)
                
            ul = 0.0_RFREAL
            vl = 0.0_RFREAL
            wl = 0.0_RFREAL    
                                             
            flx(1) = 0.0_RFREAL
            flx(2) = pl*nx*nm
            flx(3) = pl*ny*nm
            flx(4) = pl*nz*nm
            flx(5) = pl*fs*nm                                                
          END IF ! minj
        ELSE 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! gasModel

        mfMixt(indMf*ifl) = flx(1)

        pPatch%cp(ifl)          = iCpRef*(pl - pRef)
        pPatch%cmass(ifl)       = iCmassRef*flx(1)/nm
        pPatch%cmom(XCOORD,ifl) = iCpRef*0.5_RFREAL*(ul+ur)*flx(1)/nm
        pPatch%cmom(YCOORD,ifl) = iCpRef*0.5_RFREAL*(vl+vr)*flx(1)/nm
        pPatch%cmom(ZCOORD,ifl) = iCpRef*0.5_RFREAL*(wl+wr)*flx(1)/nm

        rhs(CV_MIXT_DENS,c1) = rhs(CV_MIXT_DENS,c1) + flx(1)
        rhs(CV_MIXT_XMOM,c1) = rhs(CV_MIXT_XMOM,c1) + flx(2)
        rhs(CV_MIXT_YMOM,c1) = rhs(CV_MIXT_YMOM,c1) + flx(3)
        rhs(CV_MIXT_ZMOM,c1) = rhs(CV_MIXT_ZMOM,c1) + flx(4)
        rhs(CV_MIXT_ENER,c1) = rhs(CV_MIXT_ENER,c1) + flx(5)

        sd(SD_XMOM,c1*indSd) = sd(SD_XMOM,c1*indSd) + ul*flx(1)
        sd(SD_YMOM,c1*indSd) = sd(SD_YMOM,c1*indSd) + vl*flx(1)
        sd(SD_ZMOM,c1*indSd) = sd(SD_ZMOM,c1*indSd) + wl*flx(1)

        IF ( pRegion%irkStep == 1 ) THEN
          global%massIn  = global%massIn  - MIN(flx(1),0.0_RFREAL)
          global%massOut = global%massOut + MAX(flx(1),0.0_RFREAL)
        END IF ! pRegion
      END DO ! ifl

! ==============================================================================
!   Boundaries for which fluxes must not or need not be computed
! ==============================================================================

    CASE ( BC_PERIODIC, &
           BC_SYMMETRY, & 
           BC_VIRTUAL )

! ==============================================================================
!   Default
! ==============================================================================

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pPatch%bcType

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CentralFirstPatch

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CentralFirstPatch.F90,v $
! Revision 1.5  2016/04/06 19:15:57  rahul
! Removed vFrac definitions.
!
! Revision 1.4  2016/02/08 20:04:19  fred
! Fixing non-SPEC flag compiling error
!
! Revision 1.3  2016/02/04 20:00:08  fred
! Adding iterative JWL EOS capabilities for the cylindrical detonation case
!
! Revision 1.2  2015/12/19 01:13:16  rahul
! Updated to suit Ausm+up implementation.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.6  2009/07/09 20:44:10  mparmar
! Added support for GAS_MODEL_MIXT_TCPERF in BC_INFLOW_VELTEMP
!
! Revision 1.5  2009/07/08 19:12:07  mparmar
! Removed TEMPORARY statements
!
! Revision 1.4  2008/12/06 08:43:47  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:16:59  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/06/18 18:08:04  mparmar
! Implemented curvature corrected slipwall BC
!
! Revision 1.1  2007/04/09 18:49:56  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:01  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.31  2006/10/20 21:31:38  mparmar
! Added computation of mass and moment coeffs
!
! Revision 1.30  2006/08/19 15:39:50  mparmar
! Renamed patch variables
!
! Revision 1.29  2006/04/07 15:19:22  haselbac
! Removed tabs
!
! Revision 1.28  2006/03/25 22:00:43  haselbac
! Added CASEs for sype patches
!
! Revision 1.27  2005/11/14 17:00:29  haselbac
! Added support for pseudo-gas model
!
! Revision 1.26  2005/11/10 02:31:05  haselbac
! Added support for mixture of tcperf gases
!
! Revision 1.25  2005/10/31 21:09:37  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.24  2005/10/05 14:17:38  haselbac
! Adapted to changes in no-slip wall bc defs
!
! Revision 1.23  2005/05/16 20:44:09  haselbac
! Converted to pRegion and pPatch, cosmetics
!
! Revision 1.22  2005/04/27 02:12:24  haselbac
! Added treatment for velocity-temperature inlet
!
! Revision 1.21  2005/04/20 14:43:12  haselbac
! Removed CHECK_UNIFLOW code section
!
! Revision 1.20  2005/03/31 17:05:06  haselbac
! Changed computation of sd
!
! Revision 1.19  2005/03/09 15:08:11  haselbac
! Added virtual boundary
!
! Revision 1.18  2004/12/27 23:31:30  haselbac
! Changed setting of bc farf rind state (not completed)
!
! Revision 1.17  2004/10/19 19:29:07  haselbac
! Modified fluxes for injecting boundaries, clean-up
!
! Revision 1.16  2004/07/30 22:47:36  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.15  2004/06/16 20:01:05  haselbac
! Added setting of cp, cosmetics
!
! Revision 1.14  2004/04/14 02:09:07  haselbac
! Added proper setting of rind state for slip walls
!
! Revision 1.13  2004/02/13 03:00:37  haselbac
! Fixed comment
!
! Revision 1.12  2004/01/31 03:58:36  haselbac
! Improved printing of reverse flow messages
!
! Revision 1.11  2004/01/29 22:59:16  haselbac
! Removed hardcode, improved mfMixt use, adapted bcondInflowPerf call
!
! Revision 1.10  2004/01/22 16:04:33  haselbac
! Changed declaration to eliminate warning on ALC
!
! Revision 1.9  2003/12/04 03:29:58  haselbac
! Changed logic for printing reverse flow information, clean-up
!
! Revision 1.8  2003/11/25 21:04:37  haselbac
! Added computation of mass and volume fluxes on patches
!
! Revision 1.7  2003/09/04 20:19:37  haselbac
! Removed temporary fix for Rocburn problem in GENx
!
! Revision 1.6  2003/08/23 21:44:21  haselbac
! Removed restriction to checkLevel for inflow-outflow checks
!
! Revision 1.5  2003/08/22 20:31:56  haselbac
! Added temporary fix for ignition problems
!
! Revision 1.4  2003/07/23 22:14:43  haselbac
! Fixed cut-and-paste bug from last check-in
!
! Revision 1.3  2003/07/23 15:47:15  haselbac
! Enclosed WARNINGS within ifs, otherwise too verbose
!
! Revision 1.2  2003/07/22 02:06:32  haselbac
! Cosmetics and inflow-outflow warnings
!
! Revision 1.1  2003/06/04 22:49:27  haselbac
! Renamed from RFLU_RoeCentralFirstPatch
!
! Revision 1.12  2003/01/28 14:47:47  haselbac
! Use parameters in fn
!
! Revision 1.11  2002/12/19 23:00:50  haselbac
! Rewrote injection fluxes to use minj directly
!
! Revision 1.10  2002/11/27 20:26:51  haselbac
! Cosmetic change only (use of vals pointer)
!
! Revision 1.9  2002/11/26 15:18:11  haselbac
! Cosmetic changes only
!
! Revision 1.8  2002/11/19 23:45:30  haselbac
! Corrected flux on slip wall with grid motion
!
! Revision 1.7  2002/11/08 21:33:45  haselbac
! Removed grid speed from slip-wall bc
!
! Revision 1.6  2002/10/05 19:34:03  haselbac
! GENX integration: changed injection wall, fixed bug for uniform flow check
!
! Revision 1.5  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.4  2002/09/09 15:51:56  haselbac
! global and mixtInput now under region, added grid speeds, fixed pointer bug for vals
!
! Revision 1.3  2002/06/22 16:58:07  haselbac
! Added USE ModInterface, modified argument lists for inflow and injection
!
! Revision 1.2  2002/06/10 21:34:03  haselbac
! Added boundary fluxes, calls to bcond*Perf, changed flag to CHECK_UNIFLOW
!
! Revision 1.1  2002/05/04 17:01:59  haselbac
! Initial revision
!
! ******************************************************************************

