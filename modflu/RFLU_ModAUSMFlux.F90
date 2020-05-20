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
! Purpose: Collection of AUSM flux routines.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModAUSMFlux.F90,v 1.2 2016/02/04 22:20:29 fred Exp $
!
! Copyright: (c) 2005-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModAUSMFlux

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModGrid, ONLY: t_grid
  USE ModMPI

#ifdef SPEC
  USE ModSpecies, ONLY: t_spec_input, &
                        t_spec_type
  USE ModInterfacesSpecies, ONLY: SPEC_GetSpeciesIndex
#endif

  USE RFLU_ModGridSpeedUtils, ONLY: RFLU_DescaleGridSpeed
  USE RFLU_ModRindStates

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_AUSM_ComputeFlux
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModAUSMFlux.F90,v $ $Revision: 1.2 $' 
              
! *****************************************************************************
! Routines
! *****************************************************************************

  CONTAINS






! ******************************************************************************
!
! Purpose: Wrapper function for AUSM flux functions.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_AUSM_ComputeFlux(pRegion)
                       
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

  CALL RegisterFunction(global,'RFLU_AUSM_ComputeFlux',__FILE__)

! ******************************************************************************
! Call flux functions
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%gasModel ) 

! ==============================================================================
!   Thermally and calorically perfect gas
! ==============================================================================

    CASE ( GAS_MODEL_TCPERF )
      SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
        CASE ( DISCR_ORDER_1 ) 
          CALL RFLU_AUSM_ComputeFlux1(pRegion)
        CASE ( DISCR_ORDER_2 )
          CALL RFLU_AUSM_ComputeFlux2_TCP(pRegion)               
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pRegion%mixtInput%spaceOrder    

! ==============================================================================
!   Mixture of thermally and calorically perfect gases
! ==============================================================================

    CASE ( GAS_MODEL_MIXT_TCPERF ) 
      SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
        CASE ( DISCR_ORDER_1 ) 
          CALL RFLU_AUSM_ComputeFlux1(pRegion)
        CASE ( DISCR_ORDER_2 )
          CALL RFLU_AUSM_ComputeFlux2_MTCP(pRegion)               
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pRegion%mixtInput%spaceOrder

! ==============================================================================
!   Pseudo-gas
! ==============================================================================

    CASE ( GAS_MODEL_MIXT_PSEUDO ) 
      SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
! TO DO
!        CASE ( DISCR_ORDER_1 ) 
!          CALL RFLU_AUSM_ComputeFlux1_MPSD(pRegion)
! END TO DO 
        CASE ( DISCR_ORDER_2 )
          CALL RFLU_AUSM_ComputeFlux2_MPSD(pRegion)               
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pRegion%mixtInput%spaceOrder                           

! ==============================================================================
!   Mixture of thermally and calorically perfect gas and detonation products
! ==============================================================================

    CASE ( GAS_MODEL_MIXT_JWL ) 
      SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
        CASE ( DISCR_ORDER_1 ) 
          CALL RFLU_AUSM_ComputeFlux1_MJWL(pRegion)
        CASE ( DISCR_ORDER_2 )
          CALL RFLU_AUSM_ComputeFlux2_MJWL(pRegion)               
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pRegion%mixtInput%spaceOrder

! ==============================================================================
!   Default
! ==============================================================================

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)             
  END SELECT ! pRegion%mixtInput%gasModel       

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AUSM_ComputeFlux







! ******************************************************************************
!
! Purpose: Compute convective fluxes using AUSM+ scheme.
!
! Description: None.
!
! Input: 
!   nx          x-component of face normal
!   ny          y-component of face normal
!   nz          z-component of face normal
!   nm          Magnitude of face normal
!   fs          Face speed
!   rl          Density of left state
!   ul          x-component of velocity of left state
!   vl          y-component of velocity of left state
!   wl          z-component of velocity of left state   
!   Hl		Total enthalpy of left state
!   al		Speed of sound of left state
!   pl          Pressure of left state
!   rr          Density of right state
!   ur          x-component of velocity of right state
!   vr          y-component of velocity of right state
!   wr          z-component of velocity of right state  
!   pr          Pressure of right state
!   Hr		Total enthalpy of right state
!   ar		Speed of sound of right state
!
! Output: 
!   flx         Fluxes
!   vf          Face velocities
!
! Notes: 
!   1. Liou M.-S., Progress towards an improved CFD method: AUSM+, AIAA Paper
!      95-1701, 1995
!   2. Do not use computation of face speed of sound which leads to exact 
!      capturing of isolated normal shock waves because of robustness problems
!      for unsteady flows and because that formulation is not applicable to 
!      anything but calorically and thermally perfect gases.
!
! ******************************************************************************

SUBROUTINE RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur, &
                                  vr,wr,pr,Hr,ar,flx,vf)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  REAL(RFREAL), INTENT(IN) :: al,ar,fs,Hl,Hr,nm,nx,ny,nz,pl,pr,rl,rr,ul,ur,vl, &
                              vr,wl,wr
  REAL(RFREAL), INTENT(OUT) :: flx(5),vf(3)

! ==============================================================================
! Locals
! ==============================================================================

  REAL(RFREAL) :: af,mf,mfa,mfm,mfp,ml,mla,mlp,mr,mra,mrm,pf,ql,qr,vml,vmr, &
                  wtl,wtr

! ******************************************************************************
! Start, compute face state
! ******************************************************************************

  ql = ul*nx + vl*ny + wl*nz - fs
  qr = ur*nx + vr*ny + wr*nz - fs
    
  af = 0.5_RFREAL*(al+ar) ! NOTE not using original formulation, see note

  ml  = ql/af
  mla = ABS(ml)

  mr  = qr/af
  mra = ABS(mr)    

  IF ( mla <= 1.0_RFREAL ) THEN 
    mlp = 0.25_RFREAL*(ml+1.0_RFREAL)*(ml+1.0_RFREAL) & 
        + 0.125_RFREAL*(ml*ml-1.0_RFREAL)*(ml*ml-1.0_RFREAL)
    wtl = 0.25_RFREAL*(ml+1.0_RFREAL)*(ml+1.0_RFREAL)*(2.0_RFREAL-ml) & 
        + 0.1875_RFREAL*ml*(ml*ml-1.0_RFREAL)*(ml*ml-1.0_RFREAL)
  ELSE
    mlp = 0.5_RFREAL*(ml+mla)
    wtl = 0.5_RFREAL*(1.0_RFREAL+ml/mla)
  END IF ! mla

  IF ( mra <= 1.0_RFREAL ) THEN 
    mrm = -0.25_RFREAL*(mr-1.0_RFREAL)*(mr-1.0_RFREAL) & 
          -0.125_RFREAL*(mr*mr-1.0_RFREAL)*(mr*mr-1.0_RFREAL)
    wtr = 0.25_RFREAL*(mr-1.0_RFREAL)*(mr-1.0_RFREAL)*(2.0_RFREAL+mr) & 
        - 0.1875_RFREAL*mr*(mr*mr-1.0_RFREAL)*(mr*mr-1.0_RFREAL)
  ELSE
    mrm = 0.5_RFREAL*(mr-mra)
    wtr = 0.5_RFREAL*(1.0_RFREAL-mr/mra)
  END IF ! mla

  mf  = mlp + mrm
  mfa = ABS(mf)
  mfp = 0.5_RFREAL*(mf+mfa)
  mfm = 0.5_RFREAL*(mf-mfa)

  pf = wtl*pl + wtr*pr 

! ******************************************************************************
! Compute fluxes
! ******************************************************************************

  vf(1) = mfp*ul + mfm*ur
  vf(2) = mfp*vl + mfm*vr
  vf(3) = mfp*wl + mfm*wr    

  flx(1) = (af*(mfp*rl    + mfm*rr   )        )*nm
  flx(2) = (af*(mfp*rl*ul + mfm*rr*ur) + pf*nx)*nm
  flx(3) = (af*(mfp*rl*vl + mfm*rr*vr) + pf*ny)*nm
  flx(4) = (af*(mfp*rl*wl + mfm*rr*wr) + pf*nz)*nm
  flx(5) = (af*(mfp*rl*Hl + mfm*rr*Hr) + pf*fs)*nm
  !IF (ml > 0.0_RFREAL .OR. mr > 0.0_RFREAL)  WRITE(*,*) 'ml,mr=',ml,mr 
  
! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE RFLU_AUSM_FluxFunction









! ******************************************************************************
!
! Purpose: Compute convective fluxes using AUSM+ scheme.
!
! Description: None.
!
! Input: 
!   nx          x-component of face normal
!   ny          y-component of face normal
!   nz          z-component of face normal
!   nm          Magnitude of face normal
!   fs          Face speed
!   rl          Density of left state
!   ul          x-component of velocity of left state
!   vl          y-component of velocity of left state
!   wl          z-component of velocity of left state   
!   Hl		Total enthalpy of left state
!   al		Speed of sound of left state
!   pl          Pressure of left state
!   rr          Density of right state
!   ur          x-component of velocity of right state
!   vr          y-component of velocity of right state
!   wr          z-component of velocity of right state  
!   pr          Pressure of right state
!   Hr		Total enthalpy of right state
!   ar		Speed of sound of right state
!
! Output: 
!   flx         Fluxes
!   vf          Face velocities
!
! Notes: 
!   1. Liou M.-S., Progress towards an improved CFD method: AUSM+, AIAA Paper
!      95-1701, 1995
!   2. Do not use computation of face speed of sound which leads to exact 
!      capturing of isolated normal shock waves because of robustness problems
!      for unsteady flows and because that formulation is not applicable to 
!      anything but calorically and thermally perfect gases.
!
! ******************************************************************************

SUBROUTINE RFLU_AUSM_FluxFunction1(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur, &
                                   vr,wr,pr,Hr,ar,flx,vf)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  REAL(RFREAL), INTENT(IN) :: al,ar,fs,Hl,Hr,nm,nx,ny,nz,pl,pr,rl,rr,ul,ur,vl, &
                              vr,wl,wr
  REAL(RFREAL), INTENT(OUT) :: flx(5),vf(3)

! ==============================================================================
! Locals
! ==============================================================================

  REAL(RFREAL) :: af,mf,mfa,mfm,mfp,ml,mla,mlp,mr,mra,mrm,pf,ql,qr,vml,vmr, &
                  wtl,wtr

! ******************************************************************************
! Start, compute face state
! ******************************************************************************

  ql = ul*nx + vl*ny + wl*nz - fs
  qr = ur*nx + vr*ny + wr*nz - fs
    
  af = 0.5_RFREAL*(al+ar) ! NOTE not using original formulation, see note
WRITE(*,'(A,E24.16)') "af=",af

  ml  = ql/af
  mla = ABS(ml)

  mr  = qr/af
  mra = ABS(mr)    
WRITE(*,'(2(A,E24.16,1X))') "ml=",ml," mr=",mr

  IF ( mla <= 1.0_RFREAL ) THEN 
    mlp = 0.25_RFREAL*(ml+1.0_RFREAL)*(ml+1.0_RFREAL) & 
        + 0.125_RFREAL*(ml*ml-1.0_RFREAL)*(ml*ml-1.0_RFREAL)
    wtl = 0.25_RFREAL*(ml+1.0_RFREAL)*(ml+1.0_RFREAL)*(2.0_RFREAL-ml) & 
        + 0.1875_RFREAL*ml*(ml*ml-1.0_RFREAL)*(ml*ml-1.0_RFREAL)
  ELSE
    mlp = 0.5_RFREAL*(ml+mla)
    wtl = 0.5_RFREAL*(1.0_RFREAL+ml/mla)
  END IF ! mla
WRITE(*,'(2(A,E24.16,1X))') "mlp=",mlp," wtl=",wtl

  IF ( mra <= 1.0_RFREAL ) THEN 
    mrm = -0.25_RFREAL*(mr-1.0_RFREAL)*(mr-1.0_RFREAL) & 
          -0.125_RFREAL*(mr*mr-1.0_RFREAL)*(mr*mr-1.0_RFREAL)
    wtr = 0.25_RFREAL*(mr-1.0_RFREAL)*(mr-1.0_RFREAL)*(2.0_RFREAL+mr) & 
        - 0.1875_RFREAL*mr*(mr*mr-1.0_RFREAL)*(mr*mr-1.0_RFREAL)
  ELSE
    mrm = 0.5_RFREAL*(mr-mra)
    wtr = 0.5_RFREAL*(1.0_RFREAL-mr/mra)
  END IF ! mla
WRITE(*,'(2(A,E24.16,1X))') "mrm=",mrm," wtr=",wtr

  mf  = mlp + mrm
  mfa = ABS(mf)
  mfp = 0.5_RFREAL*(mf+mfa)
  mfm = 0.5_RFREAL*(mf-mfa)
WRITE(*,'(3(A,E24.16,1X))') "mf=",mf," mfp=",mfp," mfm=",mfm

  pf = wtl*pl + wtr*pr 
WRITE(*,'(A,E24.16)') "pf=",pf

! ******************************************************************************
! Compute fluxes
! ******************************************************************************

  vf(1) = mfp*ul + mfm*ur
  vf(2) = mfp*vl + mfm*vr
  vf(3) = mfp*wl + mfm*wr    

  flx(1) = (af*(mfp*rl    + mfm*rr   )        )*nm
  flx(2) = (af*(mfp*rl*ul + mfm*rr*ur) + pf*nx)*nm
  flx(3) = (af*(mfp*rl*vl + mfm*rr*vr) + pf*ny)*nm
  flx(4) = (af*(mfp*rl*wl + mfm*rr*wr) + pf*nz)*nm
  flx(5) = (af*(mfp*rl*Hl + mfm*rr*Hr) + pf*fs)*nm

! ******************************************************************************
! Compute fluxes
! ******************************************************************************

!  vf(1) = mfp*ul + mfm*ur
!  vf(2) = mfp*vl + mfm*vr
!  vf(3) = mfp*wl + mfm*wr    

!  flx(1) = ((ql*rl   )        )*nm
!  flx(2) = ((ql*rl*ul) + pl*nx)*nm
!  flx(3) = ((ql*rl*vl) + pl*ny)*nm
!  flx(4) = ((ql*rl*wl) + pl*nz)*nm
!  flx(5) = ((ql*rl*Hl) + pl*fs)*nm

! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE RFLU_AUSM_FluxFunction1









! ******************************************************************************
!
! Purpose: Compute convective fluxes using first-order accurate AUSM+ scheme.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine can also be used for a mixture of thermally and calorically 
!      perfect gases because, the face states being identical to the cell 
!      states, the mixture properties are constant.
!
! ******************************************************************************

SUBROUTINE RFLU_AUSM_ComputeFlux1(pRegion)

  USE ModInterfaces, ONLY: MixtPerf_Ho_CpTUVW, & 
                           RFLU_CentralFirstPatch

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

  LOGICAL :: interfaceFlag,updateFlux
  INTEGER :: c1,c2,ifg,indCp,indGs,indMf,indMol,indSd,iPatch
  REAL(RFREAL) :: al,ar,cpl,cpr,fs,fsu,Hl,Hr,irl,irr,mwl,mwr,nm,nTol,nx,ny,nz, &
                  pl,pr,rl,rr,rul,rur,rvl,rvr,rwl,rwr,tl,tr,ul,ur,vl,vr,wl,wr
  REAL(RFREAL) :: flx(5),vf(3)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pRhs,pSd
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
#ifdef SPEC
  INTEGER :: iCvSpecExplosive,iCvSpecProducts,prodInterfaceType
  REAL(RFREAL) :: ud,YExpl,YExpr,YProdl,YProdr
  TYPE(t_spec_type), POINTER :: pSpecType
  TYPE(t_spec_input), POINTER :: pSpecInput
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AUSM_ComputeFlux1',__FILE__)

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::AUSM_ComputeFlux1")
#endif

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs  = pRegion%grid%indGs

  indCp  = pRegion%mixtInput%indCp
  indMf  = pRegion%mixtInput%indMfMixt  
  indMol = pRegion%mixtInput%indMol
  indSd  = pRegion%mixtInput%indSd

  pCv  => pRegion%mixt%cv
  pDv  => pRegion%mixt%dv
  pGv  => pRegion%mixt%gv
  
  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs
  pSd  => pRegion%mixt%sd  

#ifdef SPEC
  pSpecInput => pRegion%specInput
#endif

  nTol = 1.0E-14_RFREAL

! ******************************************************************************
! Compute fluxes through interior faces
! ******************************************************************************

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

! ==============================================================================
!   Get face geometry and grid speed
! ==============================================================================

    nx = pGrid%fn(XCOORD,ifg)
    ny = pGrid%fn(YCOORD,ifg)
    nz = pGrid%fn(ZCOORD,ifg)
    nm = pGrid%fn(XYZMAG,ifg)

    fs = pGrid%gs(indGs*ifg)

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! ------------------------------------------------------------------------------
!   Left state
! ------------------------------------------------------------------------------

    cpl = pGv(GV_MIXT_CP,indCp*c1)
    
    rl  = pCv(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl

    ul  = pCv(CV_MIXT_XMOM,c1)*irl
    vl  = pCv(CV_MIXT_YMOM,c1)*irl
    wl  = pCv(CV_MIXT_ZMOM,c1)*irl
    
    pl  = pDv(DV_MIXT_PRES,c1)
    tl  = pDv(DV_MIXT_TEMP,c1)
    al  = pDv(DV_MIXT_SOUN,c1)

! TEMPORARY: Manoj            
#ifdef PLAG
    IF ( global%plagUsed .AND. (pRegion%plag%nPcls > 1) ) THEN
      rl = rl/(1.0_RFREAL - pRegion%plag%vFracE(1,c1))
    END IF ! global%plagUsed
#endif
! END TEMPORARY: Manoj

! DEBUG: Manoj-PBA1D, Notes: To match Jianghui's implementation
!                         1. Use conventional definiation of enthalpy
!                    Changing it to earlier implementation of Manoj
! END DEBUG    
! DEBUG: Manoj-PBA1D    
!    Hl  = MixtPerf_Ho_CpTUVW(cpl,tl,ul,vl,wl)
    Hl  = pCv(CV_MIXT_ENER,c1)*irl + pl/rl
! END DEBUG

! TEMPORARY: Manoj            
!#ifdef PLAG
!    IF ( global%plagUsed .AND. (pRegion%plag%nPcls > 1) ) THEN
!      rl = rl*(1.0_RFREAL - pRegion%plag%vFracE(1,c1))
!    END IF ! global%plagUsed
!#endif
! END TEMPORARY: Manoj 

! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------

    cpr = pGv(GV_MIXT_CP,indCp*c2)

    rr  = pCv(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr

    ur  = pCv(CV_MIXT_XMOM,c2)*irr
    vr  = pCv(CV_MIXT_YMOM,c2)*irr
    wr  = pCv(CV_MIXT_ZMOM,c2)*irr
    
    pr  = pDv(DV_MIXT_PRES,c2)
    tr  = pDv(DV_MIXT_TEMP,c2)    
    ar  = pDv(DV_MIXT_SOUN,c2)
    
! TEMPORARY: Manoj            
#ifdef PLAG
    IF ( global%plagUsed .AND. (pRegion%plag%nPcls > 1) ) THEN
      rr = rr/(1.0_RFREAL - pRegion%plag%vFracE(1,c2))
    END IF ! global%plagUsed
#endif
! END TEMPORARY: Manoj 

! DEBUG: Manoj-PBA1D    
!    Hr  = MixtPerf_Ho_CpTUVW(cpr,tr,ur,vr,wr)
    Hr  = pCv(CV_MIXT_ENER,c2)*irr + pr/rr
! END DEBUG

! TEMPORARY: Manoj            
!#ifdef PLAG
!    IF ( global%plagUsed .AND. (pRegion%plag%nPcls > 1) ) THEN
!      rr = rr*(1.0_RFREAL - pRegion%plag%vFracE(1,c2))
!    END IF ! global%plagUsed
!#endif
! END TEMPORARY: Manoj 

! ==============================================================================
!   Compute fluxes 
! ==============================================================================

! DEBUG: Manoj-PBA1D, checking fluxes
IF (1==2) THEN
  IF ( (ifg > 499) .AND. (ifg < 501) ) THEN
    WRITE(*,'(2(A,I4,E16.8))') "c1=",c1,pGrid%cofg(1,c1)," ==  c2=",c2,pGrid%cofg(1,c2)       
    WRITE(*,'(A,5(1X,E23.16))') "Left  state:",rl,ul,pl,Hl,al
    WRITE(*,'(A,5(1X,E23.16))') "Right state:",rr,ur,pr,Hr,ar

    CALL RFLU_AUSM_FluxFunction1(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur,vr, &
                                 wr,pr,Hr,ar,flx,vf)

    WRITE(*,'(A,I4,3(1X,E23.16))') "ifg=",ifg,flx(1),flx(2),flx(5)
    WRITE(*,*) "--------------------------------------------------------------------------"
    WRITE(*,*) " "
  ELSE
    CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur,vr, &
                                wr,pr,Hr,ar,flx,vf)
  END IF ! ifg
ELSE
    CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur,vr, &
                                wr,pr,Hr,ar,flx,vf)
END IF ! 1==2    
! END DEBUG

! ==============================================================================
!   Program burn: If c1 or c2 contains YExplosive>0, then treat face as slipwall
! ==============================================================================

    updateFlux = .TRUE.

! DEBUG: Manoj-PBA, commenting out flux treatment at explosive-air, explosive-product faces 
IF (1==1) THEN
#ifdef SPEC
    IF ( (global%pbaFlag .EQV. .TRUE.) .AND. &
         (global%pbaBurnFlag .EQV. .TRUE.)) THEN      
      iCvSpecExplosive = SPEC_GetSpeciesIndex(global,pSpecInput,'EXPLOSIVE')
      iCvSpecProducts  = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')
      YExpl  = irl*pRegion%spec%cv(iCvSpecExplosive,c1)
      YExpr  = irr*pRegion%spec%cv(iCvSpecExplosive,c2)
      YProdl = irl*pRegion%spec%cv(iCvSpecProducts,c1)
      YProdr = irr*pRegion%spec%cv(iCvSpecProducts,c2)

! ------------------------------------------------------------------------------
!     Identify the interface between explosive and non-explosive  
! ------------------------------------------------------------------------------

      interfaceFlag = .FALSE.

      IF ( (ABS(YExpl) >= nTol) .AND. (ABS(YExpr) < nTol) ) THEN
        interfaceFlag = .TRUE.
      END IF ! YExpl, YExpr

      IF ( (ABS(YExpl) < nTol) .AND. (ABS(YExpr) >= nTol) ) THEN
        interfaceFlag = .TRUE.
      END IF ! YExpl, YExpr

      prodInterfaceType = 0

      IF ( ABS(1.0_RFREAL-YProdl) < nTol ) THEN
        prodInterfaceType = 1
      END IF ! YProdl

      IF ( ABS(1.0_RFREAL-YProdr) < nTol ) THEN
        prodInterfaceType = 2
      END IF ! YProdr

! ------------------------------------------------------------------------------
!     Compute fluxes if it is interface between explosive and non-explosive  
! ------------------------------------------------------------------------------

      IF ( interfaceFlag .EQV. .TRUE. ) THEN
        fsu = RFLU_DescaleGridSpeed(pRegion,fs)

! ----- Its an Air-Explosive interface -----------------------------------------        
        IF ( prodInterfaceType == 0 ) THEN
          Hl = Hl - 0.5_RFREAL*(ul*ul + vl*vl + wl*wl)
          Hr = Hr - 0.5_RFREAL*(ur*ur + vr*vr + wr*wr)

          ul = 0.0_RFREAL
          vl = 0.0_RFREAL
          wl = 0.0_RFREAL
          ur = 0.0_RFREAL
          vr = 0.0_RFREAL
          wr = 0.0_RFREAL

!          Hl  = MixtPerf_Ho_CpTUVW(cpl,tl,ul,vl,wl)
!          Hr  = MixtPerf_Ho_CpTUVW(cpr,tr,ur,vr,wr)
        END IF ! prodInterfaceType

! ----- Product is on left side ------------------------------------------------        
        IF ( prodInterfaceType == 1 ) THEN
          rr = rl
          ur = ul
          vr = vl
          wr = wl
          pr = pl
          ar = al
          Hr = Hl
        END IF ! prodInterfaceType

! ----- Product is on right side -----------------------------------------------        
        IF ( prodInterfaceType == 2 ) THEN
          rl = rr
          ul = ur
          vl = vr
          wl = wr
          pl = pr
          al = ar
          Hl = Hr
        END IF ! prodInterfaceType

! ------------------------------------------------------------------------------
!       Left state, compute flux with pr = pl
! ------------------------------------------------------------------------------

        CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al, &
                                    rr,ur,vr,wr,pl,Hr,ar,flx,vf)

! ----  Accumulate into residual of left cell ----------------------------------

        pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + flx(1)
        pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + flx(2)
        pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + flx(3)
        pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + flx(4)
        pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + flx(5)

! ----  Accumulate into substantial derivative of left cell --------------------
    
        pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + vf(1)*flx(1)
        pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vf(2)*flx(1)
        pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + vf(3)*flx(1)

! ------------------------------------------------------------------------------
!       Right state, compute flux with velocities at face=0, and pl = pr
! ------------------------------------------------------------------------------

        CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pr,Hl,al, &
                                    rr,ur,vr,wr,pr,Hr,ar,flx,vf)

! ----  Accumulate into residual of right cell ---------------------------------
    
        pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - flx(1)
        pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - flx(2)
        pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - flx(3)
        pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - flx(4)
        pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - flx(5)
    
! ----  Accumulate into substantial derivative of right cell -------------------
    
        pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - vf(1)*flx(1)
        pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vf(2)*flx(1)
        pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - vf(3)*flx(1)   
    
! ------------------------------------------------------------------------------
!       Store mass flux
! ------------------------------------------------------------------------------

! TEMPORARY: Manoj-PBA, making mfMixt=0 so that species cv is not modified
        pMf(indMf*ifg) = 0.0_RFREAL
! END TEMPORARY

! ------------------------------------------------------------------------------
!       Set flag for not updating flux later
! ------------------------------------------------------------------------------

        updateFlux = .FALSE.

      END IF ! interfaceFlag
    END IF ! global%pbaFlag
#endif
END IF ! 1==2
! END DEBUG

! ==============================================================================
!   Program burn: Update mass flux, residual, and substantial detivative only
!   when face is not in contact with explosive region with YExplosive=1
! ==============================================================================

    IF ( updateFlux ) THEN

! ==============================================================================
!     Store mass flux
! ==============================================================================

      pMf(indMf*ifg) = flx(1)

! ==============================================================================
!     Accumulate into residual
! ==============================================================================

      pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + flx(1)
      pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + flx(2)
      pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + flx(3)
      pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + flx(4)
      pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + flx(5)

      pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - flx(1)
      pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - flx(2)
      pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - flx(3)
      pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - flx(4)
      pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - flx(5)
    
! ==============================================================================
!     Accumulate into substantial derivative
! ==============================================================================

      pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + vf(1)*flx(1)
      pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vf(2)*flx(1)
      pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + vf(3)*flx(1)
    
      pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - vf(1)*flx(1)
      pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vf(2)*flx(1)
      pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - vf(3)*flx(1)   
    END IF ! updateFlux
  END DO ! ifg

! DEBUG: Manoj-PBA1D, checking fluxes
IF (1==2) THEN
  IF ( (global%pbaFlag .EQV. .TRUE.) .AND. &
       (global%pbaBurnFlag .EQV. .TRUE.)) THEN      
    c1 = 511
    WRITE(*,'(A,I3,3(2X,E24.16))') "icell=",c1, &
    pRhs(CV_MIXT_DENS,c1),pRhs(CV_MIXT_XMOM,c1),pRhs(CV_MIXT_ENER,c1)
  END IF ! glboal%pbaFlag  
END IF ! 1==2    
! END DEBUG

! DEBUG: Manoj-PBA1D
IF (1==2) THEN
WRITE(*,*) 'Stopping here...'
STOP
END IF ! 1==2
! END DEBUG
! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    CALL RFLU_CentralFirstPatch(pRegion,pPatch)
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF 
  CALL FPROFILER_ENDS("RFLU::AUSM_ComputeFlux1")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AUSM_ComputeFlux1







! ******************************************************************************
!
! Purpose: Compute convective fluxes using first-order accurate AUSM+ scheme 
!   for mixture of thermally and calorically perfect gas and detonation
!   products.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine can also be used for a mixture of thermally and calorically 
!      perfect gases because, the face states being identical to the cell 
!      states, the mixture properties are constant.
!
! ******************************************************************************

SUBROUTINE RFLU_AUSM_ComputeFlux1_MJWL(pRegion)

  USE ModInterfaces, ONLY: MixtPerf_Ho_CpTUVW, & 
                           RFLU_CentralFirstPatch

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

  INTEGER :: c1,c2,ifg,indCp,indGs,indMf,indMol,indSd,iPatch
  REAL(RFREAL) :: al,ar,cpl,cpr,fs,Hl,Hr,irl,irr,nm,nx,ny,nz,pl,pr,rl,rr, &
                  tl,tr,ul,ur,vl,vr,wl,wr
  REAL(RFREAL) :: flx(5),vf(3)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pRhs,pSd
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AUSM_ComputeFlux1_MJWL',__FILE__)

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::AUSM_ComputeFlux1_MJWL")
#endif

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs  = pRegion%grid%indGs

  indCp  = pRegion%mixtInput%indCp
  indMf  = pRegion%mixtInput%indMfMixt  
  indMol = pRegion%mixtInput%indMol
  indSd  = pRegion%mixtInput%indSd

  pCv  => pRegion%mixt%cv
  pDv  => pRegion%mixt%dv
  pGv  => pRegion%mixt%gv
  
  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs
  pSd  => pRegion%mixt%sd  

! ******************************************************************************
! Compute fluxes through interior faces
! ******************************************************************************

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

! ==============================================================================
!   Get face geometry and grid speed
! ==============================================================================

    nx = pGrid%fn(XCOORD,ifg)
    ny = pGrid%fn(YCOORD,ifg)
    nz = pGrid%fn(ZCOORD,ifg)
    nm = pGrid%fn(XYZMAG,ifg)

    fs = pGrid%gs(indGs*ifg)

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! ------------------------------------------------------------------------------
!   Left state
! ------------------------------------------------------------------------------

    cpl = pGv(GV_MIXT_CP,indCp*c1)
    
    rl  = pCv(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl

    ul  = pCv(CV_MIXT_XMOM,c1)*irl
    vl  = pCv(CV_MIXT_YMOM,c1)*irl
    wl  = pCv(CV_MIXT_ZMOM,c1)*irl
    
    pl  = pDv(DV_MIXT_PRES,c1)
    tl  = pDv(DV_MIXT_TEMP,c1)
    al  = pDv(DV_MIXT_SOUN,c1)

! TEMPORARY: Manoj-JWL, Confirm it    
!    Hl  = MixtPerf_Ho_CpTUVW(cpl,tl,ul,vl,wl)
    Hl  = pCv(CV_MIXT_ENER,c1)*irl + pl/rl
! END TEMPORARY

! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------

    cpr = pGv(GV_MIXT_CP,indCp*c2)

    rr  = pCv(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr

    ur  = pCv(CV_MIXT_XMOM,c2)*irr
    vr  = pCv(CV_MIXT_YMOM,c2)*irr
    wr  = pCv(CV_MIXT_ZMOM,c2)*irr
    
    pr  = pDv(DV_MIXT_PRES,c2)
    tr  = pDv(DV_MIXT_TEMP,c2)    
    ar  = pDv(DV_MIXT_SOUN,c2)
    
! TEMPORARY: Manoj-JWL, Confirm it    
!    Hr  = MixtPerf_Ho_CpTUVW(cpr,tr,ur,vr,wr)
    Hr  = pCv(CV_MIXT_ENER,c2)*irr + pr/rr
! END TEMPORARY

! ==============================================================================
!   Compute fluxes 
! ==============================================================================

    CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur,vr, &
                                wr,pr,Hr,ar,flx,vf)

! ==============================================================================
!   Store mass flux
! ==============================================================================

    pMf(indMf*ifg) = flx(1)

! ==============================================================================
!   Accumulate into residual
! ==============================================================================

    pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + flx(1)
    pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + flx(2)
    pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + flx(3)
    pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + flx(4)
    pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + flx(5)

    pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - flx(1)
    pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - flx(2)
    pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - flx(3)
    pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - flx(4)
    pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - flx(5)
    
! ==============================================================================
!   Accumulate into substantial derivative
! ==============================================================================

    pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + vf(1)*flx(1)
    pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vf(2)*flx(1)
    pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + vf(3)*flx(1)
    
    pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - vf(1)*flx(1)
    pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vf(2)*flx(1)
    pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - vf(3)*flx(1)   
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    CALL RFLU_CentralFirstPatch(pRegion,pPatch)
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF 
  CALL FPROFILER_ENDS("RFLU::AUSM_ComputeFlux1_MJWL")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AUSM_ComputeFlux1_MJWL








! ******************************************************************************
!
! Purpose: Compute convective fluxes using second-order accurate AUSM+ scheme 
!   for mixture of thermally and calorically perfect gas and detonation
!   products
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

SUBROUTINE RFLU_AUSM_ComputeFlux2_MJWL(pRegion)

  USE ModInterfaces, ONLY: MixtPerf_C_GRT, &
                           MixtPerf_G_CpR, &
                           MixtPerf_Ho_CpTUVW, & 
                           MixtPerf_R_M, & 
                           MixtPerf_T_DPR, & 
                           RFLU_CentralFirstPatch, & 
                           RFLU_CentralSecondPatch

  USE RFLU_ModJWL

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

  INTEGER :: c1,c2,ifg,indCp,indGs,indMf,indMol,indSd,iPatch
  REAL(RFREAL) :: al,ar,cp,cpl,cpr,dx,dy,dz,el,er,fs,g,gc,gcl,gcr,gl,gr,Hl,Hr, &
                  irl,irr,mw,nm,nx,ny,nz,pl,pr,rl,rr,tl,tr,ul,ur,vl,vr,wl,wr, &
                  xc,yc,zc
  REAL(RFREAL) :: flx(5),vf(3)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pRhs,pSd
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

#ifdef SPEC
  INTEGER :: iCvSpecAir,iCvSpecProducts,iSpec
  REAL(RFREAL) :: cpAir,gAir,gcAir,mml,mmr,mwAir,YProducts,Y1,Y2
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec 
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGcSpec 
  TYPE(t_spec_type), POINTER :: pSpecType
  TYPE(t_spec_input), POINTER :: pSpecInput
#endif  

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AUSM_ComputeFlux2_MJWL',__FILE__)

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::AUSM_ComputeFlux2_MJWL")
#endif

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs  = pGrid%indGs
  
  indCp  = pRegion%mixtInput%indCp
  indMf  = pRegion%mixtInput%indMfMixt  
  indMol = pRegion%mixtInput%indMol 
  indSd  = pRegion%mixtInput%indSd
  
  pCv  => pRegion%mixt%cv
  pDv  => pRegion%mixt%dv
  pGv  => pRegion%mixt%gv
  
  pGc  => pRegion%mixt%gradCell  
  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs
  pSd  => pRegion%mixt%sd  

#ifdef SPEC
  pSpecInput => pRegion%specInput
  pCvSpec => pRegion%spec%cv  
  pGcSpec => pRegion%spec%gradCell    
#endif

  IF ( pRegion%spec%cvState /= CV_MIXT_STATE_CONS ) THEN 
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! pRegion%spec%cvState  

  IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_JWL ) THEN
    CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__)  
  END IF ! pRegion%mixtInput%gasModel

  IF ( (indCp /= 1) .OR. (indMol /= 1) ) THEN 
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
  END IF ! indCp

! ******************************************************************************
! Compute fluxes through interior faces
! ******************************************************************************

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

! ==============================================================================
!   Get face geometry and grid speed
! ==============================================================================

    nx = pGrid%fn(XCOORD,ifg)
    ny = pGrid%fn(YCOORD,ifg)
    nz = pGrid%fn(ZCOORD,ifg)
    nm = pGrid%fn(XYZMAG,ifg)

    xc = pGrid%fc(XCOORD,ifg)
    yc = pGrid%fc(YCOORD,ifg)
    zc = pGrid%fc(ZCOORD,ifg)

    fs = pGrid%gs(indGs*ifg)

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! ------------------------------------------------------------------------------
!   Left state
! ------------------------------------------------------------------------------

    rl  = pCv(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl
    
    dx  = xc - pGrid%cofg(XCOORD,c1)
    dy  = yc - pGrid%cofg(YCOORD,c1)
    dz  = zc - pGrid%cofg(ZCOORD,c1)

    ul = pCv(CV_MIXT_XMOM,c1)*irl
    vl = pCv(CV_MIXT_YMOM,c1)*irl
    wl = pCv(CV_MIXT_ZMOM,c1)*irl
    pl = pDv(DV_MIXT_PRES,c1)

    rl = rl + pGc(XCOORD,GRC_MIXT_DENS,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_DENS,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_DENS,c1)*dz
    ul = ul + pGc(XCOORD,GRC_MIXT_XVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_XVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_XVEL,c1)*dz
    vl = vl + pGc(XCOORD,GRC_MIXT_YVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_YVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_YVEL,c1)*dz
    wl = wl + pGc(XCOORD,GRC_MIXT_ZVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_ZVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_ZVEL,c1)*dz
    pl = pl + pGc(XCOORD,GRC_MIXT_PRES,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_PRES,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_PRES,c1)*dz    

#ifdef SPEC
    iCvSpecAir = SPEC_GetSpeciesIndex(global,pSpecInput,'AIR')
    mwAir = pSpecInput%specType(iCvSpecAir)%pMaterial%molw
    cpAir = pSpecInput%specType(iCvSpecAir)%pMaterial%spht
    gcAir = MixtPerf_R_M(mwAir)
    gAir  = MixtPerf_G_CpR(cpAir,gcAir)

    iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')

    YProducts = irl*pCvSpec(iCvSpecProducts,c1) &
              + pGcSpec(XCOORD,iCvSpecProducts,c1)*dx &
              + pGcSpec(YCOORD,iCvSpecProducts,c1)*dy &
              + pGcSpec(ZCOORD,iCvSpecProducts,c1)*dz
 
    !CALL RFLU_JWL_ComputeEnergyMixt(gAir,gcAir,pl,rl,YProducts,al,el,tl)
    CALL RFLU_JWL_ComputeEnergyMixt(pRegion,c1,gAir,gcAir,pl,rl,YProducts,al,el,tl)

    Hl = el + 0.5_RFREAL*(ul*ul+vl*vl+wl*wl) + pl/rl
#endif                   

! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------

    rr  = pCv(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr
    
    dx  = xc - pGrid%cofg(XCOORD,c2)
    dy  = yc - pGrid%cofg(YCOORD,c2)
    dz  = zc - pGrid%cofg(ZCOORD,c2)

    ur = pCv(CV_MIXT_XMOM,c2)*irr
    vr = pCv(CV_MIXT_YMOM,c2)*irr
    wr = pCv(CV_MIXT_ZMOM,c2)*irr
    pr = pDv(DV_MIXT_PRES,c2)
        
    rr = rr + pGc(XCOORD,GRC_MIXT_DENS,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_DENS,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_DENS,c2)*dz
    ur = ur + pGc(XCOORD,GRC_MIXT_XVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_XVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_XVEL,c2)*dz
    vr = vr + pGc(XCOORD,GRC_MIXT_YVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_YVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_YVEL,c2)*dz
    wr = wr + pGc(XCOORD,GRC_MIXT_ZVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_ZVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_ZVEL,c2)*dz
    pr = pr + pGc(XCOORD,GRC_MIXT_PRES,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_PRES,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_PRES,c2)*dz

#ifdef SPEC
    iCvSpecAir = SPEC_GetSpeciesIndex(global,pSpecInput,'AIR')
    mwAir = pSpecInput%specType(iCvSpecAir)%pMaterial%molw
    cpAir = pSpecInput%specType(iCvSpecAir)%pMaterial%spht
    gcAir = MixtPerf_R_M(mwAir)
    gAir  = MixtPerf_G_CpR(cpAir,gcAir)

    iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')

    YProducts = irr*pCvSpec(iCvSpecProducts,c2) &
              + pGcSpec(XCOORD,iCvSpecProducts,c2)*dx &
              + pGcSpec(YCOORD,iCvSpecProducts,c2)*dy &
              + pGcSpec(ZCOORD,iCvSpecProducts,c2)*dz

    !CALL RFLU_JWL_ComputeEnergyMixt(gAir,gcAir,pr,rr,YProducts,ar,er,tr)
    CALL RFLU_JWL_ComputeEnergyMixt(pRegion,c2,gAir,gcAir,pr,rr,YProducts,ar,er,tr)

    Hr = er + 0.5_RFREAL*(ur*ur+vr*vr+wr*wr) + pr/rr
#endif                   

! ==============================================================================
!   Compute fluxes 
! ==============================================================================

    CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur,vr, &
                                wr,pr,Hr,ar,flx,vf)

! ==============================================================================
!   Store mass flux
! ==============================================================================

    pMf(indMf*ifg) = flx(1)

! ==============================================================================
!   Accumulate into residual
! ==============================================================================

    pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + flx(1)
    pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + flx(2)
    pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + flx(3)
    pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + flx(4)
    pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + flx(5)

    pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - flx(1)
    pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - flx(2)
    pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - flx(3)
    pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - flx(4)
    pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - flx(5)

! ==============================================================================
!   Accumulate into substantial derivative
! ==============================================================================

    pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + vf(1)*flx(1)
    pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vf(2)*flx(1)
    pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + vf(3)*flx(1)
    
    pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - vf(1)*flx(1)
    pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vf(2)*flx(1)
    pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - vf(3)*flx(1)   
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    SELECT CASE ( pPatch%spaceOrder )
      CASE ( 1 ) 
        CALL RFLU_CentralFirstPatch(pRegion,pPatch)      
      CASE ( 2 ) 
        CALL RFLU_CentralSecondPatch(pRegion,pPatch)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pPatch%spaceOrder
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF 
  CALL FPROFILER_ENDS("RFLU::AUSM_ComputeFlux2_MJWL")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AUSM_ComputeFlux2_MJWL








! ******************************************************************************
!
! Purpose: Compute convective fluxes using second-order accurate AUSM+ scheme 
!   for mixture of thermally and calorically perfect gaseous species and 
!   particulate species.
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

SUBROUTINE RFLU_AUSM_ComputeFlux2_MPSD(pRegion)

#ifdef SPEC
  USE ModSpecies, ONLY: t_spec_type
#endif

  USE ModInterfaces, ONLY: MixtPerf_C_GRT, &
                           MixtPerf_G_CpR, &
                           MixtPerf_Ho_CpTUVW, & 
                           MixtPerf_M_R, & 
                           MixtPerf_R_M, & 
                           MixtPerf_T_DPR, & 
                           RFLU_CentralFirstPatch, &
                           RFLU_CentralSecondPatch

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

  INTEGER :: c1,c2,ifg,indCp,indGs,indMf,indMol,indSd,iPatch
  REAL(RFREAL) :: al,ar,cpl,cpr,dx,dy,dz,fs,gcl,gcr,gl,gr,Hl,Hr,irl,irr,mml, &
                  mmr,nm,nx,ny,nz,pl,pr,rl,rr,tl,tr,ul,ur,vl,vr,wl,wr,xc,yc,zc
  REAL(RFREAL) :: flx(5),vf(3)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pRhs,pSd
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

#ifdef SPEC
  INTEGER :: iSpec
  REAL(RFREAL) :: gcg,immg,mmg,phip,Yg,Y1,Y2
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec 
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGcSpec 
  TYPE(t_spec_type), POINTER :: pSpecType
#endif  

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AUSM_ComputeFlux2_MPSD',__FILE__)

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::AUSM_ComputeFlux2_MPSD")
#endif

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs  = pGrid%indGs
  
  indCp  = pRegion%mixtInput%indCp
  indMf  = pRegion%mixtInput%indMfMixt  
  indMol = pRegion%mixtInput%indMol 
  indSd  = pRegion%mixtInput%indSd
  
  pCv  => pRegion%mixt%cv
  pDv  => pRegion%mixt%dv
  pGv  => pRegion%mixt%gv
  
  pGc  => pRegion%mixt%gradCell  
  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs
  pSd  => pRegion%mixt%sd  

#ifdef SPEC
  pCvSpec => pRegion%spec%cv  
  pGcSpec => pRegion%spec%gradCell    
#endif

  IF ( pRegion%spec%cvState /= CV_MIXT_STATE_CONS ) THEN 
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! pRegion%spec%cvState  

  IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_PSEUDO ) THEN
    CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__)  
  END IF ! pRegion%mixtInput%gasModel

  IF ( (indCp /= 1) .OR. (indMol /= 1) ) THEN 
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
  END IF ! indCp

! ******************************************************************************
! Compute fluxes through interior faces
! ******************************************************************************

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

! ==============================================================================
!   Get face geometry and grid speed
! ==============================================================================

    nx = pGrid%fn(XCOORD,ifg)
    ny = pGrid%fn(YCOORD,ifg)
    nz = pGrid%fn(ZCOORD,ifg)
    nm = pGrid%fn(XYZMAG,ifg)

    xc = pGrid%fc(XCOORD,ifg)
    yc = pGrid%fc(YCOORD,ifg)
    zc = pGrid%fc(ZCOORD,ifg)

    fs = pGrid%gs(indGs*ifg)

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! ------------------------------------------------------------------------------
!   Left state
! ------------------------------------------------------------------------------

    rl  = pCv(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl
    
    dx  = xc - pGrid%cofg(XCOORD,c1)
    dy  = yc - pGrid%cofg(YCOORD,c1)
    dz  = zc - pGrid%cofg(ZCOORD,c1)

#ifdef SPEC
    immg = 0.0_RFREAL
    cpl  = 0.0_RFREAL
    Yg   = 0.0_RFREAL
    phip = 0.0_RFREAL

    DO iSpec = 1,pRegion%specInput%nSpecies
      pSpecType => pRegion%specInput%specType(iSpec)

      Y1 = irl*pCvSpec(iSpec,c1) + pGcSpec(XCOORD,iSpec,c1)*dx & 
                                 + pGcSpec(YCOORD,iSpec,c1)*dy & 
                                 + pGcSpec(ZCOORD,iSpec,c1)*dz

      IF ( pSpecType%discreteFlag .EQV. .TRUE. ) THEN 
        cpl  = cpl  +    Y1*pSpecType%pMaterial%spht
        phip = phip + rl*Y1/pSpecType%pMaterial%dens
      ELSE 
        immg = immg + Y1/pSpecType%pMaterial%molw
        cpl  = cpl  + Y1*pSpecType%pMaterial%spht
        Yg   = Yg   + Y1                                                              
      END IF ! pSpecType%discreteFlag
    END DO ! iSpec

    mmg = Yg/immg
    gcg = MixtPerf_R_M(mmg)
    gcl = Yg/(1.0_RFREAL-phip)*gcg
    mml = MixtPerf_M_R(gcl)
#endif                    

    ul = pCv(CV_MIXT_XMOM,c1)*irl
    vl = pCv(CV_MIXT_YMOM,c1)*irl
    wl = pCv(CV_MIXT_ZMOM,c1)*irl
    pl = pDv(DV_MIXT_PRES,c1)
        
    rl = rl + pGc(XCOORD,GRC_MIXT_DENS,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_DENS,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_DENS,c1)*dz
    ul = ul + pGc(XCOORD,GRC_MIXT_XVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_XVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_XVEL,c1)*dz
    vl = vl + pGc(XCOORD,GRC_MIXT_YVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_YVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_YVEL,c1)*dz
    wl = wl + pGc(XCOORD,GRC_MIXT_ZVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_ZVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_ZVEL,c1)*dz
    pl = pl + pGc(XCOORD,GRC_MIXT_PRES,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_PRES,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_PRES,c1)*dz    

    tl = MixtPerf_T_DPR(rl,pl,gcl)
    al = MixtPerf_C_GRT(gl,gcl,tl)
    Hl = MixtPerf_Ho_CpTUVW(cpl,tl,ul,vl,wl)

! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------

    rr  = pCv(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr
    
    dx  = xc - pGrid%cofg(XCOORD,c2)
    dy  = yc - pGrid%cofg(YCOORD,c2)
    dz  = zc - pGrid%cofg(ZCOORD,c2)

#ifdef SPEC
    immg = 0.0_RFREAL
    cpr  = 0.0_RFREAL
    Yg   = 0.0_RFREAL
    phip = 0.0_RFREAL

    DO iSpec = 1,pRegion%specInput%nSpecies
      pSpecType => pRegion%specInput%specType(iSpec)

      Y2 = irr*pCvSpec(iSpec,c2) + pGcSpec(XCOORD,iSpec,c2)*dx & 
                                 + pGcSpec(YCOORD,iSpec,c2)*dy & 
                                 + pGcSpec(ZCOORD,iSpec,c2)*dz

      IF ( pSpecType%discreteFlag .EQV. .TRUE. ) THEN 
        cpr  = cpr  +    Y2*pSpecType%pMaterial%spht
        phip = phip + rr*Y2/pSpecType%pMaterial%dens
      ELSE 
        immg = immg + Y2/pSpecType%pMaterial%molw
        cpr  = cpr  + Y2*pSpecType%pMaterial%spht
        Yg   = Yg   + Y2                                                              
      END IF ! pSpecType%discreteFlag
    END DO ! iSpec

    mmg = Yg/immg
    gcg = MixtPerf_R_M(mmg)
    gcr = Yg/(1.0_RFREAL-phip)*gcg
    mmr = MixtPerf_M_R(gcr)
#endif  
    
    ur = pCv(CV_MIXT_XMOM,c2)*irr
    vr = pCv(CV_MIXT_YMOM,c2)*irr
    wr = pCv(CV_MIXT_ZMOM,c2)*irr
    pr = pDv(DV_MIXT_PRES,c2)
        
    rr = rr + pGc(XCOORD,GRC_MIXT_DENS,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_DENS,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_DENS,c2)*dz
    ur = ur + pGc(XCOORD,GRC_MIXT_XVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_XVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_XVEL,c2)*dz
    vr = vr + pGc(XCOORD,GRC_MIXT_YVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_YVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_YVEL,c2)*dz
    wr = wr + pGc(XCOORD,GRC_MIXT_ZVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_ZVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_ZVEL,c2)*dz
    pr = pr + pGc(XCOORD,GRC_MIXT_PRES,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_PRES,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_PRES,c2)*dz

    tr = MixtPerf_T_DPR(rr,pr,gcr)
    ar = MixtPerf_C_GRT(gr,gcr,tr)
    Hr = MixtPerf_Ho_CpTUVW(cpr,tr,ur,vr,wr)
    
! ==============================================================================
!   Compute fluxes 
! ==============================================================================

    gl = cpl/(cpl-gcl)
    gr = cpr/(cpr-gcr)
    
    IF ( ABS(gr-gl) < 0.05_RFREAL ) THEN 
      CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al, & 
                                  rr,ur,vr,wr,pr,Hr,ar,flx,vf)

      pMf(indMf*ifg) = flx(1)

      pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + flx(1)
      pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + flx(2)
      pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + flx(3)
      pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + flx(4)
      pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + flx(5)

      pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - flx(1)
      pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - flx(2)
      pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - flx(3)
      pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - flx(4)
      pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - flx(5)

      pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + vf(1)*flx(1)
      pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vf(2)*flx(1)
      pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + vf(3)*flx(1)

      pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - vf(1)*flx(1)
      pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vf(2)*flx(1)
      pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - vf(3)*flx(1)
    ELSE 
      CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al, & 
                                  rr,ur,vr,wr,pr,Hr,ar,flx,vf)

      pMf(indMf*ifg) = 0.5_RFREAL*flx(1)

      pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + flx(1)
      pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + flx(2)
      pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + flx(3)
      pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + flx(4)
      pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + flx(5)

      pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + vf(1)*flx(1)
      pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vf(2)*flx(1)
      pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + vf(3)*flx(1)
      
      CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al, & 
                                  rr,ur,vr,wr,pr,Hr,ar,flx,vf)

      pMf(indMf*ifg) = pMf(indMf*ifg) + 0.5_RFREAL*flx(1)

      pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - flx(1)
      pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - flx(2)
      pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - flx(3)
      pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - flx(4)
      pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - flx(5)            
      
      pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - vf(1)*flx(1)
      pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vf(2)*flx(1)
      pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - vf(3)*flx(1)      
    END IF ! ABS(gr-gl)   
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    SELECT CASE ( pPatch%spaceOrder )
      CASE ( 1 ) 
        CALL RFLU_CentralFirstPatch(pRegion,pPatch)      
      CASE ( 2 ) 
        CALL RFLU_CentralSecondPatch(pRegion,pPatch)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pPatch%spaceOrder
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF 
  CALL FPROFILER_ENDS("RFLU::AUSM_ComputeFlux2_MPSD")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AUSM_ComputeFlux2_MPSD









! ******************************************************************************
!
! Purpose: Compute convective fluxes using second-order accurate AUSM+ scheme 
!   for mixture of thermally and calorically perfect gases.
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

SUBROUTINE RFLU_AUSM_ComputeFlux2_MTCP(pRegion)

#ifdef SPEC
  USE ModSpecies, ONLY: t_spec_type
#endif

  USE ModInterfaces, ONLY: MixtPerf_C_GRT, &
                           MixtPerf_G_CpR, &
                           MixtPerf_Ho_CpTUVW, & 
                           MixtPerf_R_M, & 
                           MixtPerf_T_DPR, & 
                           RFLU_CentralFirstPatch, & 
                           RFLU_CentralSecondPatch

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

  LOGICAL :: interfaceFlag,updateFlux
  INTEGER :: c1,c2,ifg,indCp,indGs,indMf,indMol,indSd,iPatch
  REAL(RFREAL) :: al,ar,cpl,cpr,dx,dy,dz,fs,fsu,gcl,gcr,gl,gr,Hl,Hr,irl,irr, &
                  nm,nTol,nx,ny,nz,pl,pr,rl,rr,tl,tr,ul,ur,vl,vr,wl,wr,xc,yc,zc
  REAL(RFREAL) :: flx(5),vf(3)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pRhs,pSd
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

#ifdef SPEC
  INTEGER :: iSpec,iCvSpecExplosive,iCvSpecProducts,prodInterfaceType
  REAL(RFREAL) :: cpProducts,gcProducts,gProducts,mml,mmr,mwProducts,ud, &
                  YExpl,YExpr,YProdl,YProdr,Y1,Y2
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec 
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGcSpec 
  TYPE(t_spec_type), POINTER :: pSpecType
  TYPE(t_spec_input), POINTER :: pSpecInput
#endif  

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AUSM_ComputeFlux2_MTCP',__FILE__)

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::AUSM_ComputeFlux2_MTCP")
#endif

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs  = pGrid%indGs
  
  indCp  = pRegion%mixtInput%indCp
  indMf  = pRegion%mixtInput%indMfMixt  
  indMol = pRegion%mixtInput%indMol 
  indSd  = pRegion%mixtInput%indSd
  
  pCv  => pRegion%mixt%cv
  pDv  => pRegion%mixt%dv
  pGv  => pRegion%mixt%gv
  
  pGc  => pRegion%mixt%gradCell  
  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs
  pSd  => pRegion%mixt%sd  

#ifdef SPEC
  pCvSpec    => pRegion%spec%cv  
  pGcSpec    => pRegion%spec%gradCell    
  pSpecInput => pRegion%specInput
#endif

  nTol = 1.0E-14_RFREAL

  IF ( pRegion%spec%cvState /= CV_MIXT_STATE_CONS ) THEN 
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! pRegion%spec%cvState  

  IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_TCPERF ) THEN
    CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__)  
  END IF ! pRegion%mixtInput%gasModel

  IF ( (indCp /= 1) .OR. (indMol /= 1) ) THEN 
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
  END IF ! indCp

! ******************************************************************************
! Compute fluxes through interior faces
! ******************************************************************************

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

! ==============================================================================
!   Get face geometry and grid speed
! ==============================================================================

    nx = pGrid%fn(XCOORD,ifg)
    ny = pGrid%fn(YCOORD,ifg)
    nz = pGrid%fn(ZCOORD,ifg)
    nm = pGrid%fn(XYZMAG,ifg)

    xc = pGrid%fc(XCOORD,ifg)
    yc = pGrid%fc(YCOORD,ifg)
    zc = pGrid%fc(ZCOORD,ifg)

    fs = pGrid%gs(indGs*ifg)

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! ------------------------------------------------------------------------------
!   Left state
! ------------------------------------------------------------------------------

    rl  = pCv(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl
    
    dx  = xc - pGrid%cofg(XCOORD,c1)
    dy  = yc - pGrid%cofg(YCOORD,c1)
    dz  = zc - pGrid%cofg(ZCOORD,c1)

#ifdef SPEC
    mml = 0.0_RFREAL
    cpl = 0.0_RFREAL

    DO iSpec = 1,pRegion%specInput%nSpecies
      pSpecType => pRegion%specInput%specType(iSpec)

      Y1 = irl*pCvSpec(iSpec,c1) + pGcSpec(XCOORD,iSpec,c1)*dx & 
                                 + pGcSpec(YCOORD,iSpec,c1)*dy & 
                                 + pGcSpec(ZCOORD,iSpec,c1)*dz

      mml = mml + Y1/pSpecType%pMaterial%molw
      cpl = cpl + Y1*pSpecType%pMaterial%spht
    END DO ! iSpec
    
    mml = 1.0_RFREAL/mml
    gcl = MixtPerf_R_M(mml)
    gl  = MixtPerf_G_CpR(cpl,gcl)    
#endif                   

    ul = pCv(CV_MIXT_XMOM,c1)*irl
    vl = pCv(CV_MIXT_YMOM,c1)*irl
    wl = pCv(CV_MIXT_ZMOM,c1)*irl
    pl = pDv(DV_MIXT_PRES,c1)
        
    rl = rl + pGc(XCOORD,GRC_MIXT_DENS,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_DENS,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_DENS,c1)*dz
    ul = ul + pGc(XCOORD,GRC_MIXT_XVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_XVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_XVEL,c1)*dz
    vl = vl + pGc(XCOORD,GRC_MIXT_YVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_YVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_YVEL,c1)*dz
    wl = wl + pGc(XCOORD,GRC_MIXT_ZVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_ZVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_ZVEL,c1)*dz
    pl = pl + pGc(XCOORD,GRC_MIXT_PRES,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_PRES,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_PRES,c1)*dz    

    tl = MixtPerf_T_DPR(rl,pl,gcl)
    al = MixtPerf_C_GRT(gl,gcl,tl)
    Hl = MixtPerf_Ho_CpTUVW(cpl,tl,ul,vl,wl)

#ifdef SPEC
! DEBUG: Manoj-PBA1D, enthalpy for program burn
    IF ( (global%pbaFlag .EQV. .TRUE.) .AND. &
         (global%pbaBurnFlag .EQV. .TRUE.) ) THEN
      iCvSpecExplosive = SPEC_GetSpeciesIndex(global,pSpecInput,'EXPLOSIVE')
      iCvSpecProducts  = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')

      mwProducts = pSpecInput%specType(iCvSpecProducts)%pMaterial%molw
      cpProducts = pSpecInput%specType(iCvSpecProducts)%pMaterial%spht
      gcProducts = MixtPerf_R_M(mwProducts)
      gProducts  = MixtPerf_G_CpR(cpProducts,gcProducts)

      YExpl  = irl*pRegion%spec%cv(iCvSpecExplosive,c1) &
             + pGcSpec(XCOORD,iCvSpecExplosive,c1)*dx & 
             + pGcSpec(YCOORD,iCvSpecExplosive,c1)*dy & 
             + pGcSpec(ZCOORD,iCvSpecExplosive,c1)*dz

      YProdl = irl*pRegion%spec%cv(iCvSpecProducts,c1) &
             + pGcSpec(XCOORD,iCvSpecProducts,c1)*dx & 
             + pGcSpec(YCOORD,iCvSpecProducts,c1)*dy & 
             + pGcSpec(ZCOORD,iCvSpecProducts,c1)*dz

      IF ( ABS(YExpl) > nTol ) THEN
        IF ( ABS(1.0_RFREAL-YExpl) < nTol ) THEN
          al = SQRT(pl/rl)
          Hl = pCv(CV_MIXT_ENER,c1)*irl + pl/rl
        ELSE
          al = SQRT((pl/rl)*(1.0_RFREAL+(gProducts-1.0_RFREAL)*YProdl))
          Hl = (pl/rl)*(1.0_RFREAL/((gProducts-1.0_RFREAL)*YProdl)+1.0_RFREAL) &
               + 0.5_RFREAL*(ul*ul+vl*vl+wl*wl)
        END IF
      END IF ! YExpl
    END IF ! global%pbaFlag 
! END DEBUG
#endif
! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------

    rr  = pCv(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr
    
    dx  = xc - pGrid%cofg(XCOORD,c2)
    dy  = yc - pGrid%cofg(YCOORD,c2)
    dz  = zc - pGrid%cofg(ZCOORD,c2)

#ifdef SPEC
    mmr = 0.0_RFREAL
    cpr = 0.0_RFREAL


    DO iSpec = 1,pRegion%specInput%nSpecies
      pSpecType => pRegion%specInput%specType(iSpec)

      Y2 = irr*pCvSpec(iSpec,c2) + pGcSpec(XCOORD,iSpec,c2)*dx & 
                                 + pGcSpec(YCOORD,iSpec,c2)*dy & 
                                 + pGcSpec(ZCOORD,iSpec,c2)*dz


      mmr = mmr + Y2/pSpecType%pMaterial%molw
      cpr = cpr + Y2*pSpecType%pMaterial%spht
    END DO ! iSpec
    
    mmr = 1.0_RFREAL/mmr 
    gcr = MixtPerf_R_M(mmr)
    gr  = MixtPerf_G_CpR(cpr,gcr)           
#endif
              
    ur = pCv(CV_MIXT_XMOM,c2)*irr
    vr = pCv(CV_MIXT_YMOM,c2)*irr
    wr = pCv(CV_MIXT_ZMOM,c2)*irr
    pr = pDv(DV_MIXT_PRES,c2)
        
    rr = rr + pGc(XCOORD,GRC_MIXT_DENS,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_DENS,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_DENS,c2)*dz
    ur = ur + pGc(XCOORD,GRC_MIXT_XVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_XVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_XVEL,c2)*dz
    vr = vr + pGc(XCOORD,GRC_MIXT_YVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_YVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_YVEL,c2)*dz
    wr = wr + pGc(XCOORD,GRC_MIXT_ZVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_ZVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_ZVEL,c2)*dz
    pr = pr + pGc(XCOORD,GRC_MIXT_PRES,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_PRES,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_PRES,c2)*dz

    tr = MixtPerf_T_DPR(rr,pr,gcr)
    ar = MixtPerf_C_GRT(gr,gcr,tr)
    Hr = MixtPerf_Ho_CpTUVW(cpr,tr,ur,vr,wr)

#ifdef SPEC
! DEBUG: Manoj-PBA1D, enthalpy for program burn
    IF ( (global%pbaFlag .EQV. .TRUE.) .AND. &
         (global%pbaBurnFlag .EQV. .TRUE.) ) THEN
      iCvSpecExplosive = SPEC_GetSpeciesIndex(global,pSpecInput,'EXPLOSIVE')
      iCvSpecProducts  = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')

      mwProducts = pSpecInput%specType(iCvSpecProducts)%pMaterial%molw
      cpProducts = pSpecInput%specType(iCvSpecProducts)%pMaterial%spht
      gcProducts = MixtPerf_R_M(mwProducts)
      gProducts  = MixtPerf_G_CpR(cpProducts,gcProducts)

      YExpr  = irr*pRegion%spec%cv(iCvSpecExplosive,c2) &
             + pGcSpec(XCOORD,iCvSpecExplosive,c2)*dx & 
             + pGcSpec(YCOORD,iCvSpecExplosive,c2)*dy & 
             + pGcSpec(ZCOORD,iCvSpecExplosive,c2)*dz

      YProdr = irr*pRegion%spec%cv(iCvSpecProducts,c2) &
             + pGcSpec(XCOORD,iCvSpecProducts,c2)*dx & 
             + pGcSpec(YCOORD,iCvSpecProducts,c2)*dy & 
             + pGcSpec(ZCOORD,iCvSpecProducts,c2)*dz

      IF ( ABS(YExpr) > nTol ) THEN
        IF ( ABS(1.0_RFREAL-YExpr) < nTol ) THEN
          ar = SQRT(pr/rr)
          Hr = pCv(CV_MIXT_ENER,c2)*irr + pr/rr
        ELSE
          ar = SQRT((pr/rr)*(1.0_RFREAL+(gProducts-1.0_RFREAL)*YProdr))
          Hr = (pr/rr)*(1.0_RFREAL/((gProducts-1.0_RFREAL)*YProdr)+1.0_RFREAL) &
               + 0.5_RFREAL*(ur*ur+vr*vr+wr*wr)
        END IF
      END IF ! YExpl
    END IF ! global%pbaFlag 
! END DEBUG
#endif
! ==============================================================================
!   Compute fluxes 
! ==============================================================================

    CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur,vr, &
                                wr,pr,Hr,ar,flx,vf)

! ==============================================================================
!   Program burn: If c1 or c2 contains YExplosive>0, then treat face as slipwall
! ==============================================================================

    updateFlux = .TRUE.

#ifdef SPEC
! DEBUG: Manoj-PBA, commenting out flux treatment at explosive-air, explosive-product faces 
IF (1==1) THEN
    IF ( (global%pbaFlag .EQV. .TRUE.) .AND. &
         (global%pbaBurnFlag .EQV. .TRUE.)) THEN      
      iCvSpecExplosive = SPEC_GetSpeciesIndex(global,pSpecInput,'EXPLOSIVE')
      iCvSpecProducts  = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')
      YExpl  = irl*pRegion%spec%cv(iCvSpecExplosive,c1)
      YExpr  = irr*pRegion%spec%cv(iCvSpecExplosive,c2)
      YProdl = irl*pRegion%spec%cv(iCvSpecProducts,c1)
      YProdr = irr*pRegion%spec%cv(iCvSpecProducts,c2)

! ------------------------------------------------------------------------------
!     Identify the interface between explosive and non-explosive  
! ------------------------------------------------------------------------------

      interfaceFlag = .FALSE.

      IF ( (ABS(YExpl) >= nTol) .AND. (ABS(YExpr) < nTol) ) THEN
        interfaceFlag = .TRUE.
      END IF ! YExpl, YExpr

      IF ( (ABS(YExpl) < nTol) .AND. (ABS(YExpr) >= nTol) ) THEN
        interfaceFlag = .TRUE.
      END IF ! YExpl, YExpr

      prodInterfaceType = 0

      IF ( ABS(1.0_RFREAL-YProdl) < nTol ) THEN
        prodInterfaceType = 1
      END IF ! YProdl

      IF ( ABS(1.0_RFREAL-YProdr) < nTol ) THEN
        prodInterfaceType = 2
      END IF ! YProdr

! ------------------------------------------------------------------------------
!     Compute fluxes if it is interface between explosive and non-explosive  
! ------------------------------------------------------------------------------

      IF ( interfaceFlag .EQV. .TRUE. ) THEN
        fsu = RFLU_DescaleGridSpeed(pRegion,fs)

! ----- Its an Air-Explosive interface -----------------------------------------        
        IF ( prodInterfaceType == 0 ) THEN
          ul = 0.0_RFREAL
          vl = 0.0_RFREAL
          wl = 0.0_RFREAL
          ur = 0.0_RFREAL
          vr = 0.0_RFREAL
          wr = 0.0_RFREAL

          Hl  = MixtPerf_Ho_CpTUVW(cpl,tl,ul,vl,wl)
          Hr  = MixtPerf_Ho_CpTUVW(cpr,tr,ur,vr,wr)
        END IF ! prodInterfaceType

! ----- Product is on left side ------------------------------------------------        
        IF ( prodInterfaceType == 1 ) THEN
          rr = rl
          ur = ul
          vr = vl
          wr = wl
          pr = pl
          ar = al
          Hr = Hl
        END IF ! prodInterfaceType

! ----- Product is on right side -----------------------------------------------        
        IF ( prodInterfaceType == 2 ) THEN
          rl = rr
          ul = ur
          vl = vr
          wl = wr
          pl = pr
          al = ar
          Hl = Hr
        END IF ! prodInterfaceType

! ------------------------------------------------------------------------------
!       Left state, compute flux with pr = pl
! ------------------------------------------------------------------------------

        CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al, &
                                    rr,ur,vr,wr,pl,Hr,ar,flx,vf)

! ----  Accumulate into residual of left cell ----------------------------------

        pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + flx(1)
        pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + flx(2)
        pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + flx(3)
        pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + flx(4)
        pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + flx(5)

! ----  Accumulate into substantial derivative of left cell --------------------
    
        pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + vf(1)*flx(1)
        pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vf(2)*flx(1)
        pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + vf(3)*flx(1)

! ------------------------------------------------------------------------------
!       Right state, compute flux with velocities at face=0, and pl = pr
! ------------------------------------------------------------------------------

        CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pr,Hl,al, &
                                    rr,ur,vr,wr,pr,Hr,ar,flx,vf)

! ----  Accumulate into residual of right cell ---------------------------------
    
        pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - flx(1)
        pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - flx(2)
        pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - flx(3)
        pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - flx(4)
        pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - flx(5)
    
! ----  Accumulate into substantial derivative of right cell -------------------
    
        pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - vf(1)*flx(1)
        pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vf(2)*flx(1)
        pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - vf(3)*flx(1)   
    
! ------------------------------------------------------------------------------
!       Store mass flux
! ------------------------------------------------------------------------------

! TEMPORARY: Manoj-PBA, making mfMixt=0 so that species cv is not modified
        pMf(indMf*ifg) = 0.0_RFREAL
! END TEMPORARY

! ------------------------------------------------------------------------------
!       Set flag for not updating flux later
! ------------------------------------------------------------------------------

        updateFlux = .FALSE.

      END IF ! interfaceFlag
    END IF ! global%pbaFlag
END IF ! 1==2
! END DEBUG
#endif

! ==============================================================================
!   Program burn: Update mass flux, residual, and substantial detivative only
!   when face is not in contact with explosive region with YExplosive=1
! ==============================================================================

    IF ( updateFlux ) THEN

! ==============================================================================
!     Store mass flux
! ==============================================================================

      pMf(indMf*ifg) = flx(1)

! ==============================================================================
!     Accumulate into residual
! ==============================================================================

      pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + flx(1)
      pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + flx(2)
      pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + flx(3)
      pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + flx(4)
      pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + flx(5)

      pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - flx(1)
      pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - flx(2)
      pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - flx(3)
      pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - flx(4)
      pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - flx(5)
    
! ==============================================================================
!     Accumulate into substantial derivative
! ==============================================================================

      pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + vf(1)*flx(1)
      pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vf(2)*flx(1)
      pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + vf(3)*flx(1)
    
      pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - vf(1)*flx(1)
      pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vf(2)*flx(1)
      pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - vf(3)*flx(1)   
    END IF ! updateFlux
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    SELECT CASE ( pPatch%spaceOrder )
      CASE ( 1 ) 
        CALL RFLU_CentralFirstPatch(pRegion,pPatch)      
      CASE ( 2 ) 
        CALL RFLU_CentralSecondPatch(pRegion,pPatch)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pPatch%spaceOrder
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF 
  CALL FPROFILER_ENDS("RFLU::AUSM_ComputeFlux2_MTCP")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AUSM_ComputeFlux2_MTCP








! ******************************************************************************
!
! Purpose: Compute convective fluxes using second-order accurate AUSM+ scheme
!   for thermally and calorically perfect gas.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!  1. Only applicable to thermally and calorically perfect gas because would
!     need to recompute mixture properties at faces depending on face states.
!
! ******************************************************************************

SUBROUTINE RFLU_AUSM_ComputeFlux2_TCP(pRegion)

  USE ModInterfaces, ONLY: MixtPerf_C_GRT, &
                           MixtPerf_Ho_CpTUVW, & 
                           MixtPerf_R_CpG, & 
                           MixtPerf_T_DPR, & 
                           RFLU_CentralFirstPatch, & 
                           RFLU_CentralSecondPatch

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

  INTEGER :: c1,c2,ifg,indGs,indMf,indSd,iPatch
  REAL(RFREAL) :: al,ar,cp,dx,dy,dz,fs,gc,g,Hl,Hr,irl,irr,nm,nx,ny,nz,pl,pr, &
                  rl,rr,tl,tr,ul,ur,vl,vr,wl,wr,xc,yc,zc
  REAL(RFREAL) :: flx(5),vf(3)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pRhs,pSd
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AUSM_ComputeFlux2_TCP',__FILE__)

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::AUSM_ComputeFlux2_TCP")
#endif

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs = pGrid%indGs
  
  indMf = pRegion%mixtInput%indMfMixt  
  indSd = pRegion%mixtInput%indSd
  
  pCv  => pRegion%mixt%cv
  pDv  => pRegion%mixt%dv
  
  pGc  => pRegion%mixt%gradCell  
  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs
  pSd  => pRegion%mixt%sd  

  IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_TCPERF ) THEN
    CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__)  
  END IF ! pRegion%mixtInput%gasModel

  cp = global%refCp
  g  = global%refGamma  
  gc = MixtPerf_R_CpG(cp,g)    

! ******************************************************************************
! Compute fluxes through interior faces
! ******************************************************************************

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

! ==============================================================================
!   Get face geometry and grid speed
! ==============================================================================

    nx = pGrid%fn(XCOORD,ifg)
    ny = pGrid%fn(YCOORD,ifg)
    nz = pGrid%fn(ZCOORD,ifg)
    nm = pGrid%fn(XYZMAG,ifg)

    xc = pGrid%fc(XCOORD,ifg)
    yc = pGrid%fc(YCOORD,ifg)
    zc = pGrid%fc(ZCOORD,ifg)

    fs = pGrid%gs(indGs*ifg)

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! ------------------------------------------------------------------------------
!   Left state
! ------------------------------------------------------------------------------

    rl  = pCv(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl

    ul  = pCv(CV_MIXT_XMOM,c1)*irl
    vl  = pCv(CV_MIXT_YMOM,c1)*irl
    wl  = pCv(CV_MIXT_ZMOM,c1)*irl
    pl  = pDv(DV_MIXT_PRES,c1)
    
    dx  = xc - pGrid%cofg(XCOORD,c1)
    dy  = yc - pGrid%cofg(YCOORD,c1)
    dz  = zc - pGrid%cofg(ZCOORD,c1)
        
    rl = rl + pGc(XCOORD,GRC_MIXT_DENS,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_DENS,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_DENS,c1)*dz
    ul = ul + pGc(XCOORD,GRC_MIXT_XVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_XVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_XVEL,c1)*dz
    vl = vl + pGc(XCOORD,GRC_MIXT_YVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_YVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_YVEL,c1)*dz
    wl = wl + pGc(XCOORD,GRC_MIXT_ZVEL,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_ZVEL,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_ZVEL,c1)*dz
    pl = pl + pGc(XCOORD,GRC_MIXT_PRES,c1)*dx &
            + pGc(YCOORD,GRC_MIXT_PRES,c1)*dy &
            + pGc(ZCOORD,GRC_MIXT_PRES,c1)*dz

    ! Subbu - vFracE Correction
    IF ( global%plagUsed .AND. (pRegion%plag%nPcls > 1) ) THEN
      tl = MixtPerf_T_DPR( rl/(1.0_RFREAL - pRegion%plag%vFracE(1,c1)), pl,gc)
    ELSE 
      tl = MixtPerf_T_DPR(rl,pl,gc)
    END IF ! global%plagUsed
    ! Subbu - vFracE Correction

    !tl = MixtPerf_T_DPR(rl,pl,gc)
    al = MixtPerf_C_GRT(g,gc,tl)
    Hl = MixtPerf_Ho_CpTUVW(cp,tl,ul,vl,wl)
    
! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------

    rr  = pCv(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr

    ur  = pCv(CV_MIXT_XMOM,c2)*irr
    vr  = pCv(CV_MIXT_YMOM,c2)*irr
    wr  = pCv(CV_MIXT_ZMOM,c2)*irr
    pr  = pDv(DV_MIXT_PRES,c2)

    dx  = xc - pGrid%cofg(XCOORD,c2)
    dy  = yc - pGrid%cofg(YCOORD,c2)
    dz  = zc - pGrid%cofg(ZCOORD,c2)
        
    rr = rr + pGc(XCOORD,GRC_MIXT_DENS,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_DENS,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_DENS,c2)*dz
    ur = ur + pGc(XCOORD,GRC_MIXT_XVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_XVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_XVEL,c2)*dz
    vr = vr + pGc(XCOORD,GRC_MIXT_YVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_YVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_YVEL,c2)*dz
    wr = wr + pGc(XCOORD,GRC_MIXT_ZVEL,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_ZVEL,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_ZVEL,c2)*dz
    pr = pr + pGc(XCOORD,GRC_MIXT_PRES,c2)*dx &
            + pGc(YCOORD,GRC_MIXT_PRES,c2)*dy &
            + pGc(ZCOORD,GRC_MIXT_PRES,c2)*dz
    
    ! Subbu - vFracE Correction
    IF ( global%plagUsed .AND. (pRegion%plag%nPcls > 1) ) THEN
      tr = MixtPerf_T_DPR( rr/(1.0_RFREAL - pRegion%plag%vFracE(1,c2)),pr,gc)
    ELSE 
      tr = MixtPerf_T_DPR(rr,pr,gc)
    END IF ! global%plagUsed
    ! Subbu - vFracE Correction

    !tr = MixtPerf_T_DPR(rr,pr,gc)
    ar = MixtPerf_C_GRT(g,gc,tr)
    Hr = MixtPerf_Ho_CpTUVW(cp,tr,ur,vr,wr)
    
! ==============================================================================
!   Compute fluxes 
! ==============================================================================

    CALL RFLU_AUSM_FluxFunction(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur,vr, & 
                                wr,pr,Hr,ar,flx,vf)

! ==============================================================================
!   Store mass flux
! ==============================================================================

    pMf(indMf*ifg) = flx(1)

! ==============================================================================
!   Accumulate into residual
! ==============================================================================

    pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + flx(1)
    pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + flx(2)
    pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + flx(3)
    pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + flx(4)
    pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + flx(5)

    pRhs(CV_MIXT_DENS,c2) = pRhs(CV_MIXT_DENS,c2) - flx(1)
    pRhs(CV_MIXT_XMOM,c2) = pRhs(CV_MIXT_XMOM,c2) - flx(2)
    pRhs(CV_MIXT_YMOM,c2) = pRhs(CV_MIXT_YMOM,c2) - flx(3)
    pRhs(CV_MIXT_ZMOM,c2) = pRhs(CV_MIXT_ZMOM,c2) - flx(4)
    pRhs(CV_MIXT_ENER,c2) = pRhs(CV_MIXT_ENER,c2) - flx(5)

! ==============================================================================
!   Accumulate into substantial derivative
! ==============================================================================

    pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + vf(1)*flx(1)
    pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vf(2)*flx(1)
    pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + vf(3)*flx(1)
    
    pSd(SD_XMOM,c2*indSd) = pSd(SD_XMOM,c2*indSd) - vf(1)*flx(1)
    pSd(SD_YMOM,c2*indSd) = pSd(SD_YMOM,c2*indSd) - vf(2)*flx(1)
    pSd(SD_ZMOM,c2*indSd) = pSd(SD_ZMOM,c2*indSd) - vf(3)*flx(1)   
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    SELECT CASE ( pPatch%spaceOrder )
      CASE ( 1 ) 
        CALL RFLU_CentralFirstPatch(pRegion,pPatch)      
      CASE ( 2 ) 
        CALL RFLU_CentralSecondPatch(pRegion,pPatch)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pPatch%spaceOrder
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF 
  CALL FPROFILER_ENDS("RFLU::AUSM_ComputeFlux2_TCP")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AUSM_ComputeFlux2_TCP






  
! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModAUSMFlux


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModAUSMFlux.F90,v $
! Revision 1.2  2016/02/04 22:20:29  fred
! Adding JWL EOS iterative capabilities for cylindrical detonation problem
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
! Revision 1.2  2008/11/19 22:16:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:23  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:39  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.9  2006/05/01 22:20:28  haselbac
! Removed debug statements
!
! Revision 1.8  2006/05/01 21:00:47  haselbac
! Rewrite for consistency and cleanliness
!
! Revision 1.7  2006/04/15 16:58:42  haselbac
! Added capability of running 1st order boundary fluxes with 2nd order volume 
! fluxes
!
! Revision 1.6  2006/04/07 15:19:18  haselbac
! Removed tabs
!
! Revision 1.5  2005/11/15 17:25:07  haselbac
! Bug fix: Missing declarations lead to compilation failure without SPEC=1
!
! Revision 1.4  2005/11/14 16:58:05  haselbac
! Added flux for pseudo-gas, reordered routines
!
! Revision 1.3  2005/11/10 02:26:06  haselbac
! Complete rewrite to allow for computations with variable properties
!
! Revision 1.2  2005/07/19 20:06:29  haselbac
! Modified af to prevent crashes for Skews problem
!
! Revision 1.1  2005/07/14 21:39:02  haselbac
! Initial revision
!
! ******************************************************************************
  
  

