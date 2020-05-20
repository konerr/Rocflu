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
! Purpose: Collection of AUSM+up flux routines.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: INRT_DiscreteDelta.F90,v 1.3 2016/02/05 16:42:23 rahul Exp $
!
! Copyright: (c) 2005-2006 by the University of Illinois
!
! ******************************************************************************

MODULE INRT_DiscreteDelta

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
  PUBLIC :: RFLU_AUSMPlusUp_ComputeFlux
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: INRT_DiscreteDelta.F90,v $ $Revision: 1.3 $' 
              
! *****************************************************************************
! Routines
! *****************************************************************************

  CONTAINS






! ******************************************************************************
!
! Purpose: Wrapper function for AUSM+up flux functions.
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

!SUBROUTINE RFLU_AUSMPlusUp_ComputeFlux(pRegion)
!                       
!  IMPLICIT NONE
!
!! ******************************************************************************
!! Definitions and declarations
!! ******************************************************************************
!
!! ==============================================================================
!! Arguments
!! ==============================================================================
!
!  TYPE(t_region), POINTER :: pRegion
!
!! ==============================================================================
!! Locals
!! ==============================================================================
!  
!  TYPE(t_global), POINTER :: global
!
!! ******************************************************************************
!! Start
!! ******************************************************************************
!
!  global => pRegion%global
!
!  CALL RegisterFunction(global,'RFLU_AUSMPlusUp_ComputeFlux',__FILE__)
!
!! ******************************************************************************
!! Call flux functions
!! ******************************************************************************
!
!  SELECT CASE ( pRegion%mixtInput%gasModel ) 
!
!! ==============================================================================
!!   Thermally and calorically perfect gas
!! ==============================================================================
!
!    CASE ( GAS_MODEL_TCPERF )
!      SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
!        CASE ( DISCR_ORDER_1 ) 
!          CALL INRT_CalcDiscreteDelta(pRegion)
!        CASE ( DISCR_ORDER_2 )
!          CALL RFLU_AUSMPlusUp_ComputeFlux2_TCP(pRegion)               
!        CASE DEFAULT
!          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
!      END SELECT ! pRegion%mixtInput%spaceOrder    
!
!! ==============================================================================
!!   Mixture of thermally and calorically perfect gases
!! ==============================================================================
!
!    CASE ( GAS_MODEL_MIXT_TCPERF ) 
!      SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
!        CASE ( DISCR_ORDER_1 ) 
!          CALL INRT_CalcDiscreteDelta(pRegion)
!        CASE ( DISCR_ORDER_2 )
!          CALL RFLU_AUSMPlusUp_ComputeFlux2_MTCP(pRegion)               
!        CASE DEFAULT
!          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
!      END SELECT ! pRegion%mixtInput%spaceOrder
!
!! ==============================================================================
!!   Pseudo-gas
!! ==============================================================================
!
!    CASE ( GAS_MODEL_MIXT_PSEUDO ) 
!      SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
!! TO DO
!!        CASE ( DISCR_ORDER_1 ) 
!!          CALL INRT_CalcDiscreteDelta_MPSD(pRegion)
!! END TO DO 
!        CASE ( DISCR_ORDER_2 )
!          CALL RFLU_AUSMPlusUp_ComputeFlux2_MPSD(pRegion)               
!        CASE DEFAULT
!          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
!      END SELECT ! pRegion%mixtInput%spaceOrder                           
!
!! ==============================================================================
!!   Mixture of thermally and calorically perfect gas and detonation products
!! ==============================================================================
!
!    CASE ( GAS_MODEL_MIXT_JWL ) 
!      SELECT CASE ( pRegion%mixtInput%spaceOrder ) 
!        CASE ( DISCR_ORDER_1 ) 
!          CALL INRT_CalcDiscreteDelta_MJWL(pRegion)
!        CASE ( DISCR_ORDER_2 )
!          CALL RFLU_AUSMPlusUp_ComputeFlux2_MJWL(pRegion)               
!        CASE DEFAULT
!          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
!      END SELECT ! pRegion%mixtInput%spaceOrder
!
!! ==============================================================================
!!   Default
!! ==============================================================================
!
!    CASE DEFAULT
!      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)             
!  END SELECT ! pRegion%mixtInput%gasModel       
!
!! ******************************************************************************
!! End
!! ******************************************************************************
!
!  CALL DeregisterFunction(global)
!
!END SUBROUTINE RFLU_AUSMPlusUp_ComputeFlux







! ******************************************************************************
!
! Purpose: Compute convective fluxes using AUSMPlusUp+ scheme.
!
! Description: None.
!
! Input: 
!   r          x-component of face normal
!
! Output: 
!   psi         Fluxes
!
! Notes: 
!   1. Liou M.-S., Progress towards an improved CFD method: AUSMPlusUp+, AIAA Paper
!      95-1701, 1995
!
! ******************************************************************************

SUBROUTINE INRT_ComputePsi(r,psi)
                          

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  REAL(RFREAL), INTENT(IN) :: r
                              
  REAL(RFREAL), INTENT(OUT) :: psi

! ==============================================================================
! Locals
! ==============================================================================

  REAL(RFREAL) :: af,mf,mfa,mfm,mfp,ml,mla,mlp,mr,mra,mrm,ql,qr,vml,vmr, &
                  wtl,wtr,mbar,mo,minf,fa,sigma,mp,kp,pu,ku,rf,alpha,mbarsq, &
                  mosq,sigmabar,mmax

! ******************************************************************************
! Start, compute face state
! ******************************************************************************
    
  oneThird = 1.0_RFREAL/3.0_RFREAL
  oneSixth = 1.0_RFREAL/6.0_RFREAL

  IF ( (ABS(r) .LE. 1.5_RFREAL) .OR. &
       (ABS(r) .GE. 0.5_RFREAL) ) THEN 
    psi = 5.0_RFREAL-3.0_RFREAL*ABS(r)- &
          DSQRT( 1.0_RFREAL-3.0_RFREAL* &
          (1.0_RFREAL-ABS(r))*(1.0_RFREAL-ABS(r)) )
    psi = oneSixth*r
  ELSEIF ( ABS(r) .LE. 0.5_RFREAL ) THEN
    psi = 1.0_RFREAL + DSQRT(1-3.0_RFREAL*r*r)
    psi = oneThird*r
  ELSE
    psi = 0.0_RFREAL 
  END IF ! ABS(r)

! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE INRT_ComputePsi



SUBROUTINE INRT_CalcDiscreteDelta(pRegion)

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
                  pl,pr,rl,rr,rul,rur,rvl,rvr,rwl,rwr,tl,tr,ul,ur,vl,vr,wl,wr, &
                  vfracGl,vfracGr,vFracF,pf
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

  CALL RegisterFunction(global,'INRT_CalcDiscreteDelta',__FILE__)

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::AUSMPlusUp_ComputeFlux1")
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

    vfracGl = 1.0_RFREAL - pRegion%plag%vFracE(1,c1) ! rahul
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
    vfracGr = 1.0_RFREAL - pRegion%plag%vFracE(1,c2) ! rahul

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

    CALL INRT_ComputePsi1(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur,vr, &
                                 wr,pr,Hr,ar,flx,vf)

    WRITE(*,'(A,I4,3(1X,E23.16))') "ifg=",ifg,flx(1),flx(2),flx(5)
    WRITE(*,*) "--------------------------------------------------------------------------"
    WRITE(*,*) " "
  ELSE
    CALL INRT_ComputePsi2(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur,vr, &
                                wr,pr,Hr,ar,vFracF,flx,vf)
  END IF ! ifg
ELSE
    CALL INRT_ComputePsi2(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur,vr, &
                                wr,pr,Hr,ar,vFracF,flx,vf)
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

        CALL INRT_ComputePsi(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al, &
                                 rr,ur,vr,wr,pl,Hr,ar,flx,vf,vfracGl,vfracGr,pf)

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

        CALL INRT_ComputePsi(nx,ny,nz,nm,fs,rl,ul,vl,wl,pr,Hl,al, &
                                 rr,ur,vr,wr,pr,Hr,ar,flx,vf,vfracGl,vfracGr,pf)

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
  CALL FPROFILER_ENDS("RFLU::AUSMPlusUp_ComputeFlux1")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE INRT_CalcDiscreteDelta







! ******************************************************************************
!
! Purpose: Compute convective fluxes using first-order accurate AUSMPlusUp+ scheme 
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

SUBROUTINE INRT_CalcDiscreteDelta_MJWL(pRegion)

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
                  tl,tr,ul,ur,vl,vr,wl,wr,vfracGl,vfracGr,pf
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

  CALL RegisterFunction(global,'INRT_CalcDiscreteDelta_MJWL',__FILE__)

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::AUSMPlusUp_ComputeFlux1_MJWL")
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

    CALL INRT_ComputePsi(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur,vr, &
                                wr,pr,Hr,ar,flx,vf,vfracGl,vfracGr,pf)

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
  CALL FPROFILER_ENDS("RFLU::AUSMPlusUp_ComputeFlux1_MJWL")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE INRT_CalcDiscreteDelta_MJWL








! ******************************************************************************
!
! Purpose: Compute convective fluxes using second-order accurate AUSMPlusUp+ scheme 
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

SUBROUTINE RFLU_AUSMPlusUp_ComputeFlux2_MJWL(pRegion)

  USE ModInterfaces, ONLY: MixtPerf_C_GRT, &
                           MixtPerf_G_CpR, &
                           MixtPerf_Ho_CpTUVW, & 
                           MixtPerf_R_M, & 
                           MixtPerf_T_DPR, & 
                           RFLU_CentralFirstPatch, & 
                           RFLU_CentralSecondPatch

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

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: c1,c2,ifg,indCp,indGs,indMf,indMol,indSd,iPatch
  REAL(RFREAL) :: al,ar,cp,cpl,cpr,dx,dy,dz,el,er,fs,g,gc,gcl,gcr,gl,gr,Hl,Hr, &
                  irl,irr,mw,nm,nx,ny,nz,pl,pr,rl,rr,tl,tr,ul,ur,vl,vr,wl,wr, &
                  xc,yc,zc,vfracGl,vfracGr,pf,vFracF
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

#ifdef PLAG
  REAL(RFREAL) :: rcL,rcR,wtL,wtR
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AUSMPlusUp_ComputeFlux2_MJWL',__FILE__)

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::AUSMPlusUp_ComputeFlux2_MJWL")
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

  vFracF = 1.0_RFREAL

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

#ifdef PLAG
    IF ( global%plagUsed .AND. (pRegion%plag%nPcls .GE. 1) ) THEN
     rl  = pCv(CV_MIXT_DENS,c1)/(1.0_RFREAL - pRegion%plag%vFracE(1,c1))
    END IF
#endif

!    irl = 1.0_RFREAL/rl
    irl = 1.0_RFREAL/pCv(CV_MIXT_DENS,c1)
    
    dx  = xc - pGrid%cofg(XCOORD,c1)
    dy  = yc - pGrid%cofg(YCOORD,c1)
    dz  = zc - pGrid%cofg(ZCOORD,c1)

#ifdef PLAG
    rcL = (dx*dx + dy*dy + dz*dz)**0.5_RFREAL
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

#ifdef PLAG
    IF ( global%plagUsed .AND. (pRegion%plag%nPcls .GE. 1) ) THEN
     rr  = pCv(CV_MIXT_DENS,c2)/(1.0_RFREAL - pRegion%plag%vFracE(1,c2))
    END IF
#endif

!    irl = 1.0_RFREAL/rl
    irr = 1.0_RFREAL/pCv(CV_MIXT_DENS,c2)
    
    dx  = xc - pGrid%cofg(XCOORD,c2)
    dy  = yc - pGrid%cofg(YCOORD,c2)
    dz  = zc - pGrid%cofg(ZCOORD,c2)

#ifdef PLAG
    rcR = (dx*dx + dy*dy + dz*dz)**0.5_RFREAL
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

! Rahul - Compute volume fraction at the face by a weighted average. Weights are
! inversely proprtional to the distance between face and cell centroid.

#ifdef PLAG
     IF ( global%plagUsed .AND. (pRegion%plag%nPcls .GE. 1) ) THEN
      rcL = 1.0_RFREAL/rcL
      rcR = 1.0_RFREAL/rcR
 
      wtL = rcL/(rcL+rcR)
      wtR = rcR/(rcL+rcR)
 
      vFracGl = 1.0_RFREAL-pRegion%plag%vFracE(1,c1)
      vFracGr = 1.0_RFREAL-pRegion%plag%vFracE(1,c2)
!      vFracF  = 0.5_RFREAL*(vFracGl + vFacGr)
      vFracF  = wtL*vFracGl + wtR*vFracGr
     END IF

     IF ((vFracF .LT. 0) .OR. (vFracF .GT. 1)) THEN
      write(STDOUT,*) 'face #', ifg
      write(STDOUT,*) 'c1, c2', c1,c2
      write(STDOUT,*) 'phiF',vFracF
      CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__,'Invalid volume fraction')
     END IF
#endif

! Rahul - end

! ==============================================================================
!   Compute fluxes 
! ==============================================================================

!   CALL INRT_ComputePsi(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur,vr, &
!                               wr,pr,Hr,ar,flx,vf,vfracGl,vfracGr,pf)
    CALL INRT_ComputePsi2(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,&
                                        ur,vr,wr,pr,Hr,ar,vFracF,flx,vf)

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
  CALL FPROFILER_ENDS("RFLU::AUSMPlusUp_ComputeFlux2_MJWL")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AUSMPlusUp_ComputeFlux2_MJWL








! ******************************************************************************
!
! Purpose: Compute convective fluxes using second-order accurate AUSMPlusUp+ scheme 
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

SUBROUTINE RFLU_AUSMPlusUp_ComputeFlux2_MPSD(pRegion)

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
                  mmr,nm,nx,ny,nz,pl,pr,rl,rr,tl,tr,ul,ur,vl,vr,wl,wr,xc,yc,zc,&
                  vfracGl,vfracGr,pf
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

  CALL RegisterFunction(global,'RFLU_AUSMPlusUp_ComputeFlux2_MPSD',__FILE__)

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::AUSMPlusUp_ComputeFlux2_MPSD")
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
      CALL INRT_ComputePsi(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al, & 
                               rr,ur,vr,wr,pr,Hr,ar,flx,vf,vfracGl,vfracGr,pf)

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
      CALL INRT_ComputePsi(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al, & 
                               rr,ur,vr,wr,pr,Hr,ar,flx,vf,vfracGl,vfracGr,pf)

      pMf(indMf*ifg) = 0.5_RFREAL*flx(1)

      pRhs(CV_MIXT_DENS,c1) = pRhs(CV_MIXT_DENS,c1) + flx(1)
      pRhs(CV_MIXT_XMOM,c1) = pRhs(CV_MIXT_XMOM,c1) + flx(2)
      pRhs(CV_MIXT_YMOM,c1) = pRhs(CV_MIXT_YMOM,c1) + flx(3)
      pRhs(CV_MIXT_ZMOM,c1) = pRhs(CV_MIXT_ZMOM,c1) + flx(4)
      pRhs(CV_MIXT_ENER,c1) = pRhs(CV_MIXT_ENER,c1) + flx(5)

      pSd(SD_XMOM,c1*indSd) = pSd(SD_XMOM,c1*indSd) + vf(1)*flx(1)
      pSd(SD_YMOM,c1*indSd) = pSd(SD_YMOM,c1*indSd) + vf(2)*flx(1)
      pSd(SD_ZMOM,c1*indSd) = pSd(SD_ZMOM,c1*indSd) + vf(3)*flx(1)
      
      CALL INRT_ComputePsi(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al, & 
                              rr,ur,vr,wr,pr,Hr,ar,flx,vf,vfracGl,vfracGr,pf)

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
  CALL FPROFILER_ENDS("RFLU::AUSMPlusUp_ComputeFlux2_MPSD")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AUSMPlusUp_ComputeFlux2_MPSD









! ******************************************************************************
!
! Purpose: Compute convective fluxes using second-order accurate AUSMPlusUp+ scheme 
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

SUBROUTINE RFLU_AUSMPlusUp_ComputeFlux2_MTCP(pRegion)

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
                  nm,nTol,nx,ny,nz,pl,pr,rl,rr,tl,tr,ul,ur,vl,vr,wl,wr,xc,yc,zc,&
                  vfracGl,vfracGr,pf
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

  CALL RegisterFunction(global,'RFLU_AUSMPlusUp_ComputeFlux2_MTCP',__FILE__)

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::AUSMPlusUp_ComputeFlux2_MTCP")
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

    CALL INRT_ComputePsi(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,ur,vr, &
                                wr,pr,Hr,ar,flx,vf,vfracGl,vfracGr,pf)

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

        CALL INRT_ComputePsi(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al, &
                                rr,ur,vr,wr,pl,Hr,ar,flx,vf,vfracGl,vfracGr,pf)

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

        CALL INRT_ComputePsi(nx,ny,nz,nm,fs,rl,ul,vl,wl,pr,Hl,al, &
                                 rr,ur,vr,wr,pr,Hr,ar,flx,vf,vfracGl,vfracGr,pf)

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
  CALL FPROFILER_ENDS("RFLU::AUSMPlusUp_ComputeFlux2_MTCP")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AUSMPlusUp_ComputeFlux2_MTCP








! ******************************************************************************
!
! Purpose: Compute convective fluxes using second-order accurate AUSMPlusUp+ scheme
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

SUBROUTINE RFLU_AUSMPlusUp_ComputeFlux2_TCP(pRegion)

  USE ModInterfaces, ONLY: MixtPerf_C_GRT, &
                           MixtPerf_Ho_CpTUVW, & 
                           MixtPerf_R_CpG, & 
                           MixtPerf_T_DPR, & 
                           RFLU_CentralFirstPatch, & 
                           RFLU_CentralSecondPatch
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

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: c1,c2,ifg,indGs,indMf,indSd,iPatch
  REAL(RFREAL) :: al,ar,cp,dx,dy,dz,fs,gc,g,Hl,Hr,irl,irr,nm,nx,ny,nz,pl,pr, &
                  rl,rr,tl,tr,ul,ur,vl,vr,wl,wr,xc,yc,zc,vFracF,vFracGr,vFracGl
  REAL(RFREAL) :: flx(5),vf(3)
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pRhs,pSd
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

#ifdef PLAG
! Plag variabless
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGvolFracE,pGvolFracEg
  TYPE(t_plag), POINTER :: pPlag
  REAL(RFREAL) :: rcL,rcR,wtL,wtR
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AUSMPlusUp_ComputeFlux2_TCP',__FILE__)

#ifdef ROCPROF 
  CALL FPROFILER_BEGINS("RFLU::AUSMPlusUp_ComputeFlux2_TCP")
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

! rahul-pointer to gradient of Eul. vol. fraction  
#ifdef PLAG
  pPlag => pRegion%plag  ! rahul
  pGvolFracE=> pPlag%gradVFracE  
  pGvolFracEg=> pPlag%gradVFracEg  
#endif
    
  pMf  => pRegion%mixt%mfMixt  
  pRhs => pRegion%mixt%rhs
  pSd  => pRegion%mixt%sd  

  IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_TCPERF ) THEN
    CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__)  
  END IF ! pRegion%mixtInput%gasModel

  cp = global%refCp
  g  = global%refGamma  
  gc = MixtPerf_R_CpG(cp,g)    

! For gas only simulations
  vFracF  = 1.0_RFREAL

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
 
#ifdef PLAG
    IF ( global%plagUsed .AND. (pRegion%plag%nPcls .GE. 1) ) THEN
     rl  = pCv(CV_MIXT_DENS,c1)/(1.0_RFREAL - pRegion%plag%vFracE(1,c1))
    END IF
#endif

!    irl = 1.0_RFREAL/rl
    irl = 1.0_RFREAL/pCv(CV_MIXT_DENS,c1)

    dx  = xc - pGrid%cofg(XCOORD,c1)
    dy  = yc - pGrid%cofg(YCOORD,c1)
    dz  = zc - pGrid%cofg(ZCOORD,c1)

#ifdef PLAG
    rcL = (dx*dx + dy*dy + dz*dz)**0.5_RFREAL
#endif

    ul  = pCv(CV_MIXT_XMOM,c1)*irl
    vl  = pCv(CV_MIXT_YMOM,c1)*irl
    wl  = pCv(CV_MIXT_ZMOM,c1)*irl
    pl  = pDv(DV_MIXT_PRES,c1)
    
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

    tl = MixtPerf_T_DPR(rl,pl,gc)

!#ifdef PLAG
    ! Subbu - vFracE Correction
!    IF ( global%plagUsed .AND. (pRegion%plag%nPcls .GE. 1) ) THEN
!      tl = MixtPerf_T_DPR( rl/(1.0_RFREAL - pRegion%plag%vFracE(1,c1)), pl,gc)
!    ELSE 
!      tl = MixtPerf_T_DPR(rl,pl,gc)
!    END IF ! global%plagUsed
    ! Subbu - vFracE Correction
!#endif

    al = MixtPerf_C_GRT(g,gc,tl)
    Hl = MixtPerf_Ho_CpTUVW(cp,tl,ul,vl,wl)
    
! ------------------------------------------------------------------------------
!   Right state
! ------------------------------------------------------------------------------

    rr  = pCv(CV_MIXT_DENS,c2)

#ifdef PLAG
    IF ( global%plagUsed .AND. (pRegion%plag%nPcls .GE. 1) ) THEN
     rr  = pCv(CV_MIXT_DENS,c2)/(1.0_RFREAL - pRegion%plag%vFracE(1,c2))
    END IF
#endif

!    irr = 1.0_RFREAL/rr
    irr = 1.0_RFREAL/pCv(CV_MIXT_DENS,c2)

    dx  = xc - pGrid%cofg(XCOORD,c2)
    dy  = yc - pGrid%cofg(YCOORD,c2)
    dz  = zc - pGrid%cofg(ZCOORD,c2)

#ifdef PLAG
    rcR = (dx*dx + dy*dy + dz*dz)**0.5_RFREAL
#endif

    ur  = pCv(CV_MIXT_XMOM,c2)*irr
    vr  = pCv(CV_MIXT_YMOM,c2)*irr
    wr  = pCv(CV_MIXT_ZMOM,c2)*irr
    pr  = pDv(DV_MIXT_PRES,c2)
        
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

    tr = MixtPerf_T_DPR(rr,pr,gc)

!#ifdef PLAG
    ! Subbu - vFracE Correction
!    IF ( global%plagUsed .AND. (pRegion%plag%nPcls .GE. 1) ) THEN
!      tr = MixtPerf_T_DPR( rr/(1.0_RFREAL - pRegion%plag%vFracE(1,c2)),pr,gc)
!    ELSE 
!      tr = MixtPerf_T_DPR(rr,pr,gc)
!    END IF ! global%plagUsed
    ! Subbu - vFracE Correction
!#endif

    ar = MixtPerf_C_GRT(g,gc,tr)
    Hr = MixtPerf_Ho_CpTUVW(cp,tr,ur,vr,wr)
    
! ==============================================================================
!   Rahul - 11/10/15
!   Check validity of volume fraction computed from 2nd order scheme. If the
!   value becomes < 0 or > 1, then switch to first order scheme i.e left/right
!   face value = cell center value. The value at face will be the average of
!   left and right states computed in ..FluxFunction2. Otherwise flux will be
!   computed by ..FluxFunction.
!   
!   Edit: Rahul - 12/14/15
!   The second order computation of phi is abandoned in place of an weighted 
!   average. Weights are inversely proportional to the distance between the 
!   face and the cell centroid.
! ==============================================================================

#ifdef PLAG

!       vFracGl = vFracGl + pGvolFracEg(XCOORD,1,c1)*dx &
!                         + pGvolFracEg(YCOORD,1,c1)*dy &
!                         + pGvolFracEg(ZCOORD,1,c1)*dz
 
!       vFracGr = vFracGr + pGvolFracEg(XCOORD,1,c2)*dx &
!                         + pGvolFracEg(YCOORD,1,c2)*dy &
!                         + pGvolFracEg(ZCOORD,1,c2)*dz
!    IF ((vFracGl .LT. 0 ).OR. (vFracGr .LT. 0) .OR. &
!        (vFracGl .GT. 1 ).OR. (vFracGr .GT. 1)) THEN
!     write(STDOUT,*) 'face #', ifg
!     write(STDOUT,*) 'c1, c2', c1,c2
!     write(STDOUT,*) 'phiL,phiR',vFracGl,vFracGr
!    END IF
  
    IF ( global%plagUsed .AND. (pRegion%plag%nPcls .GE. 1) ) THEN
      rcL = 1.0_RFREAL/rcL
      rcR = 1.0_RFREAL/rcR
    
      wtL = rcL/(rcL+rcR)
      wtR = rcR/(rcL+rcR)

      vFracGl = 1.0_RFREAL-pPlag%vFracE(1,c1)
      vFracGr = 1.0_RFREAL-pPlag%vFracE(1,c2)
!      vFracF  = 0.5_RFREAL*(vFracGl + vFacGr)
      vFracF  = wtL*vFracGl + wtR*vFracGr
    END IF

!    IF ((vFracF .LT. 0.0_RFREAL) .OR. (vFracF .GT. 1.0_RFREAL)) THEN
!     write(STDOUT,*) 'face #', ifg
!     write(STDOUT,*) 'c1, c2', c1,c2
!     write(STDOUT,*) 'phiF',vFracF
!     CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__,'Invalid value for volume &
!                                                                      fraction')
!    END IF
#endif

! ==============================================================================
!   Compute fluxes using vFracF
! ==============================================================================
    CALL INRT_ComputePsi2(nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,Hl,al,rr,&
                                       ur,vr,wr,pr,Hr,ar,vFracF,flx,vf)

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
  CALL FPROFILER_ENDS("RFLU::AUSMPlusUp_ComputeFlux2_TCP")
#endif

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AUSMPlusUp_ComputeFlux2_TCP






  
! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE INRT_DiscreteDelta


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_DiscreteDelta.F90,v $
! Revision 1.3  2016/02/05 16:42:23  rahul
! 1. Fixed a bug in conversion of phi*rho to primitive variable for multi-
! phase simulations.
! 2. Fixed a minor bug in JWL flux routine.
!
! Revision 1.2  2016/02/04 22:21:57  fred
! Adding JWL EOS iterative capabilities for cylindrical detonation problem
!
! Revision 1.1  2015/12/18 22:40:11  rahul
! AUSM+up flux computation routines.
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
  
  

