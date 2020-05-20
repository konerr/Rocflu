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
!******************************************************************************
!
! Purpose: compute interaction source for drag forces on Lagrangian particles.
!
! Description: none.
!
! Input: region  = current region.
!
! Output: region%levels(iLev)%plag%inrtSources
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_RFLU_CalcForceInterface.F90,v 1.2 2016/05/16 15:19:38 rahul Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_CalcForceInterface( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartLag,    ONLY : t_plag
  USE ModError
  USE ModParameters
  USE INRT_ModParameters

#ifdef PLAG
  USE PLAG_ModParameters
#endif

  USE RFLU_ModDifferentiationCells, ONLY: RFLU_ComputeGradCellsGGScalar, &
                                          RFLU_ComputeGradCellsGGVector

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region

! ... loop variables
  INTEGER :: iPcls

! ... local variables
  CHARACTER(CHRLEN)  :: RCSIdentString

  INTEGER :: dragUnsteady, iCont, nCont, nPcls
  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass

  REAL(RFREAL) :: CdTotal,diamL,factor,gamma,machL,massL,mixtVolR,pi,psiL, &
                  relVelMagL,reyL,tauLR
  REAL(RFREAL),          DIMENSION(3)   :: relVel, accelL
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv, pDv, pTv
  REAL(RFREAL), POINTER, DIMENSION(:,:,:) :: dudtMixt,dudtPlag

  TYPE(t_plag)  , POINTER :: pPlag
  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! New variables for unsteady force ----------
  INTEGER :: icg,iT,nUnsteadyData
  REAL(RFREAL), DIMENSION(3) :: forceIU,forcePG,forceTotal,forceVU
  REAL(RFREAL) :: A,B,CamEff,dt,fH,kernelVU,mf,mu,nu,refArea,rhoMixt, &
                  speedSound,time,vFrac,vFracCorr,volL,volMixt
  ! Subbu - Vars cyldet case
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGc
  REAL(RFREAL) :: ConvFlux,tBeg,tEnd
  INTEGER :: Coord
  ! Subbu - End Vars cyldet case 

! New variables for augmenting rhs ----------
  REAL(RFREAL) :: coeffIU,contFac,energydotg,energydotp
  REAL(RFREAL) :: drudtMixt,drvdtMixt,drwdtMixt,drudtPlag,drvdtPlag,drwdtPlag 
  REAL(RFREAL), DIMENSION(3) :: ug,up

! Rahul - New variables for Delpforce
!  INTEGER :: icg,iT,nUnsteadyData
!  REAL(RFREAL), POINTER, DIMENSION(:)   :: pInt
  REAL(RFREAL), DIMENSION(3) :: forceDelP
  REAL(RFREAL) :: u_g,v_g,w_g,Qg,u_p,v_p,w_p,Qp,p,pint,delp,rhog,CpUp,eps
! Rahul - end

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_CalcForceInterface.F90,v $'

  global => region%global
  pRegion => region

  CALL RegisterFunction( global,'PLAG_RFLU_CalcForceInterface',__FILE__ )

#ifdef PLAG
! Check if there are any particles

  nPcls = 0
  IF (global%plagUsed) nPcls = region%plag%nPcls

  IF (nPcls < 1) GO TO 999

! Get dimensions --------------------------------------------------------------

  pi = global%pi
  dt = global%dtMin

  nCont        = region%plagInput%nCont

    pPlag     => region%plag

    pCv       => pPlag%cv
    pDv       => pPlag%dv
    pTv       => pPlag%tv

!    pInt      => pPlag%pInt

    dudtMixt  => pPlag%dudtMixt
    dudtPlag  => pPlag%dudtPlag

    pCvPlagMass => pPlag%cvPlagMass

    nUnsteadyData = region%plagInput%nUnsteadyData
IF(1==1) THEN
! =============================================================================
! Rahul: Hyperbolicity pressure correction to paricle phase 
! =============================================================================

    DO iPcls = 1,nPcls
      icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcls)
      diamL   = pDv(DV_PLAG_DIAM,iPcls)
      volL    = (pi*(diamL**3.0_RFREAL)/6.0_RFREAL)
      volMixt = region%grid%vol(icg) 

      vFrac = pPlag%vFracL(1,iPcls)

      massL = SUM( pCv(pCvPlagMass(:),iPcls) )
      CpUp  = 2.0_RFREAL
      eps   = 0.01_RFREAL

      u_g    = pDv(DV_PLAG_UVELMIXT,iPcls)
      v_g    = pDv(DV_PLAG_VVELMIXT,iPcls)
      w_g    = pDv(DV_PLAG_WVELMIXT,iPcls)
      Qg    = (u_g*u_g + v_g*v_g + w_g*w_g)**0.5_RFREAL

      u_p    = pDv(DV_PLAG_UVEL,iPcls)
      v_p    = pDv(DV_PLAG_VVEL,iPcls)
      w_p    = pDv(DV_PLAG_WVEL,iPcls)
      Qp    = (u_p*u_p + v_p*v_p + w_p*w_p)**0.5_RFREAL

      p = pDv(DV_PLAG_PRESMIXT,iPcls)
      rhog  = pDv(DV_PLAG_DENSMIXT,iPcls)/(1.0_RFREAL-vFrac)
      delp  = CpUp*vFrac*rhog*(Qg-Qp)**2.0_RFREAL
      delp  = MIN(delp,eps*p)

!      pint = p - delp

      forceDelp(XCOORD) = volL*delp*pPlag%gradVFracL(XCOORD,1,iPcls)
      forceDelp(YCOORD) = volL*delp*pPlag%gradVFracL(YCOORD,1,iPcls)
      forceDelp(ZCOORD) = volL*delp*pPlag%gradVFracL(ZCOORD,1,iPcls)

!      contFac = pPlag%arv(ARV_PLAG_SPLOAD,iPcls)

! Rahul - Consistent energy and momentum formulation for AUSM+up 
! NOTE: forceIU being inviscid should not contribute to fluid phase internal energy.
      IF ( region%dummyStep .EQV. .FALSE. ) THEN
!        up = region%plag%dv(DV_PLAG_UVEL:DV_PLAG_WVEL,iPcls)
        energydotp = DOT_PRODUCT(forceDelp,up)
     
! ----- Augment Particle Sources

         region%plag%rhs(CV_PLAG_XMOM:CV_PLAG_ZMOM,iPcls) &
                         = region%plag%rhs(CV_PLAG_XMOM:CV_PLAG_ZMOM,iPcls) &
                           - forceDelp(XCOORD:ZCOORD)
         region%plag%rhs(CV_PLAG_ENER,iPcls) &
                         = region%plag%rhs(CV_PLAG_ENER,iPcls) &
                           - energydotp
      END IF ! region%dummyStep
    END DO ! iPcls
END IF ! 1==2

IF (1==2) THEN
 IF (region%irkStep == 1 .AND. region%iRegionGlobal == 1 ) THEN
  IF(global%currentTime .EQ. 0.00000E+00)THEN
    OPEN(1113,file='forcesPart.dat',form='formatted',status='unknown')
  ELSE
    OPEN(1113,file='forcesPart.dat',form='formatted',status='old' &
        ,position='append')
  END IF !currentTime
!  DO iPcls = 1,4
!   WRITE(1113,'(6(E23.16,1X))')global%currentTime,&
!                               pPlag%inrtSources(1,iPcls),& 
!                               pPlag%inrtPG(1,iPcls),&
!                               pPlag%inrtIU(1,iPcls),&
!                               pPlag%inrtVU(1,iPcls),&
!                               pPlag%inrtTOT(1,iPcls)
!  END DO
   CLOSE(1113)
 END IF !irkStep
END IF

999  CONTINUE
#endif
  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLU_CalcForceInterface

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_CalcForceInterface.F90,v $
! Revision 1.2  2016/05/16 15:19:38  rahul
! Minor bug fix. Commented out WRITE statement that is preventing compilation.
!
! Revision 1.1  2016/05/06 00:06:39  rahul
! Computes force and work done on the particle due to interface pressure to
! preserve hyperbolicity of governing equations. This is a consequence of
! E-L AUSM+up scheme.
!
! Revision 1.3  2016/02/08 22:26:29  rahul
! Subtracted work done contribution from F_pg to fluid phase energy equation.
!
! Revision 1.2  2015/12/18 23:34:41  rahul
! Suppressed the computation of pressure gradient force. This is a
! consequence of governing equations' formulation of multiphase AUSM+up
! scheme.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
!******************************************************************************

