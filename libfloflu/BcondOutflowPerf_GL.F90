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
! *****************************************************************************
!
! Purpose: Set outflow boundary condition for one cell for gas-liquid model.
!
! Description: 
!
! Input:
!   bcOpt    boundary treatment: subsonic, supersonic, or mixed
!   betaP    compressibility at const. pressure
!   betaT    compressibility at const. temperature
!   cvg      specific heat of gas at constant volume (boundary cell) 
!   cvl      specific heat of liquid at constant volume (boundary cell)
!   cvv      specific heat of vapor at constant volume (boundary cell)
!   pin      static pressure at boundary cell
!   Po       reference pressure
!   ro       reference density
!   To       reference Temperature
!   pout     given static outlet pressure
!   Rg       gas constant
!   Rv       vapor constant
!   rho      density at boundary cell       
!   rhoe     density * total energy at boundary cell
!   rhou/v/w density * velocity components at boundary cell   
!   rhogpg   density of gas * volume fraction of gas at boundary cell
!   rhovpv   density of vapor * volume fraction of vapor at boundary cell
!   sx/y/zn  components of ortho-normalized face vector (outward facing)
!
! Output:
!   rhob      density at boundary
!   rhou/v/wb density * velocity components at boundary
!   rhoeb     density * total energy at boundary
!   rhogpgb   density of gas * volume fraction of gas at boundary
!   rhovpvb   density of vapor * volume fraction of vapor at boundary
!         
! Notes: None. 
!
! *****************************************************************************
!
! $Id: BcondOutflowPerf_GL.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $ 
!
! Copyright: (c) 2006 by the University of Illinois
!
! *****************************************************************************

SUBROUTINE BcondOutflowPerf_GL(bcOpt,ro,Po,To,betaP,betaT,cvl,cvv,cvg,Rg,Rv, &
                                pout,sxn,syn,szn,rho,rhou,rhov,rhow,rhoe, &
                                rhogpg,rhovpv,pin,rhob,rhoub,rhovb,rhowb, &
                                rhoeb,rhogpgb,rhovpvb)

  USE ModDataTypes
  USE ModParameters
  USE ModInterfaces, ONLY: MixtGasLiq_C, &
                           MixtLiq_C2_Bp, &
                           MixtPerf_C2_GRT, &
                           MixtLiq_D_DoBpPPoBtTTo, &
                           MixtPerf_D_PRT, &
                           MixtPerf_T_CvEoVm2

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================
  
  INTEGER :: bcOpt
  REAL(RFREAL), INTENT(IN) :: betaP,betaT,cvg,cvl,cvv,pin,Po,pout,Rg,rho,rhoe, &
                              rhogpg,rhou,rhov,rhovpv,rhow,ro,Rv,sxn,syn,szn,To
  REAL(RFREAL), INTENT(OUT) :: rhob,rhoeb,rhogpgb,rhoub,rhovb,rhovpvb,rhowb 

! ==============================================================================
! Locals
! ==============================================================================  

  REAL(RFREAL) :: Bg2,Bl2,Bv2,Cg2,Cl2,cm,ct2,Cv2,cvm,deltp,e,fg,fv,ic2p,ic2pb, &
                  mach,rhg,rhgb,rhl,rhlb,rholpl,rhv,rhvb,rrhoc,t,tb,u,ub,v,vb, &
                  Vel2,vfg,vfgb,vfl,vflb,vfv,vfvb,vnd,w,wb                  

! ******************************************************************************
! Start
! ******************************************************************************

  u    = rhou/rho
  v    = rhov/rho
  w    = rhow/rho
  Vel2 = (u*u + v*v + w*w)

  rholpl = rho - rhovpv - rhogpg 
  cvm    = (rholpl*cvl + rhovpv*cvv + rhogpg*cvg)/rho
  e      = rhoe/rho
  t      = MixtPerf_T_CvEoVm2(cvm,e,Vel2)

  Cl2 = MixtLiq_C2_Bp(betaP)
  Cv2 = MixtPerf_C2_GRT(1.0_RFREAL,Rv,t)
  Cg2 = MixtPerf_C2_GRT(1.0_RFREAL,Rg,t)

  rhl = MixtLiq_D_DoBpPPoBtTTo(ro,betaP,betaT,pin,Po,t,To)
  rhv = MixtPerf_D_PRT(pin,Rv,t)
  rhg = MixtPerf_D_PRT(pin,Rg,t)

  Bl2 = -betaT/betaP
  Bv2 = rhv*Rv
  Bg2 = rhg*Rg

  vfl = rholpl/rhl
  vfv = rhovpv/rhv
  vfg = rhogpg/rhg
  cm  = MixtGasLiq_C(cvm,rho,pin,rhl,rhv,rhg,vfl,vfv,vfg,Cl2,Cv2,Cg2, &
                     Bl2,Bv2,Bg2)
  mach = SQRT(Vel2)/cm

! *****************************************************************************
! Subsonic outflow
! *****************************************************************************

  IF ( mach < 1.0_RFREAL ) THEN
    rrhoc = 1.0_RFREAL/(rho*cm)
    deltp = pin - pout  
    ub    = u + sxn*deltp*rrhoc
    vb    = v + syn*deltp*rrhoc
    wb    = w + szn*deltp*rrhoc

    tb = t - (deltp*pin)/(rho*rho*cvm*cm*cm)
    fg = (vfg*(rhg*Cg2 - rho*cm*cm +(Bg2*pin)/(rho*cvm)))/(rhg*Cg2*rho*cm*cm)
    fv = (vfv*(rhv*Cv2 - rho*cm*cm +(Bv2*pin)/(rho*cvm)))/(rhv*Cv2*rho*cm*cm)

    vfgb = vfg -  deltp*fg
    vfvb = vfv -  deltp*fv
    vflb = 1.0_RFREAL - vfgb - vfvb

    ic2p  = vfl/Cl2 + vfg/Cg2 + vfv/Cv2
    ic2pb = (Bl2*vfl)/Cl2 + (Bg2*vfg)/Cg2 + (Bv2*vfv)/Cv2

    rhob = rho - (rhv-rhl)*(vfl-vflb) - (rhg-rhl)*(vfg-vfgb) &
            - deltp*ic2p + ic2pb*(t-tb)
    rhlb = rhl - deltp/Cl2 + (Bl2*(t-tb))/Cl2
    rhvb = rhv - deltp/Cv2 + (Bv2*(t-tb))/Cv2
    rhgb = rhg - deltp/Cg2 + (Bg2*(t-tb))/Cg2 

! =============================================================================
!   Special treatment to prevent "deltp" from changing the sign of velocity
!   components. This may happen for very small u, v, w.
! =============================================================================

    vnd = ub*sxn + vb*syn + wb*szn
    
    IF ( vnd < 0.0_RFREAL ) THEN ! inflow at outflow boundary
      ub = SIGN(1.0_RFREAL,u)*MAX(ABS(ub),ABS(u))
      vb = SIGN(1.0_RFREAL,v)*MAX(ABS(vb),ABS(v))
      wb = SIGN(1.0_RFREAL,w)*MAX(ABS(wb),ABS(w))
    END IF ! vnd

    rhoub = rhob*ub
    rhovb = rhob*vb
    rhowb = rhob*wb
    rhoeb = (rhlb*vflb*cvl +rhgb*vfgb*cvg +rhvb*vfvb*cvv)*tb + &
              0.5_RFREAL*rhob*(ub*ub + vb*vb + wb*wb)
    rhogpgb = rhgb*vfgb
    rhovpvb = rhvb*vfvb

! *****************************************************************************
! Supersonic flow
! *****************************************************************************

  ELSE
    rhob    = rho
    rhoub   = rhou
    rhovb   = rhov
    rhowb   = rhow
    rhoeb   = rhoe
    rhogpgb = rhogpg
    rhovpvb = rhovpv
  END IF ! mach

! *****************************************************************************
! End
! *****************************************************************************

END SUBROUTINE BcondOutflowPerf_GL

! *****************************************************************************
!
! RCS Revision history:
!
! $Log: BcondOutflowPerf_GL.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:24  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2006/03/26 20:20:45  haselbac
! Initial revision
! 
! *****************************************************************************

