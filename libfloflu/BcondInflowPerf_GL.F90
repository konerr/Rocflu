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
! Purpose: Set inflow boundary condition for one cell for gas-liquid model.
!
! Input: 
!   bcOptType  boundary treatment: subsonic, supersonic, or mixed  
!   Bp         compressibility at const. pressure
!   Bt         compressibility at const. temperature
!   cvg        specific heat of gas at constant volume (boundary cell)
!   cvl        specific heat of liquid at constant volume (boundary cell)
!   cvv        specific heat of vapor at constant volume (boundary cell)
!   po         reference pressure
!   ro         reference density
!   to         reference Temperature
!   rl         given density
!   rel        given energy
!   ru/v/wl    given velocity components
!   nx/y/z     components of normalized face vector (outward facing)
!   Rg         gas constant
!   Rv         vapor constant
!   rgpgl      given density of gas * volume fraction of gas 
!   rvpvl      given density of vapor * volume fraction of vapor 
!   pl         given pressure
!   press      given total pressure 
!   temp       given total temperature
!
! Output:
!   rr        density at boundary
!   ru/v/wr   density * velocity components at boundary
!   rer       density * total internal energy at boundary
!   rgpgr     density of gas * volume fraction of gas at boundary
!   rvpvr     density of vapor * volume fraction of vapor at boundary
!   pr        pressure at boundary
!
! Notes: None.
!
! *****************************************************************************
!
! $Id: 
!
! Copyright: (c) 2006 by the University of Illinois
!
! *****************************************************************************

SUBROUTINE BcondInflowPerf_GL(bcOptType,ro,po,to,Bp,Bt,cvl,cvv,cvg,Rg,Rv,ur,&  
                              vr,wr,vfgr,vfvr,vflr,temp,press,nx,ny,nz,rl, & 
                              rul,rvl,rwl,rel,rgpgl,rvpvl,pl,rr,rur,rvr, &
                              rwr,rer,rgpgr,rvpvr,pr)
  
  USE ModDataTypes
  USE ModGlobal , ONLY : t_global
  USE ModParameters
  USE ModTools, ONLY : MakeNonZero

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
  
  INTEGER, INTENT(IN) :: bcOptType
  REAL(RFREAL), INTENT(IN) :: Bp,Bt,cvg,cvl,cvv,nx,ny,nz,pl,po,press,rel, & 
                              Rg,rgpgl,rl,ro,rul,Rv,rvl,rvpvl,rwl,temp,to, &
                              ur,vfgr,vflr,vfvr,vr,wr
  REAL(RFREAL), INTENT(OUT):: pr,rer,rgpgr,rr,rur,rvr,rvpvr,rwr
  TYPE(t_global), POINTER :: global
   
! ==============================================================================
! Locals
! ==============================================================================
  
  REAL(RFREAL) :: Bg2i,Bg2l,Bl2i,Bl2l,Bv2i,Bv2l,Cg2i,Cg2l,Cl2i,Cl2l,cmi,cmi2, &
                  cml,Cv2i,Cv2l,cvmi,cvml,el,rgi,rhogl,rhogr,rholl,rholr, & 
                  rhovl,rhovr,rli,rlpll,rtcvtr,rti,rvi,tl,tr,ul,vfgi,vfgl, &
                  vfli,vfll,vfvi,vfvl,vl,vml2,vmr2,wl 

! ******************************************************************************
! Start
! ******************************************************************************

  rli = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,press,po,temp,to)
  rvi = MixtPerf_D_PRT(press,Rv,temp)
  rgi = MixtPerf_D_PRT(press,Rg,temp)

  vfli = vflr 
  vfvi = vfvr 
  vfgi = vfgr 

  Cl2i = MixtLiq_C2_Bp(Bp)
  Cv2i = MixtPerf_C2_GRT(1.0_RFREAL,Rv,temp)
  Cg2i = MixtPerf_C2_GRT(1.0_RFREAL,Rg,temp)

  Bl2i = -Bt/Bp
  Bv2i = rvi*Rv
  Bg2i = rgi*Rg

  rti  = rli*vfli + rvi*vfvi + rgi*vfgi
  cvmi = (rli*vfli*cvl + rvi*vfvi*cvv + rgi*vfgi*cvg)/rti
  cmi  = MixtGasLiq_C(cvmi,rti,press,rli,rvi,rgi,vfli,vfvi,vfgi,Cl2i, &
                      Cv2i,Cg2i,Bl2i,Bv2i,Bg2i)
  cmi2 = cmi*cmi
  vmr2 = ur*ur + vr*vr + wr*wr
 
! ==============================================================================
! Supersonic inflow
! ==============================================================================
           
  IF ( vmr2 > cmi2 )  THEN  
    rr    = rti
    rur   = rr*ur
    rvr   = rr*vr
    rwr   = rr*wr
    rer   = rti*cvmi*temp + 0.5_RFREAL*rr*vmr2
    rgpgr = rgi*vfgi
    rvpvr = rvi*vfvi
    pr    = press

! ==============================================================================
! Subsonic outflow
! ==============================================================================
     
  ELSE   
    ul = rul/rl
    vl = rvl/rl
    wl = rwl/rl

    vml2  = ul*ul + vl*vl + wl*wl
    rlpll = rl - rvpvl - rgpgl
    cvml  = (rlpll*cvl + rvpvl*cvv + rgpgl*cvg)/rl          
    el    = rel/rl
    tl    = MixtPerf_T_CvEoVm2(cvml,el,vml2)

    rholl = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,pl,po,tl,to)
    rhovl = MixtPerf_D_PRT(pl,Rv,tl)
    rhogl = MixtPerf_D_PRT(pl,Rg,tl)

    vfll = rlpll/rholl
    vfvl = rvpvl/rhovl
    vfgl = rgpgl/rhogl

    Cl2l = MixtLiq_C2_Bp(Bp)
    Cv2l = MixtPerf_C2_GRT(1.0_RFREAL,Rv,tl)
    Cg2l = MixtPerf_C2_GRT(1.0_RFREAL,Rg,tl)

    Bl2l = -Bt/Bp
    Bv2l = rhovl*Rv
    Bg2l = rhogl*Rg
    cml  = MixtGasLiq_C(cvml,rl,pl,rholl,rhovl,rhogl,vfll,vfvl,vfgl,Cl2l, &
                        Cv2l,Cg2l,Bl2l,Bv2l,Bg2l) 
    pr   = pl + (SQRT(vmr2) - SQRT(vml2))*rl*cml

    rholr = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,pr,po,temp,to)
    rhovr = MixtPerf_D_PRT(pr,Rv,temp)
    rhogr = MixtPerf_D_PRT(pr,Rg,temp)

    rr     = rholr*vflr + rhovr*vfvr + rhogr*vfgr
    rtcvtr = rholr*vflr*cvl + rhovr*vfvr*cvv + rhogr*vfgr*cvg
    rer    = rtcvtr*temp + 0.5_RFREAL*rr*vmr2
    rur    = rr*ur
    rvr    = rr*vr
    rwr    = rr*wr
    rgpgr  = rhogr*vfgr
    rvpvr  = rhovr*vfvr
  END IF ! vmr2

! ******************************************************************************
! End
! ******************************************************************************
 
END SUBROUTINE BcondInflowPerf_GL

! *****************************************************************************
!
! RCS Revision history:
!
! $Log: BcondInflowPerf_GL.F90,v $
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
! Revision 1.1  2006/03/26 20:20:42  haselbac
! Initial revision
!
! *****************************************************************************

