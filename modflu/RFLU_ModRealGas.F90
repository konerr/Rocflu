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
! Purpose: Collect relations for real gas.
!
! Description: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_ModRealGas.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

! Copied a sample function from perfectgas routines
FUNCTION MixtPerf_P_DRT(D,R,T)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: D,R,T
  REAL(RFREAL) :: MixtPerf_P_DRT
   
  MixtPerf_P_DRT = D*R*T

END FUNCTION MixtPerf_P_DRT

! -----------------------------------------------------------------------------

! Coding not completed yet.....
! ST-JWL----------------------
           tgas1, &
           tgas3, &
           tgas4, &
           tgas5, &
           ugas2, &
           ugas3, &
           ugas4, &
           JWL1, &
           JWL2, &
 
    CALL tgas5(p,r,e)
    CALL tgas1(e,r,dummyReal,a,T,3) 


  SUBROUTINE tgas5(p,r,e)

   DOUBLE PRECISION, INTENT(IN) :: p, r
   DOUBLE PRECISION, INTENT(OUT):: e

   INTEGER :: iter, iterMax
   DOUBLE PRECISION :: h,dummyReal
   DOUBLE PRECISION :: e0,e01,p0,dfde0,f0,a0,T0,p1,e1,f1,dfde1,s, a1,T1

!ST-------temporary-use perfect gas law
   go to 88
 
   CALL tgas4(p,r,h,e0)
   CALL tgas1(e0,r,p0,a0,T0,3)
   CALL tgas4(p0,r,h,e01)
   dfde0 = (p-p0)/(e0-e01)
   f0 = p0-p

   iterMax = 100
   DO iter = 1,iterMax
     s  = - f0/dfde0
     e1 = e0 + s 
     CALL tgas1(e1,r,p1,a1,T1,3)
     f1 = p1 - p
     dfde1 = (f1-f0)/s

     IF ( ABS(s) < 1.0D-6) EXIT
     e0 = e1 
     f0 = f1
     dfde0 = dfde1
   END DO ! iter
   e = e1

88 continue
! TEMPORARY: Manoj, changing gamma to 1.3988
!   e=p/r/(1.4D0-1.0D0)
   e=p/r/(1.3988D0-1.0D0)
! END TEMPORARY
   
  END SUBROUTINE tgas5 

 

      subroutine tgas1(e,r,p,a,t,mflag)
!
!     Last update:   Mar 17 0 1993
!     Last edited:   Mar 17 0 1993
!***********************************************************************
!***********                                                 ***********
!***********    pressure, temperature, and speed of sound    ***********
!***********            returned as functions of             ***********
!***********          internal energy and density            ***********
!***********                                                 ***********
!***********************************************************************
!                  new routines by srinivasan et.al.
!     ************************************************************
!
!     inputs:
!
!     e =  internal energy, j/kg ;or; (m/sec)**2 
!     r =  density, kg/m**3
!
!     outputs:
!
!     p =  pressure, pa. ;or; in newtons/m**2
!     a =  speed of sound, m/sec ;
!     t =  temperature, kelvin .
!
!     if mflag = 0, return p ;
!     if mflag = 1, return p and a ;
!     if mflag = 2, return p and t ;
!     if mflag = 3, return p, a and t .
!
!     ************************************************************
!
!***********************************************************************
!     new update by Stanley: Apr 22, 2010
!     Purpose: define variables used in the subroutines, and change data
!              type to double precision 
!     begin update

      implicit none

      integer :: mflag
      double precision :: e, r
      double precision, intent (out):: p, a, t

      integer :: lflag, kflag, iflag, jflag
      double precision :: gammr, gamme, gamm
      double precision :: e0,r0,p0,t0,gascon
      double precision :: rratio, eratio, x, y, z, z1, tnon
      double precision :: rsave, ym, yhigh, ylow, phigh, plow, ahigh, alow, deno, asq
      double precision :: gas1 ,gas2 ,gas3 ,gas4 ,gas5 ,gas6 ,gas7 ,gas8 ,gas9
      double precision :: gas1r,gas2r,gas3r,gas4r,gas5r,gas6r,gas7r,gas8r,gas9r
      double precision :: gas1e,gas2e,gas3e,gas4e,gas5e,gas6e,gas7e,gas8e,gas9e

!     end update
!***********************************************************************
      common /vinkur/ gammr,gamme

      data e0,r0,p0,t0,gascon/78408.4D00,1.292D00,1.0133D05,            &
     &273.15D00,287.06D00/

!ST----------------------
      go to 77
!ST----------------------

      rratio=r/r0
      eratio=e/e0
      y=log10(rratio)
      z=log10(eratio)
      lflag=0
      kflag=0
      if (mflag.gt.1) lflag=1
      if ((mflag.eq.1).or.(mflag.eq.3)) kflag=1
      if(abs(y+4.5D00).lt.2.5D-02) go to 20
      if(abs(y+0.5D00).lt.5.0D-03) go to 50
      iflag=-1
      go to 90
10    if(lflag.eq.1) go to 300
      return
20    iflag=0
      rsave=r
      ym=y
      y=-4.5D00+2.5D-02
      yhigh=y
      r=(10.D0**y)*r0
      jflag=-1
      go to 90
30    phigh=p
      ahigh=a
      y=-4.5D00-2.5D-02
      ylow=y
      r=(10.D0**y)*r0
      jflag=0
      go to 90
40    plow=p
      alow=a
      go to 80
50    iflag=1
      rsave=r
      ym=y
      y=-0.5D00+0.5e-02
      yhigh=y
      r=(10.D0**y)*r0
      jflag=-1
      go to 90
60    phigh=p
      ahigh=a
      y=-0.5D00-0.5D-02
      ylow=y
      r=(10.D0**y)*r0
      jflag=0
      go to 90
70    plow=p
      alow=a
80    p=plow+(phigh-plow)/(yhigh-ylow)*(ym-ylow)
      a=alow+(ahigh-alow)/(yhigh-ylow)*(ym-ylow)
      r=rsave
      if(lflag.eq.1) go to 300
      return
90    continue
      if(y.gt.-0.5D00) go to 200
      if(y.gt.-4.5D00) go to 150
      if(z.gt.0.65D00) go to 100
      gamm=1.3965D00
      go to 250
100   if(z.gt.1.5D00) go to 110
      gas1=1.52792D00-1.26953D-02*y
      gas2=(-6.13514D-01-5.08262D-02*y)*z
      gas3=(-5.49384D-03+4.75120D-05*z-3.18468D-04*y)*y*y
      gas4=(6.31835D-01+3.34012D-02*y-2.19921D-01*z)*z*z
      gas5=-4.96286D01-1.17932D+01*y
      gas6=(6.91028D01+4.40405D+01*y)*z
      gas7=(5.09249D00-1.40326D00*z+2.08988D-01*y)*y*y
      gas8=(1.37308D01-1.78726D01*y-1.86943D01*z)*z*z
      gas9=exp(24.60452D00-2.D00*y-2.093022D01*z)
      deno=1.-gas9
      gamm=gas1+gas2+gas3+gas4+(gas5+gas6+gas7+gas8)/deno
      if(kflag.eq.0) go to 260
      gas1r=-1.26953D-02
      gas2r=-5.08262D-02*z
      gas3r=(-1.098768D-02-9.50240D-05*z-9.554040D-04*y)*y
      gas4r=3.34012D-02*z*z
      gas5r=-1.17932D01
      gas6r=4.40405D01*z
      gas7r=(1.018498D01-2.80652D00*z+6.269641D-01*y)*y
      gas8r=-1.78726D01*z*z
      gas9r=-2.0D00
      gas2e=gas2/z
      gas3e=4.75120D-05*y*y
      gas4e=(1.26367D00+6.68024D-02*y-6.59763D-01*z)*z
      gas6e=gas6/z
      gas7e=-1.40326D00*y*y
      gas8e=(2.74616D01-3.57452D01*y-5.60829D01*z)*z
      gas9e=-2.093022D01
      gammr=gas1r+gas2r+gas3r+gas4r+(gas5r+gas6r+gas7r+gas8r)/deno      &
     &+(gas5+gas6+gas7+gas8)*gas9r*gas9/(deno**2)
      gamme=gas2e+gas3e+gas4e+(gas6e+gas7e+gas8e)/deno                  &
     &+(gas5+gas6+gas7+gas8)*gas9e*gas9/(deno**2)
      go to 260
110   if(z.gt.2.2D00) go to 120
      gas1=-1.70333D01-5.08545D-01*y
      gas2=(2.46299D01+4.45617D-01*y)*z
      gas3=(-8.95298D-03+2.29618D-03*z-2.89186D-04*y)*y*y
      gas4=(-1.10204D01-9.89727D-02*y+1.62903D00*z)*z*z
      gas5=1.86797D01+5.19662D-01*y
      gas6=(-2.41338D01-4.34837D-01*y)*z
      gas7=(9.16089D-03-1.52082D-03*z+3.46482D-04*y)*y*y
      gas8=(1.02035D01+9.70762D-02*y-1.39460D00*z)*z*z
      gas9=(-1.42762D02-1.647088D00*y+7.660312D01*z                     &
     &+8.259346D-01*y*z)
      if(kflag.eq.0) go to 240
      gas1r=-5.08545D-01
      gas2r=4.45617D-01*z
      gas3r=(-1.790596D-02+4.59236D-03*z-8.67558D-04*y)*y
      gas4r=-9.89727D-02*z*z
      gas5r=5.19662D-01
      gas6r=-4.34837D-01*z
      gas7r=(1.832178d-02-3.04164d-03*z+1.039446d-03*y)*y
      gas8r=9.70762d-02*z*z
      gas9r=-1.647088d00+8.259346d-01*z
      gas2e=gas2/z
      gas3e=2.29618d-03*y*y
      gas4e=(-2.20408d01-1.979454d-01*y+4.88709d00*z)*z
      gas6e=gas6/z
      gas7e=-1.52082e-03*y*y
      gas8e=(2.0407d01+1.941524d-01*y-4.1838d00*z)*z
      gas9e=7.660312d01+8.259346d-01*y
      go to 240
120   if (z.gt.3.05e00) go to 130
      gas1=2.24374d00+1.03073d-01*y
      gas2=(-5.32238d-01-5.59852d-02*y)*z
      gas3=(3.56484d-03-1.01359d-04*z+1.59127d-04*y)*y*y
      gas4=(-4.80156d-02+1.06794d-02*y+3.66035d-02*z)*z*z
      gas5=-5.70378d00-3.10056d-01*y
      gas6=(5.01094d00+1.80411d-01*y)*z
      gas7=(-9.49361d-03+1.94839d-03*z-2.24908d-04*y)*y*y
      gas8=(-1.40331d00-2.79718d-02*y+1.20278d-01*z)*z*z
      gas9=(1.139755d02-4.985467d00*y-4.223833d01*z                     &
     &+2.009706d00*y*z)
      if(kflag.eq.0) go to 240
      gas1r=1.03073d-01
      gas2r=-5.59852d-02*z
      gas3r=(7.12968d-03-2.0218d-04*z+4.77381d-04*y)*y
      gas4r=1.06794d-02*z*z
      gas5r=-3.10056d-01
      gas6r=1.80411d-01*z
      gas7r=(-1.898722d-02+3.89678d-03*z-6.74724d-04*y)*y
      gas8r=-2.79718d-02*z*z
      gas9r=-4.985467d00+2.009706d00*z
      gas2e=gas2/z
      gas3e=-1.01359d-04*y*y
      gas4e=(-9.60312d-02+2.13588d-02*y+1.098105d-01*z)*z
      gas6e=gas6/z
      gas7e=1.94839d-03*y*y
      gas8e=(-2.80662d00-5.59436d-02*y+3.60834d-01*z)*z
      gas9e=-4.223833d01+2.009706d00*y
      go to 240
130   if(z.gt.3.4d00) go to 140
      gas1=-0.20807d02+0.40197d00*y
      gas2=(0.22591d02-0.25660d00*y)*z
      gas3=(-0.95833d-03+0.23966d-02*z+0.33671d-03*y)*y*y
      gas4=(-0.77174d01+0.4606d-01*y+0.878d00*z)*z*z
      gas5=-0.21737d03-0.46927d01*y
      gas6=(0.18101d03+0.26621d01*y)*z
      gas7=(-0.34759d-01+0.64681d-02*z-0.70391d-03*y)*y*y
      gas8=(-0.50019d02-0.38381d00*y+0.45795d01*z)*z*z
      gas9=(0.4544373d03+0.1250133d02*y-0.1376001d03*z                  &
     &-0.3641774d01*y*z)
      if(kflag.eq.0) go to 240
      gas1r=0.40197d00
      gas2r=-0.25660d00*z
      gas3r=(-1.91666d-03+4.7932d-03*z+1.01013d-03*y)*y
      gas4r=0.4606d-01*z*z
      gas5r=-0.46927d01
      gas6r=0.26621d01*z
      gas7r=(-6.9518d-02+1.29362d-02*z-2.11173d-03*y)*y
      gas8r=-0.38381d00*z*z
      gas9r=0.1250133d02-0.3641774d01*z
      gas2e=gas2/z
      gas3e=0.23966d-02*y*y
      gas4e=(-1.54348d01+9.212d-02*y+2.634d00*z)*z
      gas6e=gas6/z
      gas7e=0.64681d-02*y*y
      gas8e=(-1.00038d02-7.6762d-01*y+1.37385d01*z)*z
      gas9e=-0.1376001d03-0.3641774d01*y
      go to 240
140   continue
!140   if(z.gt.3.69e00) write(6,1000) r,e
      gas1=-5.22951d01-4.00011d-01*y
      gas2=(4.56439d01+2.24484d-01*y)*z
      gas3=(-3.73775d-03+2.43161d-03*z+2.24755d-04*y)*y*y
      gas4=(-1.29756d01-2.79517d-02*y+1.22998d00*z)*z*z
      gamm=gas1+gas2+gas3+gas4
      if(kflag.eq.0) go to 260
      gas1r=-4.00011d-01
      gas2r=2.24484d-01*z
      gas3r=(-7.4755d-03+4.86322d-03*z+6.74265d-04*y)*y
      gas4r=-2.79517d-02*z*z
      gas2e=gas2/z
      gas3e=2.43161d-03*y*y
      gas4e=(-2.59512d01-5.59034d-02*y+3.68994d00*z)*z
      gammr=gas1r+gas2r+gas3r+gas4r
      gamme=gas2e+gas3e+gas4e
      go to 260
150   if(z.gt.0.65e00) go to 160
      gamm=1.398d00
      go to 250
160   if (z.gt.1.5d00) go to 170
      gas1=1.39123d00-4.08321d-03*y
      gas2=(1.42545d-02+1.41769d-02*y)*z
      gas3=(2.57225d-04+6.52912d-04*z+8.46912d-05*y)*y*y
      gas4=(6.2555d-02-7.83637d-03*y-9.78720d-02*z)*z*z
      gas5=5.80955-1.82302d-01*y
      gas6=(-9.62396d00+1.79619d-01*y)*z
      gas7=(-2.30518d-02+1.18720d-02*z-3.35499d-04*y)*y*y
      gas8=(5.27047d00-3.65507d-02*y-9.19897d-01*z)*z*z
      gas9=(-10.0d00*z+14.2d00)
      if(kflag.eq.0) go to 240
      gas1r=-4.08321d-03
      gas2r=1.41769d-02*z
      gas3r=(5.1445d-04+1.305824d-03*z+2.540736d-04*y)*y
      gas4r=-7.83637d-03*z*z
      gas5r=-1.82302d-01
      gas6r=1.79619d-01*z
      gas7r=(-4.61036d-02+2.3744d-02*z-1.006497d-03*y)*y
      gas8r=-3.65507d-02*z*z
      gas9r=0.0d00
      gas2e=gas2/z
      gas3e=6.52912d-04*y*y
      gas4e=(1.2511d-01-1.567274d-02*y-2.93616d-01*z)*z
      gas6e=gas6/z
      gas7e=1.1872d-02*y*y
      gas8e=(1.054094d01-7.31014d-02*y-2.759691d00*z)*z
      gas9e=-10.0d00
      go to 240
170   if (z.gt.2.22d00) go to 180
      gas1=-1.20784d00-2.57909d-01*y
      gas2=(5.02307d00+2.87201d-01*y)*z
      gas3=(-9.95577d-03+5.23524d-03*z-1.45574d-04*y)*y*y
      gas4=(-3.20619d00-7.50405d-02*y+6.51564d-01*z)*z*z
      gas5=-6.62841d00+2.77112d-02*y
      gas6=(7.30762d00-7.68230d-02*y)*z
      gas7=(7.19421d-03-3.62463d-03*z+1.62777d-04*y)*y*y
      gas8=(-2.33161d00+3.04767d-02*y+1.66856d-01*z)*z*z
      gas9=(1.255324d02+2.015335d00*y-6.390747d01*z-                    &
     &6.515225d-01*y*z)
      if(kflag.eq.0) go to 240
      gas1r=-2.57909d-01
      gas2r=2.87201d-01*z
      gas3r=(-1.991154d-02+1.047048d-02*z-4.36722d-04*y)*y
      gas4r=-7.50405d-02*z*z
      gas5r=2.77112d-02
      gas6r=-7.6823d-02*z
      gas7r=(1.438842d-02-7.24926d-03*z+4.88331d-04*y)*y
      gas8r=3.04767d-02*z*z
      gas9r=2.015335d00-6.515225d-01*z
      gas2e=gas2/z
      gas3e=5.23524d-03*y*y
      gas4e=(-6.41238d00-1.50081d-01*y+1.954692d00*z)*z
      gas6e=gas6/z
      gas7e=-3.62463d-03*y*y
      gas8e=(-4.66322d00+6.09534d-02*y+5.00568d-01*z)*z
      gas9e=-6.390747d01-6.515225d-01*y
      go to 240
180   if (z.gt.2.95) go to 190
      gas1=-2.26460d00-7.82263d-02*y
      gas2=(4.90497d00+7.18096d-02*y)*z
      gas3=(-3.06443d-03+1.74209d-03*z+2.84214d-05*y)*y*y
      gas4=(-2.24750d00-1.31641d-02*y+3.33658d-01*z)*z*z
      gas5=-1.47904d01-1.76627d-01*y
      gas6=(1.35036d01+8.77280d-02*y)*z
      gas7=(-2.13327d-03+7.15487d-04*z+7.30928d-05*y)*y*y
      gas8=(-3.95372d00-8.96151d-03*y+3.63229d-01*z)*z*z
      gas9=(1.788542d02+6.317894d00*y-6.756741d01*z-                    &
     &2.460060d00*y*z)
      if(kflag.eq.0) go to 240
      gas1r=-7.82263d-02
      gas2r=7.18096d-02*z
      gas3r=(-6.12886d-03+3.48418d-03*z+8.52642d-05*y)*y
      gas4r=-1.31641d-02*z*z
      gas5r=-1.76627d-01
      gas6r=8.7728d-02*z
      gas7r=(-4.26654d-03+1.430974d-03*z+2.192784d-04*y)*y
      gas8r=-8.96151d-03*z*z
      gas9r=6.317894d00-2.46006d00*z
      gas2e=gas2/z
      gas3e=1.74209d-03*y*y
      gas4e=(-4.495d00-2.63282d-02*y+1.000974d00*z)*z
      gas6e=gas6/z
      gas7e=7.15487d-04*y*y
      gas8e=(-7.90744d00-1.792302d-02*y+1.089687d00*z)*z
      gas9e=-6.756741d01-2.46006d00*y
      go to 240
 190  continue
!190   if(z.gt.3.4e00) write(6,1000) r,e
      gas1=-1.66904d01-2.58318d-01*y
      gas2=(1.78350d01+1.54898d-01*y)*z
      gas3=(-9.71263d-03+3.97740d-03*z+9.04300d-05*y)*y*y
      gas4=(-5.94108d00-2.01335d-02*y+6.60432d-01*z)*z*z
      gas5=8.54690d01+1.17554d01*y
      gas6=(-7.21760d01-7.15723d00*y)*z
      gas7=(-4.16150d-02+1.38147d-02*z+5.45184d-04*y)*y*y
      gas8=(2.01758d01+1.08990d00*y-1.86438d00*z)*z*z
      gas9=(2.883262d02+1.248536d01*y-8.816985d01*z-                    &
     &3.720309d00*y*z)
      if(kflag.eq.0) go to 240
      gas1r=-2.58318d-01
      gas2r=1.54898d-01*z
      gas3r=(-1.942526d-02+7.9548d-03*z+2.7129d-04*y)*y
      gas4r=-2.01335d-02*z*z
      gas5r=1.17554d01
      gas6r=-7.15723d00*z
      gas7r=(-8.323d-02+2.76294d-02*z+1.635552d-03*y)*y
      gas8r=1.0899d00*z*z
      gas9r=1.248536d01-3.720309d00*z
      gas2e=gas2/z
      gas3e=3.9774d-03*y*y
      gas4e=(-1.188216d01-4.0267d-02*y+1.981296d00*z)*z
      gas6e=gas6/z
      gas7e=1.38147d-02*y*y
      gas8e=(4.03516d01+2.1798d00*y-5.59314d00*z)*z
      gas9e=-8.816985d01-3.720309d00*y
      go to 240
200   if(z.gt.0.65d00) go to 210
      gamm=1.3988d00
      go to 250
210   if(z.gt.1.7d00) go to 220
      gas1=1.37062d00+1.29673d-02*y
      gas2=(1.11418d-01-3.26912d-02*y)*z
      gas3=(1.06869d-03-2.00286d-03*z+2.38305d-04*y)*y*y
      gas4=(-1.06133d-01+1.90251d-02*y+3.02210d-03*z)*z*z
      gamm=gas1+gas2+gas3+gas4
      if(kflag.eq.0) go to 260
      gas1r=1.29673d-02
      gas2r=-3.26912d-02*z
      gas3r=(2.13738d-03-4.00572d-03*z+7.14915d-04*y)*y
      gas4r=1.90251d-02*z*z
      gas2e=gas2/z
      gas3e=-2.00286d-03*y*y
      gas4e=(-2.12266d-01+3.80502d-02*y+9.0663d-03*z)*z
      gammr=gas1r+gas2r+gas3r+gas4r
      gamme=gas2e+gas3e+gas4e
      go to 260
220   if(z.gt.2.35) go to 230
      gas1=3.43846d-02-2.33584d-01*y
      gas2=(2.85574d00+2.59787d-01*y)*z
      gas3=(-10.89927d-03+4.23659d-03*z+3.85712d-04*y)*y*y
      gas4=(-1.94785d00-6.73865d-02*y+4.08518d-01*z)*z*z
      gas5=-4.20569d00+1.33139d-01*y
      gas6=(4.51236d00-1.66341d-01*y)*z
      gas7=(1.67787d-03-1.10022d-03*z+3.06676d-04*y)*y*y
      gas8=(-1.35516d00+4.91716d-02*y+7.52509d-02*z)*z*z
      gas9=(1.757042d02-2.163278d00*y-8.833702d01*z+                    &
     &1.897543d00*y*z)
      if(kflag.eq.0) go to 240
      gas1r=-2.33584d-01
      gas2r=2.59787d-01*z
      gas3r=(-2.179854d-02+8.47318d-03*z+1.157136d-03*y)*y
      gas4r=-6.73865d-02*z*z
      gas5r=1.33139d-01
      gas6r=-1.66341d-01*z
      gas7r=(3.35574d-03-2.20044d-03*z+9.20028d-04*y)*y
      gas8r=4.91716d-02*z*z
      gas9r=-2.163278d00+1.897543d00*z
      gas2e=gas2/z
      gas3e=4.23659d-03*y*y
      gas4e=(-3.8957d00-1.34773d-01*y+1.225554d00*z)*z
      gas6e=gas6/z
      gas7e=-1.10022d-03*y*y
      gas8e=(-2.71032d00+9.83432d-02*y+2.257527d-01*z)*z
      gas9e=-8.833702d01+1.897543d00*y
      go to 240
 230  continue
!230   if(z.gt.2.9e00) write(6,1000) r,e
      gas1=-1.70633d00-1.48403d-01*y
      gas2=(4.23104d00+1.37290d-01*y)*z
      gas3=(-9.10934d-03+3.85707d-03*z+2.69026d-04*y)*y*y
      gas4=(-1.97292d00-2.81830d-02*y+2.95882d-01*z)*z*z
      gas5=3.41580d01-1.89972d01*y
      gas6=(-4.0858d01+1.30321d01*y)*z
      gas7=(-8.01272d-01+2.75121d-01*z-1.77969d-04*y)*y*y
      gas8=(1.60826d01-2.23386d00*y-2.08853d00*z)*z*z
      gas9=(2.561323d02+1.737089d02*y-9.058890d01*z-                    &
     &5.838803d01*y*z)
      if(gas9.gt.30.d00) gas9=30.d00
      if(gas9.lt.-30.d00) gas9=-30.d00
      if(kflag.eq.0) go to 240
      gas1r=-1.48403d-01
      gas2r=1.3729d-01*z
      gas3r=(-1.821868d-02+7.71414d-03*z+8.07078d-04*y)*y
      gas4r=-2.8183d-02*z*z
      gas5r=-1.89972d01
      gas6r=1.30321d01*z
      gas7r=(-1.602544d00+5.50242d-01*z-5.33907d-04*y)*y
      gas8r=-2.23386d00*z*z
      gas9r=1.737089d02-5.838803d01*z
      gas2e=gas2/z
      gas3e=3.85707d-03*y*y
      gas4e=(-3.94584d00-5.6366d-02*y+8.87646d-01*z)*z
      gas6e=gas6/z
      gas7e=2.75121d-01*y*y
      gas8e=(3.21652d01-4.46772d00*y-6.26559d00*z)*z
      gas9e=-9.05889d01-5.838803d01*y
240   gas9=exp(gas9)
      gamm=gas1+gas2+gas3+gas4+(gas5+gas6+gas7+gas8)/(1.+gas9)
      if(kflag.eq.0) go to 260
      gammr=gas1r+gas2r+gas3r+gas4r+(gas5r+gas6r+gas7r+gas8r)/          &
     &(1.+gas9)-(gas5+gas6+gas7+gas8)*gas9r*gas9/((1.+gas9)**2)
      gamme=gas2e+gas3e+gas4e+(gas6e+gas7e+gas8e)/(1.+gas9)             &
     &-(gas5+gas6+gas7+gas8)*gas9e*gas9/((1.+gas9)**2)
      go to 260
250   continue
      if(kflag.eq.0) go to 260
      gammr=0.0d00
      gamme=0.0d00
260   p=(gamm-1)*e*r
      if(kflag.eq.0) go to 270
      gammr=gammr/2.302585d00
      gamme=gamme/2.302585d00
      asq=e*((gamm-1.d00)*(gamm+gamme)+gammr)
      a=sqrt(asq)
270   if(iflag) 10,280,290
280   if(jflag) 30,40,10
290   if(jflag) 60,70,10
300   x=log10(p/p0)
      y=log10(r/r0)
      z1=x-y
      if (y.gt.-0.5d00) go to 400
      if(y.gt.-4.5d00) go to 350
      if(z1.gt.0.25d00) go to 310
      t=p/(gascon*r)
      return
310   if(z1.gt.0.95d00) go to 320
      gas1=1.44824d-01+1.36744d-02*y
      gas2=(1.17099d-01-8.22299d-02*y)*z1
      gas3=(-6.75303d-04-1.47314d-03*z1-7.90851d-05*y)*y*y
      gas4=(1.3937d00+6.83066d-02*y-6.65673d-01*z1)*z1*z1
      tnon=gas1+gas2+gas3+gas4
      go to 450
320   if(z1.gt.1.4d00) go to 330
      gas1=-9.325d00-9.32017d-01*y
      gas2=(2.57176d01+1.61292d00*y)*z1
      gas3=(-3.00242d-02+2.62959d-02*z1-2.77651d-04*y)*y*y
      gas4=(-2.1662d01-6.81431d-01*y+6.26962d00*z1)*z1*z1
      gas5=-3.38534d00+1.82594d-01*y
      gas6=(1.84928d-01-7.01109d-01*y)*z1
      gas7=(1.10150d-02-1.60570d-02*z1+1.57701d-05*y)*y*y
      gas8=(5.4702d00+4.11624d-01*y-2.81498d00*z1)*z1*z1
      gas9=exp(-3.887015d01-2.908228d01*y+4.070557d01*z1                &
     &+2.682347d01*y*z1)
      go to 440
330   if( z1.gt.1.95d00) go to 340
      gas1=-1.93082d01-1.54557d00*y
      gas2=(3.69035d01+1.92214d00*y)*z1
      gas3=(-3.59027d-02+2.31827d-02*z1-2.01327d-04*y)*y*y
      gas4=(-2.20440d01-5.80935d-01*y+4.43367d00*z1)*z1*z1
      gas5=-3.83069d00+1.32864d-01*y
      gas6=(-3.91902d00-6.79564d-01*y)*z1
      gas7=(6.06341d-04-8.12997d-03*z1-1.61012d-04*y)*y*y
      gas8=(7.24632d00+3.15461d-01*y-2.17879d00*z1)*z1*z1
      gas9=exp(2.08d01-2.56d01*y+1.0d00*z1+1.80d01*y*z1)
      go to 440
340   gas1=-2.59721d01-1.77419d00*y
      gas2=(3.62495d01+1.55383d00*y)*z1
      gas3=(-4.51359d-02+2.43648d-02*z1+1.2804d-04*y)*y*y
      gas4=(-1.59988d01-3.17807d-01*y+2.40584d00*z1)*z1*z1
      gas5=-1.81433d01+1.54896d-01*y
      gas6=(1.26582d01-3.66275d-01*y)*z1
      gas7=(3.24496d-02-1.66385d-02*z1+3.02177d-04*y)*y*y
      gas8=(-1.41759d00+1.11241d-01*y-3.10983d-01*z1)*z1*z1
      gas9=exp(1.115884d02-6.452606d00*y-5.337863d01*z1                 &
     &+2.026986d00*y*z1)
      go to 440
350   if(z1.gt.0.25d00) go to 360
      t=p/(gascon*r)
      return
360   if(z1.gt.0.95d00) go to 370
      gas1=2.94996d-02+7.24997d-03*y
      gas2=(7.81783d-01-3.27402d-02*y)*z1
      gas3=(3.23357d-04-9.69989d-04*z1-8.93240d-06*y)*y*y
      gas4=(3.95198d-01+2.92926d-02*y-2.12182d-01*z1)*z1*z1
      tnon=gas1+gas2+gas3+gas4
      go to 450
370   if (z1.gt.1.4e00) go to 380
      gas1=-5.53324d00-3.53749d-01*y
      gas2=(1.63638d01+5.87547d-01*y)*z1
      gas3=(-1.16081d-02+7.99571d-03*z1-2.79316d-04*y)*y*y
      gas4=(-1.41239d01-2.35146d-01*y+4.28891d00*z1)*z1*z1
      gas5=9.07979d00+1.01308d00*y
      gas6=(-2.29428d01-1.52122d00*y)*z1
      gas7=(3.78390d-02-2.63115d-02*z1+5.46402d-04*y)*y*y
      gas8=(1.95657d01+5.73839d-01*y-5.63057d00*z1)*z1*z1
      gas9=exp(7.619803d01-1.501155d01*y-6.770845d01*z1                 &
     &+1.273147d01*y*z1)
      go to 440
380   if (z1.gt.2.00d00) go to 390
      gas1=-1.13598d01-1.02049d00*y
      gas2=(2.22793d01+1.24038d00*y)*z1
      gas3=(-3.10771d-02+1.92551d-02*z1-2.69140d-04*y)*y*y
      gas4=(-1.31512d01-3.62875d-01*y+2.64544d00*z1)*z1*z1
      gas5=8.72852d00+1.27564d00*y
      gas6=(-1.79172d01-1.52051d00*y)*z1
      gas7=(4.91264d-02-2.81731d-02*z1+5.23383d-04*y)*y*y
      gas8=(1.16719d01+4.45413d-01*y-2.45584d00*z1)*z1*z1
      gas9=exp(1.84792d02+9.583443d00*y-1.020835d02*z1                  &
     &-4.166727d00*y*z1)
      go to 440
390   gas1=-1.76079d01-1.26579d00*y
      gas2=(2.48544d01+1.09442d00*y)*z1
      gas3=(-3.65534d-02+1.54346d-02*z1-4.59822d-04*y)*y*y
      gas4=(-1.08166d01-2.27803d-01*y+1.60641d00*z1)*z1*z1
      gas5=2.60669d01+2.31791d00*y
      gas6=(-3.22433d01-1.82645d00*y)*z1
      gas7=(4.94621d-02-1.85542d-02*z1+5.04815d-04*y)*y*y
      gas8=(1.33829d01+3.59744d-01*y-1.86517d00*z1)*z1*z1
      gas9=exp(3.093755d02+1.875018d01*y-1.375004d02*z1                 &
     &-8.333418d00*y*z1)
      go to 440
400   if (z1.gt.0.25d00) go to 410
      t=p/(gascon*r)
      return
410   if (z1.gt.0.95d00) go to 420
      gas1=-2.94081d-03+5.73915d-04*y
      gas2=(9.88883d-01-3.71241d-03*y)*z1
      gas3=(1.12387d-04-3.76528d-04*z1+1.76192d-05*y)*y*y
      gas4=(2.86656d-02+4.56059d-03*y-1.99498d-02*z1)*z1*z1
      tnon=gas1+gas2+gas3+gas4
      go to 450
420   if (z1.gt.1.45d00) go to 430
      gas1=1.32396d00+8.52771d-02*y
      gas2=(-3.24257d00-2.00937d-01*y)*z1
      gas3=(5.68146d-03-6.85856d-03*z1+1.98366d-04*y)*y*y
      gas4=(4.53823d00+1.18123d-01*y-1.6246d00*z1)*z1*z1
      gas5=-5.26673d-01-1.58691d-01*y
      gas6=(2.61600d00+3.16356d-01*y)*z1
      gas7=(-1.90755d-02+1.70124d-02*z1-5.58398d-04*y)*y*y
      gas8=(-3.3793d00-1.52212d-01*y+1.30757d00*z1)*z1*z1
      gas9=exp(1.442206d02-2.544727d01*y-1.277055d02*z1                 &
     &+2.236647d01*y*z1)
      go to 440
430   gas1=-1.60643d00-5.07368d-02*y
      gas2=(3.95872d00+3.69383d-02*y)*z1
      gas3=(-1.59378d-03+1.06057d-03*z1+6.53278d-05*y)*y*y
      gas4=(-1.71201d00+9.25124d-03*y+2.71039d-01*z1)*z1*z1
      gas5=1.80476d01+1.62964d00*y
      gas6=(-2.73124d01-1.57430d00*y)*z1
      gas7=(5.85277d-02-2.77313d-02*z1+1.16146d-03*y)*y*y
      gas8=(1.36342d01+3.70714d-01*y-2.23787d00*z1)*z1*z1
      gas9=exp(1.292515d02+1.360552d00*y-7.07482d01*z1                  &
     &+1.360532d00*y*z1)
440   tnon=gas1+gas2+gas3+gas4+(gas5+gas6+gas7+gas8)/(1.d0+gas9)
450   t=(10.e0**tnon)*t0
1000  format(/20x,48hwarning!  outside of validity range of curve fit   &
     &,/,20x,5hrho =,1pd15.8,5x,3he =,1pd15.8,/)

!ST----------------------
77    continue
! DEBUG - ideal gas eos
! TEMPORARY: Manoj, changing gamma=1.3988
!      p=r*e*(1.4D0-1.0D0)
!      a=sqrt(1.4D0*p/r)
!      T=p/r/287.06D0
      p=r*e*(1.3988D0-1.0D0)
      a=sqrt(1.3988D0*p/r)
      T=p/r/287.06D0
! END TEMPORARY
! END DEBUG
!ST----------------------

      return
      end subroutine tgas1



      subroutine tgas3(p,rho,t)
!
!     Last update:   Mar 17 0 1993
!     Last edited:   Mar 17 0 1993
!***********************************************************************
!**********                                                 ************
!**********                  temperature                    ************
!**********           returned as a function of             ************
!**********             pressure and density                ************
!**********                                                 ************
!***********************************************************************
!                  new routine by srinivasan et.al.
!     *********************************************************
!
!     inputs:
!
!     p = pressure, in newtons/m**2.
!     rho =  density, kg/m**3.
!
!     outputs:
!
!     t = temperature, in kelvin.
!
!     *********************************************************
!
!***********************************************************************
!     new update by Stanley: Apr 22, 2010
!     Purpose: define variables used in the subroutines
!     begin update

      implicit none

      double precision :: p, rho
      double precision, intent (out):: t

      integer :: lflag, kflag, iflag, jflag
      double precision :: gammr, gamme, gamm
      double precision :: e0,r0,p0,t0,gascon
      double precision :: rratio, eratio, x, y, z, z1, tnon
      double precision :: rsave, ym, yhigh, ylow, phigh, plow, thigh, tlow
      double precision :: gas1 ,gas2 ,gas3 ,gas4 ,gas5 ,gas6 ,gas7 ,gas8 ,gas9
      double precision :: gas1r,gas2r,gas3r,gas4r,gas5r,gas6r,gas7r,gas8r,gas9r
      double precision :: gas1e,gas2e,gas3e,gas4e,gas5e,gas6e,gas7e,gas8e,gas9e

!     end update
!***********************************************************************

      data r0,p0,t0,gascon/1.292e00,1.0133e05,273.15e00,287.06e00/
      y=log10(rho/r0)
      x=log10(p/p0)
      if(abs(y+4.5e00).lt.2.5e-02) go to 20
      if(abs(y+0.5).lt.5.0e-03) go to 50
      iflag=-1
      go to 90
10    return
20    iflag=0
      rsave=rho
      ym=y
      y=-4.5e00+2.5e-02
      yhigh=y
      rho=(10.e0**y)*r0
      jflag=-1
      go to 90
30    thigh=t
      y=-4.5e00-2.5e-02
      ylow=y
      rho=(10.e0**y)*r0
      jflag=0
      go to 90
40    tlow=t
      go to 80
50    iflag=1
      rsave=rho
      ym=y
      y=-0.5+0.5e-02
      yhigh=y
      rho=(10.e0**y)*r0
      jflag=-1
      go to 90
60    thigh=t
      y=-0.5-0.5e-02
      ylow=y
      rho=(10.e0**y)*r0
      jflag=0
      go to 90
70    tlow=t
80    t=tlow+(thigh-tlow)/(yhigh-ylow)*(ym-ylow)
      rho=rsave
      return
90    z1=x-y
      if (y.gt.-0.5) go to 190
      if (y.gt.-4.5e00) go to 140
      if (z1.gt.0.25e00) go to 100
      t=p/(rho*gascon)
      go to 250
100   if (z1.gt.0.95e00) go to 110
      gas1=1.23718e-01+1.08623e-02*y
      gas2=(2.24239e-01-8.24608e-02*y)*z1
      gas3=(-1.17615e-03-1.87566e-03*z1-1.19155e-04*y)*y*y
      gas4=(1.18397e00+6.48520e-02*y-5.52634e-01*z1)*z1*z1
      tnon=gas1+gas2+gas3+gas4
      go to 240
110   if (z1.gt.1.4e00) go to 120
      gas1=-8.12952e00-8.28637e-01*y
      gas2=(2.26904e01+1.41132e00*y)*z1
      gas3=(-2.98633e-02+2.70066e-02*z1-2.28103e-04*y)*y*y
      gas4=(-1.91806e01-5.78875e-01*y+5.62580e00*z1)*z1*z1
      gas5=-3.99845e00+2.26369e-01*y
      gas6=(2.52876e00-7.28448e-01*y)*z1
      gas7=(1.09769e-02-1.83819e-02*z1-1.51380e-04*y)*y*y
      gas8=(2.99238e00+3.91440e-01*y-2.04463e00*z1)*z1*z1
      gas9=exp(-3.887015e01-2.908228e01*y+4.070557e01*z1                &
     &+2.682347e01*y*z1)
      go to 230
 120  if (z1.gt.1.95e00) go to 130
      gas1=-1.98573e01-1.67225e00*y
      gas2=(3.76159e01+2.10964e00*y)*z1
      gas3=(-3.40174e-02+2.31712e-02*z1-9.80275e-05*y)*y*y
      gas4=(-2.22215e01-6.44596e-01*y+4.40486e00*z1)*z1*z1
      gas5=-5.36809e00+2.41201e-01*y
      gas6=(-1.25881e00-8.62744e-01*y)*z1
      gas7=(-3.79774e-03-7.81335e-03*z1-3.80005e-04*y)*y*y
      gas8=(5.58609e00+3.78963e-01*y-1.81566e00*z1)*z1*z1
      gas9=exp(2.08e01-2.56e01*y+1.0e00*z1+1.80e01*y*z1)
      go to 230
 130  continue
!130   if (z1.gt.2.60e00) write (6,1000) rho,p
      gas1=-2.33271e01-1.89958e00*y
      gas2=(3.21440e01+1.68622e00*y)*z1
      gas3=(-4.42123e-02+2.82629e-02*z1+6.63272e-04*y)*y*y
      gas4=(-1.38645e01-3.40976e-01*y+2.04466e00*z1)*z1*z1
      gas5=8.35474e00+1.71347e00*y
      gas6=(-1.60715e01-1.63139e00*y)*z1
      gas7=(4.14641e-02-2.30068e-02*z1+1.53246e-05*y)*y*y
      gas8=(8.70275e00+3.60966e-01*y-1.46166e00*z1)*z1*z1
      gas9=exp(1.115884e02-6.452606e00*y-5.337863e01*z1                 &
     &+2.026986e00*y*z1)
      go to 230
140   if (z1.gt.0.25e00) go to 150
      t=p/(rho*gascon)
      go to 250
150   if (z1.gt.0.95e00) go to 160
      gas1=2.03910e-02+7.67310e-03*y
      gas2=(8.48581e-01-2.93086e-02*y)*z1
      gas3=(8.40269e-04-1.47701e-03*z1+3.13687e-05*y)*y*y
      gas4=(2.67251e-01+2.37262e-02*y-1.41973e-01*z1)*z1*z1
      tnon=gas1+gas2+gas3+gas4
      go to 240
160   if (z1.gt.1.45e00) go to 170
      gas1=-5.12404e00-2.84740e-01*y
      gas2=(1.54532e01+4.52475e-01*y)*z1
      gas3=(-1.22881e-02+8.56845e-03*z1-3.25256e-04*y)*y*y
      gas4=(-1.35181e01-1.68725e-01*y+4.18451e00*z1)*z1*z1
      gas5=7.52564e00+8.35238e-01*y
      gas6=(-1.95558e01-1.23393e00*y)*z1
      gas7=(3.34510e-02-2.34269e-02*z1+4.81788e-04*y)*y*y
      gas8=(1.71779e01+4.54628e-01*y-5.09936e00*z1)*z1*z1
      gas9=exp(6.148442e01-1.828123e01*y-5.468755e01*z1                 &
     &+1.562500e01*y*z1)
      go to 230
170   if (z1.gt.2.05e00) go to 180
      gas1=-1.23779e01-1.14728e00*y
      gas2=(2.41382e01+1.38957e00*y)*z1
      gas3=(-3.63693e-02+2.24265e-02*z1-3.23888e-04*y)*y*y
      gas4=(-1.42844e01-4.06553e-01*y+2.87620e00*z1)*z1*z1
      gas5=4.40782e00+1.33046e00*y
      gas6=(-1.15405e01-1.59892e00*y)*z1
      gas7=(5.30580e-02-3.10376e-02*z1+4.77650e-04*y)*y*y
      gas8=(8.57309e00+4.71274e-01*y-1.96233e00*z1)*z1*z1
      gas9=exp(1.4075e02-6.499992e00*y-7.75e01*z1+5.0e00*y*z1)
      go to 230
 180  continue
!180   if (z1.gt.2.50e00) write (6,1000) rho,p
      gas1=-1.27244e01-1.66684e00*y
      gas2=(1.72708e01+1.45307e00*y)*z1
      gas3=(-3.64515e-02+1.90463e-02*z1+4.80787e-04*y)*y*y
      gas4=(-6.97208e00-3.04323e-01*y+9.67524e-01*z1)*z1*z1
      gas5=7.71330e00+5.08340e-01*y
      gas6=(-9.82110e00-4.49138e-01*y)*z1
      gas7=(-9.41787e-04-2.40293e-03*z1-8.28450e-04*y)*y*y
      gas8=(4.16530e00+9.63923e-02*y-5.88807e-01*z1)*z1*z1
      gas9=exp(-1.092654e03-3.05312e02*y+4.656243e02*z1+                &
     &1.312498e02*y*z1)
      go to 230
190   if (z1.gt.0.25e00) go to 200
      t=p/(rho*gascon)
      go to 250
200   if (z1.gt.1.00e00) go to 210
      gas1=-1.54141e-03+6.58337e-04*y
      gas2=(9.82201e-01-3.85028e-03*y)*z1
      gas3=(1.23111e-04-4.08210e-04*z1+2.13592e-05*y)*y*y
      gas4=(3.77441e-02+4.56963e-03*y-2.35172e-02*z1)*z1*z1
      tnon=gas1+gas2+gas3+gas4
      go to 240
210   if (z1.gt.1.45e00) go to 220
      gas1=8.06492e-01+9.91293e-02*y
      gas2=(-1.70742e00-2.28264e-01*y)*z1
      gas3=(5.03500e-03-6.13927e-03*z1+1.69824e-04*y)*y*y
      gas4=(3.02351e00+1.31574e-01*y-1.12755e00*z1)*z1*z1
      gas5=-1.17930e-01-2.12207e-01*y
      gas6=(1.36524e00+4.05886e-01*y)*z1
      gas7=(-1.88260e-02+1.65486e-02*z1-5.11400e-04*y)*y*y
      gas8=(-2.10926e00-1.89881e-01*y+8.79806e-01*z1)*z1*z1
      gas9=exp(1.959604e02-4.269391e01*y-1.734931e02*z1+                &
     &3.762898e01*y*z1)
      go to 230
 220  continue
!220   if (z1.gt.2.3e00) write (6,1000) rho,p
      gas1=-1.66249e00-8.91113e-02*y
      gas2=(4.11648e00+8.78093e-02*y)*z1
      gas3=(-3.09742e-03+1.99879e-03*z1+6.85472e-05*y)*y*y
      gas4=(-1.84445e00-7.50324e-03*y+3.05784e-01*z1)*z1*z1
      gas5=1.11555e01+1.32100e00*y
      gas6=(-1.71236e01-1.29190e00*y)*z1
      gas7=(6.28124e-02-3.07949e-02*z1+1.57743e-03*y)*y*y
      gas8=(8.63804e00+3.07809e-01*y-1.42634e00*z1)*z1*z1
      gas9=exp(1.330611e02+8.979635e00*y-7.265298e01*z1                 &
     &-2.449009e00*y*z1)
230   tnon=gas1+gas2+gas3+gas4+(gas5+gas6+gas7+gas8)/(1.e0+gas9)
240   t=(10.e0**tnon)*t0
250   if (iflag) 10,260,270
260   if (jflag) 30,40,10
270   if (jflag) 60,70,10
1000  format(/20x,48hwarning!  outside of validity range of curve fit   &
     &,/,20x,5hrho =,1pd15.8,5x,3hp =,1pd15.8,/)
      end subroutine tgas3

      subroutine tgas4(p,rho,h,e)
!
!     Last update:   Mar 17 0 1993
!     Last edited:   Mar 17 0 1993
!***********************************************************************
!**********                                                 ************
!**********                    enthalpy                     ************
!**********            returned as a function of            ************
!**********              pressure and density               ************
!**********                                                 ************
!***********************************************************************
!                  new routine by srinivasan et.al.
!     *********************************************************
!
!     inputs:
!
!     p = pressure, in newtons/m**2.
!     rho =  density, kg/m**3.
!
!     outputs:
!
!     h = specific enthalpy, in (m/sec)**2.
!
!     *********************************************************
!
!***********************************************************************
!     update by Stanley: Apr 22, 2010
!     Purpose: define variables used in the subroutines
! ************************************************************************
!     update by Stanley: Apr 23,2010
!     Purpose: Compute internal energy also, e=h/gamm
!     Based on the NASA report by Srinivasan, Tannehill and Weilmuenster
!     begin update

      implicit none

      double precision :: p, rho
      double precision, intent (out):: h,e 

      integer :: lflag, kflag, iflag, jflag
      double precision :: gammr, gamme, gamm
      double precision :: e0,r0,p0,t0,gascon
      double precision :: rratio, eratio, x, y, z, z1, tnon
      double precision :: rsave, ym, yhigh, ylow, phigh, plow, hhigh, hlow
      double precision :: gas1 ,gas2 ,gas3 ,gas4 ,gas5 ,gas6 ,gas7 ,gas8 ,gas9
      double precision :: gas1r,gas2r,gas3r,gas4r,gas5r,gas6r,gas7r,gas8r,gas9r
      double precision :: gas1e,gas2e,gas3e,gas4e,gas5e,gas6e,gas7e,gas8e,gas9e

!     end update
!***********************************************************************
      data r0,p0/1.292d00,1.0133d05/
!
      y=log10(rho/r0)
      x=log10(p/p0)
      if(abs(y+4.5e00).lt.2.5e-02) go to 20
      if(abs(y+0.5e00).lt.5.0e-03) go to 50
      iflag=-1
      go to 90
10    return
20    iflag=0
      rsave=rho
      ym=y
      y=-4.5e00+2.5e-02
      yhigh=y
      rho=(10.**y)*r0
      jflag=-1
      go to 90
30    hhigh=h
      y=-4.5e00-2.5e-02
      ylow=y
      rho=(10.**y)*r0
      jflag=0
      go to 90
40    hlow=h
      go to 80
50    iflag=1
      rsave=rho
      ym=y
      y=-0.5e00+0.5e-02
      yhigh=y
      rho=(10.**y)*r0
      jflag=-1
      go to 90
60    hhigh=h
      y=-0.5e00-0.5e-02
      ylow=y
      rho=(10.**y)*r0
      jflag=0
      go to 90
70    hlow=h
80    h=hlow+(hhigh-hlow)/(yhigh-ylow)*(ym-ylow)
      rho=rsave
      go to 10
90    z1=x-y
      if (y.gt.-0.5e00) go to 190
      if (y.gt.-4.5e00) go to 140
      if (z1.gt.0.10e00) go to 100
      gamm=1.3986e00
      go to 240
100   if (z1.gt.0.85e00) go to 110
      gas1=2.53908e02+1.01491e02*y
      gas2=(-3.87199e02-1.54304e02*y)*z1
      gas3=(7.28532e00-8.04378e00*z1-1.82577e-03*y)*y*y
      gas4=(9.86233e01+4.63763e01*y+2.18994e01*z1)*z1*z1
      gas5=-2.52423e02-1.01445e02*y
      gas6=(3.87210e02+1.54298e02*y)*z1
      gas7=(-7.2773e00+8.04277e00*z1+2.28399e-03*y)*y*y
      gas8=(-9.87576e01-4.63883e01*y-2.19438e01*z1)*z1*z1
      gas9=exp(-11.e00+2.e00*y+11.e00*z1-2.e00*y*z1)
      gamm=gas1+gas2+gas3+gas4+(gas5+gas6+gas7+gas8)/(1.-gas9)
      go to 240
110   if (z1.gt.1.30e00) go to 120
      gas1=-1.05745e01-1.93693e00*y
      gas2=(3.07202e01+3.35578e00*y)*z1
      gas3=(-7.79965e-02+6.68790e-02*z1-9.86882e-04*y)*y*y
      gas4=(-2.60637e01-1.42391e00*y+7.23223e00*z1)*z1*z1
      gas5=-1.86342e01+2.41997e-02*y
      gas6=(3.20880e01-7.46914e-01*y)*z1
      gas7=(3.75161e-02-4.10125e-02*z1+5.74637e-04*y)*y*y
      gas8=(-1.69985e01+5.39041e-01*y+2.56253e00*z1)*z1*z1
      gas9=exp(2.768567e02+2.152383e01*y-2.164837e02*z1                 &
     &-1.394837e01*y*z1)
      go to 230
120   if (z1.gt.1.95e00) go to 130
      gas1=6.17584e-01-2.40690e-01*y
      gas2=(1.95904e00+3.41644e-01*y)*z1
      gas3=(-1.01073e-02+6.77631e-03*z1-1.15922e-04*y)*y*y
      gas4=(-1.68951e00-1.10932e-01*y+4.26058e-01*z1)*z1*z1
      gas5=-1.34222e01-5.43713e-01*y
      gas6=(1.81528e01+3.95928e-01*y)*z1
      gas7=(-7.41105e-03+1.67768e-03*z1-3.32714e-06*y)*y*y
      gas8=(-7.97425e00-5.80593e-02*y+1.12448e00*z1)*z1*z1
      gas9=exp(8.677803e01-8.370349e00*y-4.074084e01*z1                 &
     &+7.407405e00*y*z1)
      go to 230
 130  continue
!130   if (z1.gt.2.60e00) write (6,1000) rho,p
      gas1=-8.32595e00-3.50219e-01*y
      gas2=(1.36455e01+3.59350e-01*y)*z1
      gas3=(-3.70109e-03+3.30836e-03*z1+1.10018e-04*y)*y*y
      gas4=(-6.49007e00-8.38594e-02*y+1.02443e00*z1)*z1*z1
      gas5=-3.08441e01-1.49510e00*y
      gas6=(3.00585e01+9.19650e-01*y)*z1
      gas7=(-3.60024e-02+1.02522e-02*z1-4.68760e-04*y)*y*y
      gas8=(-9.33522e00-1.35228e-01*y+8.92634e-01*z1)*z1*z1
      gas9=exp(8.800047e01-1.679356e01*y-3.333353e01*z1                 &
     &+8.465574e00*y*z1)
      go to 230
140   if (z1.gt.0.1e00) go to 150
      gamm=1.399e00
      go to 240
150   if (z1.gt.0.95e00) go to 160
      gas1=-1.33083e02-9.98707e00*y
      gas2=(3.94734e02+2.35810e01*y)*z1
      gas3=(1.43957e00-1.43175e00*z1+1.77068e-05*y)*y*y
      gas4=(-3.84712e02-1.36367e01*y+1.24325e02*z1)*z1*z1
      gas5=1.34486e02+9.99122e00*y
      gas6=(-3.94719e02-2.35853e01*y)*z1
      gas7=(-1.43799e00+1.43039e00*z1+1.44367e-04*y)*y*y
      gas8=(3.84616e02+1.36318e01*y-1.24348e02*z1)*z1*z1
      gas9=exp(-2.141444e01+1.381584e00*y+2.039473e01*z1                &
     &-1.315789e00*y*z1)
      gamm=gas1+gas2+gas3+gas4+(gas5+gas6+gas7+gas8)/(1.-gas9)
      go to 240
160   if (z1.gt.1.50e00) go to 170
      gas1=-7.36684e00-1.13247e00*y
      gas2=(2.47879e01+1.99625e00*y)*z1
      gas3=(-4.91630e-02+4.16673e-02*z1-6.58149e-04*y)*y*y
      gas4=(-2.32990e01-8.59418e-01*y+7.19016e00*z1)*z1*z1
      gas5=-2.42647e00+5.57912e-01*y
      gas6=(-2.03055e00-1.22031e00*y)*z1
      gas7=(3.74866e-02-3.39278e-02*z1+5.21042e-04*y)*y*y
      gas8=(7.75414e00+6.08488e-01*y-3.68326e00*z1)*z1*z1
      gas9=exp(8.077385e01-1.273807e01*y-6.547623e01*z1                 &
     &+1.190475e01*y*z1)
      go to 230
170   if (z1.gt.2.00e00) go to 180
      gas1=4.31520e-01-2.83857e-01*y
      gas2=(2.27791e00+3.99159e-01*y)*z1
      gas3=(-1.29444e-02+8.78724e-03*z1-1.60583e-04*y)*y*y
      gas4=(-1.84314e00-1.28136e-01*y+4.45362e-01*z1)*z1*z1
      gas5=-1.03883e01-3.58718e-01*y
      gas6=(1.35068e01+1.87268e-01*y)*z1
      gas7=(-4.28184e-03-9.52016e-04*z1-4.10506e-05*y)*y*y
      gas8=(-5.63894e00-1.45626e-03*y+7.39915e-01*z1)*z1*z1
      gas9=exp(2.949221e02+1.368660e01*y-1.559335e02*z1                 &
     &-3.787766e00*y*z1)
      go to 230
 180  continue
!180   if (z1.gt.2.5e00) write (6,1000) rho,p
      gas1=-3.77766e00-5.53738e-01*y
      gas2=(6.60834e00+4.87181e-01*y)*z1
      gas3=(-2.11045e-02+9.67277e-03*z1-2.19420e-04*y)*y*y
      gas4=(-2.94754e00-1.02365e-01*y+4.39620e-01*z1)*z1*z1
      gas5=4.05813e01+3.25692e00*y
      gas6=(-4.79583e01-2.53660e00*y)*z1
      gas7=(9.06436e-02-3.47578e-02*z1+1.00077e-03*y)*y*y
      gas8=(1.89040e01+4.94114e-01*y-2.48554e00*z1)*z1*z1
      gas9=exp(5.34718e02+7.495657e01*y-2.219822e02*z1                  &
     &-3.017229e01*y*z1)
      go to 230
190   if (z1.gt.0.1e00) go to 200
      gamm=1.4017e00
      go to 240
200   if (z1.gt.1.05e00) go to 210
      gas1=-9.67488e01+2.05296e-01*y
      gas2=(2.69927e02-1.92887e00*y)*z1
      gas3=(3.78392e-01-3.24965e-01*z1-3.61036e-03*y)*y*y
      gas4=(-2.46711e02+1.54416e00*y+7.48760e01*z1)*z1*z1
      gas5=9.81502e01-2.05448e-01*y
      gas6=(-2.69913e02+1.93052e00*y)*z1
      gas7=(-3.78527e-01+3.24832e-01*z1+3.66182e-03*y)*y*y
      gas8=(2.46630e02-1.54646e00*y-7.48980e01*z1)*z1*z1
      gas9=exp(-2.659865e01+1.564631e00*y+2.312926e01*z1                &
     &-1.360543e00*y*z1)
      gamm=gas1+gas2+gas3+gas4+(gas5+gas6+gas7+gas8)/(1.-gas9)
      go to 240
210   if (z1.gt.1.60e00) go to 220
      gas1=-2.67593e-01-1.87457e-01*y
      gas2=(5.07693e00+2.72286e-01*y)*z1
      gas3=(1.04541e-02-1.42211e-02*z1+6.38962e-04*y)*y*y
      gas4=(-5.08520e00-7.81935e-02*y+1.58711e00*z1)*z1*z1
      gas5=2.87969e00+3.9009e-01*y
      gas6=(-8.06179e00-5.51250e-01*y)*z1
      gas7=(-1.01903e-02+1.35906e-02*z1-8.97772e-04*y)*y*y
      gas8=(7.29592e00+1.83861e-01*y-2.15153e00*z1)*z1*z1
      gas9=exp(1.828573e02-3.428596e01*y-1.51786e02*z1                  &
     &+2.976212e01*y*z1)
      go to 230
 220  continue
!220   if (z1.gt.2.30e00) write (6,1000) rho,p
      gas1=9.21537e-01-2.39670e-01*y
      gas2=(1.30714e00+3.42990e-01*y)*z1
      gas3=(-2.18847e-02+1.36691e-02*z1-4.90274e-04*y)*y*y
      gas4=(-1.20916e00-1.10206e-01*y+3.087920e-01*z1)*z1*z1
      gas5=-6.77089e00-6.90476e-02*y
      gas6=(8.18168e00-9.52708e-02*y)*z1
      gas7=(2.98487e-02-1.78706e-02*z1+6.28419e-04*y)*y*y
      gas8=(-3.07662e00+6.60408e-02*y+3.38590e-01*z1)*z1*z1
      gas9=exp(1.5916669e02+3.976192e01*y-7.966199e01*z1                &
     &-1.66667e01*y*z1)
230   gamm=gas1+gas2+gas3+gas4+(gas5+gas6+gas7+gas8)/(1.+gas9)
240   h=gamm/(gamm-1.0e00)*p/rho
! ************************************************************************
! Compute internal energy, added by Stanley Apr23,2010
! Based on the NASA report by Srinivasan, Tannehill and Weilmuenster
! (NASA Reference Publication 1181, 1987)
      e=h/gamm
! ************************************************************************
      if (iflag) 10,250,260
250   if (jflag) 30,40,10
260   if (jflag) 60,70,10
1000  format(/20x,48hwarning!  outside of validity range of curve fit   &
     &,/,20x,5hrho =,1pe15.8,5x,3hp =,1pe15.8,/)
      end subroutine tgas4


      subroutine ugas3(e,rho,mu)
!
!     Last update:   Mar 17 0 1993
!     Last edited:   Mar 17 0 1993
!***********************************************************************
!**********                                                 ************
!**********                  viscosity                      ************
!**********            returned as a function of            ************
!**********          internal energy and density            ************
!**********                                                 ************
!***********************************************************************
!                  new routine by srinivasan et.al.
!     *********************************************************
!
!     inputs:
!
!     e =  internal energy, j/kg ;or; (m/s)**2
!     r =  density, kg/(m**3)
!
!     output:
!
!     mu =  coefficient of dynamic viscosity, kg/m-s ;
!
!     *********************************************************
!
!***********************************************************************
!     update by Stanley: Mar 02, 2010
!     Purpose: define variables used in the subroutines
! ************************************************************************

      implicit none

      double precision :: e, rho
      double precision, intent (out):: mu

      double precision :: x,y,z,rho0,e0,t,f
      double precision :: gas1 ,gas2 ,gas3 ,gas4 ,gas5 ,gas6 ,gas7 ,gas8 ,gas9

!     end update
!***********************************************************************

      data rho0,e0/1.243e0,78408.4e0/
!
      z    = log10(e/e0)
      y    = log10(rho/rho0)
!
!***  new modifications (11-6-87)
!
!      if ((y .lt. -5.0e00) .or. (y .gt. 1.0e00))
!     >   write (6,1000) rho,e
      if (z .gt. 0.44e00) go to 5
      t    = 0.4e00*e/287.06e00
      mu   = 1.462e-06*sqrt(t) / (1.0e00+112.0e00/t)
      return
!***
 5    if (z .gt. 0.67e00) go to 10
      gas1 = 4.84547e-01+4.67135e-01*z
      gas2 = (5.71205e-04-1.43629e-03*z)*y
      gas3 = (2.55110e00-2.33472e-04*y-1.44102e00*z)*z*z
      gas4 = (2.53416e-04-4.72375e-04*z+1.86899e-05*y)*y*y
      f    = gas1+gas2+gas3+gas4
      go to 90
 10   if (z.gt.1.75e00) go to 20
      gas1 = -3.71666e01+6.67883e01*z
      gas2 = (-2.43998e00+2.12309e00*z)*y
      gas3 = (-3.69259e01-3.08426e-01*y+7.36486e00*z)*z*z
      gas4 = (-1.46446e-01+7.54423e-02*z-2.91464e-03*y)*y*y
      gas5 = 3.61757e01-6.11102e01*z
      gas6 = (2.40531e00-2.05914e00*z)*y
      gas7 = (3.23911e01+2.79149e-01*y-5.07640e00*z)*z*z
      gas8 = (1.37916e-01-6.72041e-02*z+2.61987e-03*y)*y*y
      gas9 = exp(-3.433e01-1.823e00*y+2.499e01*z+6.503e-01*z*y)
      go to 80
 20   if (z.gt.2.50e00) go to 30
      gas1 = -1.65147e02+2.11028e02*z
      gas2 = (-4.70948e00+2.78258e00*z)*y
      gas3 = (-8.78308e01-1.28671e-01*y+1.27639e01*z)*z*z
      gas4 = (-3.19867e-01+1.73179e-01*z+3.86106e-03*y)*y*y
      gas5 = 2.30407e02-2.98055e02*z
      gas6 = (-6.18307e00+8.44595e00*z)*y
      gas7 = (1.26933e02-2.61671e00*y-1.77257e01*z)*z*z
      gas8 = (-2.30229e-02+2.25458e-02*z-4.41072e-03*y)*y*y
      gas9 = exp(-6.882e01+8.824e00*y+3.203e01*z-5.359e00*z*y)
      go to 80
 30   if (z.gt.2.85e00) go to 40
      gas1 = -7.09274e03+7.13648e03*z
      gas2 = (-2.46014e02+1.65826e02*z)*y
      gas3 = (-2.37952e03-2.75487e01*y+2.63465e02*z)*z*z
      gas4 = (-3.49744e00+1.28641e00*z-3.13711e-03*y)*y*y
      gas5 = 5.26158e03-4.96701e03*z
      gas6 = (2.03138e02-1.32984e02*z)*y
      gas7 = (1.52424e03+2.15081e01*y-1.50450e02*z)*z*z
      gas8 = (3.32432e00-1.15997e00*z+1.14862e-02*y)*y*y
      gas9 = exp(-3.594e02-3.763e01*y+1.319e02*z+1.348e01*z*y)
      go to 80
 40   if (z.gt.3.15e00) go to 50
      gas1 = -1.27748e03+1.29400e03*z
      gas2 = (-3.60724e01+2.63194e01*z)*y
      gas3 = (-4.22958e02-4.38228e00*y+4.50571e01*z)*z*z
      gas4 = (-4.74425e-01+2.89684e-01*z+1.64048e-02*y)*y*y
      f    = gas1+gas2+gas3+gas4
      go to 90
 50   if (y.gt.-3.80e00) go to 70
      if (z.gt.3.19e00) go to 60
      gas1 = 4.55919e03-4.21057e03*z
      gas2 = (1.03001e01-2.63478e01*z)*y
      gas3 = (1.29069e03+6.59587e00*y-1.31413e02*z)*z*z
      gas4 = (-8.28137e00+1.9827e00*z-1.7287e-01*y)*y*y
      f    = gas1+gas2+gas3+gas4
      go to 90
 60   z    = e/e0
      gas1 = -4.41792e02+9.7986e-02*z
      gas2 = (-3.03148e02+7.6065e-03*z)*y
      gas3 = (-5.5711e-05-3.52836e-06*y+8.86148e-09*z)*z*z
      gas4 = (-7.561e01-4.76816e-04*z-6.48859e00*y)*y*y
      gas5 = 6.72387e04+3.28398e00*z
      gas6 = (3.55009e04+2.72616e00*z)*y
      gas7 = (2.13714e-03+3.42377e-04*y-6.84897e-08*z)*z*z
      gas8 = (6.50886e03+3.8056e-01*z+4.14116e02*y)*y*y
      gas9 = exp(2.978e01+5.415e00*y+1.713e-03*z+3.115e-04*y*z)
      f    = gas1+gas2+gas3+gas4+(gas5+gas6+gas7+gas8)/(1.e00-gas9)
      go to 90
 70   gas1 = -6.4029e03+6.24254e03*z
      gas2 = (1.03279e02-8.73181e01*z)*y
      gas3 = (-2.02865e03+1.71878e01*y+2.19907e02*z)*z*z
      gas4 = (-1.22397e01+3.57830e00*z-1.27953e-01*y)*y*y
      f    = gas1+gas2+gas3+gas4
      go to 90
 80   f    = gas1+gas2+gas3+gas4+(gas5+gas6+gas7+gas8)/(1.e00+gas9)
 90   mu   = 1.748583e-05*f
      go to 110
 1000 format(/20x,48hwarning! outside of validity range of curve fit    &
     &,/,20x,5hrho =,1pd15.8,5x,3he =,1pd15.8,/)
 110  return
      end subroutine ugas3

      subroutine ugas4(e,rho,k)
!
!     Last update:   Mar 17 0 1993
!     Last edited:   Mar 17 0 1993
!***********************************************************************
!**********                                                 ************
!**********              thermal conductivity               ************
!**********            returned as a function of            ************
!**********          internal energy and density            ************
!**********                                                 ************
!***********************************************************************
!                  new routine by srinivasan et.al.
!     *********************************************************
!
!     inputs:
!
!     e   =  internal energy, j/kg ;or; (m/s)**2
!     rho =  density, kg/(m**3)
!
!     output:
!
!     k = coefficient of thermal conductivity, in j/(kelvin*m*s)
!
!     *********************************************************
!
!***********************************************************************
!     update by Stanley: Mar 02, 2010
!     Purpose: define variables used in the subroutines
! ************************************************************************

      implicit none

      double precision :: e, rho
      double precision, intent (out):: k 

      double precision :: x,y,z,rho0,e0,t,f
      double precision :: gas1 ,gas2 ,gas3 ,gas4 ,gas5 ,gas6 ,gas7 ,gas8 ,gas9

!     end update
!***********************************************************************

      data rho0,e0/1.243e00,78408.4e00/
      z     = log10(e/e0)
      y     = log10(rho/rho0)
!
!***  new modifications (11-6-87)
!
!      if ((y .lt. -5.0e00) .or. (y .gt. 1.0e00))
!     >   write (6,1000) rho,e
      if (z .gt. 0.44e00) go to 5
      t     = 0.4e00*e/287.06e00
      k     = 1.994e-03*sqrt(t) / (1.0e00+112.0e00/t)
      return
!***
 5    if (z .gt. 0.65e00) go to 10
      gas1  = 1.8100369e-01+4.8126802e00*z
      gas2  = (-2.7231116e-02+1.2691337e-01*z)*y
      gas3  = (-8.9913034e00-1.2624085e-01*y+8.9649105e00*z)*z*z
      gas4  = (-4.7198236e-03+9.2328079e-03*z-2.9488327e-04*y)*y*y
      f     = gas1+gas2+gas3+gas4
      go to 200
 10   if (y .gt. -1.00e00) go to 130
      if (y .gt. -3.00e00) go to 70
      if (z .gt. 1.25e00) go to 20
      gas1  = -1.05935e04+2.31470e04*z
      gas2  = (-7.41294e02+1.21724e03*z)*y
      gas3  = (-1.67601e04-4.43184e02*y+4.06631e03*z)*z*z
      gas4  = (1.35105e01+4.94914e00*z+1.55386e00*y)*y*y
      gas5  = 1.06032e04-2.31560e04*z
      gas6  = (7.46951e02-1.22465e03*z)*y
      gas7  = (1.67604e04+4.45919e02*y-4.06258e03*z)*z*z
      gas8  = (-1.28615e01-5.32398e00*z-1.52956e00*y)*y*y
      gas9  = exp(-4.219e01-4.687e00*y+2.812e01*z+3.125e00*y*z)
      f     = gas1+gas2+gas3+gas4+(gas5+gas6+gas7+gas8)/(1.0-gas9)
      go to 200
 20   if (z .gt. 1.775e00) go to 30
      gas1  = 3.79375e03-7.40351e03*z
      gas2  = (3.29698e02-3.55916e02*z)*y
      gas3  = (4.77122e03+1.00241e02*y-1.00740e03*z)*z*z
      gas4  = (1.97061e01-8.42554e00*z+4.80494e-01*y)*y*y
      gas5  = -4.53603e03+9.05605e03*z
      gas6  = (-4.95870e02+6.33563e02*z)*y
      gas7  = (-5.95317e03-2.05442e02*y+1.28945e03*z)*z*z
      gas8  = (-2.00087e01+1.18851e01*z-1.71735e-01*y)*y*y
      gas9  = exp(-3.318e01+3.158e-01*y+1.863e01*z-1.035e00*y*z)
      go to 190
 30   if (z .gt. 1.93e00) go to 40
      gas1  = 2.06651875e05-3.165645e05*z
      gas2  = (-3.07322021e02+4.57036377e02*z)*y
      gas3  = (1.61824937e05-1.55508453e02*y-2.7603957e04*z)*z*z
      gas4  = (1.92260265e00-2.24788094e00*z-3.06226015e-01*y)*y*y
      gas5  = -2.06564312e05+3.18191312e05*z
      gas6  = (2.17542285e03-2.46670776e03*z)*y
      gas7  = (-1.63597062e05+7.16753174e02*y+2.80926367e04*z)*z*z
      gas8  = (3.39526825e01-7.53846645e00*z+1.91214371e00*y)*y*y
      gas9  = exp(-3.924e02-5.206e01*y+2.054e02*z+2.679e01*y*z)
      go to 190
 40   if (z .gt. 2.60e00) go to 50
      gas1  = 7.1572625e04-9.2471625e04*z
      gas2  = (1.9646323e03-2.0280527e03*z)*y
      gas3  = (3.9446105e04+4.5673853e02*y-5.5728672e03*z)*z*z
      gas4  = (-9.2131958e01+1.2724541e01*z-5.0568476e00*y)*y*y
      gas5  = -3.2910781e04+4.2551211e04*z
      gas6  = (1.4566331e03-2.2653745e03*z)*y
      gas7  = (-1.9476277e04+8.4370288e02*y+3.2389702e03*z)*z*z
      gas8  = (-1.3324594e02+1.0591533e02*z+5.8639469e00*y)*y*y
      gas9  = exp(4.917e01+2.415e01*y-2.455e01*z-1.181e01*y*z)
      go to 190
 50   if (z .gt. 2.69e00) go to 60
      gas1  = 1.145683e06-1.237525e06*z
      gas2  = (1.4024508e04-9.3467227e03*z)*y
      gas3  = (4.4593056e05+1.533074e03*y-5.3608352e04*z)*z*z
      gas4  = (2.8485107e02-1.0968916e02*z-1.0955791e00*y)*y*y
      gas5  = -1.752087e06+1.79675e06*z
      gas6  = (-1.3278737e05+9.8215562e04*z)*y
      gas7  = (-6.0791744e05-1.811943e04*y+6.7709875e04*z)*z*z
      gas8  = (-1.3384084e03+5.2707324e02*z+2.5904894e00*y)*y*y
      gas9  = exp(-1.798e02+7.371e00*y+6.731e01*z-3.205e00*y*z)
      go to 190
 60   gas1  = -8.5499625e04+1.1739656e05*z
      gas2  = (6.4563168e04-3.9551203e04*z)*y
      gas3  = (-4.8170254e04+6.0816055e03*y+6.2052031e03*z)*z*z
      gas4  = (2.3473167e-01+1.8871567e01*z+4.0757723e00*y)*y*y
      gas5  = 5.8546883e04-9.4634875e04*z
      gas6  = (-6.6513812e04+4.0899945e04*z)*y
      gas7  = (4.2127227e04-6.3717305e03*y-5.7495195e03*z)*z*z
      gas8  = (-1.0260344e00-5.343277e01*z-1.1017392e01*y)*y*y
      gas9  = exp(5.411e00+1.162e01*y-1.082e00*z-3.391e00*y*z)
      f     = gas1+gas2+gas3+gas4+(gas5+gas6+gas7+gas8)/(1.0-gas9)
      go to 200
 70   if (z .gt. 1.29e00) go to 80
      gas1  = -1.22493e04+2.41071e04*z
      gas2  = (-1.61829e03+2.22535e03*z)*y
      gas3  = (-1.59261e04-7.53213e02*y+3.53376e03*z)*z*z
      gas4  = (1.98026e00+5.18483e00*z+1.47851e00*y)*y*y
      gas5  = 1.22486e04-2.41023e04*z
      gas6  = (1.61810e03-2.22571e03*z)*y
      gas7  = (1.59235e04+7.53746e02*y-3.53168e03*z)*z*z
      gas8  = (-2.15482e00-5.05115e00*z-1.48795e00*y)*y*y
      gas9  = exp(-3.111e01-4.444e00*y+1.944e01*z+2.778e00*y*z)
      f     = gas1+gas2+gas3+gas4+(gas5+gas6+gas7+gas8)/(1.0-gas9)
      go to 200
 80   if (z .gt. 1.85e00) go to 90
      gas1  = 3.18060e03-6.69664e03*z
      gas2  = (4.33382e01-2.14649e02*z)*y
      gas3  = (4.41377e03+9.41359e01*y-9.29758e02*z)*z*z
      gas4  = (-3.62190e01+1.15538e01*z-2.14621e00*y)*y*y
      gas5  = -5.98764e03+1.29243e04*z
      gas6  = (-2.72261e02+5.42378e02*z)*y
      gas7  = (-9.03293e03-2.11787e02*y+2.07831e03*z)*z*z
      gas8  = (2.74179e01-5.68578e00*z+1.91217e00*y)*y*y
      gas9  = exp(-1.854e01+7.11e00*y+1.068e01*z-5.449e00*y*z)
      go to 190
 90   if (z .gt. 2.0e00) go to 100
      gas1  = 5.14024e04-7.52733e04*z
      gas2  = (-3.30889e02+3.11550e02*z)*y
      gas3  = (3.66539e04-7.41227e01*y-5.93015e03*z)*z*z
      gas4  = (-4.84164e01+2.23133e01*z-9.19118e-01*y)*y*y
      gas5  = -1.80898e05+2.82532e05*z
      gas6  = (-1.01053e03+9.75576e02*z)*y
      gas7  = (-1.47220e05-2.33631e02*y+2.55940e04*z)*z*z
      gas8  = (3.28681e00-1.76588e00*z-1.54962e-01*y)*y*y
      gas9  = exp(-4.104e01+6.507e01*y+2.083e01*z-3.472e01*z*y)
      go to 190
 100  if (z .gt. 2.58e00) go to 110
      gas1  = 5.1131824e04-6.664875e04*z
      gas2  = (2.02171e03-1.9306292e03*z)*y
      gas3  = (2.8762395e04+4.3353467e02*y-4.1064609e03*z)*z*z
      gas4  = (-8.4970047e01+1.7925919e01*z-6.2576542e00*y)*y*y
      gas5  = -6.2768156e04+8.6015875e04*z
      gas6  = (-1.0002036e03+6.2537280e02*z)*y
      gas7  = (-3.957827e04-3.8467377e01*y+6.12953e03*z)*z*z
      gas8  = (-1.0591702e02+7.636142e01*z+5.938859e00*y)*y*y
      gas9  = exp(-3.901e00+2.418e01*y+1.374e00*z-1.145e01*y*z)
      go to 190
 110  if (z .gt. 2.73e00) go to 120
      gas1  = 1.0088046e06-1.086321e06*z
      gas2  = (1.3844801e04-9.7268516e03*z)*y
      gas3  = (3.8985325e05+1.7091665e03*y-4.6621066e04*z)*z*z
      gas4  = (1.4840726e02-5.2645004e01*z-1.5477133e-01*y)*y*y
      gas5  = -1.073351e06+1.14571e06*z
      gas6  = (-1.9343957e04+1.3366211e04*z)*y
      gas7  = (-4.0670987e05-2.2955198e03*y+4.7999871e04*z)*z*z
      gas8  = (-4.1016724e02+1.4994148e02*z-1.9779787e00*y)*y*y
      gas9  = exp(-1.026e02+6.302e01*y+3.819e01*z-2.431e01*y*z)
      go to 190
 120  gas1  = -9.6638500e04+1.3206488e04*z
      gas2  = (-4.7458105e04+2.3596875e04*z)*y
      gas3  = (1.8602773e04-2.306802e03*y-4.0413552e03*z)*z*z
      gas4  = (-5.3564258e03+2.2433904e03*z+2.5188145e02*y)*y*y
      gas5  = 1.0962581e05-2.990116e04*z
      gas6  = (4.7883496e04-2.3785383e04*z)*y
      gas7  = (-1.1753969e04+2.2905522e03*y+3.1304399e03*z)*z*z
      gas8  = (5.473418e03-2.3208018e03*z-2.6570068e02*y)*y*y
      gas9  = exp(-3.107e01+1.082e01*y+1.047e01*z-3.047e00*y*z)
      f     = gas1+gas2+gas3+gas4+(gas5+gas6+gas7+gas8)/(1.0-gas9)
      go to 200
 130  if (z .gt. 1.40e00) go to 140
      gas1  = -1.58386e03+3.49223e03*z
      gas2  = (-8.39834e02+1.09565e03*z)*y
      gas3  = (-2.56175e03-3.56197e02*y+6.25145e02*z)*z*z
      gas4  = (-1.22407e01+7.65634e00*z+2.58235e-01*y)*y*y
      gas5  = 1.58025e03-3.47664e03*z
      gas6  = (8.39588e02-1.09490e03*z)*y
      gas7  = (2.54682e03+3.55674e02*y-6.18504e02*z)*z*z
      gas8  = (1.20843e01-7.44857e00*z-2.91202e-01*y)*y*y
      gas9  = exp(-2.171e01-4.342e00*y+1.316e01*z+2.632e00*y*z)
      f     = gas1+gas2+gas3+gas4+(gas5+gas6+gas7+gas8)/(1.0-gas9)
      go to 200
 140  if (z .gt. 1.91e00) go to 150
      gas1  = 7.89255e02-1.91743e03*z
      gas2  = (3.59227e02-4.44070e02*z)*y
      gas3  = (1.39463e03+1.34083e02*y-3.13446e02*z)*z*z
      gas4  = (1.90681e01-1.09285e01*z+4.24933e-02*y)*y*y
      gas5  = -1.31401e03+3.13134e03*z
      gas6  = (-5.18755e02+6.80268e02*z)*y
      gas7  = (-2.32493e03-2.21393e02*y+5.52563e02*z)*z*z
      gas8  = (-3.32001e01+2.11819e01*z-4.75163e-01*y)*y*y
      gas9  = exp(-5.025e01-8.412e00*y+2.982e01*z+3.509e00*y*z)
      go to 190
 150  if (z .gt. 2.05e00) go to 160
      gas1  = 3.58691e04-5.16852e04*z
      gas2  = (-6.30189e02+6.63314e02*z)*y
      gas3  = (2.47471e04-1.73538e02*y-3.93167e03*z)*z*z
      gas4  = (-4.23871e01+2.08048e01*z-1.05512e00*y)*y*y
      gas5  = -1.10522e05+1.67591e05*z
      gas6  = (4.61877e03-4.94930e03*z)*y
      gas7  = (-8.46558e04+1.32441e03*y+1.42438e04*z)*z*z
      gas8  = (2.25065e01-1.10316e01*z+9.62887e-01*y)*y*y
      gas9  = exp(-1.681e02+7.063e01*y+8.75e01*z-3.75e01*y*z)
      go to 190
 160  if (z .gt. 2.57e00) go to 170
      gas1  = 3.1899562e04-4.2186664e04*z
      gas2  = (2.3055603e03-1.9897017e03*z)*y
      gas3  = (1.849998e04+4.2561816e02*y-2.6808696e03*z)*z*z
      gas4  = (-1.6195114e01+5.8640623e00*z-3.6172504e00*y)*y*y
      gas5  = -5.7594039e04+7.9328437e04*z
      gas6  = (-1.9275989e03+1.6730544e03*z)*y
      gas7  = (-3.6473008e04-3.6100732e02*y+5.597543e03*z)*z*z
      gas8  = (-7.920808e01+4.0542084e01*z+2.1495867e00*y)*y*y
      gas9  = exp(-5.733e01+2.088e01*y+2.592e01*z-9.793e00*y*z)
      go to 190
 170  if (z .gt. 2.75e00) go to 180
      gas1  = 7.0838087e05-7.5619919e05*z
      gas2  = (3.9503091e03-2.7381802e03*z)*y
      gas3  = (2.6888181e05+4.7728687e02*y-3.183816e04*z)*z*z
      gas4  = (-1.2532251e02+4.7734787e01*z-4.0148029e00*y)*y*y
      gas5  = -2.5216325e05+2.1727769e05*z
      gas6  = (9.2882383e03-7.780918e03*z)*y
      gas7  = (-5.6539297e04+1.6120212e03*y+3.9419248e03*z)*z*z
      gas8  = (1.8537296e02-7.1010757e01*z+1.1307096e00*y)*y*y
      gas9  = exp(-1.786e02+2.18e-01*y+6.714e01*z-4.739e-01*y*z)
      go to 190
 180  gas1  = 3.1855037e05-3.3041156e05*z
      gas2  = (2.2983352e04-1.6623461e04*z)*y
      gas3  = (1.13848e05+3.0098223e03*y-1.3020133e04*z)*z*z
      gas4  = (-1.8599039e02+6.9840683e01*z-7.7371645e00*y)*y*y
      f     = gas1+gas2+gas3+gas4
      go to 200
 190  f     = gas1+gas2+gas3+gas4+(gas5+gas6+gas7+gas8)/(1.e00+gas9)
 200  k     = 1.87915e-02*f
 1000 format(/20x,48hwarning! outside of validity range of curve fit    &
     &,/,20x,5hrho =,1pd15.8,5x,3he =,1pd15.8,/)
      return
      end subroutine ugas4


      SUBROUTINE UGAS2 (T,RHO,PR)
!
!     first edited:   May 03 0 2010
!***********************************************************************
!**********                                                 ************
!**********              thermal conductivity               ************
!**********            returned as a function of            ************
!**********          internal energy and density            ************
!**********                                                 ************
!***********************************************************************
!            routine by srinivasan et.al., edited by Stanley Y. Ling
!     *********************************************************
!
!     inputs:
!
!     T = TEMPERATURE, IN KELVIN
!     RHO = DENSITY, IN KG/M**3
!
!     output:
!
!     PR = PRANDTL NUMBER
!
!     *********************************************************
!
!***********************************************************************
!     update by Stanley: Mar 03, 2010
!     Purpose: define variables used in the subroutines
! ************************************************************************

      implicit none

      double precision :: T, RHO
      double precision, intent (out):: PR

      double precision :: Z,X,Y, RHO0, GAS1,GAS2,GAS3,GAS4 
      double precision :: GAS5,GAS6,GAS7,GAS8,GAS9

!     end update
!***********************************************************************

    DATA RHO0/1.243E0/
    X=T/1000.0E00
    Y=LOG10(RHO/RHO0)
    IF ((Y.LT.-5.0E00).OR.(Y.GT.1.0E00).OR.(X.GT.1.5E01)) &
    WRITE (6,1000) RHO,T
    IF (T.GT.500.E00) GO TO 10
    GAS1=7.16321E-01+1.1135E00*X
    GAS2=(5.58243E-06-7.16815E-05*X)*Y
    GAS3=(-7.72911E00+2.25827E-04*Y+1.44166E01*X)*X*X
    GAS4=(-1.47156E-07-2.28926E-07*X-2.88338E-08*Y)*Y*Y
    GAS5=-1.4099E-01-3.35055E-01*X
    GAS6=(-2.55975E-05+1.5853E-04*X)*Y
    GAS7=(6.09194E00-3.18345E-04*Y-1.32747E01*X)*X*X
    GAS8=(1.3742E-06-1.29479E-06*X+1.48302E-07*Y)*Y*Y
    GAS9=EXP(8.636-3.03E01*X)
    GO TO 100
10  IF (T.GT.2.0E03) GO TO 20
    GAS1=6.766E-01+5.33391E-02*X
    GAS2=(-2.01021E-02+4.04905E-03*X)*X*X
    Z=GAS1+GAS2
    GO TO 110
20  IF (T.GT.4.0E03) GO TO 30
    GAS1=5.35204E-01+1.64262E-01*X
    GAS2=(-6.72637E-02+3.42314E-02*X)*Y
    GAS3=(-3.88497E-02-3.16248E-03*Y+3.05280E-03*X)*X*X
    GAS4=(-7.81832E-03+1.84389E-03*X-3.46855E-04*Y)*Y*Y
    Z=GAS1+GAS2+GAS3+GAS4
    GO TO 110
30  IF (T.GT.6.5E03) GO TO 40
    GAS1=-2.39283E00+1.28399E00*X
    GAS2=(-7.675E-01+1.89502E-01*X)*Y
    GAS3=(-1.79581E-01-1.20286E-02*Y+8.30322E-03*X)*X*X
    GAS4=(-7.09301E-02+7.19471E-03*X-2.78371E-03*Y)*Y*Y
    GAS5=3.06018E00-1.20461E00*X
    GAS6=(6.77677E-01-1.43868E-01*X)*Y
    GAS7=(1.62407E-01+7.85746E-03*Y-7.39086E-03*X)*X*X
    GAS8=(4.14157E-02-8.3635E-04*X-3.70369E-04*Y)*Y*Y
    GAS9=EXP(-26.39E00+2.969E00*X-5.042E00*Y-0.112E00*X*Y)
    GO TO 100
40  IF (T.GE.9.40E03) GO TO 50
    GAS1=6.13473E00-1.54169E00*X
    GAS2=(1.08128E00-2.04154E-01*X)*Y
    GAS3=(1.43737E-01+9.91640E-03*Y-4.54467E-03*X)*X*X
    GAS4=(6.19987E-02-5.05808E-03*X+1.56791E-03*Y)*Y*Y
    GAS5=-5.44445E00+1.58459E00*X
    GAS6=(-1.10792E00+2.13203E-01*X)*Y
    GAS7=(-1.51000E-01-1.00257E-02*Y+4.72964E-03*X)*X*X
    GAS8=(-7.80793E-02+7.29918E-03*X-2.29357E-03*Y)*Y*Y
    GAS9=EXP(13.39E00-4.258E00*X+2.298E00*Y-1.233E00*X*Y)
    GO TO 100
50  IF (T.GT.11.5E03) GO TO 70
    IF (Y.GT.-2.5E00) GO TO 60
    GAS1=-3.24776E01+8.72772E00*X
    GAS2=(-2.58872E00+3.59002E-01*X)*Y
    GAS3=(-7.61542E-01+2.10923E-02*Y+2.67953E-02*X)*X*X
    GAS4=(-1.9321E-01+8.94685E-02*X+5.64303E-02*Y)*Y*Y
    GAS5=3.99935E01-9.68334E00*X
    GAS6=(6.78337E00-9.07345E-01*X)*Y
    GAS7=(7.87932E-01-3.99108E-03*Y-2.64764E-02*X)*X*X
    GAS8=(6.97742E-01-1.25709E-01*X-4.08833E-02*Y)*Y*Y
    GAS9=EXP(105.8E00-11.67E00*X+31.67E00*Y-3.33E00*X*Y)
    GO TO 100
60  GAS1=-2.80755E01+6.80406E00*X
    GAS2=(-2.63243E00+4.0185E-01*X)*Y
    GAS3=(-5.45283E-01-1.63614E-02*Y+1.45424E-02*X)*X*X
    GAS4=(2.12026E-02-3.62386E-03*X+3.50018E-03*Y)*Y*Y
    GAS5=2.82604E01-6.62279E00*X
    GAS6=(2.06694E00-2.89135E-01*X)*Y
    GAS7=(5.26582E-01+1.13732E-02*Y-1.40944E-02*X)*X*X
    GAS8=(-1.31445E-01+1.50468E-02*X-9.13033E-03*Y)*Y*Y
    GAS9=EXP(-35.41E00+2.148E00*X-1.481E00*Y-0.3704E00*X*Y)
    GO TO 100
70  IF (Y.GT.-2.5E00) GO TO 90
    IF (T.GT.13.5E03) GO TO 80
    GAS1=6.08811E01-9.88231E00*X
    GAS2=(9.51872E00-9.95583E-01*X)*Y
    GAS3=(5.48699E-01+2.67619E-02*Y-1.03794E-02*X)*X*X
    GAS4=(5.0199E-01-2.55257E-02*X+8.45834E-03*Y)*Y*Y
    GAS5=-7.1667E01+1.2239E01*X
    GAS6=(-1.09959E01+1.23865E00*X)*Y
    GAS7=(-7.0869E-01-3.61494E-02*Y+1.38445E-02*X)*X*X
    GAS8=(-3.64368E-01+1.89712E-02*X+4.04684E-03*Y)*Y*Y
    GAS9=EXP(-71.89E00+2.248E00*X+0.4746E00*Y-0.9469E00*X*Y)
    GO TO 100
80  GAS1=2.99485E01-4.12112E00*X
    GAS2=(5.92879E00-5.27093E-01*X)*Y
    GAS3=(1.93397E-01+1.19371E-02*Y-3.08939E-03*X)*X*X
    GAS4=(4.09472E-01-1.78772E-02*X+9.49505E-03*Y)*Y*Y
    GAS5=-2.66557E01+3.05342E00*X
    GAS6=(-9.53775E00+8.98359E-01*X)*Y
    GAS7=(-9.53141E-02-1.98247E-02*Y+4.69853E-04*X)*X*X
    GAS8=(-7.62232E-01+4.34126E-02*X-5.10053E-03*Y)*Y*Y
    GAS9=EXP(-540.2E00+34.3E00*X-146.4E00*Y+9.148E00*X*Y)
    GO TO 100
90  GAS1=-3.18666E00+8.08818E-01*X
    GAS2=(-4.00164E-01+3.59959E-02*X)*Y
    GAS3=(-6.06519E-02-1.04205E-03*Y+1.45243E-03*X)*X*X
    GAS4=(1.6658E-02-4.36487E-03*X-1.86593E-03*Y)*Y*Y
    GAS5=2.68501E00-4.32123E-01*X
    GAS6=(1.36103E-01+2.5886E-02*X)*Y
    GAS7=(2.32842E-02-1.82391E-03*Y-4.09433E-04*X)*X*X
    GAS8=(-5.37705E-02+1.0074E-02*X+2.05852E-03*Y)*Y*Y
    GAS9=EXP(-31.16E00+1.633E00*X+2.395E00*Y-0.8707E00*X*Y)
100 Z=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/(1.0E00+GAS9)
110 PR=Z
1000 FORMAT(/20X,48HWARNING! OUTSIDE OF VALIDITY RANGE OF CURVE FIT    &
    &,/,20X,5HRHO =,1PE15.8,5X,3HT =,1PE15.8,/)
    RETURN
    END SUBROUTINE UGAS2

! ST-JWL----------------------

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModRealGas.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
!
!******************************************************************************
