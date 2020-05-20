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
! Purpose: set outflow boundary condition for one cell.
!
! Description: the subsonic boundary condition is based on non-reflecting,
!              characteristics method of Whitfield and Janus: Three-Dimensional
!              Unsteady Euler Equations Solution Using Flux Vector Splitting.
!              AIAA Paper 84-1552, 1984. The supersonic boundary condition
!              consists of simple extrapolation.
!
! Input: bcOpt    = boundary treatment: subsonic, supersonic, or mixed
!        pout     = given static outlet pressure
!        sx/y/zn  = components of ortho-normalized face vector (outward facing)
!        cpgas    = specific heat at constant pressure (boundary cell)
!        mol      = molecular mass at boundary cell
!        rho      = density at boundary cell
!        rhou/v/w = density * velocity components at boundary cell
!        rhoe     = density * total energy at boundary cell
!        press    = static pressure at boundary cell
!
! Output: rhob      = density at boundary
!         rhou/v/wb = density * velocity components at boundary
!         rhoeb     = density * total energy at boundary
!
! Notes: this condition is valid only for thermally and calorically
!        perfect gas (supersonic outflow valid for all gases).
!
!******************************************************************************
!
! $Id: BcondOutflowPerf.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE BcondOutflowPerf( bcOpt,pout,sxn,syn,szn,cpgas,mol, & 
                             rho,rhou,rhov,rhow,rhoe,press, &
                             rhob,rhoub,rhovb,rhowb,rhoeb )

  USE ModDataTypes
  USE ModParameters
  USE ModInterfaces, ONLY: MixtPerf_C_DGP, MixtPerf_Eo_DGPUVW, & 
                           MixtPerf_G_CpR, MixtPerf_R_M, MixtPerf_P_DEoGVm2

  IMPLICIT NONE

! ... parameters
  INTEGER :: bcOpt

  REAL(RFREAL) :: pout
  REAL(RFREAL) :: rho, rhou, rhov, rhow, rhoe, press
  REAL(RFREAL) :: sxn, syn, szn, cpgas, mol
  REAL(RFREAL) :: rhob, rhoub, rhovb, rhowb, rhoeb

! ... local variables
  REAL(RFREAL) :: csound, rgas, gamma, gam1, u, v, w, mach, rrhoc, deltp, &
                  ub, vb, wb, vnd

!******************************************************************************
! gas properties; velocity components; Mach number

  rgas  = MixtPerf_R_M( mol )
  gamma = MixtPerf_G_CpR( cpgas,rgas )
  gam1  = gamma - 1.0_RFREAL

  u      = rhou/rho
  v      = rhov/rho
  w      = rhow/rho
  csound = MixtPerf_C_DGP( rho,gamma,press )
  mach   = SQRT(u*u+v*v+w*w)/csound

! subsonic flow ---------------------------------------------------------------

  IF (mach < 1.0_RFREAL .AND. &
      (bcOpt == BCOPT_SUBSONIC .OR. bcOpt == BCOPT_MIXED)) THEN
    rrhoc = 1.0_RFREAL/(rho*csound)
    deltp = press - pout
    rhob  = rho - deltp/(csound*csound)
    ub    = u + sxn*deltp*rrhoc
    vb    = v + syn*deltp*rrhoc
    wb    = w + szn*deltp*rrhoc

! - special treatment to prevent "deltp" from changing the sign
!   of velocity components. This may happen for very small u, v, w.

    vnd = ub*sxn + vb*syn + wb*szn
    IF ( vnd < 0.0_RFREAL ) THEN ! inflow at outflow boundary
      ub = SIGN(1.0_RFREAL,u)*MAX(ABS(ub),ABS(u))
      vb = SIGN(1.0_RFREAL,v)*MAX(ABS(vb),ABS(v))
      wb = SIGN(1.0_RFREAL,w)*MAX(ABS(wb),ABS(w))
    END IF ! vnd

    rhoub = rhob*ub
    rhovb = rhob*vb
    rhowb = rhob*wb
    rhoeb = rhob*MixtPerf_Eo_DGPUVW( rhob,gamma,pout,ub,vb,wb )

! supersonic flow -------------------------------------------------------------

  ELSE
    rhob  = rho
    rhoub = rhou
    rhovb = rhov
    rhowb = rhow
    rhoeb = rhoe
  END IF ! mach

END SUBROUTINE BcondOutflowPerf

!******************************************************************************
!
! RCS Revision history:
!
! $Log: BcondOutflowPerf.F90,v $
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
! Revision 1.2  2006/03/26 20:21:09  haselbac
! Fix mistake in declarations
!
! Revision 1.1  2004/12/01 16:48:04  haselbac
! Initial revision after changing case
!
! Revision 1.5  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.2  2002/06/22 00:49:50  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.1  2002/06/10 21:19:34  haselbac
! Initial revision
!
!******************************************************************************

