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
! Purpose: compute dependent variables for a thermally and calorically
!          perfect gas.
!
! Description: none.
!
! Input: inBeg  = first value to update
!        inEnd  = last value to update
!        indCp  = indicates if cp varies over cells (=1) or is constant (=0)
!        indMol = indicates if the mol mass varies over cells (=1) or is
!                 constant (=0)
!        cv     = conservative variables
!        gv     = gas variables (cp, Mol)
!
! Output: dv = dependent variables (p, T, c and u, v, w for RocfloMP)
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PerfgasDependentVars.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PerfgasDependentVars( inBeg,inEnd,indCp,indMol,cv,gv,dv )

  USE ModDataTypes
  USE ModParameters

  USE ModInterfaces, ONLY: MixtPerf_C_GRT,MixtPerf_G_CpR, & 
                           MixtPerf_P_DEoGVm2,MixtPerf_R_M,MixtPerf_T_DPR

  IMPLICIT NONE

! ... parameters
  INTEGER :: inBeg, inEnd, indCp, indMol

  REAL(RFREAL), POINTER :: cv(:,:), gv(:,:), dv(:,:)

! ... loop variables
  INTEGER :: ic

! ... local variables
  REAL(RFREAL) :: rgas,rrho,Vm2
  REAL(RFREAL) :: Eo,gamma,rho

!******************************************************************************

  DO ic=inBeg,inEnd
    rgas  = MixtPerf_R_M(gv(GV_MIXT_MOL,ic*indMol))
    gamma = MixtPerf_G_CpR(gv(GV_MIXT_CP,ic*indCp),rgas)

    rho  = cv(CV_MIXT_DENS,ic)
    rrho = 1.0_RFREAL/rho
    Eo   = cv(CV_MIXT_ENER,ic)*rrho

    Vm2 = (cv(CV_MIXT_XMOM,ic)*cv(CV_MIXT_XMOM,ic) + &
           cv(CV_MIXT_YMOM,ic)*cv(CV_MIXT_YMOM,ic) + &
           cv(CV_MIXT_ZMOM,ic)*cv(CV_MIXT_ZMOM,ic))*rrho*rrho

    dv(DV_MIXT_PRES,ic) = MixtPerf_P_DEoGVm2(rho,Eo,gamma,Vm2)
    dv(DV_MIXT_TEMP,ic) = MixtPerf_T_DPR(rho,dv(DV_MIXT_PRES,ic),rgas)
    dv(DV_MIXT_SOUN,ic) = MixtPerf_C_GRT(gamma,rgas,dv(DV_MIXT_TEMP,ic))
  ENDDO

END SUBROUTINE PerfgasDependentVars

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PerfgasDependentVars.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 16:49:56  haselbac
! Initial revision after changing case
!
! Revision 1.12  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.9  2003/05/13 23:46:52  haselbac
! Reverted to use of MixtPerf routines for RFLU
!
! Revision 1.8  2003/05/06 20:05:39  jblazek
! Corrected bug in grid motion (corner "averaging").
!
! Revision 1.7  2002/07/05 23:20:46  jblazek
! Corrected bug in perfgasDependentVars.F90; did some cosmetics.
!
! Revision 1.6  2002/06/05 18:33:00  haselbac
! Converted to use of mixtPerf routines
!
! Revision 1.5  2002/03/18 22:25:45  jblazek
! Finished multiblock and MPI.
!
! Revision 1.4  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.3  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.2  2002/01/23 03:51:24  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.1  2002/01/10 00:02:06  jblazek
! Added calculation of mixture properties.
!
!******************************************************************************

