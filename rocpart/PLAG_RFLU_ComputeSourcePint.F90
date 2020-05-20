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
! Purpose: Add source terms to the residual.
!
! Description: None.
!
! Input: 
!   region      data of current region.
!
! Output: None.
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: PLAG_RFLU_ComputeSourcePint.F90,v 1.4 2016/05/16 21:01:57 rahul Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_ComputeSourcePint(region)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region,t_grid,t_mixt,t_mixt_input
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError

  USE ModPartLag, ONLY: t_plag, t_plag_input
  USE PLAG_ModParameters
  USE ModInterfaces, ONLY: MixtPerf_P_DEoGVm2

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region) :: region
  TYPE(t_grid)   :: grid
  TYPE(t_mixt)   :: mixt
  TYPE(t_mixt_input) :: mixtInput
  TYPE(t_plag), POINTER   :: pPlag

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: ic,errorflag
  REAL(RFREAL) :: rhoVol
                  
  REAL(RFREAL), POINTER :: cv(:,:),rhs(:,:),vol(:)
  TYPE(t_global), POINTER :: global

  INTEGER :: cID,ipcl,indx
  REAL(RFREAL) :: g,r,cpup,eps,up,vp,wp,ug,vg,wg,Qg,Qp,p,term1,&
                  term2,delp,pint,vFracG,vFracGOld,pcl
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv
  REAL(RFREAL), POINTER :: Cvp(:,:),Dvp(:,:)
  REAL(RFREAL), ALLOCATABLE :: upE(:)
  INTEGER, ALLOCATABLE :: npcl(:)

! ******************************************************************************
! Start
! ******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'PLAG_RFLU_ComputeSourcePint',__FILE__ )

! ******************************************************************************
! Get dimensions and pointers
! ******************************************************************************

  cv  => region%mixt%cv
  rhs => region%mixt%rhs
  vol => region%grid%vol 

  pAiv=> region%plag%aiv
  Dvp => region%plag%dv  ! rahul - particle derived quantities
  Cvp => region%plag%cv  ! rahul - particle conserved quantities

  ALLOCATE(upE(region%grid%nCells), STAT=errorflag)
  IF (errorflag /= 0) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'upE')
  END IF

  ALLOCATE(npcl(region%grid%nCells), STAT=errorflag)
  IF (errorflag /= 0) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'npcl')
  END IF

! ******************************************************************************
! Non-conservative source terms, ref. Chang et al. 2007, ASME 
! ******************************************************************************
  g  = global%refGamma
 
  DO ic = 1,region%grid%nCells
   upE(ic)  = 0.0_RFREAL
   npcl(ic) = 0.0_RFREAL
  END DO

  Qp   = 0.0_RFREAL
  cpup = 2.0_RFREAL
  eps  = 0.01_RFREAL

! *****************************************************************************
! Compute particle velocity in a cell. This is the particle number averaged
! velocity per cell
! *****************************************************************************
!write(*,*) 'in sourcePint'

  DO ipcl = 1,region%plag%npcls
    cID = pAiv(AIV_PLAG_ICELLS,ipcl)  
    up = Dvp(DV_PLAG_UVEL,ipcl)
    vp = Dvp(DV_PLAG_VVEL,ipcl)
    wp = Dvp(DV_PLAG_WVEL,ipcl)
    Qp = Qp + (up*up + vp*vp + wp*wp)**0.5_RFREAL

    upE(cID)  = upE(cID) + Qp
    npcl(cID) = npcl(cID) + 1 
  END DO

! *****************************************************************************
! End particle velocity computation
! *****************************************************************************

  DO ic = 1,region%grid%nCells
   IF (npcl(ic) .NE. 0) THEN
     vFracG    = 1.0_RFREAL - region%plag%vFracE(1,ic)   ! Volume fraction of gas
     vFracGOld = 1.0_RFREAL - region%plag%vFracEOld(1,ic)

     upE(ic) = upE(ic)/npcl(ic)
 
     r      = cv(CV_MIXT_DENS,ic) ! phig*rhog
     rhovol = vol(ic)*r

     ug = cv(CV_MIXT_XMOM,ic)/r 
     vg = cv(CV_MIXT_YMOM,ic)/r 
     wg = cv(CV_MIXT_ZMOM,ic)/r  

     Qg = (ug*ug + vg*vg + wg*wg)**0.5_RFREAL
        
!     p  = MixtPerf_P_DEoGVm2(r/vFracG,cv(CV_MIXT_ENER,ic)/r,g,Qg*Qg)
     p = region%mixt%dv(DV_MIXT_PRES,ic)
     term1 = cpup*(1.0_RFREAL-vFracG)*r/vFracG
     term2 = (Qg-upE(ic))**2.0_RFREAL

     delp  = term1*term2
     delp  = MIN(delp,eps*p)

     pint = p - delp

!write(007,*) ,ic,npcl(ic),upE(ic),Qg,delp,pint
IF(1==2) THEN
!! gradVFracEg is the Eulerian cell-centered gradient of gas volume fraction
     rhs(CV_MIXT_XMOM,ic) = rhs(CV_MIXT_XMOM,ic) - &
                            pint*region%plag%gradVFracEg(XCOORD,1,ic)*vol(ic)  

     rhs(CV_MIXT_YMOM,ic) = rhs(CV_MIXT_YMOM,ic) - &
                            pint*region%plag%gradVFracEg(YCOORD,1,ic)*vol(ic)  

     rhs(CV_MIXT_ZMOM,ic) = rhs(CV_MIXT_ZMOM,ic) - &
                            pint*region%plag%gradVFracEg(ZCOORD,1,ic)*vol(ic)  

     rhs(CV_MIXT_ENER,ic) = rhs(CV_MIXT_ENER,ic) + &
                            delp*( ug*region%plag%gradVFracEg(XCOORD,1,ic) + &
                                   vg*region%plag%gradVFracEg(YCOORD,1,ic) + &
                                   wg*region%plag%gradVFracEg(ZCOORD,1,ic) )*vol(ic) + &
                            pint*(vFracG - vFracGOld)*vol(ic)/global%dtMin
ELSE ! 1==2

     rhs(CV_MIXT_XMOM,ic) = rhs(CV_MIXT_XMOM,ic) + &
                            pint*region%plag%gradvFracE(XCOORD,1,ic)*vol(ic)  

     rhs(CV_MIXT_YMOM,ic) = rhs(CV_MIXT_YMOM,ic) + &
                            pint*region%plag%gradvFracE(YCOORD,1,ic)*vol(ic)  

     rhs(CV_MIXT_ZMOM,ic) = rhs(CV_MIXT_ZMOM,ic) + &
                            pint*region%plag%gradvFracE(ZCOORD,1,ic)*vol(ic)  

     rhs(CV_MIXT_ENER,ic) = rhs(CV_MIXT_ENER,ic) - &
                            delp*( ug*region%plag%gradvFracE(XCOORD,1,ic) + &
                                   vg*region%plag%gradvFracE(YCOORD,1,ic) + &
                                   wg*region%plag%gradvFracE(ZCOORD,1,ic) )*vol(ic) - &
                            !pint*(vFracGOld - vFracG)
                            pint*(region%plag%vFracE(1,ic) - &
                                  region%plag%vFracEOld(1,ic))*vol(ic)/global%dtMin
END IF
   END IF ! upE .NE. 0
  END DO ! ic

 DEALLOCATE(upE, STAT=errorflag)
 IF (errorflag /= 0) THEN
   CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'upE')
 END IF
 
 DEALLOCATE(npcl, STAT=errorflag)
 IF (errorflag /= 0) THEN
   CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'npcl')
 END IF


! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_ComputeSourcePint

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_ComputeSourcePint.F90,v $
! Revision 1.4  2016/05/16 21:01:57  rahul
! Minor fix in energy equation.
!
! Revision 1.3  2016/05/16 15:22:13  rahul
! Bug fix. Correct formulation for time integration of energy equation.
!
! Revision 1.2  2016/05/06 00:43:08  rahul
! Cosmetic changes.
!
! Revision 1.1  2016/05/06 00:39:39  rahul
! 1. Compute non-conservative terms in fluid phase momentum and energy for
!    E-L AUSM+up scheme.
! 2. Major bug fix: Fixed a bug in non-conservative term in energy rhs.
!
! *****************************************************************************

