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
! Purpose: Correct the mixture cv varaibles due to particle volume fraction.
!
! Description: none.
!
! Input: pRegion = current region.
!
! Output: region%levels%mixt%cv 
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_RFLU_ModifyFlowFieldVolFrac.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_ModifyFlowFieldVolFrac(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModMixture, ONLY    : t_mixt
  USE ModError  
  USE ModParameters
  USE PLAG_ModParameters

  USE ModInterfaces, ONLY: MixtPerf_Eo_DGPUVW, &
                           MixtPerf_R_M, & 
                           MixtPerf_G_CpR, &
                           MixtPerf_P_DEoGVm2

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: icg,indCp,indMol,iPcl,nPcls
  REAL(RFREAL) :: e,Eo,ir,u,v,w,vFrac,mw,cp,gc,g,Vm2,p,r
  REAL(RFREAL), POINTER, DIMENSION(:,:)    :: pCv,pGv
  
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag),   POINTER :: pPlag 
  TYPE(t_mixt),   POINTER :: pMixt 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_ModifyFlowFieldVolFrac.F90,v $'

  global => pRegion%global
  
  CALL RegisterFunction( global, 'PLAG_RFLU_ModifyFlowFieldVolFrac',__FILE__ )

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pMixt => pRegion%mixt 
  pGrid => pRegion%grid
  pPlag => pRegion%plag

  pCv => pRegion%mixt%cv
!  pGv => pRegion%mixt%gv

!  indCp  = pRegion%mixtInput%indCp
!  indMol = pRegion%mixtInput%indMol

! ******************************************************************************        
! Correct mixture cv due to particle volume fraction
! ******************************************************************************

  DO icg = 1,pGrid%nCellsTot    
!    mw = pGv(GV_MIXT_MOL,indMol*icg)
!    cp = pGv(GV_MIXT_CP ,indCp *icg)
!        
!    gc = MixtPerf_R_M(mw)
!    g  = MixtPerf_G_CpR(cp,gc)  
!
!    r  = pCv(CV_MIXT_DENS,icg)
!    ir = 1.0_RFREAL/r
!    u  = ir*pCv(CV_MIXT_XMOM,icg)
!    v  = ir*pCv(CV_MIXT_YMOM,icg)
!    w  = ir*pCv(CV_MIXT_ZMOM,icg)
!    Eo = ir*pCv(CV_MIXT_ENER,icg)
!
!    Vm2 = ir*ir*(pCv(CV_MIXT_XMOM,icg)*pCv(CV_MIXT_XMOM,icg) + &
!                 pCv(CV_MIXT_YMOM,icg)*pCv(CV_MIXT_YMOM,icg) + &
!                 pCv(CV_MIXT_ZMOM,icg)*pCv(CV_MIXT_ZMOM,icg))
!
!    p = MixtPerf_P_DEoGVm2(r,Eo,g,Vm2)

    vFrac = 1.0_RFREAL - pPlag%vFracE(1,icg)

!    pCv(CV_MIXT_DENS,icg) = vFrac*r
!    pCv(CV_MIXT_XMOM,icg) = vFrac*r*u
!    pCv(CV_MIXT_YMOM,icg) = vFrac*r*v
!    pCv(CV_MIXT_ZMOM,icg) = vFrac*r*w
!    pCv(CV_MIXT_ENER,icg) = vFrac*r*MixtPerf_Eo_DGPUVW(r,g,p,u,v,w)
    
    pCv(CV_MIXT_DENS,icg) = vFrac*pCv(CV_MIXT_DENS,icg)
    pCv(CV_MIXT_XMOM,icg) = vFrac*pCv(CV_MIXT_XMOM,icg)
    pCv(CV_MIXT_YMOM,icg) = vFrac*pCv(CV_MIXT_YMOM,icg)
    pCv(CV_MIXT_ZMOM,icg) = vFrac*pCv(CV_MIXT_ZMOM,icg)
    pCv(CV_MIXT_ENER,icg) = vFrac*pCv(CV_MIXT_ENER,icg)

  END DO ! icg
      
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLU_ModifyFlowFieldVolFrac

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_ModifyFlowFieldVolFrac.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
!******************************************************************************

