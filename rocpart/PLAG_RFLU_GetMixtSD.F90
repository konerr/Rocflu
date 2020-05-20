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
! Purpose: Correct the mixture properties using higher-order
!          interpolation schemed onto particle locations
!
! Description: none.
!
! Input: pRegion = current region.
!
! Output: region%levels%plag%dv 
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_RFLU_GetMixtSD.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_GetMixtSD(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModMixture, ONLY    : t_mixt
  USE ModError  
  USE ModParameters
  USE PLAG_ModParameters

  USE ModInterfaces, ONLY: MixtPerf_R_M, & 
                           MixtPerf_T_DPR

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
  INTEGER :: icg,iPcl,nPcls
  INTEGER, POINTER, DIMENSION(:,:)    :: pAiv
  
  REAL(RFREAL) :: dx,dy,dz
  REAL(RFREAL), POINTER, DIMENSION(:,:)    :: pCv,pEv,pPlagRhs
  REAL(RFREAL), POINTER, DIMENSION(:,:,:) :: pGrad
  
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag),   POINTER :: pPlag 
  TYPE(t_mixt),   POINTER :: pMixt 

! TEMPORARY: Manoj: Adding interparticle force  
  REAL(RFREAL) :: beta,diamL,factor,phi,pi,phiCP,Ps,volL
! END TEMPORARY

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_GetMixtSD.F90,v $'

  global => pRegion%global
  
  CALL RegisterFunction( global, 'PLAG_RFLU_GetMixtSD',__FILE__ )

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pMixt => pRegion%mixt 
  pGrid => pRegion%grid
  pPlag => pRegion%plag

  pi = global%pi

  nPcls = pRegion%plag%nPcls

! ******************************************************************************        
! Correct discrete particle dv
! ******************************************************************************
    
  DO iPcl = 1, nPcls
    icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcl)

    pPlag%dudtMixt(XCOORD,1,iPcl) = pRegion%mixt%sd(SD_XMOM,icg)
    pPlag%dudtMixt(YCOORD,1,iPcl) = pRegion%mixt%sd(SD_YMOM,icg)
    pPlag%dudtMixt(ZCOORD,1,iPcl) = pRegion%mixt%sd(SD_ZMOM,icg)
! TEMPORARY: Manoj        
!    pPlag%dudtMixt(XCOORD,1,iPcl) = -pPlag%pgMixt(XCOORD,iPcl)/pMixt%cv(CV_MIXT_DENS,icg)
!    pPlag%dudtMixt(YCOORD,1,iPcl) = -pPlag%pgMixt(YCOORD,iPcl)/pMixt%cv(CV_MIXT_DENS,icg)
!    pPlag%dudtMixt(ZCOORD,1,iPcl) = -pPlag%pgMixt(ZCOORD,iPcl)/pMixt%cv(CV_MIXT_DENS,icg)
! END TEMPORARY: Manoj    
  END DO ! iPcl
! TEMPORARY: Manoj        
!WRITE(*,*) "Setting DuDt = -dpdx/rho -----------------------"
! END TEMPORARY: Manoj    
      
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLU_GetMixtSD

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_GetMixtSD.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
!******************************************************************************

