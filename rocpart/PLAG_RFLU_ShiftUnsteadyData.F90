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
! $Id: PLAG_RFLU_ShiftUnsteadyData.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_ShiftUnsteadyData(pRegion)

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
  INTEGER :: icg,iPcl,iT,nPcls,nUnsteadyData
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

  RCSIdentString = '$RCSfile: PLAG_RFLU_ShiftUnsteadyData.F90,v $'

  global => pRegion%global
  
  CALL RegisterFunction( global, 'PLAG_RFLU_ShiftUnsteadyData',__FILE__ )

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pMixt => pRegion%mixt 
  pGrid => pRegion%grid
  pPlag => pRegion%plag

  pi = global%pi

  nPcls = pRegion%plag%nPcls

  nUnsteadyData = pRegion%plagInput%nUnsteadyData

! ******************************************************************************        
! Correct discrete particle dv
! ******************************************************************************

  DO iPcl = 1, nPcls
    DO iT = nUnsteadyData,2,-1
      pPlag%dudtMixt(XCOORD,iT,iPcl) = pPlag%dudtMixt(XCOORD,iT-1,iPcl)
      pPlag%dudtMixt(YCOORD,iT,iPcl) = pPlag%dudtMixt(YCOORD,iT-1,iPcl)
      pPlag%dudtMixt(ZCOORD,iT,iPcl) = pPlag%dudtMixt(ZCOORD,iT-1,iPcl)

      pPlag%dudtPlag(XCOORD,iT,iPcl) = pPlag%dudtPlag(XCOORD,iT-1,iPcl)
      pPlag%dudtPlag(YCOORD,iT,iPcl) = pPlag%dudtPlag(YCOORD,iT-1,iPcl)
      pPlag%dudtPlag(ZCOORD,iT,iPcl) = pPlag%dudtPlag(ZCOORD,iT-1,iPcl)
    END DO ! iT
  END DO ! iPcl

  IF (pRegion%plagInput%nTimeBH < nUnsteadyData) THEN
    pRegion%plagInput%nTimeBH = pRegion%plagInput%nTimeBH+1
  END IF ! pRegion%plagInput%nTimeBH

  DO iT = pRegion%plagInput%nTimeBH,2,-1
    pPlag%timeBH(iT) = pPlag%timeBH(iT-1) + global%dtMin
  END DO ! iT
 
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLU_ShiftUnsteadyData

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_ShiftUnsteadyData.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
!******************************************************************************

