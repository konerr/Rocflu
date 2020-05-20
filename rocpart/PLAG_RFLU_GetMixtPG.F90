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
! $Id: PLAG_RFLU_GetMixtPG.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_GetMixtPG(pRegion)

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
  REAL(RFREAL) :: term 
  REAL(RFREAL), DIMENSION(:), POINTER :: vol
  REAL(RFREAL), DIMENSION(:,:), POINTER :: sd,rhs,cv
  TYPE(t_global), POINTER :: global
  TYPE(t_plag),   POINTER :: pPlag 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_GetMixtPG.F90,v $'

  global => pRegion%global
  
  CALL RegisterFunction( global, 'PLAG_RFLU_GetMixtPG',__FILE__ )

! ******************************************************************************
! Checks: Defensive coding, should never occur
! ******************************************************************************

  IF ( pRegion%mixtInput%moveGrid .EQV. .TRUE. ) THEN
    CALL ErrorStop( global,ERR_UNKNOWN_OPTION,__LINE__, &
      'Equilibrium Eulerian method not yet implemented with moving grid')
  END IF ! region%mixtInput%moveGrid

  IF ( pRegion%mixtInput%indSd /= 1 ) THEN
    CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__)
  END IF ! region%mixtInput%indSd

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pPlag => pRegion%plag

  vol => pRegion%grid%vol

  sd  => pRegion%mixt%sd
  rhs => pRegion%mixt%rhs
  cv  => pRegion%mixt%cv

  nPcls = pRegion%plag%nPcls

! ******************************************************************************        
! Get mixt stress gradient (pressure gradient + viscous gradient) at particle
! ******************************************************************************
    
  DO iPcl = 1, nPcls
    icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcl)

    term = 1.0_RFREAL/(vol(icg))

    pPlag%pgMixt(XCOORD,iPcl) = -term*(sd(SD_XMOM,icg) - rhs(CV_MIXT_XMOM,icg))
    pPlag%pgMixt(YCOORD,iPcl) = -term*(sd(SD_YMOM,icg) - rhs(CV_MIXT_YMOM,icg))
    pPlag%pgMixt(ZCOORD,iPcl) = -term*(sd(SD_ZMOM,icg) - rhs(CV_MIXT_ZMOM,icg))
  END DO ! iPcl
      
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLU_GetMixtPG

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_GetMixtPG.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
!******************************************************************************

