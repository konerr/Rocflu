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
! Purpose: Compute particle volume fraction for Eulerian and Lagrangian fields.
!
! Description: none.
!
! Input: pRegion = current region.
!
! Output: plag%vFracE and plag%vFracL
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_RFLU_ComputeMaxImpulse.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_ComputeMaxImpulse(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModMixture, ONLY    : t_mixt
  USE ModError  
  USE ModParameters
  USE PLAG_ModParameters

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
  INTEGER :: errorFlag,icg,iCont,iPcl,iVar,nCont,nPcls 
  INTEGER, POINTER, DIMENSION(:) :: cvMass    

  REAL(RFREAL) :: ir,massE,massL,uE,uL,vE,vL,wE,wL 
                    
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_ComputeMaxImpulse.F90,v $'

  global => pRegion%global
  
  CALL RegisterFunction( global, 'PLAG_RFLU_ComputeMaxImpulse',__FILE__ )

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

  nCont = pRegion%plagInput%nCont

  cvMass => pPlag%cvPlagMass

! ******************************************************************************
! Initialize Eulerian field
! ******************************************************************************

  DO icg = 1,pGrid%nCellsTot
    DO iVar=1,3
      pPlag%forceTotal(iVar,icg) = 0.0_RFREAL
      pPlag%dImpulseMax(iVar,icg) = 0.0_RFREAL
    END DO ! iVar
  END DO ! icg

! ==============================================================================
!   Compute volume fraction if there are particles in the region
! ==============================================================================

  IF ( pPlag%nPcls > 0 ) THEN

! ==============================================================================
!   Compute Eulerian field for mass and momentum transfer of particles
!   Use forceTotal(1,:) to store total mass of all particles in that cell
!   Use dImpulseMax(i,:) to store mp*(uf-upi)
! ==============================================================================

    DO iPcl = 1,pPlag%nPcls
      icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcl)

      ir = 1.0_RFREAL/pRegion%mixt%cv(CV_MIXT_DENS,icg)
      uE = ir*pRegion%mixt%cv(CV_MIXT_XMOM,icg)
      vE = ir*pRegion%mixt%cv(CV_MIXT_YMOM,icg)
      wE = ir*pRegion%mixt%cv(CV_MIXT_ZMOM,icg)

      massE = pRegion%grid%vol(icg)*pRegion%mixt%cv(CV_MIXT_DENS,icg)

      uL = pPlag%dv(DV_PLAG_UVEL,iPcl)
      vL = pPlag%dv(DV_PLAG_VVEL,iPcl)
      wL = pPlag%dv(DV_PLAG_WVEL,iPcl)

      massL = SUM(pPlag%cv(cvMass(:),iPcl))

      pPlag%forceTotal(1,icg) = pPlag%forceTotal(1,icg) + massL

      pPlag%dImpulseMax(XCOORD,icg) = pPlag%dImpulseMax(XCOORD,icg)+massL*ABS(uE-uL)
      pPlag%dImpulseMax(YCOORD,icg) = pPlag%dImpulseMax(YCOORD,icg)+massL*ABS(vE-vL)
      pPlag%dImpulseMax(ZCOORD,icg) = pPlag%dImpulseMax(ZCOORD,icg)+massL*ABS(wE-wL)
    END DO ! iPcl

! ==============================================================================
!   Finalize computation of max Impulse transfer possible
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot
      massE = pRegion%grid%vol(icg)*pRegion%mixt%cv(CV_MIXT_DENS,icg)

      pPlag%dImpulseMax(XCOORD,icg) = massE*pPlag%dImpulseMax(XCOORD,icg) &
                                     /(massE+pPlag%forceTotal(1,icg))
      pPlag%dImpulseMax(YCOORD,icg) = massE*pPlag%dImpulseMax(YCOORD,icg) &
                                     /(massE+pPlag%forceTotal(1,icg))
      pPlag%dImpulseMax(ZCOORD,icg) = massE*pPlag%dImpulseMax(ZCOORD,icg) &
                                     /(massE+pPlag%forceTotal(1,icg))
!IF (1==2) THEN
!  IF (ABS(pPlag%forceTotal(1,icg)) > 0.0_RFREAL) THEN
!    WRITE(*,'(2(A,1X,I7,1X),3(E12.6,2X))') &
!                                   "iReg=",pRegion%iRegionGlobal," icg=",icg, &
!                                   pPlag%dImpulseMax(XCOORD,icg), &
!                                   pPlag%dImpulseMax(YCOORD,icg), &
!                                   pPlag%dImpulseMax(ZCOORD,icg)
!  END IF
!END IF ! 1==1
    END DO ! icg  
  END IF ! pPlag%nPcls

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLU_ComputeMaxImpulse

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_ComputeMaxImpulse.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
!
!******************************************************************************

