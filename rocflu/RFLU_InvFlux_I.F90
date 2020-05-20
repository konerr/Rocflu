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
! Purpose: Compute convective fluxes using centered scheme.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_InvFlux_I.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_InvFlux_I(pRegion)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters

  USE ModInterfaces, ONLY: RFLU_CentralFirstPatch, & 
                           RFLU_GetCvLoc

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: cvMixtXVel,cvMixtYVel,cvMixtZVel,c1,c2,ifg,iPatch
  REAL(RFREAL) :: term
  REAL(RFREAL) :: flx(3)
  REAL(RFREAL), DIMENSION(:), POINTER :: pVf
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pRhs
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_InvFlux_I.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_InvFlux_I',__FILE__)

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  pCv  => pRegion%mixt%cv
  pRhs => pRegion%mixt%rhs
  pVf  => pRegion%mixt%vfMixt    

  cvMixtXVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_XVEL)
  cvMixtYVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_YVEL)
  cvMixtZVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_ZVEL)

! ******************************************************************************
! Compute fluxes through interior faces
! ******************************************************************************

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

    term = 0.5_RFREAL*pVf(ifg)*pGrid%fn(XYZMAG,ifg)        

    flx(1) = pCv(cvMixtXVel,c2)*term
    flx(2) = pCv(cvMixtYVel,c2)*term
    flx(3) = pCv(cvMixtZVel,c2)*term

    pRhs(cvMixtXVel,c1) = pRhs(cvMixtXVel,c1) - flx(1)
    pRhs(cvMixtYVel,c1) = pRhs(cvMixtYVel,c1) - flx(2)
    pRhs(cvMixtZVel,c1) = pRhs(cvMixtZVel,c1) - flx(3)

    flx(1) = -pCv(cvMixtXVel,c1)*term
    flx(2) = -pCv(cvMixtYVel,c1)*term
    flx(3) = -pCv(cvMixtZVel,c1)*term

    pRhs(cvMixtXVel,c2) = pRhs(cvMixtXVel,c2) - flx(1)
    pRhs(cvMixtYVel,c2) = pRhs(cvMixtYVel,c2) - flx(2)
    pRhs(cvMixtZVel,c2) = pRhs(cvMixtZVel,c2) - flx(3)
  END DO ! ifg

! ******************************************************************************
! Compute fluxes through boundary faces
! ******************************************************************************

! TO DO 
!  DO iPatch = 1,region%grid%nPatches
!    CALL RFLU_CentralFirstPatch(region,region%patches(iPatch))
!  END DO ! iPatch
! END TO DO 

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InvFlux_I

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_InvFlux_I.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:48  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:57  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:01  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/19 15:41:11  haselbac
! Initial revision
!
! ******************************************************************************

