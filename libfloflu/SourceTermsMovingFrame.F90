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
! *****************************************************************************
!
! Purpose: Add source terms due to moving frame to the residual.
!
! Description: None.
!
! Input: 
!   region                 data of current region.
!
! Output: 
!   region%levels%mixt%rhs complete right-hand side (residual).
!
! Notes: None.
!
! *****************************************************************************
!
! $Id: SourceTermsMovingFrame.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! *****************************************************************************

SUBROUTINE SourceTermsMovingFrame(region)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModParameters

  USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, &
                               RFLU_ConvertCvPrim2Cons

  USE ModInterfaces, ONLY: MixtPerf_G_CpR, MixtPerf_R_M

  IMPLICIT NONE

! *****************************************************************************
! Arguments
! *****************************************************************************

  TYPE(t_region), TARGET :: region

! *****************************************************************************
! Locals
! *****************************************************************************

  INTEGER :: errorFlag,icg
  REAL(RFREAL) :: r,sx,sy,sz,u,v,vol,w
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pRhs
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_region), POINTER :: pRegion

! *****************************************************************************
! Start, set pointers
! *****************************************************************************

  pRegion => region
  global  => pRegion%global

  CALL RegisterFunction(global,'SourceTermsMovingFrame',__FILE__)

  pGrid => pRegion%grid
  pCv   => pRegion%mixt%cv
  pRhs  => pRegion%mixt%rhs

! *****************************************************************************
! Convert conservative variables to primary variables 
! *****************************************************************************

  CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)

! *****************************************************************************
! Loop over cells and compute source terms for moving reference frame
! *****************************************************************************

  DO icg = 1,region%grid%nCellsTot
    r = pCv(CV_MIXT_DENS,icg)
    u = pCv(CV_MIXT_XMOM,icg)
    v = pCv(CV_MIXT_YMOM,icg)
    w = pCv(CV_MIXT_ZMOM,icg)
    
    vol = pGrid%vol(icg)

    sx = -r*pRegion%mvfAcc(XCOORD)*vol
    sy = -r*pRegion%mvfAcc(YCOORD)*vol
    sz = -r*pRegion%mvfAcc(ZCOORD)*vol

    pRhs(CV_MIXT_XMOM,icg) = pRhs(CV_MIXT_XMOM,icg) - sx
    pRhs(CV_MIXT_YMOM,icg) = pRhs(CV_MIXT_YMOM,icg) - sy
    pRhs(CV_MIXT_ZMOM,icg) = pRhs(CV_MIXT_ZMOM,icg) - sz
    pRhs(CV_MIXT_ENER,icg) = pRhs(CV_MIXT_ENER,icg) - (sx*u + sy*v + sz*w)
  END DO ! ic

! *****************************************************************************
! Convert back primary variables to conservative variables 
! *****************************************************************************

  CALL RFLU_ConvertCvPrim2Cons(pRegion,CV_MIXT_STATE_CONS)

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE SourceTermsMovingFrame

! *****************************************************************************
!
! RCS Revision history:
!
! $Log: SourceTermsMovingFrame.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/06/18 17:30:25  mparmar
! Initial revision
!
! *****************************************************************************

