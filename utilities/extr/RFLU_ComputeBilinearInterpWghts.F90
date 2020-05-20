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
! Purpose: Compute weights for bilinear interpolation.
!
! Description: None.
!
! Input: 
!   pRegion	Pointer to region
!   xLoc	x-location of intersection
!   yLoc	y-location of intersection
!   zLoc	z-location of intersection
!   cs		Cell stencil
!
! Output: 
!   wghts	Interpolation weights
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ComputeBilinearInterpWghts.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ComputeBilinearInterpWghts(pRegion,xLoc,yLoc,zLoc,cs,wghts)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  USE ModSortSearch
      
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Parameters
! ==============================================================================

  INTEGER, DIMENSION(4), INTENT(IN) :: cs
  REAL(RFREAL), INTENT(IN) :: xLoc,yLoc,zLoc
  REAL(RFREAL), DIMENSION(4), INTENT(OUT) :: wghts
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: icl
  REAL(RFREAL) :: xLocNorm,xMax,xMin,yLocNorm,yMax,yMin
  REAL(RFREAL), DIMENSION(4) :: xc,yc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: RFLU_ComputeBilinearInterpWghts.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeBilinearInterpWghts', &
                        __FILE__)

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  DO icl = 1,4

    xc(icl)  = pGrid%cofg(XCOORD,cs(icl))
    yc(icl)  = pGrid%cofg(YCOORD,cs(icl))
  END DO ! icl

  xMin = MINVAL(xc(1:4))
  xMax = MAXVAL(xc(1:4))
  yMin = MINVAL(yc(1:4))
  yMax = MAXVAL(yc(1:4))

  xLocNorm = (xLoc-xMin)/(xMax-xMin)
  yLocNorm = (yLoc-yMin)/(yMax-yMin)

  wghts(1) = (1.0_RFREAL-xLocNorm)*(1.0_RFREAL-yLocNorm)
  wghts(2) =             xLocNorm *(1.0_RFREAL-yLocNorm)
  wghts(3) = (1.0_RFREAL-xLocNorm)*            yLocNorm 
  wghts(4) =             xLocNorm *            yLocNorm 

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputeBilinearInterpWghts

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeBilinearInterpWghts.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:55  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:07  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/11/27 13:17:26  haselbac
! Initial revision
!
! ******************************************************************************

