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
! Purpose: Calculate the face centroid of a triangular or quadrilateral face.
!
! Description: File contains the following subroutines:
!  - FaceCentroidTria = Centroid of a triangle
!  - FaceCentroidQuad = Centroid of a quadrilateral
!
! Input: xyzNodes = coordinates (1st index) of the nodes (2nd index) of the
!                   polygon (clockwise ordered).
!
! Output: fCenX, fCenY, fCenZ = x-,y-,z-component of the face centroid.
!
! Notes:
!   1. Quadrilateral face centroid is only approximate for distorted faces.
!
!******************************************************************************
!
! $Id: FaceCentroid.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE FaceCentroidTria( xyz,fCenX,fCenY,fCenZ )

  USE ModDataTypes

  IMPLICIT NONE

! ... parameters

  REAL(RFREAL), INTENT(IN) :: xyz(3,3)
  REAL(RFREAL), INTENT(OUT) :: fCenX, fCenY, fCenZ

! ... local variables

  REAL(RFREAL), PARAMETER :: thrd = 1.0_RFREAL/3.0_RFREAL 
  
!******************************************************************************

  fCenX = thrd*(xyz(1,1) + xyz(1,2) + xyz(1,3))
  fCenY = thrd*(xyz(2,1) + xyz(2,2) + xyz(2,3))
  fCenZ = thrd*(xyz(3,1) + xyz(3,2) + xyz(3,3))

END SUBROUTINE FaceCentroidTria

! #############################################################################
! #############################################################################

SUBROUTINE FaceCentroidQuad( xyz,fCenX,fCenY,fCenZ )

  USE ModDataTypes
  
  IMPLICIT NONE

! ... parameters

  REAL(RFREAL), INTENT(IN) :: xyz(3,4)
  REAL(RFREAL), INTENT(OUT) :: fCenX, fCenY, fCenZ

! ... local variables

!******************************************************************************

  fCenX = 0.25_RFREAL*(xyz(1,1) + xyz(1,2) + xyz(1,3) + xyz(1,4))
  fCenY = 0.25_RFREAL*(xyz(2,1) + xyz(2,2) + xyz(2,3) + xyz(2,4))
  fCenZ = 0.25_RFREAL*(xyz(3,1) + xyz(3,2) + xyz(3,3) + xyz(3,4))

END SUBROUTINE FaceCentroidQuad

!******************************************************************************
!
! RCS Revision history:
!
! $Log: FaceCentroid.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:24  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 16:48:32  haselbac
! Initial revision after changing case
!
! Revision 1.1  2002/03/14 19:01:10  haselbac
! Initial revision
!
!******************************************************************************

