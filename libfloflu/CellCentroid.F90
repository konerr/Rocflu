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
! Purpose: calculate the centroid of a cell.
!
! Description: file contains the following subroutines:
!
!  - CentroidTetra   = tetrahedron
!  - CentroidPyramid = pyramid
!  - CentroidPrism   = prism
!  - CentroidHexa    = hexahedra
!
! Input: xyzNodes = coordinates (1st index) of the nodes (2nd index) of the
!                   control volume (must be ordered!).
!
! Output: cofgX,cofgY,cofgZ = x-,y-,z-coordinate of the centroid.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: CellCentroid.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE CentroidTetra( xyzNodes,cofgX,cofgY,cofgZ )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: xyzNodes(3,4)
  REAL(RFREAL) :: cofgX, cofgY, cofgZ

! ... local variables


!******************************************************************************

END SUBROUTINE CentroidTetra

! #############################################################################
! #############################################################################

SUBROUTINE CentroidPyramid( xyzNodes,cofgX,cofgY,cofgZ )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: xyzNodes(3,5)
  REAL(RFREAL) :: cofgX, cofgY, cofgZ

! ... local variables


!******************************************************************************

END SUBROUTINE CentroidPyramid

! #############################################################################
! #############################################################################

SUBROUTINE CentroidPrism( xyzNodes,cofgX,cofgY,cofgZ )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: xyzNodes(3,6)
  REAL(RFREAL) :: cofgX, cofgY, cofgZ

! ... local variables


!******************************************************************************

END SUBROUTINE CentroidPrism

! #############################################################################
! #############################################################################

SUBROUTINE CentroidHexa( xyzNodes,cofgX,cofgY,cofgZ )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: xyzNodes(3,8)
  REAL(RFREAL) :: cofgX, cofgY, cofgZ

!******************************************************************************

  cofgX = 0.125_RFREAL*(xyzNodes(1,1)+xyzNodes(1,2)+xyzNodes(1,3)+ &
                        xyzNodes(1,4)+xyzNodes(1,5)+xyzNodes(1,6)+ &
                        xyzNodes(1,7)+xyzNodes(1,8))
  cofgY = 0.125_RFREAL*(xyzNodes(2,1)+xyzNodes(2,2)+xyzNodes(2,3)+ &
                        xyzNodes(2,4)+xyzNodes(2,5)+xyzNodes(2,6)+ &
                        xyzNodes(2,7)+xyzNodes(2,8))
  cofgZ = 0.125_RFREAL*(xyzNodes(3,1)+xyzNodes(3,2)+xyzNodes(3,3)+ &
                        xyzNodes(3,4)+xyzNodes(3,5)+xyzNodes(3,6)+ &
                        xyzNodes(3,7)+xyzNodes(3,8))

END SUBROUTINE CentroidHexa

!******************************************************************************
!
! RCS Revision history:
!
! $Log: CellCentroid.F90,v $
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
! Revision 1.1  2004/12/01 16:48:18  haselbac
! Initial revision after changing case
!
! Revision 1.1  2002/01/08 22:09:16  jblazek
! Added calculation of face vectors and volumes.
!
!******************************************************************************

