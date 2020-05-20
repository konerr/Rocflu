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
! Purpose: Scale and rotate vector.
!
! Description: None.
!
! Input: vect = vector to scale and rotate.
!
! Output: vect = scaled and rotated.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: ScaleRotateVector.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ScaleRotateVector( global,vect )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  USE ModParameters

  IMPLICIT NONE
    
! ... parameters
  REAL(RFREAL), POINTER :: vect(:,:)

  TYPE(t_global), POINTER :: global
  
! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: ic, nComp

  REAL(RFREAL) :: angleX, angleY, angleZ, cx, cy, cz, sx, sy, sz, x, y, z

!******************************************************************************

  RCSIdentString = '$RCSfile: ScaleRotateVector.F90,v $ $Revision: 1.1.1.1 $'

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Scaling and rotating vector...'
    WRITE(STDOUT,'(A,3X,A,1X,E15.6)') SOLVER_NAME,'X-scale:',global%scaleX
    WRITE(STDOUT,'(A,3X,A,1X,E15.6)') SOLVER_NAME,'Y-scale:',global%scaleY
    WRITE(STDOUT,'(A,3X,A,1X,E15.6)') SOLVER_NAME,'Z-scale:',global%scaleZ
    WRITE(STDOUT,'(A,3X,A,1X,E15.6)') SOLVER_NAME,'X-angle:',global%angleX
    WRITE(STDOUT,'(A,3X,A,1X,E15.6)') SOLVER_NAME,'Y-angle:',global%angleY
    WRITE(STDOUT,'(A,3X,A,1X,E15.6)') SOLVER_NAME,'Z-angle:',global%angleZ
  END IF ! global%verbLevel

  nComp = SIZE(vect,DIM=2)

! scale vector

  DO ic = 1,nComp
    vect(1,ic) = global%scaleX*vect(1,ic)
    vect(2,ic) = global%scaleY*vect(2,ic)
    vect(3,ic) = global%scaleZ*vect(3,ic)
  END DO ! ic

! rotate vector

  angleX = global%angleX*global%rad
  angleY = global%angleY*global%rad
  angleZ = global%angleZ*global%rad   

  cx = COS(angleX)
  sx = SIN(angleX)
  cy = COS(angleY)
  sy = SIN(angleY)
  cz = COS(angleZ)
  sz = SIN(angleZ)
  
  DO ic = 1,nComp
    x = vect(1,ic)
    y = vect(2,ic)
    z = vect(3,ic)

    vect(1,ic) = cy*cz*x - (sx*sy*cz + cx*sz)*y - (cx*sy*cz - sx*sz)*z
    vect(2,ic) = cy*sz*x - (sx*sy*sz - cx*cz)*y - (cx*sy*sz + sx*cz)*z
    vect(3,ic) =    sy*x + (           sx*cy)*y + (           cx*cy)*z
  END DO ! ic
  
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Scaling and rotating vector done.'
  END IF ! global%verbLevel
  
END SUBROUTINE ScaleRotateVector

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ScaleRotateVector.F90,v $
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
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 16:51:21  haselbac
! Initial revision after changing case
!
! Revision 1.5  2004/01/29 22:52:48  haselbac
! Clean up
!
! Revision 1.4  2003/02/01 00:28:20  haselbac
! Some clean up, added diagnostic output
!
! Revision 1.3  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/07/16 21:34:37  jblazek
! Prefixed screen output with SOLVER_NAME.
!
! Revision 1.1  2002/05/07 18:51:27  haselbac
! Initial revision
!
!******************************************************************************

