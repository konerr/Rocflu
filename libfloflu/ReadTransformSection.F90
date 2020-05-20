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
! Purpose: Read in user input related to geometric transformation of grid and 
!  solution.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: ReadTransformSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadTransformSection( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  INTEGER :: nVals
  INTEGER, PARAMETER :: NVALS_MAX = 12

  CHARACTER(CHRLEN) :: RCSIdentString
  CHARACTER(10)     :: keys(NVALS_MAX)

  LOGICAL :: defined(NVALS_MAX)

  REAL(RFREAL) :: vals(NVALS_MAX)

!******************************************************************************

  RCSIdentString = '$RCSfile: ReadTransformSection.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction( global,'ReadTranformSection',__FILE__ )

! specify keywords and search for them

  nVals = NVALS_MAX
  
  keys( 1) = 'FLAG'
  keys( 2) = 'SCALE_X'
  keys( 3) = 'SCALE_Y'
  keys( 4) = 'SCALE_Z'
  keys( 5) = 'ANGLE_X'
  keys( 6) = 'ANGLE_Y'
  keys( 7) = 'ANGLE_Z'
  keys( 8) = 'DIST_FLAG'
  keys( 9) = 'DIST_X'
  keys(10) = 'DIST_Y'
  keys(11) = 'DIST_Z'
  keys(12) = 'ENFORCE'
  
  CALL ReadSection( global,IF_INPUT,nVals,keys(1:nVals),vals(1:nVals), & 
                    defined(1:nVals) )

! transformation flag 

  IF (defined(1)) THEN
    IF ( NINT(vals(1)) == 1 ) THEN 
      global%transformFlag = .TRUE.
    ELSE 
      global%transformFlag = .FALSE.
    END IF ! NINT
  ELSE 
    global%transformFlag = .FALSE.
  ENDIF ! defined(1)

! scaling factors 

  IF (defined(2)) THEN
    global%scaleX = vals(2)
  ELSE 
    global%scaleX = 1.0_RFREAL
  ENDIF ! defined(2)

  IF (defined(3)) THEN
    global%scaleY = vals(3)
  ELSE 
    global%scaleY = 1.0_RFREAL
  ENDIF ! defined(3)

  IF (defined(4)) THEN
    global%scaleZ = vals(4)
  ELSE 
    global%scaleZ = 1.0_RFREAL
  ENDIF ! defined(4)

! rotation angles 

  IF (defined(5)) THEN
    global%angleX = vals(5)
  ELSE 
    global%angleX = 0.0_RFREAL
  ENDIF ! defined(5)

  IF (defined(6)) THEN
    global%angleY = vals(6)
  ELSE 
    global%angleY = 0.0_RFREAL
  ENDIF ! defined(6)

  IF (defined(7)) THEN
    global%angleZ = vals(7)
  ELSE 
    global%angleZ = 0.0_RFREAL
  ENDIF ! defined(7)

! distortion

  IF (defined(8)) THEN
    IF ( NINT(vals(8)) == 1 ) THEN 
      global%distortFlag = .TRUE.
    ELSE 
      global%distortFlag = .FALSE.
    END IF ! NINT
  ELSE 
    global%distortFlag = .FALSE.
  ENDIF ! defined(8)
  
  IF (defined(9)) THEN
    global%distortX = vals(9)
  ELSE 
    global%distortX = 0.0_RFREAL
  ENDIF ! defined(9)

  IF (defined(10)) THEN
    global%distortY = vals(10)
  ELSE 
    global%distortY = 0.0_RFREAL
  ENDIF ! defined(10)

  IF (defined(11)) THEN
    global%distortZ = vals(11)
  ELSE 
    global%distortZ = 0.0_RFREAL
  ENDIF ! defined(11)

! boundary patch coordinate enforcement

  IF (defined(12)) THEN 
    IF ( NINT(vals(12)) == 1 ) THEN 
      global%enforceFlag = .TRUE.
    ELSE 
      global%enforceFlag = .FALSE.
    END IF ! NINT
  ELSE 
    global%enforceFlag = .FALSE.
  END IF ! defined(12)

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadTransformSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadTransformSection.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/08/04 03:00:28  haselbac
! Added reading of distortion parameters
!
! Revision 1.1  2004/12/01 16:50:57  haselbac
! Initial revision after changing case
!
! Revision 1.6  2003/05/16 02:27:43  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.5  2003/03/15 16:32:50  haselbac
! Added KIND qualifyers
!
! Revision 1.4  2003/02/01 00:28:56  haselbac
! Added enforceFlag, transformFlag now LOGICAL
!
! Revision 1.3  2002/11/08 21:21:40  haselbac
! Bug fix for transformFlag
!
! Revision 1.2  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/05/07 18:52:26  haselbac
! Initial revision
!
!******************************************************************************

