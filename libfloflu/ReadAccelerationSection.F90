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
! Purpose: Read in user input related to acceleration terms.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ReadAccelerationSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReadAccelerationSection(global)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  
  USE ModInterfaces, ONLY: ReadSection  
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Local variables
! ==============================================================================

  LOGICAL :: defined(5)
  CHARACTER(10) :: keys(5)
  REAL(RFREAL) :: vals(5)

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction( global,'ReadAccelerationSection',__FILE__ )

! ******************************************************************************
! Specify keywords and search for them
! ******************************************************************************

  keys(1) = 'TYPE'
  keys(2) = 'ACCELX'
  keys(3) = 'ACCELY'
  keys(4) = 'ACCELZ'
  keys(5) = 'GRAVITY'
  
  CALL ReadSection(global,IF_INPUT,5,keys,vals,defined)

! ******************************************************************************
! Set values
! ******************************************************************************

  IF ( defined(1) .EQV. .TRUE. ) THEN
    IF ( NINT(vals(1)) == 1 ) THEN 
      global%accelOn = .TRUE.
    ELSE 
      global%accelOn = .FALSE.
    END IF ! NINT(vals)
  ELSE 
    global%accelOn = .FALSE.
  END IF ! defined

  IF ( defined(2) .EQV. .TRUE. ) THEN
    global%accelX = vals(2)
  END IF ! defined
  
  IF ( defined(3) .EQV. .TRUE. ) THEN
    global%accelY = vals(3)
  END IF ! defined
  
  IF ( defined(4) .EQV. .TRUE. ) THEN
    global%accelZ = vals(4)
  END IF ! defined

  IF ( defined(5) .EQV. .TRUE. ) THEN
    global%gravity = vals(5)
  END IF ! defined

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE ReadAccelerationSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadAccelerationSection.F90,v $
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
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/10/20 21:20:19  mparmar
! Added reading of gravity
!
! Revision 1.3  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.2  2005/11/10 01:55:42  haselbac
! Added Rocflu support, clean-up
!
! Revision 1.1  2004/12/01 16:50:09  haselbac
! Initial revision after changing case
!
! Revision 1.4  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.1  2003/08/28 20:05:38  jblazek
! Added acceleration terms.
!
! ******************************************************************************

