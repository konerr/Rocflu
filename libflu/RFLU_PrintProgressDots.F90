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
! Purpose: Print dots on a line to indicate progress of a length task.
!
! Description: None.
!
! Input: 
!   nowValue    Current value
!   endValye    End value
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_PrintProgressDots.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_PrintProgressDots(nowValue,endValue)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: global 
  USE ModParameters
  
  IMPLICIT NONE

! ... parameters
  INTEGER, INTENT(IN) :: endValue,nowValue

! ... loop variables


! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  REAL(RFREAL) :: pDone

!******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PrintProgressDots.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction('RFLU_PrintProgressDots',__FILE__)

! start -----------------------------------------------------------------------

  pDone = 100.0_RFREAL*nowValue/REAL(endValue,KIND=RFREAL)
    
  IF ( INT(pDone) >= 10*global%progressCounter ) THEN 
    WRITE(STDOUT,'(A)',ADVANCE="NO") '.'
    global%progressCounter = global%progressCounter + 1
  END IF ! NINT
  
  IF ( nowValue == endValue ) THEN 
    global%progressCounter = 1
    WRITE(STDOUT,'(A)') ' '
  END IF ! nowValue

! end -------------------------------------------------------------------------

  CALL DeregisterFunction()

END SUBROUTINE 

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PrintProgressDots.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:35  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:53  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2003/05/16 02:27:44  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.3  2003/03/15 17:06:38  haselbac
! Added KIND qualifyer
!
! Revision 1.2  2002/10/27 18:54:05  haselbac
! Removed tabs
!
! Revision 1.1  2002/07/25 14:34:59  haselbac
! Initial revision
!
!******************************************************************************

