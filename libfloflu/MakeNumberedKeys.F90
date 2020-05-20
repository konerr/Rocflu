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
! Purpose: Fills in a section of an array of strings with an input string
!          appended by a sequence of numbers.
!
! Description:
!
! Input:  indBegin: first element of keys to write to
!         string:   string to append numbers to
!         numBegin, numEnd, numSkip: DO loop specification of numbers to append
!
! Output: keys: the array of strings with numbers appended
!
! Notes: indBegin addresses keys as if it started at 1.
!        Does not check to see if strings written to are long enough.
!        Numbers appended must be between 0 and 10^7 - 1, inclusive.
!
!******************************************************************************
!
! $Id: MakeNumberedKeys.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE MakeNumberedKeys(keys,indBegin,string,numBegin,numEnd,numSkip)

  IMPLICIT NONE

! ... parameters
  CHARACTER(*)             :: keys(:)
  CHARACTER(*), INTENT(in) :: string
  INTEGER,      INTENT(in) :: indBegin,numBegin,numEnd,numSkip

! ... loop variables
  INTEGER :: num

! ... local variables
  INTEGER :: nKeys, iKeys, skip

!******************************************************************************

  nKeys = UBOUND(keys,1)

  iKeys = indBegin

  skip = numSkip
  IF (skip == 0) skip = 1

  DO num = numBegin,numEnd,skip

    IF (iKeys > nKeys) EXIT

    IF (iKeys > 0) THEN

      SELECT CASE (num)

      CASE (      0:      9)
        WRITE(keys(iKeys),'(A,I1)') TRIM(string), num

      CASE (     10:     99)
        WRITE(keys(iKeys),'(A,I2)') TRIM(string), num

      CASE (    100:    999)
        WRITE(keys(iKeys),'(A,I3)') TRIM(string), num

      CASE (   1000:   9999)
        WRITE(keys(iKeys),'(A,I4)') TRIM(string), num

      CASE (  10000:  99999)
        WRITE(keys(iKeys),'(A,I5)') TRIM(string), num

      CASE ( 100000: 999999)
        WRITE(keys(iKeys),'(A,I6)') TRIM(string), num

      CASE (1000000:9999999)
        WRITE(keys(iKeys),'(A,I7)') TRIM(string), num

      CASE DEFAULT
        WRITE(keys(iKeys),'(A)') TRIM(string)

      END SELECT ! num

    END IF ! iKeys

    iKeys = iKeys + 1

  END DO ! num

END SUBROUTINE MakeNumberedKeys

!******************************************************************************
!
! RCS Revision history:
!
! $Log: MakeNumberedKeys.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:31  mtcampbe
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
! Revision 1.1  2004/12/01 16:48:45  haselbac
! Initial revision after changing case
!
! Revision 1.1  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
!
!******************************************************************************

