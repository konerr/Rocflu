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
! Purpose: read in user input related to calculation of thrust.
!
! Description: none.
!
! Input: user input file.
!
! Output: global = parameters for calculation of thrust.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ReadThrustSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadThrustSection( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  CHARACTER(10) :: keys(7)

  LOGICAL :: defined(7)

  REAL(RFREAL) :: vals(7)

!******************************************************************************

  CALL RegisterFunction( global,'ReadThrustSection',__FILE__ )

! specify keywords and search for them

  keys(1) = 'TYPE'
  keys(2) = 'PLANE'
  keys(3) = 'COORD'
  keys(4) = 'WRITIME'
  keys(5) = 'WRIITER'
  keys(6) = 'OPENCLOSE'
  keys(7) = 'PAMB'
  
  CALL ReadSection( global,IF_INPUT,7,keys,vals,defined )

  IF (defined(1)) THEN
                                       global%thrustType = THRUST_NONE
    IF (vals(1)>0.9 .AND. vals(1)<1.1) global%thrustType = THRUST_MOM
    IF (vals(1) > 1.9)                 global%thrustType = THRUST_MOMP
  ENDIF
  IF (defined(2)) THEN
    IF (vals(2) < 1.1)                  global%thrustPlane = XCOORD
    IF (vals(2)>=1.1 .AND. vals(2)<2.1) global%thrustPlane = YCOORD
    IF (vals(2) >= 2.1)                 global%thrustPlane = ZCOORD
  ENDIF
  IF (defined(3)) THEN
    global%thrustCoord = vals(3)
  ENDIF
  IF (defined(4)) global%thrustSaveTime = ABS(vals(4))
  IF (defined(5)) THEN
    global%thrustSaveIter = INT(ABS(vals(5))+0.5_RFREAL)
    global%thrustSaveIter = MAX(1,global%thrustSaveIter)
  ENDIF
  IF (defined(6)) THEN
    IF (vals(6) < 0.5_RFREAL) THEN
      global%thrustOpenClose = .false.
    ELSE
      global%thrustOpenClose = .true.
    ENDIF
  ENDIF
  IF (defined(7)) THEN
    global%thrustPamb = ABS(vals(7))
  ENDIF

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadThrustSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadThrustSection.F90,v $
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
! Revision 1.1  2004/12/01 16:50:52  haselbac
! Initial revision after changing case
!
! Revision 1.4  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.1  2003/06/02 17:12:00  jblazek
! Added computation of thrust.
!
!******************************************************************************

