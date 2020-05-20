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
! Purpose: Read in user input for various specific tasks.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ReadMiscSection.F90,v 1.2 2015/02/10 18:24:22 brollin Exp $
!
! Copyright: (c) 2007-2008 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReadMiscSection(global)

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

  ! BBR - Modified to accomodate 10 Misc qty
 
  LOGICAL :: defined(10)
  CHARACTER(10) :: keys(10)
  REAL(RFREAL) :: vals(10)

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction( global,'ReadMiscSection',__FILE__ )

! ******************************************************************************
! Specify keywords and search for them
! ******************************************************************************

  keys(1) = 'NBVECFACT'
  ! BBR - begin
  keys(2) = 'CYCLETIME'
  keys(3) = 'TEST_BBR'
  keys(4) = 'TEST_SUBBU'
  keys(5) = 'TEST_CHRIS'
  keys(6) = 'TEST_YASH'
  keys(7) = 'TEST_FRED'
  keys(8) = 'TEST_SAPPY'
  keys(9) = 'TEST_RAHUL'
  keys(10) = 'TEST_BRAD'
  ! BBR -end
  
  CALL ReadSection(global,IF_INPUT,10,keys,vals,defined)

! ******************************************************************************
! Set values
! ******************************************************************************

  IF ( defined(1) .EQV. .TRUE. ) THEN
    global%nBVertEstCorrFactor = vals(1)
  END IF ! defined

  ! BBR - Begin
  IF ( defined(2) .EQV. .TRUE. ) THEN
    IF (vals(2) .EQ. 1 ) THEN
      global%TimeCycle = .TRUE.
    ELSE
      global%TimeCycle = .FALSE.
    END IF ! vals
  END IF ! defined

  IF ( defined(3) .EQV. .TRUE. ) THEN
    IF (vals(3) .EQ. 1 ) THEN
      global%Test_BBR = .TRUE.
    ELSE
      global%Test_BBR = .FALSE.
    END IF ! vals
  END IF ! defined

  IF ( defined(4) .EQV. .TRUE. ) THEN
    IF (vals(4) .EQ. 1 ) THEN
      global%Test_Subbu = .TRUE.
    ELSE
      global%Test_Subbu = .FALSE.
    END IF ! vals
  END IF ! defined

  IF ( defined(5) .EQV. .TRUE. ) THEN
    IF (vals(5) .EQ. 1 ) THEN
      global%Test_Chris = .TRUE.
    ELSE
      global%Test_Chris = .FALSE.
    END IF ! vals
  END IF ! defined

  IF ( defined(6) .EQV. .TRUE. ) THEN
    IF (vals(6) .EQ. 1 ) THEN
      global%Test_Yash = .TRUE.
    ELSE
      global%Test_Yash = .FALSE.
    END IF ! vals
  END IF ! defined

  IF ( defined(7) .EQV. .TRUE. ) THEN
    IF (vals(7) .EQ. 1 ) THEN
      global%Test_Fred = .TRUE.
    ELSE
      global%Test_Fred = .FALSE.
    END IF ! vals
  END IF ! defined

  IF ( defined(8) .EQV. .TRUE. ) THEN
    IF (vals(8) .EQ. 1 ) THEN
      global%Test_Sappy = .TRUE.
    ELSE
      global%Test_Sappy = .FALSE.
    END IF ! vals
  END IF ! defined

  IF ( defined(9) .EQV. .TRUE. ) THEN
    IF (vals(9) .EQ. 1 ) THEN
      global%Test_Rahul = .TRUE.
    ELSE
      global%Test_Rahul = .FALSE.
    END IF ! vals
  END IF ! defined

  IF ( defined(10) .EQV. .TRUE. ) THEN
    IF (vals(10) .EQ. 1 ) THEN
      global%Test_Brad = .TRUE.
    ELSE
      global%Test_Brad = .FALSE.
    END IF ! vals
  END IF ! defined

  !BBR - end

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE ReadMiscSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadMiscSection.F90,v $
! Revision 1.2  2015/02/10 18:24:22  brollin
! adding missing Metis library files
!
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
! Revision 1.1  2007/11/28 23:04:28  mparmar
! Initial revision
!
! ******************************************************************************

