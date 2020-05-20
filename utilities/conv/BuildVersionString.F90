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
! Purpose: Build version string for printing in header.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   versionString       String containing version number and date
!
! Notes: 
!   1. The strings are NOT to be edited by anyone except the main code 
!      developer. 
!   2. Marks Rocbuild program will edit the build string to insert a 
!      build number or identifier
!
! ******************************************************************************
!
! $Id: BuildVersionString.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002-2008 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE BuildVersionString(versionString)

  USE ModDataTypes
  USE ModError
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
!*******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(CHRLEN), INTENT(OUT) :: versionString

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(LEN=2) :: major,minor,patch
  CHARACTER(LEN=4) :: build
  CHARACTER(CHRLEN) :: date,RCSIdentString

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: BuildVersionString.F90,v $ $Revision: 1.1.1.1 $'

! ==============================================================================
! Set strings: DO NOT EDIT UNLESS YOU ARE MAIN DEVELOPER
! ==============================================================================

  major = '7'
  minor = '1'
  patch = '0'
  build = '0' ! to be edited by Rocbuild

  date  = '05/21/10'

! ==============================================================================
! Write into version string
! ==============================================================================

  WRITE(versionString,'(A)') TRIM(major)//'.'//TRIM(minor)//'.'//TRIM(patch)
  WRITE(versionString,'(A)') 'Version: '//TRIM(versionString)//'-'//TRIM(build)  
  WRITE(versionString,'(A)') TRIM(versionString)//', Date: '//TRIM(date)

! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE BuildVersionString

