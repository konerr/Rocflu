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
! Purpose: Read in user input related to the absorbing boundary conditions.
!
! Description: None.
!
! Input: from file.
!
! Output:
!   global  = global parameters related to grid motion (RFLO)
!   regions = region data (RFLU).
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ReadABCSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! ******************************************************************************

SUBROUTINE ReadABCSection(global)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModInterfaces, ONLY: ReadSection
  USE ModError
  USE ModParameters

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER, PARAMETER :: NVALS_MAX = 21

  LOGICAL :: defined(NVALS_MAX)
  CHARACTER(10) :: keys(NVALS_MAX)
  INTEGER :: iReg,nVals
  REAL(RFREAL) :: vals(NVALS_MAX)

! ******************************************************************************
! Start, specify keywords and search for them
! ******************************************************************************

  keys(1)  = 'FLAG'
  keys(2)  = 'ORDER'
  keys(3)  = 'COEFF0'
  keys(4)  = 'COEFF1'
  keys(5)  = 'COEFF2'
  keys(6)  = 'DOMAIN'
  keys(7)  = 'LEFT'
  keys(8)  = 'RIGHT'
  keys(9)  = 'BOTTOM'
  keys(10) = 'TOP'
  keys(11) = 'CENTERX'
  keys(12) = 'CENTERY'
  keys(13) = 'RADIUS'
  keys(14) = 'KIND'
  keys(15) = 'SPONGE'
  keys(16) = 'DENS'
  keys(17) = 'VELX'
  keys(18) = 'VELY'
  keys(19) = 'VELZ'
  keys(20) = 'PRESS'
  keys(21) = 'DISTRIB'

  nVals = NVALS_MAX

  CALL RegisterFunction(global,'ReadABCSection',__FILE__)

  CALL ReadSection(global,IF_INPUT,nVals,keys,vals,defined) 
  
  IF ( defined(1) .EQV. .TRUE. ) THEN
    IF ( NINT(vals(1)) == 1 ) THEN 
      global%abcFlag = .TRUE.
    END IF ! NINT
  END IF ! defined  
  
  IF ( defined(2) .EQV. .TRUE. ) THEN
    global%abcOrder = NINT(vals(2))
  END IF ! defined  
  
  IF ( defined(3) .EQV. .TRUE. ) THEN
    global%abcCoeff0 = vals(3)
  END IF ! defined  
  
  IF ( defined(4) .EQV. .TRUE. ) THEN
    global%abcCoeff1 = vals(4)
  END IF ! defined  
  
  IF ( defined(5) .EQV. .TRUE. ) THEN
    global%abcCoeff2 = vals(5)
  END IF ! defined  
  
  IF ( defined(6) .EQV. .TRUE. ) THEN
    global%abcDomain = NINT(vals(6))
  END IF ! defined  
  
  IF ( defined(7) .EQV. .TRUE. ) THEN
    global%abcLeft = vals(7)
  END IF ! defined  
  
  IF ( defined(8) .EQV. .TRUE. ) THEN
    global%abcRight = vals(8)
  END IF ! defined  
  
  IF ( defined(9) .EQV. .TRUE. ) THEN
    global%abcBottom = vals(9)
  END IF ! defined  
  
  IF ( defined(10) .EQV. .TRUE. ) THEN
    global%abcTop = vals(10)
  END IF ! defined  
  
  IF ( defined(11) .EQV. .TRUE. ) THEN
    global%abcCenterX = vals(11)
  END IF ! defined  
  
  IF ( defined(12) .EQV. .TRUE. ) THEN
    global%abcCenterY = vals(12)
  END IF ! defined  
  
  IF ( defined(13) .EQV. .TRUE. ) THEN
    global%abcRadius = vals(13)
  END IF ! defined  
  
  IF ( defined(14) .EQV. .TRUE. ) THEN
    global%abcKind = NINT(vals(14))
  END IF ! defined  
  
  IF ( defined(15) .EQV. .TRUE. ) THEN
    global%abcSponge = vals(15)
  END IF ! defined  
  
  IF ( defined(16) .EQV. .TRUE. ) THEN
    global%abcDens = vals(16)
  END IF ! defined  
  
  IF ( defined(17) .EQV. .TRUE. ) THEN
    global%abcVelX = vals(17)
  END IF ! defined  
  
  IF ( defined(18) .EQV. .TRUE. ) THEN
    global%abcVelY = vals(18)
  END IF ! defined  
  
  IF ( defined(19) .EQV. .TRUE. ) THEN
    global%abcVelZ = vals(19)
  END IF ! defined  
  
  IF ( defined(20) .EQV. .TRUE. ) THEN
    global%abcPress = vals(20)
  END IF ! defined  
  
  IF ( defined(21) .EQV. .TRUE. ) THEN
    IF ( NINT(vals(21)) == 0 ) THEN
      global%abcDistrib = 0
    ELSE
      global%abcDistrib = 1
    END IF ! NINT(vals(21))
  END IF ! defined  
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE ReadABCSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadABCSection.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.1  2009/07/08 19:11:15  mparmar
! Initial revision
!
!
! ******************************************************************************
