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
! Purpose: Read in user input related to the Moving Frame scheme.
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
! $Id: ReadMvFrameSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReadMvFrameSection(global)

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

  INTEGER, PARAMETER :: NVALS_MAX = 18

  LOGICAL :: defined(NVALS_MAX)
  CHARACTER(10) :: keys(NVALS_MAX)
  INTEGER :: iReg,nVals
  REAL(RFREAL) :: vals(NVALS_MAX)

! ******************************************************************************
! Start, specify keywords and search for them
! ******************************************************************************

  keys(1)  = 'FLAG'
  keys(2)  = 'IPATCH'
  keys(3)  = 'VELX'
  keys(4)  = 'VELY'
  keys(5)  = 'VELZ'
  keys(6)  = 'MASS'
  keys(7)  = 'ACCTS'
  keys(8)  = 'ACCTE'
  keys(9)  = 'ACCX'
  keys(10) = 'ACCY'
  keys(11) = 'ACCZ'
  keys(12) = 'ACCFLAG'
  keys(13) = 'LOCX'
  keys(14) = 'LOCY'
  keys(15) = 'LOCZ'
  keys(16) = 'ACCTYPE'
  keys(17) = 'OMEGA'
  keys(18) = 'PHASE'

  nVals = NVALS_MAX

  CALL RegisterFunction(global,'ReadMvFrameSection',__FILE__)

  CALL ReadSection(global,IF_INPUT,nVals,keys,vals,defined) 
  
  IF ( defined(1) .EQV. .TRUE. ) THEN
    IF ( NINT(vals(1)) == 1 ) THEN 
      global%mvFrameFlag = .TRUE.
    ELSE 
      global%mvFrameFlag = .FALSE.
    END IF ! NINT
  END IF ! defined  
  
  IF ( defined(2) .EQV. .TRUE. ) THEN
    global%iPatchGlobalMvFrame = NINT(vals(2))
  END IF ! defined  
  
  IF ( defined(3) .EQV. .TRUE. ) THEN
    global%mvfVelInitX = vals(3)
  END IF ! defined  
  
  IF ( defined(4) .EQV. .TRUE. ) THEN
    global%mvfVelInitY = vals(4)
  END IF ! defined  
  
  IF ( defined(5) .EQV. .TRUE. ) THEN
    global%mvfVelInitZ = vals(5)
  END IF ! defined  
  
  IF ( defined(6) .EQV. .TRUE. ) THEN
    global%mvfMass = vals(6)
  END IF ! defined  
  
  IF ( defined(7) .EQV. .TRUE. ) THEN
    global%mvfAccTs = vals(7)
  END IF ! defined  
  
  IF ( defined(8) .EQV. .TRUE. ) THEN
    global%mvfAccTe = vals(8)
  END IF ! defined  
  
  IF ( defined(9) .EQV. .TRUE. ) THEN
    global%mvfAccX = vals(9)
  END IF ! defined  
  
  IF ( defined(10) .EQV. .TRUE. ) THEN
    global%mvfAccY = vals(10)
  END IF ! defined  
  
  IF ( defined(11) .EQV. .TRUE. ) THEN
    global%mvfAccZ = vals(11)
  END IF ! defined  
  
  IF ( defined(12) .EQV. .TRUE. ) THEN
    IF ( NINT(vals(12)) == 1 ) THEN 
      global%mvfAccFlag = .TRUE.
    ELSE 
      global%mvfAccFlag = .FALSE.
    END IF ! NINT
  END IF ! defined  
  
  IF ( defined(13) .EQV. .TRUE. ) THEN
    global%mvfLocInitX = vals(13)
  END IF ! defined  
  
  IF ( defined(14) .EQV. .TRUE. ) THEN
    global%mvfLocInitY = vals(14)
  END IF ! defined  
  
  IF ( defined(15) .EQV. .TRUE. ) THEN
    global%mvfLocInitZ = vals(15)
  END IF ! defined  
  
  IF ( defined(16) .EQV. .TRUE. ) THEN
    global%mvfAccType = NINT(vals(16))
  END IF ! defined  
  
  IF ( defined(17) .EQV. .TRUE. ) THEN
    global%mvfOmega = vals(17)
  END IF ! defined  
  
  IF ( defined(18) .EQV. .TRUE. ) THEN
    global%mvfPhase = vals(18)
  END IF ! defined  
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE ReadMvFrameSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadMvFrameSection.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.5  2009/07/08 19:11:24  mparmar
! Removed REAL type conversion of input data
!
! Revision 1.4  2008/12/06 08:43:32  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:16:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2008/05/29 01:35:08  mparmar
! Added 3 keywords to allow different types of imposed acceleration
!
! Revision 1.1  2007/06/18 17:31:22  mparmar
! Initial revision
!
! ******************************************************************************

