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
! $Id: ReadGFMSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! ******************************************************************************

SUBROUTINE ReadGFMSection(global)

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

  INTEGER, PARAMETER :: NVALS_MAX = 12

  LOGICAL :: defined(NVALS_MAX)
  CHARACTER(10) :: keys(NVALS_MAX)
  INTEGER :: iReg,nVals
  REAL(RFREAL) :: vals(NVALS_MAX)

! ******************************************************************************
! Start, specify keywords and search for them
! ******************************************************************************

  keys(1)  = 'FLAG'
  keys(2)  = 'NRIGIDS'
  keys(3)  = 'NFLUIDS'
  keys(4)  = 'X1'
  keys(5)  = 'X2'
  keys(6)  = 'Y1'
  keys(7)  = 'Y2'
  keys(8)  = 'Z1'
  keys(9)  = 'Z2'
  keys(10) = 'NX'
  keys(11) = 'NY'
  keys(12) = 'NZ'

  nVals = NVALS_MAX

  CALL RegisterFunction(global,'ReadGFMSection',__FILE__)

  CALL ReadSection(global,IF_INPUT,nVals,keys,vals,defined) 
  
  IF ( defined(1) .EQV. .TRUE. ) THEN
    IF ( NINT(vals(1)) == 1 ) THEN 
      global%gfmFlag = .TRUE.
    END IF ! NINT
  END IF ! defined  
  
  IF ( defined(2) .EQV. .TRUE. ) THEN
    global%gfmNRigids = NINT(vals(2))
  END IF ! defined  
  
  IF ( defined(3) .EQV. .TRUE. ) THEN
    global%gfmNFluids = NINT(vals(3))
  END IF ! defined  
  
  IF ( defined(4) .EQV. .TRUE. ) THEN
    global%gfmGridX1 = vals(4)
  END IF ! defined  
  
  IF ( defined(5) .EQV. .TRUE. ) THEN
    global%gfmGridX2 = vals(5)
  END IF ! defined  
  
  IF ( defined(6) .EQV. .TRUE. ) THEN
    global%gfmGridY1 = vals(6)
  END IF ! defined  
  
  IF ( defined(7) .EQV. .TRUE. ) THEN
    global%gfmGridY2 = vals(7)
  END IF ! defined  
  
  IF ( defined(8) .EQV. .TRUE. ) THEN
    global%gfmGridZ1 = vals(8)
  END IF ! defined  
  
  IF ( defined(9) .EQV. .TRUE. ) THEN
    global%gfmGridZ2 = vals(9)
  END IF ! defined  
  
  IF ( defined(10) .EQV. .TRUE. ) THEN
    global%gfmGridNx = NINT(vals(10))
  END IF ! defined  
  
  IF ( defined(11) .EQV. .TRUE. ) THEN
    global%gfmGridNy = NINT(vals(11))
  END IF ! defined  
  
  IF ( defined(12) .EQV. .TRUE. ) THEN
    global%gfmGridNz = NINT(vals(12))
  END IF ! defined  
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE ReadGFMSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadGFMSection.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
!
!
! ******************************************************************************
