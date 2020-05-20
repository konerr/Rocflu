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
! Purpose: Compute integrals 1,2,4,5 of optimal LES approach.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to current region
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_ComputeIntegrals1245OLES.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_ComputeIntegrals1245OLES(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModMPI

  USE ModInterfaces, ONLY: RFLU_ComputeIntegral1OLES, & 
                           RFLU_ComputeIntegral2OLES, & 
                           RFLU_ComputeIntegral4OLES, & 
                           RFLU_ComputeIntegral5OLES

  IMPLICIT NONE

! --- parameters      

  TYPE(t_region), POINTER :: pRegion

! --- local variables

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Start
! ==============================================================================

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeIntegrals1245OLES',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Computing optimal LES '// & 
                             'integrals...'
  END IF ! global%verbLevel

! ==============================================================================
! Compute integrals
! ==============================================================================

  CALL RFLU_ComputeIntegral1OLES(pRegion)
  CALL RFLU_ComputeIntegral2OLES(pRegion)
  CALL RFLU_ComputeIntegral4OLES(pRegion)
  CALL RFLU_ComputeIntegral5OLES(pRegion)

! ==============================================================================
! End
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Computing optimal LES '// & 
                             'integrals done.'
  END IF ! global%verbLevel  

  CALL DeregisterFunction(global)     

END SUBROUTINE RFLU_ComputeIntegrals1245OLES

