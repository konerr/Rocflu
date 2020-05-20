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
! Purpose: Dummy flow solver for GenX.
!
! Description: None.
!
! Input: 
!   globalGenx   	Pointer to global data
!   timeSystem		System time
!   dTimeSystem		System time step
!   genxHandleBc	Handle for BC update
!   genxHandleGm	Handle for geometry update.
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_FlowSolverDummy.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_FlowSolverDummy(globalGenx,timeSystem,dTimeSystem, & 
                                genxHandleBc,genxHandleGm)

  USE ModDataTypes
  USE ModGenx, ONLY: t_globalGenx
  USE ModDataStruct, ONLY: t_level
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  
  IMPLICIT NONE

! ... parameters
  INTEGER, INTENT(IN) :: genxHandleBc,genxHandleGm
  DOUBLE PRECISION, INTENT(IN) :: dTimeSystem,timeSystem
  TYPE(t_globalGenx), POINTER :: globalGenx

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_level), POINTER :: pLevel

!******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_FlowSolverDummy.F90,v $ $Revision: 1.1.1.1 $'

! initialize some global variables

  global => globalGenx%global
  pLevel => globalGenx%levels(1)

  global%genxHandleBc = genxHandleBc
  global%genxHandleGm = genxHandleGm
  global%dTimeSystem  = dTimeSystem

! start time stepping

  CALL RegisterFunction(global,'RFLU_FlowSolverDummy',__FILE__)

  CALL ErrorStop(global,ERR_EXTERNAL_FUNCT,__LINE__)

! finalize

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_FlowSolverDummy

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_FlowSolverDummy.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:30  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:47:51  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:58:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2002/10/17 19:55:55  haselbac
! Added timeSystem as second argument
!
! Revision 1.1  2002/10/05 18:28:18  haselbac
! Initial revision
!
!******************************************************************************

