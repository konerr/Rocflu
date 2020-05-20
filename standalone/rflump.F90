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
! Purpose: Main driver of ROCFLU-MP.
!
! Description: None.
!
! Input: 
!   casename    Case name
!   verbLevel   Verbosity
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: rflump.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE rflump(caseString,verbLevel)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModParameters

  USE ModInterfaces, ONLY: RFLU_EndFlowSolver, &
                           RFLU_FlowSolver, & 
                           RFLU_InitFlowSolver

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(*) :: caseString
  INTEGER, INTENT(IN) :: verbLevel
  
! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: casename
  INTEGER :: dIterSystem,errorFlag
  REAL(RFREAL) :: dTimeSystem
  TYPE(t_level), POINTER :: levels(:)  
  TYPE(t_global), POINTER :: global
 
! ******************************************************************************
! Start, allocate global pointer
! ******************************************************************************

  ALLOCATE(global,STAT=errorFlag)
  IF ( errorFlag /= ERR_NONE ) THEN 
    WRITE(STDERR,'(A)') SOLVER_NAME//' ERROR - pointer allocation failed.'
    STOP
  END IF ! global

  casename = caseString(1:LEN(caseString))

! ******************************************************************************
! Call initialization, solver, and  finalization
! ******************************************************************************
  CALL RFLU_InitFlowSolver(casename,verbLevel,global,levels)

  dTimeSystem = MIN(global%maxTime - global%currentTime,global%runTime)
  dIterSystem = global%MaxIter - global%currentIter

  CALL RFLU_FlowSolver(dTimeSystem,dIterSystem,levels)
 
  CALL RFLU_EndFlowSolver(levels)
! ******************************************************************************
! Deallocate global pointer
! ******************************************************************************

  DEALLOCATE(global,STAT=errorFlag)
  IF ( errorFlag /= ERR_NONE ) THEN 
    WRITE(STDERR,'(A)') SOLVER_NAME//' ERROR - pointer deallocation failed.'
    STOP
  END IF ! global

! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE rflump

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: rflump.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:54  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:05  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/06/14 01:50:20  haselbac
! Now take into account global%runTime when calling RFLU_TimeStepping
!
! Revision 1.1  2007/04/09 18:51:36  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:02:03  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.5  2006/04/07 15:19:25  haselbac
! Removed tabs
!
! Revision 1.4  2005/09/14 15:59:33  haselbac
! Minor clean-up
!
! Revision 1.3  2005/09/13 20:46:01  mtcampbe
! Moved profiling calls to (Init)(End)Flowsolver
!
! Revision 1.2  2005/07/07 22:45:14  haselbac
! Added profiling calls
!
! Revision 1.1  2005/05/03 02:55:45  haselbac
! Initial revision
!
! ******************************************************************************

