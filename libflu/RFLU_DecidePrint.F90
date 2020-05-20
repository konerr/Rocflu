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
! Purpose: Determine whether to print data to screen.
!
! Description: None.
!
! Input:
!   global                      Pointer to global data
!
! Output: 
!   RFLU_DecidePrint = .TRUE.   If should print to screen
!   RFLU_DecidePrint = .FALSE.  If should not print to screen
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_DecidePrint.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
!******************************************************************************

LOGICAL FUNCTION RFLU_DecidePrint(global)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
      
  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_global), POINTER :: global

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  LOGICAL :: logical1,logical2,logical3

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_DecidePrint.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction(global,'RFLU_DecidePrint',__FILE__)

! *****************************************************************************
! Initialize
! *****************************************************************************

  RFLU_DecidePrint = .FALSE.
  
! *****************************************************************************
! Determine whether should print to screen
! *****************************************************************************

! =============================================================================
! Unsteady flow
! =============================================================================

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    logical1 = ABS(global%timeSincePrint-global%printTime) & 
             < 0.1_RFREAL*global%dtMin             
    logical2 = (global%timeSincePrint > global%printTime)    
    logical3 = (global%iterSinceRestart == 1)
  
    IF ( logical1 .OR. logical2 .OR. logical3 ) THEN
      RFLU_DecidePrint = .TRUE.
    END IF ! logical1 

! =============================================================================
! Steady flow. NOTE must have RFLU_DecidePrint be TRUE if currentIter == 1 
! because only then will RFLU_ResidualNorm be called from RFLU_TimeStepping, 
! and at the first iteration, this is needed to get a residual for correct 
! normalization.
! =============================================================================

  ELSE    
    RFLU_DecidePrint = (MOD(global%currentIter,global%printIter) == 0)

    IF ( global%currentIter == 1 ) THEN 
      RFLU_DecidePrint = .TRUE. 
    END IF ! global%currentIter
  END IF ! global%flowType

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END FUNCTION RFLU_DecidePrint

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_DecidePrint.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:35  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:52  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.3  2004/01/31 03:56:54  haselbac
! Added explanatory comment
!
! Revision 1.2  2004/01/29 22:56:28  haselbac
! Removed setting of timeSincePrint to make routine more usable
!
! Revision 1.1  2003/10/29 21:36:58  haselbac
! Initial revision
!
!******************************************************************************

