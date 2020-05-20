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
! Purpose: Determine whether to build stencils and gradient weights.
!
! Description: None.
!
! Input:
!   global                              Pointer to global data
!
! Output: 
!   RFLU_DecideBuildStencilsWeights     TRUE or FALSE
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_DecideBuildStencilsWeights.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

LOGICAL FUNCTION RFLU_DecideBuildStencilsWeights(global)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
      
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

  CHARACTER(CHRLEN) :: RCSIdentString

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_DecideBuildStencilsWeights.F90,v $ $Revision: 1.1.1.1 $'

! ******************************************************************************
! Initialize
! ******************************************************************************

  RFLU_DecideBuildStencilsWeights = .FALSE.
  
! ******************************************************************************
! Determine whether should build stencils. NOTE need initFlowFlag here also 
! because this signals that errors can be computed...
! ******************************************************************************

  IF ( (global%postDiscFlag .EQV. .TRUE.) .OR. & 
       (global%postGradFlag .EQV. .TRUE.) .OR. & 
       (global%postVortFlag .EQV. .TRUE.) .OR. & 
       (global%postVortCoreFlag .EQV. .TRUE.) ) THEN
    RFLU_DecideBuildStencilsWeights = .TRUE.
  END IF ! global%postDiscFlag

! ******************************************************************************
! End
! ******************************************************************************

END FUNCTION RFLU_DecideBuildStencilsWeights

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_DecideBuildStencilsWeights.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:55  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:07  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/12/05 13:18:30  haselbac
! Initial revision
!
! Revision 1.1  2007/04/09 18:58:09  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.1  2006/01/06 22:02:58  haselbac
! Initial revision
!
! ******************************************************************************

