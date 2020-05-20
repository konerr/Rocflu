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
! Purpose: Update boundary values.
!
! Description: None.
!
! Input: 
!   region     Data of current region
!   istage     Current RK stage
!
! Output: 
!   None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_UpdateBoundaryValues.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_UpdateBoundaryValues(region,istage)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI
  USE ModParameters

#ifdef GENX
  USE ModInterfaces, ONLY: RFLU_GetBoundaryValues, &
                           RFLU_PutBoundaryValuesAlpha, & 
                           UpdateTbc
#endif

  IMPLICIT NONE

#ifdef GENX
  INCLUDE 'roccomf90.h'
#endif

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER :: istage
  TYPE(t_region) :: region  

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  REAL(RFREAL) :: subdt,time
  REAL(RFREAL) :: trk(5)
  TYPE(t_global), POINTER :: global

#ifdef GENX
  DOUBLE PRECISION :: alpha
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_UpdateBoundaryValues.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction(global,'RFLU_UpdateBoundaryValues',__FILE__)

  trk(:) = region%mixtInput%trk(:)

! ******************************************************************************
! Compute time
! ******************************************************************************

  time = global%currentTime + global%dtMin*trk(istage)

#ifdef GENX 
! ******************************************************************************
! Fill incoming buffers
! ******************************************************************************
   
  alpha = (time-global%timeStamp)/global%dTimeSystem

  CALL COM_call_function(global%genxHandleBc,2,alpha,1)
  CALL RFLU_PutBoundaryValuesAlpha(region)
  CALL COM_call_function(global%genxHandleBc,2,alpha,2)
  CALL RFLU_GetBoundaryValues(region)
#endif

! ******************************************************************************
! Fill in time-dependent boundary condition data 
! ******************************************************************************

  IF ( istage > 1 ) THEN 
    subdt = global%dtMin*(trk(istage) - trk(istage-1))
  ELSE
    subdt = 0.0_RFREAL
  END IF ! istage

  CALL UpdateTbc(region,time,subdt,istage==global%nrkSteps)  

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_UpdateBoundaryValues

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_UpdateBoundaryValues.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:48  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:57  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:02  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2005/04/15 15:07:28  haselbac
! Removed Charm/FEM stuff, cosmetics
!
! Revision 1.3  2004/01/22 16:04:34  haselbac
! Changed declaration to eliminate warning on ALC
!
! Revision 1.2  2003/11/25 21:07:44  haselbac
! Temporarily disable UpdateTbc - leads to core dump
!
! Revision 1.1  2003/10/03 20:48:52  haselbac
! Initial revision
!
! ******************************************************************************

