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
! Purpose: Wrapper routine to get deformation.
!
! Description: None.
!
! Input: 
!   region     Region data
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_GetDeformationWrapper.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_GetDeformationWrapper(regions)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI
  USE ModParameters

#ifdef GENX
  USE ModInterfaces, ONLY: RFLU_GetDeformation
#else
  USE ModInterfaces, ONLY: RFLU_USER_GetDeformation
#endif

  IMPLICIT NONE

#ifdef GENX
  INCLUDE 'roccomf90.h'
#endif

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iReg
  TYPE(t_global), POINTER :: global

#ifdef GENX
  DOUBLE PRECISION :: dAlpha
#endif

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_GetDeformationWrapper.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_GetDeformationWrapper',__FILE__)

! *****************************************************************************
! Get deformation
! *****************************************************************************
  
#ifdef GENX 
  dAlpha = global%dtMin/global%dTimeSystem
  CALL COM_call_function(global%genxHandleGm,1,dAlpha)
  
  DO iReg = 1,global%nRegionsLocal
    CALL RFLU_GetDeformation(regions(iReg))
  END DO ! iReg  
#else
  DO iReg = 1,global%nRegionsLocal
    CALL RFLU_USER_GetDeformation(regions(iReg))
  END DO ! iReg
#endif     

! *****************************************************************************
! End
! *****************************************************************************
 
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_GetDeformationWrapper

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_GetDeformationWrapper.F90,v $
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
! Revision 1.1  2007/04/09 18:01:01  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2003/05/09 17:01:04  jiao
! Renamed the COM_call_function_handlers to COM_call_function.
!
! Revision 1.1  2003/03/15 19:06:12  haselbac
! Initial revision
!
!******************************************************************************

