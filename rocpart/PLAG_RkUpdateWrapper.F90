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
! Purpose: wrapper for update step in PLAG module.
!
! Description: none.
!
! Input: pRegion = data of current region,
!        iReg    = current region,
!        istage  = current RK stage
!
! Output: regions%levels%mixt = new solution after one time step.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_RkUpdateWrapper.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RkUpdateWrapper( region, iReg, istage )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  USE PLAG_ModInterfaces, ONLY: PLAG_InjcTileUpdate

  USE PLAG_ModInterfaces, ONLY: PLAG_RFLU_Update

  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  
  INTEGER, INTENT(IN) :: iReg, istage

! ... loop variables

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  TYPE(t_global), POINTER :: global

  TYPE(t_region), POINTER :: pRegion


!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RkUpdateWrapper.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction( global,'PLAG_rkUpdateWrapper',__FILE__ )

! update tiles ----------------------------------------------------------------
  CALL PLAG_InjcTileUpdate( region, iReg, iStage )

! update evolution equations --------------------------------------------------

!  WRITE(*,*) 'Entering PLAG_RFLU_Update: iReg, iStage = ',iReg, iStage
  pRegion => region%pRegion
  CALL PLAG_RFLU_Update( pRegion, iStage )

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RkUpdateWrapper

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RkUpdateWrapper.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:52  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/04/16 23:22:44  fnajjar
! Deleted all RFLO-pertinent calls and statements
!
! Revision 1.1  2007/04/09 18:50:27  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:34  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 20:58:17  fnajjar
! Initial revision after changing case
!
! Revision 1.4  2004/03/26 21:32:36  fnajjar
! Cleaned up routine for separate RFLO and RFLU calls within ifdef constructs and added PLAG_RFLU_Update call
!
! Revision 1.3  2004/03/08 22:26:29  fnajjar
! Removed ifdef RFLO around PLAG_InjcTileUpdate since injection is active with RFLU
!
! Revision 1.2  2004/02/26 21:02:22  haselbac
! Commented out RFLO-specific tile update routines
!
! Revision 1.1  2003/03/28 19:53:10  fnajjar
! Initial import of wrapper routines for RocfluidMP
!
!******************************************************************************

