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
! Purpose: update step for injection tiles.
!
! Description: none.
!
! Input: region = current region
!        iReg   = current region number
!        iStage = current RK stage.
!
! Output: region%plag = plag variables
!
! Notes: This corresponds to Part I Step 2 in RocfluidMP framework.
!
!******************************************************************************
!
! $Id: PLAG_InjcTileUpdate.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_InjcTileUpdate( region, iReg, iStage )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMixture, ONLY    : t_mixt_input
  USE ModPartLag, ONLY    : t_plag
  USE ModError
  USE ModParameters

  USE PLAG_ModParameters 
  
  USE PLAG_ModInterfaces, ONLY: PLAG_InjcTileRKUpdate        
  USE PLAG_ModInterfaces, ONLY: PLAG_RFLU_InjcTileCalcRhs

  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  
  INTEGER :: iReg, iStage
  
! ... loop variables

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InjcTileUpdate.F90,v $ $Revision: 1.1.1.1 $'
  
  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_InjcTileUpdate',__FILE__ )

! - Calculate rhs for tile ----------------------------------------------------

  pRegion => region%pRegion
!  WRITE(STDOUT,'(A)') '    Entering PLAG_RFLU_InjcTileCalcRhs'
  CALL PLAG_RFLU_InjcTileCalcRhs( pRegion )

! - Invoke RK update ----------------------------------------------------------

!  WRITE(STDOUT,'(A)') '    Entering PLAG_InjcTileRKUpdate'
  CALL PLAG_InjcTileRKUpdate( region, iStage )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_InjcTileUpdate

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InjcTileUpdate.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:50  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/04/16 23:20:36  fnajjar
! Deleted all RFLO-pertinent calls and statements
!
! Revision 1.1  2007/04/09 18:50:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 20:57:48  fnajjar
! Initial revision after changing case
!
! Revision 1.6  2004/03/08 22:24:39  fnajjar
! Modified routine to be RFLU-aware and added PLAG_RFLU_InjcTileCalcRhs call
!
! Revision 1.5  2004/02/25 21:56:15  fnajjar
! Moved tile pointers outside do-loop
!
! Revision 1.4  2003/05/01 22:51:12  fnajjar
! Removed PLAG_CalcFaceCentroids in PLAG_ModInterfaces list
!
! Revision 1.3  2003/04/15 23:02:38  fnajjar
! Removed dead code section
!
! Revision 1.2  2003/03/28 19:52:22  fnajjar
! Removed initialization step for tiles
!
! Revision 1.1  2003/02/04 19:09:46  f-najjar
! Initial Import
!
!******************************************************************************

