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
! Purpose: initializes Runge-Kutta data for particles
!
! Description: none.
!
! Input: none.
!
! Output: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_RkInit.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RkInit( region,iStage )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartLag,    ONLY : t_plag
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  USE PLAG_ModInterfaces, ONLY: PLAG_InjcTileZeroRhs, &
                                PLAG_ZeroRhs
  USE PLAG_ModRkInit
  
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region

  INTEGER, INTENT(IN) :: iStage

! ... loop variables

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString


  TYPE(t_plag) ,  POINTER :: pPlag
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RkInit.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction( global,'PLAG_RkInit',__FILE__ )

  IF (global%plagUsed) THEN

    pPlag => region%plag

    IF ( pPlag%nPcls > 0 ) THEN
      CALL PLAG_ZeroRhs( region )
      CALL PLAG_RkInitPrimary( region, iStage )
    END IF ! nPcls

  END IF ! plagUsed

! Main algorithm for Tile evolution -------------------------------------------

! - Zero out residuals --------------------------------------------------------

  CALL PLAG_InjcTileZeroRhs( region )

! - Load old values -----------------------------------------------------------

  CALL PLAG_InjcTileRkInit( region, iStage )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RkInit

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RkInit.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
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
! Revision 1.3  2005/05/19 16:02:20  fnajjar
! Added calls to generic routines for initialization
!
! Revision 1.2  2005/04/22 18:34:32  fnajjar
! Removed IF statements to properly handle RK3 scheme
!
! Revision 1.1  2004/12/01 20:58:16  fnajjar
! Initial revision after changing case
!
! Revision 1.8  2004/03/22 23:49:28  fnajjar
! Removed RFLO-based ifdef for tile kernel
!
! Revision 1.7  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.6  2004/02/26 21:02:21  haselbac
! Commented out RFLO-specific tile init routines
!
! Revision 1.5  2004/02/06 21:18:20  fnajjar
! Initial Integration of Rocpart with Rocflu
!
! Revision 1.4  2003/11/12 21:33:07  fnajjar
! Moved FC computations to PLAG_RFLO_SetMetrics
!
! Revision 1.3  2003/11/03 22:39:25  fnajjar
! Added call to copy face vectors
!
! Revision 1.2  2003/03/28 19:51:50  fnajjar
! Included initialization step for tiles
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!
!******************************************************************************

