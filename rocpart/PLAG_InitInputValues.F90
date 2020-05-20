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
! Purpose: initialize user input parameters for Lagrangians particles
!          to default values.
!
! Description: none.
!
! Input: none.
!
! Output: regions = initial input values.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_InitInputValues.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_InitInputValues( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartLag,    ONLY : t_plag_input
  USE ModError
  USE ModParameters
  USE ModMaterials, ONLY  : t_material
  USE PLAG_ModParameters
  USE PLAG_ModBcData, ONLY : PLAG_InitBcData
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  TYPE(t_plag_input), POINTER :: pPlagInput
  TYPE(t_region), POINTER :: pRegion
  TYPE(t_global),     POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InitInputValues.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PLAG_InitInputValues',__FILE__ )

! global values ---------------------------------------------------------------

! none currently

! region related values -------------------------------------------------------

  DO iReg=LBOUND(regions,1),UBOUND(regions,1)

    pRegion    => regions(iReg)
    pPlagInput => regions(iReg)%plagInput

! plagInput quantities

    !pPlagInput%nPclsMax          = 1000
    ! Subbu - Debug
    pPlagInput%nPclsMax          = 50000
    ! Subbu - End Debug
    pPlagInput%intrplMixtModel   = ZEROTH_ORDER
    pPlagInput%nCont             = 1
    pPlagInput%breakupModel      = PLAG_BREAKUP_NOMODEL
    pPlagInput%breakupFac        = 1.0_RFREAL
    pPlagInput%breakupWebSwi     = PLAG_BREAKUP_NOWEBSWI
    pPlagInput%readStatus        = -1 ! not read
    pPlagInput%findPclMethod     = FIND_PCL_METHOD_TRAJ_SAFE
    pPlagInput%initDimens        = 3
    pPlagInput%nUnsteadyData     = 50
    pPlagInput%nTimeBH           = 1 ! By default there is no usable history data

! index for volume fraction array 
    pPlagInput%indVFracE = 0
    
! inter-particle collision

    pPlagInput%flagCollision = .FALSE.
    pPlagInput%collisionPs   = 1.0E5

! continuous random walk for particle velocity fluctuation

    pPlagInput%flagCRW = .FALSE.

! stability criteria based on fluid/particle momentum

    pPlagInput%flagStability = .FALSE.
    pPlagInput%limitForce    = .FALSE.
    pPlagInput%cfl           = 1.0

! initial quantities

    pPlagInput%nPclsIni          = 0

! initialize patch data for boundary conditions

    CALL PLAG_InitBcData(pRegion)

  ENDDO  ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_InitInputValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InitInputValues.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.5  2008/12/06 08:43:50  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/05/16 22:27:22  fnajjar
! Deleted initialization for old datastructure and added call for PLAG_InitBcData
!
! Revision 1.2  2007/04/15 02:35:45  haselbac
! Added setting of initDimens
!
! Revision 1.1  2007/04/09 18:50:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2007/03/06 23:13:13  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.2  2005/03/11 02:22:33  haselbac
! Changed default for tracking to safe method
!
! Revision 1.1  2004/12/01 20:57:36  fnajjar
! Initial revision after changing case
!
! Revision 1.14  2004/10/11 22:12:23  haselbac
! Bug fix
!
! Revision 1.13  2004/10/09 16:37:37  fnajjar
! Removed initialization of initFlag
!
! Revision 1.12  2004/10/08 22:12:43  haselbac
! Added initialization of findPclMethod
!
! Revision 1.11  2004/08/20 23:27:13  fnajjar
! Added Infrastructure for Plag prep tool
!
! Revision 1.10  2004/06/17 15:19:03  fnajjar
! Added infrastructure for ejection model
!
! Revision 1.9  2004/06/16 23:03:49  fnajjar
! Renamed variabled for CRE kernel
!
! Revision 1.8  2004/03/10 23:09:50  fnajjar
! Added maximum buffer size for corner-edge cells
!
! Revision 1.7  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.6  2004/02/26 21:02:14  haselbac
! Changed loop limits for generality
!
! Revision 1.5  2004/02/16 23:31:11  fnajjar
! Included default values for injcDiamMin and injcDiamMax
!
! Revision 1.4  2003/11/21 22:42:16  fnajjar
! Added plagActive
!
! Revision 1.3  2003/09/13 20:14:21  fnajjar
! Added infrastructure for Breakup model
!
! Revision 1.2  2003/09/10 23:35:50  fnajjar
! Removed flags that are subsumed with Rocinteract
!
! Revision 1.1  2003/04/14 14:33:15  fnajjar
! Initial import for proper input initialization
!
!******************************************************************************

