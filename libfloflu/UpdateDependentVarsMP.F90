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
! Purpose: Update dependent variables.
!
! Description: None.
!
! Input:
!   region        Region data
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: UpdateDependentVarsMP.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE UpdateDependentVarsMP(region)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, &
                               RFLU_ConvertCvPrim2Cons

  USE ModInterfaces, ONLY: MixtureProperties

#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY: PLAG_NonCvUpdate
#endif

#ifdef SPEC
  USE ModInterfacesSpecies, ONLY: SPEC_UpdateDependentVars
#endif

  IMPLICIT NONE

! *****************************************************************************
! Declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region) :: region

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: ibc,iec
  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: UpdateDependentVarsMP.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction(global,'UpdateDependentVarsMP',__FILE__)

! *****************************************************************************
! Set variables
! *****************************************************************************

  ibc = 1
  iec = region%grid%nCellsTot

! *****************************************************************************
! Update dependent variables
! *****************************************************************************

! =============================================================================
! Mixture. NOTE here the last parameter MUST be TRUE, otherwise the gas
! properties do not get initialized correctly.
! =============================================================================

  CALL MixtureProperties(region,ibc,iec,.TRUE.)

#ifdef PLAG
! =============================================================================
! Particles
! =============================================================================

  IF ( global%plagUsed .EQV. .TRUE. ) THEN
    pRegion => region%pRegion
    CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)
    CALL PLAG_NonCvUpdate(pRegion)
    CALL RFLU_ConvertCvPrim2Cons(pRegion,CV_MIXT_STATE_CONS)
  END IF ! plagUsed
#endif

#ifdef SPEC
! =============================================================================
! Species
! =============================================================================

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    pRegion => region%pRegion
    CALL SPEC_UpdateDependentVars(pRegion,ibc,iec)
  END IF ! global%specUsed
#endif

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE UpdateDependentVarsMP

!******************************************************************************
!
! RCS Revision history:
!
! $Log: UpdateDependentVarsMP.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2005/04/15 15:06:04  haselbac
! Adapted interface to SPEC_UpdateDependentVars
!
! Revision 1.1  2004/12/01 16:51:59  haselbac
! Initial revision after changing case
!
! Revision 1.6  2004/11/29 17:15:20  wasistho
! use ModInterfacesSpecies
!
! Revision 1.5  2004/11/14 19:37:18  haselbac
! Changed interfaces to PLAG_nonCvUpdate and SPEC_UpdateDependentVars
!
! Revision 1.4  2004/07/28 15:29:18  jferry
! created global variable for spec use
!
! Revision 1.3  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.2  2004/02/26 21:01:50  haselbac
! Added PLAG support
!
! Revision 1.1  2004/01/29 22:52:33  haselbac
! Initial revision
!
!******************************************************************************

