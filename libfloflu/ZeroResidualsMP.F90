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
! Purpose: Set residuals to zero.
!
! Description: None.
!
! Input:
!   region        Data of current region.
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: ZeroResidualsMP.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ZeroResidualsMP(region)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global

#ifdef SPEC
  USE SPEC_RFLU_ModPBA, ONLY: SPEC_RFLU_PBA_ZeroResiduals
#endif

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), TARGET :: region

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: iVar
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => region%global

  CALL RegisterFunction(global,'ZeroResidualsMP',__FILE__)

! ******************************************************************************
! Zero residuals
! ******************************************************************************

  IF ( region%mixtInput%frozenFlag .EQV. .TRUE. ) THEN
    DO iVar = 1,region%mixtInput%nCv
      region%mixt%rhs(iVar,:) = 0.0_RFREAL
    END DO ! iVar
  END IF ! region%mixtInput%frozenFlag

#ifdef SPEC
! ==============================================================================
! Species
! ==============================================================================

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    DO iVar = 1,region%specInput%nSpecies
      IF ( region%specInput%specType(iVar)%frozenFlag .EQV. .TRUE. ) THEN
        region%spec%rhs(iVar,:) = 0.0_RFREAL
      END IF ! region%specInput%specType
    END DO ! iVar
  END IF ! global%specUsed
#endif

#ifdef SPEC
! DEBUG: Manoj-PBA1D, Notes: To match Jianghui's implementation
!                         1. Remove call to ZeroResiduals
!                         2. Not have frozenFlag = TRUE
! END DEBUG
! ==============================================================================
! Program burn
! ==============================================================================

  IF ( (global%specUsed .EQV. .TRUE.) .AND. (global%pbaFlag .EQV. .TRUE.) ) THEN
! DEBUG: Manoj-PBA1D
    CALL SPEC_RFLU_PBA_ZeroResiduals(region)
! END DEBUG
  END IF ! global%specUsed
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE ZeroResidualsMP

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ZeroResidualsMP.F90,v $
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
! Revision 1.1  2004/12/01 16:52:35  haselbac
! Initial revision after changing case
!
! Revision 1.3  2004/12/01 00:06:31  wasistho
! added StatBuildVersionString
!
! Revision 1.2  2004/07/28 15:29:18  jferry
! created global variable for spec use
!
! Revision 1.1  2003/11/25 21:01:50  haselbac
! Initial revision
!
!******************************************************************************

