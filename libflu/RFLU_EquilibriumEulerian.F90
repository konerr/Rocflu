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
! Purpose: Add Equilibrium Eulerian corrections
!
! Description: None.
!
! Input:
!   pRegion        Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_EquilibriumEulerian.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_EquilibriumEulerian(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
#ifdef SPEC
  USE ModSpecies, ONLY: t_spec_type
#endif

  USE ModInterfaces, ONLY: RFLU_FinishSD

#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY: PLAG_RFLU_GetMixtSD
#endif

#ifdef SPEC
  USE ModInterfacesSpecies, ONLY: SPEC_EqEulCorr, &
                                  SPEC_RFLU_SetEEv
#endif

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: iSpec
  TYPE(t_global), POINTER :: global
#ifdef SPEC
  TYPE(t_spec_type), POINTER :: pSpecType
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_EquilibriumEulerian',__FILE__)

! ******************************************************************************
! Call various subroutines to finish computation, if necessary
! ******************************************************************************

  IF ( pRegion%mixtInput%indSd == 1 ) THEN
    CALL RFLU_FinishSD(pRegion)

#ifdef PLAG
    IF ( global%plagUsed ) THEN
      IF ( pRegion%dummyStep .EQV. .FALSE. ) THEN
        CALL PLAG_RFLU_GetMixtSD(pRegion)
      END IF ! pRegion%dummyStep
    END IF ! global%plagUsed
#endif

#ifdef SPEC
    IF ( global%specUsed ) THEN
      DO iSpec = 1,pRegion%specInput%nSpecies
        pSpecType => pRegion%specInput%specType(iSpec)    
    
        IF ( pSpecType%velocityMethod == SPEC_METHV_EQEUL ) THEN
          CALL SPEC_RFLU_SetEEv(pRegion,iSpec)
          CALL SPEC_EqEulCorr(pRegion,iSpec)
        END IF ! pRegion%specInput
      END DO ! iSpec
    END IF ! global%specUsed
#endif
  END IF ! pRegion%mixtInput%indSd

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_EquilibriumEulerian

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_EquilibriumEulerian.F90,v $
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
! Revision 1.3  2005/11/27 01:48:31  haselbac
! Adapted to changes in EE implementation
!
! Revision 1.2  2004/08/02 13:57:28  haselbac
! Bug fix: Missing ifdefs, cosmetics
!
! Revision 1.1  2004/07/30 22:47:35  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! ******************************************************************************

