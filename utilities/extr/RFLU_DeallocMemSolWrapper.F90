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
! Purpose: Deallocate memory wrapper for rflupost.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.

! ******************************************************************************
!
! $Id: RFLU_DeallocMemSolWrapper.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_DeallocMemSolWrapper(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModParameters
  USE ModMPI

#ifdef PLAG
  USE ModPartLag, ONLY: t_plag
#endif

  USE RFLU_ModDeallocateMemory, ONLY: RFLU_DeallocateMemoryAuxVars, &
                                      RFLU_DeallocateMemorySolCv, & 
                                      RFLU_DeallocateMemorySolDv, &
                                      RFLU_DeallocateMemorySolGv, & 
                                      RFLU_DeallocateMemorySolTv

#ifdef PLAG
  USE RFLU_ModDeallocateMemory, ONLY: RFLU_DeallocateMemorySolTv
  USE ModInterfacesLagrangian, ONLY: PLAG_RFLU_DeallocMemSol, &
                                     PLAG_RFLU_DeallocMemSolTile
#endif

#ifdef SPEC
  USE ModInterfacesSpecies, ONLY: SPEC_RFLU_DeallocateMemoryEEv, & 
                                  SPEC_RFLU_DeallocateMemorySol
#endif

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Parameters
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

#ifdef PLAG
  TYPE(t_plag), POINTER :: pPlag
#endif

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = &
    '$RCSfile: RFLU_DeallocMemSolWrapper.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocMemSolWrapper', &
                        __FILE__)

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Deallocating memory...'
  END IF ! global%verbLevel

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  CALL RFLU_DeallocateMemorySolCv(pRegion)

  IF ( global%solverType == SOLV_IMPLICIT_HM .AND. &
       pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN
    CALL RFLU_DeallocateMemoryAuxVars(pRegion)
  END IF ! global%solverType

  CALL RFLU_DeallocateMemorySolDv(pRegion)
  CALL RFLU_DeallocateMemorySolGv(pRegion)  

  IF ( pRegion%mixtInput%computeTv .EQV. .TRUE. ) THEN 
    CALL RFLU_DeallocateMemorySolTv(pRegion)
  END IF ! pRegion%mixtInput%computeTv

#ifdef PLAG
  IF ( global%plagUsed .EQV. .TRUE. ) THEN
    pPlag => pRegion%plag

    CALL PLAG_RFLU_DeallocMemSol(pRegion,pPlag)
    CALL PLAG_RFLU_DeallocMemSolTile(pRegion)
  END IF ! plagUsed
#endif

#ifdef SPEC
  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL SPEC_RFLU_DeallocateMemorySol(pRegion)

    IF ( pRegion%specInput%nSpeciesEE > 0 ) THEN 
      CALL SPEC_RFLU_DeallocateMemoryEEv(pRegion)
    END IF ! pRegion%specInput%nSpeciesEE
  END IF ! global%specUsed
#endif

! *****************************************************************************
! End
! *****************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Deallocating memory done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocMemSolWrapper

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_DeallocMemSolWrapper.F90,v $
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
! Revision 1.2  2007/11/28 23:05:47  mparmar
! Added deallocation of aux vars
!
! Revision 1.1  2007/04/09 18:58:08  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.8  2005/11/27 01:58:43  haselbac
! Added deallocation of EEv
!
! Revision 1.7  2005/11/17 14:49:39  haselbac
! Bug fix: Tv needed by Rocspecies with EE
!
! Revision 1.6  2005/11/11 23:45:59  fnajjar
! Bug fix: Added dealloc of tv for PLAG
!
! Revision 1.5  2005/11/10 02:45:54  haselbac
! Change logic bcos of variable properties
!
! Revision 1.4  2004/07/28 15:29:21  jferry
! created global variable for spec use
!
! Revision 1.3  2004/03/19 21:32:14  haselbac
! Changed memory deallocation logic
!
! Revision 1.2  2004/03/05 22:09:05  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.1  2004/02/26 21:01:22  haselbac
! Initial revision
!
! ******************************************************************************

