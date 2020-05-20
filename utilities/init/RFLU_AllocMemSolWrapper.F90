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
! Purpose: Allocate memory wrapper.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: 
!   1. Do not allocate memory for particles because treated separately.
!
! ******************************************************************************
!
! $Id: RFLU_AllocMemSolWrapper.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_AllocMemSolWrapper(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  USE RFLU_ModAllocateMemory, ONLY: RFLU_AllocateMemoryAuxVars, &
                                    RFLU_AllocateMemorySolCv, & 
                                    RFLU_AllocateMemorySolDv, & 
                                    RFLU_AllocateMemorySolGv

#ifdef SPEC
  USE ModInterfacesSpecies, ONLY: SPEC_RFLU_AllocateMemoryEEv, &
                                  SPEC_RFLU_AllocateMemorySol
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

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_AllocMemSolWrapper.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AllocMemSolWrapper', &
                        __FILE__)

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Allocating memory...'
  END IF ! global%verbLevel

! ******************************************************************************
! Allocate memory
! ******************************************************************************

! ==============================================================================
! Mixture
! ==============================================================================

  CALL RFLU_AllocateMemorySolCv(pRegion)
  CALL RFLU_AllocateMemorySolGv(pRegion)  
  CALL RFLU_AllocateMemorySolDv(pRegion)

  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
    CALL RFLU_AllocateMemoryAuxVars(pRegion)
  END IF ! solverType

#ifdef SPEC
! ==============================================================================
! Species
! ==============================================================================

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL SPEC_RFLU_AllocateMemorySol(pRegion)
    
    IF ( pRegion%specInput%nSpeciesEE > 0 ) THEN 
      CALL SPEC_RFLU_AllocateMemoryEEv(pRegion)
    END IF ! pRegion%specInput%nSpeciesEE    
  END IF ! global%specUsed
#endif

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Allocating memory done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AllocMemSolWrapper

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_AllocMemSolWrapper.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:55  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:07  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/11/28 23:05:36  mparmar
! Allocating for auxilliary vars
!
! Revision 1.1  2007/04/09 18:55:41  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/08/24 02:14:27  haselbac
! Bug fix: Allocate dv bcos of NSCBC code
!
! Revision 1.3  2005/11/27 01:57:47  haselbac
! Added allocation of EEv
!
! Revision 1.2  2005/11/10 02:41:02  haselbac
! Now alloc gv outside GENX bcos of variable properties
!
! Revision 1.1  2005/04/15 15:08:11  haselbac
! Initial revision
!
! ******************************************************************************

