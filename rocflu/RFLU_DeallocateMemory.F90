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
! Purpose: Deallocate memory for mixture.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_DeallocateMemory.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_DeallocateMemory(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  USE RFLU_ModDeallocateMemory, ONLY: RFLU_DeallocateMemoryGSpeeds, & 
                                      RFLU_DeallocateMemorySol, & 
                                      RFLU_DeallocateMemoryTStep, & 
                                      RFLU_DeallocateMemoryTStep_C, & 
                                      RFLU_DeallocateMemoryTStep_I 

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
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

  RCSIdentString = &
    '$RCSfile: RFLU_DeallocateMemory.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_DeallocateMemory',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Deallocating memory for mixture...'
  END IF ! global%verbLevel

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%fluidModel )

! ==============================================================================
!   Compressible fluid 
! ==============================================================================

    CASE ( FLUID_MODEL_COMP ) 
      CALL RFLU_DeallocateMemorySol(pRegion)
      CALL RFLU_DeallocateMemoryTStep(pRegion)
      CALL RFLU_DeallocateMemoryTStep_C(pRegion)
      CALL RFLU_DeallocateMemoryGSpeeds(pRegion)

! ==============================================================================
!   Incompressible fluid 
! ==============================================================================

    CASE ( FLUID_MODEL_INCOMP )
      CALL RFLU_DeallocateMemorySol(pRegion)
      CALL RFLU_DeallocateMemoryTStep(pRegion)      
      CALL RFLU_DeallocateMemoryTStep_I(pRegion)
      CALL RFLU_DeallocateMemoryGSpeeds(pRegion)

! ==============================================================================
!   Default
! ==============================================================================

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixtInput%fluidModel
  
! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Deallocating memory for mixture done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_DeallocateMemory

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_DeallocateMemory.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:47  mtcampbe
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
! Revision 1.12  2006/03/26 20:22:14  haselbac
! Removed error trap for GL model
!
! Revision 1.11  2004/12/19 15:49:11  haselbac
! Modified so can select different fluid models
!
! Revision 1.10  2004/03/19 21:21:06  haselbac
! Complete rewrite
!
! Revision 1.9  2004/03/03 23:55:40  jferry
! Allowed particles to be run with Euler case
!
! Revision 1.8  2004/01/29 23:04:11  haselbac
! Added/deleted deallocation for mfMixt/vfMixt arrays, clean-up
!
! Revision 1.7  2003/12/04 03:30:00  haselbac
! Added memory deallocation for gradients, cleaned up
!
! Revision 1.6  2003/11/25 21:04:38  haselbac
! Added deallocation for mass and volume fluxes on patches
!
! Revision 1.5  2003/11/03 03:50:33  haselbac
! Removed deallocation of bf2bg list
!
! Revision 1.4  2003/03/31 16:16:49  haselbac
! Added disp array, some cosmetics
!
! Revision 1.3  2003/03/25 19:16:09  haselbac
! Fixed typo in output string
!
! Revision 1.2  2003/03/15 18:34:31  haselbac
! Adaptation for parallel computations
!
! Revision 1.1  2003/01/28 14:59:34  haselbac
! Initial revision
!
! ******************************************************************************

