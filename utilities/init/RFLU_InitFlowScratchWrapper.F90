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
! Purpose: Wrapper for initializing solution from scratch.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_InitFlowScratchWrapper.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_InitFlowScratchWrapper(pRegion)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  USE RFLU_ModConvertCv, ONLY: RFLU_ScalarConvertCvPrim2Cons

  USE ModInterfaces, ONLY: RFLU_InitFlowScratch, & 
                           RFLU_SetGasVars

#ifdef SPEC
  USE ModInterfacesSpecies, ONLY: SPEC_RFLU_InitFlowScratch
#endif

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Local variables
! ==============================================================================

  TYPE(t_global), POINTER :: global
  CHARACTER(CHRLEN) :: RCSIdentString

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_InitFlowScratchWrapper.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_InitFlowScratchWrapper',__FILE__)

#ifdef SPEC
! ******************************************************************************
! Species. NOTE needs to be done first so can define gas properties when these
! are influenced by species. NOTE initialize species state vector in primitive 
! state
! ******************************************************************************

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL SPEC_RFLU_InitFlowScratch(pRegion)
  END IF ! global%specUsed
#endif

! ******************************************************************************
! Mixture. NOTE if use pseudo-gas model, need to know mixture density so can 
! determine volume fraction of discrete species; but mixture density is only 
! known once have already initialized mixture. So initialize mixture first with
! perfect gas model, then initialize again with pseudo-gas model. 
! ******************************************************************************

  IF ( pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_PSEUDO ) THEN 
    pRegion%mixtInput%gasModel = GAS_MODEL_TCPERF
    
    CALL RFLU_SetGasVars(pRegion,1,pRegion%grid%nCellsTot)
    CALL RFLU_InitFlowScratch(pRegion) 
    
    pRegion%mixtInput%gasModel = GAS_MODEL_MIXT_PSEUDO    
  END IF ! pRegion%mixtInput%gasModel

  CALL RFLU_SetGasVars(pRegion,1,pRegion%grid%nCellsTot)
  CALL RFLU_InitFlowScratch(pRegion)

! ******************************************************************************
! Physical modules
! ******************************************************************************

! TO BE DONE

! END TO BE DONE

#ifdef SPEC
! ******************************************************************************
! Convert species state vector to conservative state
! ******************************************************************************

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL RFLU_ScalarConvertCvPrim2Cons(pRegion,pRegion%spec%cv, &
                                       pRegion%spec%cvState)
  END IF ! global%specUsed
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InitFlowScratchWrapper

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_InitFlowScratchWrapper.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:55  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:08  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:55:41  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2005/11/14 17:04:03  haselbac
! Added support for pseudo-gas model
!
! Revision 1.2  2005/11/10 02:43:30  haselbac
! Added support for variable properties
!
! Revision 1.1  2005/04/15 15:08:20  haselbac
! Initial revision
!
! ******************************************************************************

