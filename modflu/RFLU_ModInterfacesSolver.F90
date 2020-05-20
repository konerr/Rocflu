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
! Purpose: Set explicit interfaces to subroutines and functions for Rocflu
!   solver.
!
! Description: None
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModInterfacesSolver.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************
  
MODULE RFLU_ModInterfacesSolver

  IMPLICIT NONE

  INTERFACE

  SUBROUTINE RFLU_CentralFirstPatch(pRegion,pPatch)
    USE ModBndPatch, ONLY: t_patch
    USE ModDataStruct, ONLY: t_region
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion   
  END SUBROUTINE RFLU_CentralFirstPatch

  SUBROUTINE RFLU_CentralFirstPatch_GL(pRegion,pPatch)
    USE ModBndPatch, ONLY: t_patch
    USE ModDataStruct, ONLY: t_region
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_CentralFirstPatch_GL

  SUBROUTINE RFLU_CentralSecondPatch(pRegion,pPatch)
    USE ModBndPatch, ONLY: t_patch
    USE ModDataStruct, ONLY: t_region
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion 
  END SUBROUTINE RFLU_CentralSecondPatch

  SUBROUTINE RFLU_CentralSecondPatch_GL(pRegion,pPatch)
    USE ModBndPatch, ONLY: t_patch
    USE ModDataStruct, ONLY: t_region
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_CentralSecondPatch_GL

  SUBROUTINE RFLU_CheckGridSpeeds(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_CheckGridSpeeds
  
  SUBROUTINE RFLU_ComputeEnerDissOLES(region)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLU_ComputeEnerDissOLES

  SUBROUTINE RFLU_ComputeFluxInv(pRegion,fluxPart)
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: fluxPart    
    TYPE(t_region), POINTER :: pRegion     
  END SUBROUTINE RFLU_ComputeFluxInv

  SUBROUTINE RFLU_ComputeGridSpeeds(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_ComputeGridSpeeds

  SUBROUTINE RFLU_ComputeIntegral1OLES(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion     
  END SUBROUTINE RFLU_ComputeIntegral1OLES

  SUBROUTINE RFLU_ComputeIntegral2OLES(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion     
  END SUBROUTINE RFLU_ComputeIntegral2OLES

  SUBROUTINE RFLU_ComputeIntegral3OLES(region)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLU_ComputeIntegral3OLES

  SUBROUTINE RFLU_ComputeIntegral4OLES(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion     
  END SUBROUTINE RFLU_ComputeIntegral4OLES

  SUBROUTINE RFLU_ComputeIntegral5OLES(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion     
  END SUBROUTINE RFLU_ComputeIntegral5OLES

#ifndef GENX
  SUBROUTINE RFLU_ComputeIntegralValues(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLU_ComputeIntegralValues
#else
  SUBROUTINE RFLU_ComputeIntegralValues(regions,integ)
    USE ModGenx ! To access MAN_INTEG_SIZE
    USE ModDataStruct, ONLY: t_region
    DOUBLE PRECISION, DIMENSION(MAN_INTEG_SIZE) :: integ
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLU_ComputeIntegralValues
#endif

  SUBROUTINE RFLU_ComputeIntegrals1245OLES(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion     
  END SUBROUTINE RFLU_ComputeIntegrals1245OLES

  SUBROUTINE RFLU_ComputeWeightsOLES(region)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLU_ComputeWeightsOLES

  SUBROUTINE RFLU_ConvFluxOLES(region)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region     
  END SUBROUTINE RFLU_ConvFluxOLES

  SUBROUTINE RFLU_CreateFields(global,pLevel)
    USE ModGlobal, ONLY: t_global
    USE ModDataStruct, ONLY: t_level
    TYPE(t_global), POINTER :: global
    TYPE(t_level), POINTER :: pLevel
  END SUBROUTINE RFLU_CreateFields

  SUBROUTINE RFLU_EndFlowSolver(levels)
    USE ModDataStruct, ONLY: t_level
    TYPE(t_level), POINTER :: levels(:)
  END SUBROUTINE RFLU_EndFlowSolver

  SUBROUTINE RFLU_EquilibriumEulerian(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_EquilibriumEulerian

  SUBROUTINE RFLU_ExplicitMultiStage(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLU_ExplicitMultiStage

  SUBROUTINE RFLU_FinishSD(region)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region     
  END SUBROUTINE RFLU_FinishSD

#ifndef GENX
  SUBROUTINE RFLU_FlowSolver(dTimeSystem,dIterSystem,levels)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_level
    INTEGER, INTENT(IN) :: dIterSystem  
    REAL(RFREAL), INTENT(IN) :: dTimeSystem
    TYPE(t_level), POINTER :: levels(:)
  END SUBROUTINE RFLU_FlowSolver
#else
  SUBROUTINE RFLU_FlowSolver(globalGenx,timeSystem,dTimeSystem,genxHandleBc, & 
                             genxHandleGm)
    USE ModDataTypes
    USE ModGenx, ONLY: t_globalGenx 
    INTEGER, INTENT(IN) :: genxHandleBc,genxHandleGm
    DOUBLE PRECISION, INTENT(IN) :: dTimeSystem,timeSystem
    TYPE(t_globalGenx), POINTER :: globalGenx         
  END SUBROUTINE RFLU_FlowSolver
#endif

  SUBROUTINE RFLU_GetDeformationWrapper(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), DIMENSION(:), POINTER :: regions   
  END SUBROUTINE RFLU_GetDeformationWrapper

#ifndef GENX
  SUBROUTINE RFLU_InitFlowSolver(casename,verbLevel,global,levels)
    USE ModDataTypes  
    USE ModGlobal, ONLY: t_global
    USE ModDataStruct, ONLY: t_level
    CHARACTER(CHRLEN), INTENT(IN) :: casename
    INTEGER, INTENT(IN) :: verbLevel
    TYPE(t_global), POINTER :: global
    TYPE(t_level), POINTER :: levels(:)
  END SUBROUTINE RFLU_InitFlowSolver  
#else
  SUBROUTINE RFLU_InitFlowSolver(globalGenx,initialTime,communicator,genxHandle)
    USE ModGenx, ONLY: t_globalGenx    
    INTEGER, INTENT(IN) :: communicator,genxHandle
    DOUBLE PRECISION, INTENT(IN) :: initialTime
    TYPE(t_globalGenx), POINTER :: globalGenx
  END SUBROUTINE RFLU_InitFlowSolver  
#endif

  SUBROUTINE RFLU_InvFlux_I(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE RFLU_InvFlux_I 

  SUBROUTINE RFLU_MinimumTimeStep(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLU_MinimumTimeStep
  
  SUBROUTINE RFLU_MinimumTimeStepImpulse(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLU_MinimumTimeStepImpulse
  
  SUBROUTINE RFLU_MoveGridDisp(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLU_MoveGridDisp  

  SUBROUTINE RFLU_MoveGridWrapper(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLU_MoveGridWrapper

  SUBROUTINE RFLU_MoveGridXyz(regions,context)
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: context
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLU_MoveGridXyz

  SUBROUTINE RFLU_OpenConverFile(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLU_OpenConverFile

! begin BBR
  SUBROUTINE RFLU_OpenPMFile(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLU_OpenPMFile

  SUBROUTINE RFLU_OpenIntegFile(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLU_OpenIntegFile
! end BBR

  SUBROUTINE RFLU_OpenStatsFileOLES(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLU_OpenStatsFileOLES

  SUBROUTINE RFLU_OpenTotalMassFile(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLU_OpenTotalMassFile

  SUBROUTINE RFLU_PrintWriteConvergence(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLU_PrintWriteConvergence

! begin BBR
  SUBROUTINE RFLU_WritePM(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLU_WritePM

  SUBROUTINE RFLU_WriteInteg(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLU_WriteInteg
! end BBR

  SUBROUTINE RFLU_ReadIntegrals1245OLES(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_ReadIntegrals1245OLES

  SUBROUTINE RFLU_ResidualNorm(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLU_ResidualNorm
  
  SUBROUTINE RFLU_SetMoveGridOptions(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_SetMoveGridOptions

  SUBROUTINE RFLU_TimeStepImpulse(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_TimeStepImpulse

  SUBROUTINE RFLU_TimeStepInviscid(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_TimeStepInviscid

  SUBROUTINE RFLU_TimeStepViscous(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_TimeStepViscous

  SUBROUTINE RFLU_TimeStepping(dTimeSystem,dIterSystem,regions)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: dIterSystem  
    REAL(RFREAL), INTENT(IN) :: dTimeSystem
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLU_TimeStepping

  SUBROUTINE RFLU_USER_GetDeformation(region)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region   
  END SUBROUTINE RFLU_USER_GetDeformation
  
  SUBROUTINE RFLU_UpdateBoundaryValues(region,istage)
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: istage
    TYPE(t_region) :: region
  END SUBROUTINE RFLU_UpdateBoundaryValues
  
  SUBROUTINE RFLU_WriteIntegrals1245OLES(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_WriteIntegrals1245OLES

  SUBROUTINE RFLU_WriteStatsFileOLES(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLU_WriteStatsFileOLES

  END INTERFACE

END MODULE RFLU_ModInterfacesSolver

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModInterfacesSolver.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:41  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:41  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.19  2006/05/01 21:01:52  haselbac
! Added if for RFLU_ComputeFluxInv
!
! Revision 1.18  2006/03/26 20:22:05  haselbac
! Added ifs for GL routines
!
! Revision 1.17  2005/05/16 20:43:11  haselbac
! Removed interfaces for flux routines, adapted interfaces for patch flux routines
!
! Revision 1.16  2005/04/29 12:55:54  haselbac
! Removed interfaces for probe routines
!
! Revision 1.15  2005/04/15 16:30:47  haselbac
! Added interface for RFLU_PrintWriteConvergence
!
! Revision 1.14  2004/12/19 15:47:59  haselbac
! Added interface for RFLU_InvFlux_I
!
! Revision 1.13  2004/11/03 17:03:21  haselbac
! Removed interface for RFLU_HACK_UpdateVirtualCells
!
! Revision 1.12  2004/10/19 19:28:10  haselbac
! Modified interfaces for time-step routines
!
! Revision 1.11  2004/07/30 22:47:36  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.10  2004/07/06 15:14:44  haselbac
! Cosmetics only
!
! Revision 1.9  2004/02/26 21:02:04  haselbac
! Cosmetics only
!
! Revision 1.8  2004/02/13 02:58:03  haselbac
! Added interface for fast Roe routines
!
! Revision 1.7  2004/01/29 22:57:37  haselbac
! Added ifs for RFLU_RoeFirst and RFLU_RoeSecond.F90
!
! Revision 1.6  2003/12/04 03:28:50  haselbac
! Added entries for new routines, removed RFLU_UpdateCommLists (bug fix)
!
! Revision 1.5  2003/11/25 21:03:26  haselbac
! Added interface for RFLU_UpdateDummyCells
!
! Revision 1.4  2003/10/03 20:45:47  haselbac
! Added interface for RFLU_UpdateBoundaryValues
!
! Revision 1.3  2003/06/04 22:11:13  haselbac
! Added RFLU_CentralFirstPatch, rm RFLU_RoeCentralFirstPatch
!
! Revision 1.2  2003/05/16 22:09:29  haselbac
! Added RFLU_HLLCFirst
!
! Revision 1.1  2003/04/10 14:37:10  haselbac
! Initial revision
!
! ******************************************************************************

