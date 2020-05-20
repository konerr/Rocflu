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
! Purpose: set explicit interfaces to subroutines and functions.
!
! Description: none
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ModInterfacesLagrangian.F90,v 1.2 2015/07/27 04:45:42 brollin Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModInterfacesLagrangian

  IMPLICIT NONE

  INTERFACE

! =============================================================================
! Common routines
! =============================================================================

  SUBROUTINE PLAG_BuildVersionString( versionString )
    CHARACTER(*) :: versionString
  END SUBROUTINE PLAG_BuildVersionString

  SUBROUTINE PLAG_InitPatchData(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_InitPatchData
 
  SUBROUTINE PLAG_InitSolution( iReg,region )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iReg
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_InitSolution

  SUBROUTINE PLAG_NonCvUpdate( region )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region), POINTER :: region
  END SUBROUTINE PLAG_NonCvUpdate
  
  SUBROUTINE PLAG_PrintUserInput( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_PrintUserInput

  SUBROUTINE PLAG_StatMapping( global )
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE PLAG_StatMapping

  SUBROUTINE PLAG_UserInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PLAG_UserInput

  SUBROUTINE PLAG_RkInit( region,iStage )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region), INTENT(INOUT) :: region
    INTEGER,        INTENT(IN)    :: iStage
  END SUBROUTINE PLAG_RkInit
  
  SUBROUTINE PLAG_RkUpdateWrapper( region, iReg, iStage )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region)      :: region
    INTEGER, INTENT(IN) :: iReg, iStage
  END SUBROUTINE PLAG_RkUpdateWrapper

! =============================================================================
! Rocflu-specific routines
! =============================================================================

  SUBROUTINE PLAG_INRT_AllocMemTStep(pRegion,pPlag)
    USE ModDataStruct, ONLY: t_region
    USE ModPartLag, ONLY: t_plag
    TYPE(t_region), POINTER :: pRegion 
    TYPE(t_plag), POINTER :: pPlag     
  END SUBROUTINE PLAG_INRT_AllocMemTStep

  SUBROUTINE PLAG_INRT_DeallocMemTStep(pRegion,pPlag)
    USE ModDataStruct, ONLY: t_region
    USE ModPartLag, ONLY: t_plag
    TYPE(t_region), POINTER :: pRegion 
    TYPE(t_plag), POINTER :: pPlag     
  END SUBROUTINE PLAG_INRT_DeallocMemTStep

  SUBROUTINE PLAG_RFLU_AllocMemSol(pRegion,pPlag)
    USE ModDataStruct, ONLY: t_region
    USE ModPartLag, ONLY: t_plag
    TYPE(t_region), POINTER :: pRegion 
    TYPE(t_plag), POINTER :: pPlag     
  END SUBROUTINE PLAG_RFLU_AllocMemSol

  SUBROUTINE PLAG_RFLU_AllocMemSolTile(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_AllocMemSolTile

  SUBROUTINE PLAG_RFLU_AllocMemTStep(pRegion,pPlag)
    USE ModDataStruct, ONLY: t_region
    USE ModPartLag, ONLY: t_plag
    TYPE(t_region), POINTER :: pRegion 
    TYPE(t_plag), POINTER :: pPlag    
  END SUBROUTINE PLAG_RFLU_AllocMemTStep

  SUBROUTINE PLAG_RFLU_AllocMemTStepTile(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion 
  END SUBROUTINE PLAG_RFLU_AllocMemTStepTile  

  ! Subbu - Add routine to compute ncells that can
  ! be seeded with pcls
  SUBROUTINE PLAG_RFLU_ComputeCellsContainingPcls(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE PLAG_RFLU_ComputeCellsContainingPcls
  ! Subbu - End add routine to compute ncells that
  ! can be seeded with pcls

  ! Subbu - Add routine to compute number of nPcls per Reg
  SUBROUTINE PLAG_RFLU_ComputeLagPclsPerReg(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE PLAG_RFLU_ComputeLagPclsPerReg
  ! Subbu - End add routine to compute nPcls per Region
  
  SUBROUTINE PLAG_RFLU_ComputeMaxImpulse(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_ComputeMaxImpulse
  
  SUBROUTINE PLAG_RFLU_ComputeVolFrac(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_ComputeVolFrac
  
  SUBROUTINE PLAG_RFLU_ComputeVolFracGradL(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_ComputeVolFracGradL
  
  SUBROUTINE PLAG_RFLU_CorrectMixtProperties(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_CorrectMixtProperties
  
  SUBROUTINE PLAG_RFLU_DeallocMemSol(pRegion,pPlag)
    USE ModDataStruct, ONLY: t_region
    USE ModPartLag, ONLY: t_plag
    TYPE(t_region), POINTER :: pRegion 
    TYPE(t_plag), POINTER :: pPlag     
  END SUBROUTINE PLAG_RFLU_DeallocMemSol

  SUBROUTINE PLAG_RFLU_DeallocMemSolTile(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_DeallocMemSolTile

  SUBROUTINE PLAG_RFLU_DeallocMemTStep(pRegion,pPlag)
    USE ModDataStruct, ONLY: t_region
    USE ModPartLag, ONLY: t_plag
    TYPE(t_region), POINTER :: pRegion 
    TYPE(t_plag), POINTER :: pPlag    
  END SUBROUTINE PLAG_RFLU_DeallocMemTStep
  
  SUBROUTINE PLAG_RFLU_DeallocMemTStepTile(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion 
  END SUBROUTINE PLAG_RFLU_DeallocMemTStepTile  
  
  SUBROUTINE PLAG_RFLU_GetMixtPG(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_GetMixtPG
  
  SUBROUTINE PLAG_RFLU_GetMixtSD(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_GetMixtSD

  ! Subbu - Add initialization routine specific to CCMT problem
  SUBROUTINE PLAG_RFLU_InitSolutionCyldet(pRegion,nPclsSumReg)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
    INTEGER :: nPclsSumReg
  END SUBROUTINE PLAG_RFLU_InitSolutionCyldet
  ! Subbu - End add initialization routine specific to CCMT problem

  ! BBR - Add initialization routine specific to Shktb problems
  SUBROUTINE PLAG_RFLU_InitSolutionShktb(pRegion,nPclsSumReg)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
    INTEGER :: nPclsSumReg
  END SUBROUTINE PLAG_RFLU_InitSolutionShktb  
  ! BBR - end

  SUBROUTINE PLAG_RFLU_InitSolutionFile(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion 
  END SUBROUTINE PLAG_RFLU_InitSolutionFile

  SUBROUTINE PLAG_RFLU_InitSolutionHardcode(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion 
  END SUBROUTINE PLAG_RFLU_InitSolutionHardcode
    
  SUBROUTINE PLAG_RFLU_InitSolutionScratch(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion 
  END SUBROUTINE PLAG_RFLU_InitSolutionScratch
  
  SUBROUTINE PLAG_RFLU_InitSolutionRandom(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion 
  END SUBROUTINE PLAG_RFLU_InitSolutionRandom
    
  SUBROUTINE PLAG_RFLU_InitSolFromSerial(pRegion,pRegionSerial)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion,pRegionSerial
  END SUBROUTINE PLAG_RFLU_InitSolFromSerial
  
  SUBROUTINE PLAG_RFLU_InitSolFromSerialCopy(pRegion,pRegionSerial)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion,pRegionSerial
  END SUBROUTINE PLAG_RFLU_InitSolFromSerialCopy
  
  SUBROUTINE PLAG_RFLU_InitSolSerialWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE PLAG_RFLU_InitSolSerialWrapper  
    
  SUBROUTINE PLAG_RFLU_InterParticleForce(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_InterParticleForce
  
  SUBROUTINE PLAG_RFLU_ModifyFlowFieldVolFrac(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_ModifyFlowFieldVolFrac
  
  SUBROUTINE PLAG_RFLU_ReadSolutionASCII(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_ReadSolutionASCII
  
  SUBROUTINE PLAG_RFLU_ReadSolutionBinary(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_ReadSolutionBinary

  SUBROUTINE PLAG_RFLU_ReadUnsteadyDataASCII(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_ReadUnsteadyDataASCII
  
  SUBROUTINE PLAG_RFLU_ReadUnsteadyDataBinary(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_ReadUnsteadyDataBinary

  SUBROUTINE PLAG_RFLU_ShiftUnsteadyData(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_ShiftUnsteadyData
  
  SUBROUTINE PLAG_RFLU_WriteSolutionASCII(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_WriteSolutionASCII

  SUBROUTINE PLAG_RFLU_WriteSolutionBinary(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_WriteSolutionBinary
     
  SUBROUTINE PLAG_RFLU_WriteUnsteadyDataASCII(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_WriteUnsteadyDataASCII

  SUBROUTINE PLAG_RFLU_WriteUnsteadyDataBinary(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_WriteUnsteadyDataBinary
     
  END INTERFACE

END MODULE ModInterfacesLagrangian

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModInterfacesLagrangian.F90,v $
! Revision 1.2  2015/07/27 04:45:42  brollin
! 1) Corrected bug in RFLUCONV where global%gridFormat was used instead of global%gridSrcFormat
! 2) Implemented new subroutine for shock tube problems (Shktb)
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.2  2014/07/21 16:40:21  subbu
! Added PLAG_RFLU_InitSolutionFile subroutine declaration
!
! Revision 1.1.1.1  2014/05/05 21:47:46  tmish
! Initial checkin for rocflu macro.
!
! Revision 1.4  2008/12/06 08:43:37  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:16:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/04/16 23:40:16  fnajjar
! Deleted RFLO-pertinent interface calls
!
! Revision 1.1  2007/04/09 18:49:10  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:17  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.28  2007/03/20 17:32:33  fnajjar
! Deleted obsolete entries with new module PLAG_ModDimensions
!
! Revision 1.27  2007/03/15 21:59:32  haselbac
! Renamed IF for PLAG_RFLU_InitSolSerial
!
! Revision 1.26  2006/05/05 17:19:34  haselbac
! Added if for PLAG_RFLU_InitSolSerial
!
! Revision 1.25  2005/11/30 22:15:03  fnajjar
! Added IF for PLAG_RFLU_CorrectMixtProperties
!
! Revision 1.24  2005/05/18 22:07:16  fnajjar
! Added interface for new init routine
!
! Revision 1.23  2004/12/29 23:26:58  wasistho
! prepared statistics for PLAG and PEUL
!
! Revision 1.22  2004/12/01 00:09:07  wasistho
! added BuildVersionString
!
! Revision 1.21  2004/11/14 19:44:30  haselbac
! Changed interface
!
! Revision 1.20  2004/11/04 16:43:07  fnajjar
! Added Interface call to PLAG_SetDimensions
!
! Revision 1.19  2004/10/10 20:04:35  fnajjar
! Included interface for solution generated by random state
!
! Revision 1.18  2004/08/23 23:08:09  fnajjar
! Added interface calls for binary IO
!
! Revision 1.17  2004/07/28 18:55:32  fnajjar
! Included interface for PLAG_SetMaxDimensions
!
! Revision 1.16  2004/07/26 19:03:07  fnajjar
! Included interface call to PLAG_INRT_DeallocMemTStep routine
!
! Revision 1.15  2004/07/26 17:05:50  fnajjar
! moved allocation of inrtSources into Rocpart
!
! Revision 1.14  2004/03/05 23:21:50  haselbac
! Added interface for PLAG_InitPatchData
!
! Revision 1.13  2004/02/26 21:01:58  haselbac
! Added RFLU routines, modified PLAG_PatchUpdateWrapper entry
!
! Revision 1.12  2004/02/06 21:27:24  fnajjar
! Included proper INTENT to Interfaces
!
! Revision 1.11  2003/11/12 21:18:17  fnajjar
! Added Corner-Edge cells calls
!
! Revision 1.10  2003/04/14 21:12:25  fnajjar
! Bug fix to include POINTER attribute appropriately for PLAG_patchUpdateWrapper
!
! Revision 1.9  2003/04/14 18:56:07  fnajjar
! Included POINTER attribute to regions for PLAG_PatchUpdateWrapper Interface
!
! Revision 1.8  2003/04/14 18:13:28  fnajjar
! Removed iReg from Interface sequence for PLAG_PatchUpdateWrapper
!
! Revision 1.7  2003/03/28 19:39:31  fnajjar
! Aligned with routines pertinent to RocfluidMP
!
! Revision 1.6  2003/03/04 22:12:34  jferry
! Initial import of Rocinteract
!
! Revision 1.5  2003/02/04 19:32:57  f-najjar
! Added Interfaces to PLAG_InjcTileUpdate PLAG_NonCvUpdate
!
! Revision 1.4  2003/01/23 17:06:27  f-najjar
! Included Interface call to PLAG_PatchBufferSendRecv
!
! Revision 1.3  2003/01/13 18:59:19  f-najjar
! Added PLAG_allocateDataBuffers
!
! Revision 1.2  2003/01/10 19:07:55  f-najjar
! Added call to PLAG_PatchUpdate
!
! Revision 1.1  2002/12/27 22:07:14  jblazek
! Splitted up RFLO_ModInterfaces and ModInterfaces.
!
!******************************************************************************

