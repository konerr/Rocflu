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
!   utilities.
!
! Description: None
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModInterfacesUtilities.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2007 by the University of Illinois
!
! ******************************************************************************
  
MODULE RFLU_ModInterfacesUtilities

  IMPLICIT NONE

  INTERFACE

! ==============================================================================
! Common (same procedure names in different directories)
! ==============================================================================

  SUBROUTINE RFLU_AllocMemSolWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_AllocMemSolWrapper
  
  SUBROUTINE RFLU_DeallocMemSolWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_DeallocMemSolWrapper  

  LOGICAL FUNCTION RFLU_DecideBuildGeometry(global)
    USE ModGlobal, ONLY: t_global 
    TYPE(t_global), POINTER :: global  
  END FUNCTION RFLU_DecideBuildGeometry  

  LOGICAL FUNCTION RFLU_DecideBuildStencilsWeights(global)
    USE ModGlobal, ONLY: t_global 
    TYPE(t_global), POINTER :: global  
  END FUNCTION RFLU_DecideBuildStencilsWeights

! ==============================================================================
! Rocflu cloner
! ==============================================================================

  SUBROUTINE RFLU_CheckClonability(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE RFLU_CheckClonability

  SUBROUTINE RFLU_CloneCommLists(pRegionSerial,pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion,pRegionSerial    
  END SUBROUTINE RFLU_CloneCommLists
  
  SUBROUTINE RFLU_CloneGrid(pRegionSerial,pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion,pRegionSerial    
  END SUBROUTINE RFLU_CloneGrid

! ==============================================================================
! Rocflu partitioner
! ==============================================================================

  SUBROUTINE RFLU_ReadConvGridWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_ReadConvGridWrapper

  SUBROUTINE RFLU_ReadFormatsSection
  END SUBROUTINE RFLU_ReadFormatsSection

  SUBROUTINE RFLU_ReadInputFile(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), DIMENSION(:), POINTER :: regions       
  END SUBROUTINE RFLU_ReadInputFile

  SUBROUTINE RFLU_USER_EnforcePatchCoords(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_USER_EnforcePatchCoords

! ==============================================================================
! Rocflu post-processor
! ==============================================================================

  SUBROUTINE RFLU_AllocateMemoryVert(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_AllocateMemoryVert

  SUBROUTINE RFLU_AllocMemVertWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_AllocMemVertWrapper

  SUBROUTINE RFLU_ComputeExactFlowError(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_ComputeExactFlowError
  
  SUBROUTINE RFLU_ComputeExactFlowProbeError(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_ComputeExactFlowProbeError  
  
  SUBROUTINE RFLU_ComputeVertexVariables(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_ComputeVertexVariables
  
  SUBROUTINE RFLU_DeallocateMemoryVert(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_DeallocateMemoryVert
  
  SUBROUTINE RFLU_DeallocMemVertWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_DeallocMemVertWrapper

  SUBROUTINE RFLU_InterpolateWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_InterpolateWrapper
  
  SUBROUTINE RFLU_MergePostProcessRegions(levels)
    USE ModDataStruct, ONLY: t_level
    TYPE(t_level), POINTER :: levels(:)
  END SUBROUTINE RFLU_MergePostProcessRegions    
  
  SUBROUTINE RFLU_PostProcessRegions(levels)
    USE ModDataStruct, ONLY: t_level
    TYPE(t_level), POINTER :: levels(:)
  END SUBROUTINE RFLU_PostProcessRegions 
  
  SUBROUTINE RFLU_PostProcessRegionsCommon1(pRegion,postInfoFileExists)
    USE ModDataStruct, ONLY: t_region 
    LOGICAL, INTENT(IN) :: postInfoFileExists 
    TYPE(t_region), POINTER :: pRegion   
  END SUBROUTINE RFLU_PostProcessRegionsCommon1 
  
  SUBROUTINE RFLU_PostProcessRegionsCommon2(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion   
  END SUBROUTINE RFLU_PostProcessRegionsCommon2
     
  SUBROUTINE RFLU_PostProcessRegions_ENS(levels)
    USE ModDataStruct, ONLY: t_level
    TYPE(t_level), POINTER :: levels(:)
  END SUBROUTINE RFLU_PostProcessRegions_ENS       
  
  SUBROUTINE RFLU_ReadPostInfo(pRegion,readMode)
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: readMode
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_ReadPostInfo

  SUBROUTINE RFLU_SetPatchPlotFlags(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_SetPatchPlotFlags

! ==============================================================================
! Rocflu initializor
! ==============================================================================

  SUBROUTINE RFLU_InitAuxVars(pRegion)
   USE ModDataStruct, ONLY: t_region 
   TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_InitAuxVars

  SUBROUTINE RFLU_InitBcDataHardCode(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_InitBcDataHardCode
  
  SUBROUTINE RFLU_InitFlowHardCode(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_InitFlowHardCode

  SUBROUTINE RFLU_InitFlowHardCodeLim(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_InitFlowHardCodeLim
  
  SUBROUTINE RFLU_InitFlowHardCodeLimWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_InitFlowHardCodeLimWrapper 

  SUBROUTINE RFLU_InitFlowHardCodeWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_InitFlowHardCodeWrapper 

  SUBROUTINE RFLU_InitFlowScratch(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_InitFlowScratch

  SUBROUTINE RFLU_InitFlowScratchWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_InitFlowScratchWrapper 

  SUBROUTINE RFLU_InitFlowSerialWrapper(pRegion,pRegionSerial)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion,pRegionSerial  
  END SUBROUTINE RFLU_InitFlowSerialWrapper 

! ==============================================================================
! Rocflu picker
! ==============================================================================

  SUBROUTINE RFLU_PickRegionsCoord(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE RFLU_PickRegionsCoord

  SUBROUTINE RFLU_PickRegionsManual(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), DIMENSION(:), POINTER :: regions   
  END SUBROUTINE RFLU_PickRegionsManual

  SUBROUTINE RFLU_PickSpecialCells(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_PickSpecialCells
  
  SUBROUTINE RFLU_PickSpecialFaces(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_PickSpecialFaces  

  SUBROUTINE RFLU_WritePostInfo(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_WritePostInfo

! ==============================================================================
! Rocflu extractor
! ==============================================================================

  SUBROUTINE RFLU_AllocateMemoryXSect(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_AllocateMemoryXSect
  
  SUBROUTINE RFLU_AllocateReadComputeVars(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_AllocateReadComputeVars  

  SUBROUTINE RFLU_BuildStencilInterp(pRegion,ifg,xLoc,yLoc,zLoc,cs)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region  
    INTEGER, INTENT(IN) :: ifg
    INTEGER, DIMENSION(4), INTENT(OUT) :: cs
    REAL(RFREAL), INTENT(IN) :: xLoc,yLoc,zLoc
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_BuildStencilInterp

  SUBROUTINE RFLU_ComputeBilinearInterpWghts(pRegion,xLoc,yLoc,zLoc,cs,wghts)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region  
    INTEGER, DIMENSION(4), INTENT(OUT) :: cs
    REAL(RFREAL), INTENT(IN) :: xLoc,yLoc,zLoc
    REAL(RFREAL), DIMENSION(4), INTENT(OUT) :: wghts
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_ComputeBilinearInterpWghts

  SUBROUTINE RFLU_DeallocateMemoryXSect(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_DeallocateMemoryXSect

  SUBROUTINE RFLU_DeallocateReadComputeVars(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_DeallocateReadComputeVars  

  SUBROUTINE RFLU_ExtractLineDataQuad2D(regions,iRegStart,icgStart)
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: icgStart,iRegStart
    TYPE(t_region), DIMENSION(:), POINTER :: regions 
  END SUBROUTINE RFLU_ExtractLineDataQuad2D

  SUBROUTINE RFLU_WriteLineData(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), DIMENSION(:), POINTER :: regions 
  END SUBROUTINE RFLU_WriteLineData

  END INTERFACE

END MODULE RFLU_ModInterfacesUtilities

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModInterfacesUtilities.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.7  2008/12/06 08:43:41  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:16:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2007/12/05 13:21:14  haselbac
! Added common IF section, added to extr section
!
! Revision 1.4  2007/11/28 23:05:23  mparmar
! Added RFLU_InitAuxVars
!
! Revision 1.3  2007/11/27 13:15:10  haselbac
! Added IFs for rfluextr
!
! Revision 1.2  2007/08/07 17:19:23  haselbac
! Added interfaces for rflucone routines
!
! Revision 1.1  2007/04/09 18:49:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:41  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.23  2006/01/06 22:22:11  haselbac
! Added if for RFLU_DecideBuildStencilsWeights
!
! Revision 1.22  2005/11/10 02:27:36  haselbac
! Modified interfaces for limited hard-coded init following name change
!
! Revision 1.21  2005/10/05 20:07:16  haselbac
! Added interfaces, removed some old unnecessary ones
!
! Revision 1.20  2005/08/09 00:58:43  haselbac
! Added entry for RFLU_SetPatchPlotFlags
!
! Revision 1.19  2005/04/29 12:56:19  haselbac
! Added interface for RFLU_ComputeExactFlowProbeError
!
! Revision 1.18  2005/04/15 15:07:00  haselbac
! Added section for rfluinit interfaces, changed section from rfluprep to rflupart
!
! Revision 1.17  2005/03/29 22:30:12  haselbac
! Added interface for RFLU_InitFlowHardCodeLimited
!
! Revision 1.16  2004/12/29 21:07:53  haselbac
! Added entries for new procedures
!
! Revision 1.15  2004/10/19 19:28:13  haselbac
! Added interface for RFLU_ConvGridWrapper
!
! Revision 1.14  2004/09/27 01:44:42  haselbac
! Clean-up and added interf for RFLU_PickSpecialFaces
!
! Revision 1.13  2004/07/21 14:59:28  haselbac
! Added interface for RFLU_DecideBuildGeometry
!
! Revision 1.12  2004/07/06 15:14:45  haselbac
! Removed many subroutines and functions bcos moved into modules
!
! Revision 1.11  2004/03/19 21:20:26  haselbac
! Removed RFLU_InitFlowScratch from prep, added RFLU_ReInitFlowWrapper for 
! rinit
!
! Revision 1.10  2004/02/26 21:02:05  haselbac
! Added interfaces for alloc/dealloc routines
!
! Revision 1.9  2004/01/29 22:57:39  haselbac
! Added ifs for RFLU_InitBcDataHardCode and RFLU_ComputeExactFlowError
!
! Revision 1.8  2003/11/25 21:03:27  haselbac
! Added interfaces for new routines
!
! Revision 1.7  2003/09/15 00:37:43  haselbac
! Added RFLU_InitFlowHardCode/Scratch
!
! Revision 1.6  2003/08/19 22:49:32  haselbac
! Added interfaces for COBALT conversion routines
!
! Revision 1.5  2003/08/07 15:33:03  haselbac
! Added RFLU_PickRegionsCoord and RFLU_PickRegionsManual
!
! Revision 1.4  2003/06/04 22:11:41  haselbac
! Added, deleted, and modified some interfaces
!
! Revision 1.3  2003/05/05 18:40:59  haselbac
! Added new merge routines
!
! Revision 1.2  2003/04/10 18:48:06  haselbac
! Added and deleted interf for reading CENTAUR grids
!
! Revision 1.1  2003/04/10 14:37:10  haselbac
! Initial revision
!
! ******************************************************************************

