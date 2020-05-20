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
! Purpose: 
!
! Description: None.
!
! Input: 
!   pRegion		Pointer to region
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_AllocateReadComputeVars.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_AllocateReadComputeVars(pRegion)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModMPI
  USE ModParameters

  USE RFLU_ModAllocateMemory
  USE RFLU_ModBoundLists
  USE RFLU_ModCellMapping
  USE RFLU_ModDimensions
  USE RFLU_ModFaceList
  USE RFLU_ModGeometry
  USE RFLU_ModPlottingVars
  USE RFLU_ModReadBcInputFile
  USE RFLU_ModReadWriteFlow, ONLY: RFLU_ReadFlowWrapper
  USE RFLU_ModReadWriteGrid, ONLY: RFLU_ReadGridWrapper  
  USE RFLU_ModStencilsCells
  USE RFLU_ModVertexLists
  USE RFLU_ModWeights
    
  USE ModInterfaces, ONLY: RFLU_AllocMemSolWrapper, & 
                           RFLU_CreateGrid, &
                           RFLU_DecideBuildGeometry, &
                           RFLU_DecideBuildStencilsWeights, &    
                           RFLU_PrintGridInfo, & 
                           RFLU_SetVarsWrapper
  
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

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_AllocateReadComputeVars.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AllocReadComputeVars',__FILE__)

! ******************************************************************************
! Solution variables. NOTE create grid and read bc file here because allocation
! of particles includes allocation of tiles which requires patch dimensions and
! bc type. NOTE in principle would not need to allocate particle solution here
! but easier to do rather than having an IF and only doing so if plotting vars
! and conversion from PLAG to PEUL desired by user.
! ******************************************************************************

  CALL RFLU_ReadDimensionsWrapper(pRegion)

  CALL RFLU_CreateGrid(pRegion) 

  IF ( pRegion%grid%nPatches > 0 ) THEN      
    CALL RFLU_ReadBCInputFileWrapper(pRegion)    
  END IF ! pRegion%grid%nPatches

  CALL RFLU_AllocMemSolWrapper(pRegion)  
  CALL RFLU_ReadFlowWrapper(pRegion)  
  CALL RFLU_SetVarsWrapper(pRegion,1,pRegion%grid%nCellsTot)

! ******************************************************************************
! Plotting variables 
! ******************************************************************************

  IF ( RFLU_DecideComputePlottingVars(pRegion) .EQV. .TRUE. ) THEN

! ==============================================================================
!   Read grid
! ==============================================================================

    CALL RFLU_ReadGridWrapper(pRegion)  

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      CALL RFLU_PrintGridInfo(pRegion)
    END IF ! global%verbLevel

! ==============================================================================
!   Allocate memory for data structure and grid
! ==============================================================================

    CALL RFLU_CreateCellMapping(pRegion)
    CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
    CALL RFLU_BuildGlob2LocCellMapping(pRegion)

    CALL RFLU_CreateBVertexLists(pRegion)
    CALL RFLU_BuildBVertexLists(pRegion)              

    CALL RFLU_CreateFaceList(pRegion)
    CALL RFLU_BuildFaceList(pRegion) 
    CALL RFLU_RenumberBFaceLists(pRegion)      

    IF ( RFLU_DecideBuildGeometry(global) .EQV. .TRUE. ) THEN
      CALL RFLU_CreateGeometry(pRegion)       
      CALL RFLU_BuildGeometry(pRegion)      
    END IF ! RFLU_DecideBuildGeometry

! ==============================================================================
!   Build stencils and weights
! ==============================================================================

    IF ( RFLU_DecideBuildStencilsWeights(global) .EQV. .TRUE. ) THEN
      CALL RFLU_CreateVert2CellList(pRegion)
      CALL RFLU_BuildVert2CellList(pRegion)

      IF ( pRegion%mixtInput%stencilDimensCells == 1 ) THEN 
        CALL RFLU_CreateCell2FaceList(pRegion)
        CALL RFLU_BuildCell2FaceList(pRegion)
      END IF ! pRegion%mixtInput%stencilDimensCells

      CALL RFLU_SetInfoC2CStencilWrapper(pRegion,1)
      CALL RFLU_CreateC2CStencilWrapper(pRegion)
      CALL RFLU_BuildC2CStencilWrapper(pRegion,constrInput=CONSTR_NONE)

      IF ( pRegion%mixtInput%stencilDimensCells == 1 ) THEN 
        CALL RFLU_DestroyCell2FaceList(pRegion)
      END IF ! pRegion%mixtInput%stencilDimensCells

      CALL RFLU_DestroyVert2CellList(pRegion) 

      CALL RFLU_CreateWtsC2CWrapper(pRegion,1)
      CALL RFLU_ComputeWtsC2CWrapper(pRegion,1)
    END IF ! RFLU_DecideBuildStencilsWeights    
    
! ==============================================================================
!   Compute plotting variables
! ==============================================================================

    CALL RFLU_CountPlottingVars(pRegion)
    CALL RFLU_CreatePlottingVarMaps(pRegion)
    CALL RFLU_BuildPlottingVarMaps(pRegion)
    CALL RFLU_PrintPlottingVarsInfo(pRegion)
    CALL RFLU_CreatePlottingVars(pRegion)         
    CALL RFLU_ComputePlottingVarsWrapper(pRegion)  

! ==============================================================================
!   Destroy stencils and weights  
! ==============================================================================
  
    IF ( RFLU_DecideBuildStencilsWeights(global) .EQV. .TRUE. ) THEN
      CALL RFLU_DestroyWtsC2CWrapper(pRegion)   
      CALL RFLU_DestroyC2CStencilWrapper(pRegion)         
    END IF ! RFLU_DecideBuildStencilsWeights   

! ==============================================================================
!   Destroy data structures
! ==============================================================================
    
    IF ( RFLU_DecideBuildGeometry(pRegion%global) .EQV. .TRUE. ) THEN
      CALL RFLU_DestroyGeometry(pRegion)                            
    END IF ! RFLU_DecideBuildGeometry 

    CALL RFLU_DestroyFaceList(pRegion) 
    CALL RFLU_DestroyBVertexLists(pRegion)                   
    CALL RFLU_DestroyCellMapping(pRegion)       
  END IF ! RFLU_DecideComputePlottingVars

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AllocateReadComputeVars

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_AllocateReadComputeVars.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:55  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:06  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/12/05 13:18:30  haselbac
! Initial revision
!
!******************************************************************************

