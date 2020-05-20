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
! Purpose: Driver routine for rfluconv.
!
! Description: None.
!
! Input: 
!   caseString  String with casename
!   stampString String with iteration or time stamp
!   verbLevel   Verbosity level
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: rfluconv.F90,v 1.2 2015/07/23 23:11:19 brollin Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE rfluconv(caseString,stampString,verbLevel)

  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModDataTypes
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModParameters

  USE RFLU_ModAllocateMemory, ONLY: RFLU_AllocateMemoryAuxVars, & 
                                    RFLU_AllocateMemorySolDv, &
                                    RFLU_AllocateMemorySolGv
  USE RFLU_ModBoundLists             
  USE RFLU_ModCellMapping
  USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, &
                               RFLU_ConvertCvPrim2Cons
  USE RFLU_ModDeallocateMemory, ONLY: RFLU_DeallocateMemoryAuxVars, &
                                      RFLU_DeallocateMemorySolDv, &
                                      RFLU_DeallocateMemorySolGv
  USE RFLU_ModDimensions
  USE RFLU_ModFaceList 
  USE RFLU_ModGeometry 
  USE RFLU_ModHouMahesh, ONLY: RFLU_HM_ConvCvD2ND,RFLU_HM_ConvCvND2D
  USE RFLU_ModPatchCoeffs
  USE RFLU_ModReadBcInputFile 
  USE RFLU_ModReadWriteAuxVars
  USE RFLU_ModReadWriteFlow 
  USE RFLU_ModReadWriteGrid 
  USE RFLU_ModRegionMapping        
  USE RFLU_ModSTL
  USE RFLU_ModTETMESH
               
  USE ModInterfaces, ONLY: RFLU_AllocMemSolWrapper, & 
                           RFLU_BuildDataStruct, &
                           RFLU_CreateGrid, &
                           RFLU_DeallocMemSolWrapper, &
                           RFLU_DestroyGrid, &
                           RFLU_GetUserInput, &
                           RFLU_InitAuxVars, &
                           RFLU_InitGlobal, &
                           RFLU_PrintFlowInfo, &
                           RFLU_PrintGridInfo, &
                           RFLU_PrintHeader, &
                           RFLU_PrintWarnInfo, &
                           RFLU_SetDependentVars, &
                           RFLU_SetGasVars, &
                           RFLU_SetVarInfoWrapper, &   
                           RFLU_WriteVersionString
                              
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(*) :: caseString,stampString
  INTEGER, INTENT(IN) :: verbLevel

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: casename,nRegions,RCSIdentString,stamp
!  INTEGER, PARAMETER :: CONV_RFLU2RFLU_ASC2BIN_GRID          = 10, & 
!                        CONV_RFLU2RFLU_ASC2BIN_GRIDFLOW      = 11, & 
!                        CONV_RFLU2RFLU_ASC2BIN_GRIDFLOWPATCH = 12, & 
!                        CONV_RFLU2RFLU_BIN2ASC_GRID          = 20, & 
!                        CONV_RFLU2RFLU_BIN2ASC_GRIDFLOW      = 21, & 
!                        CONV_RFLU2RFLU_BIN2ASC_GRIDFLOWPATCH = 22, & 
!                        CONV_RFLU2RFLU_DISS2NONDISS          = 30, & 
!                        CONV_RFLU2RFLU_GRID2GRID             = 70, & 
!                        CONV_RFLU2RFLU_STEADY2UNSTEADY       = 31, & 
!                        CONV_RFLU2RFLU_STEADY2STEADY         = 32, & 
!                        CONV_RFLU2RFLU_NONDISS2DISS          = 40, & 
!                        CONV_RFLU2RFLU_UNSTEADY2STEADY       = 41, & 
!                        CONV_RFLU2RFLU_UNSTEADY2UNSTEADY     = 42, & 
!                        CONV_RFLU2TETM_GRID                  = 50, & 
!                        CONV_RFLU2STL_GRID                   = 60 

! BBR - begin
  INTEGER, PARAMETER :: CONV_RFLU2RFLU_ASC2BINL_GRID          = 10, &
                        CONV_RFLU2RFLU_ASC2BINL_GRIDFLOW      = 11, &
                        CONV_RFLU2RFLU_ASC2BINL_GRIDFLOWPATCH = 12, &
                        CONV_RFLU2RFLU_ASC2BINB_GRID          = 13, &
                        CONV_RFLU2RFLU_ASC2BINB_GRIDFLOW      = 14, &
                        CONV_RFLU2RFLU_ASC2BINB_GRIDFLOWPATCH = 15, &
                        CONV_RFLU2RFLU_BINL2ASC_GRID          = 20, &
                        CONV_RFLU2RFLU_BINL2ASC_GRIDFLOW      = 21, &
                        CONV_RFLU2RFLU_BINL2ASC_GRIDFLOWPATCH = 22, &
                        CONV_RFLU2RFLU_BINB2ASC_GRID          = 23, &
                        CONV_RFLU2RFLU_BINB2ASC_GRIDFLOW      = 24, &
                        CONV_RFLU2RFLU_BINB2ASC_GRIDFLOWPATCH = 25, &
                        CONV_RFLU2RFLU_BINL2BINB_GRID          = 80, &
                        CONV_RFLU2RFLU_BINL2BINB_GRIDFLOW      = 81, &
                        CONV_RFLU2RFLU_BINL2BINB_GRIDFLOWPATCH = 82, &
                        CONV_RFLU2RFLU_BINB2BINL_GRID          = 83, &
                        CONV_RFLU2RFLU_BINB2BINL_GRIDFLOW      = 84, &
                        CONV_RFLU2RFLU_BINB2BINL_GRIDFLOWPATCH = 85, &
                        CONV_RFLU2RFLU_DISS2NONDISS          = 30, &
                        CONV_RFLU2RFLU_GRID2GRID             = 70, &
                        CONV_RFLU2RFLU_STEADY2UNSTEADY       = 31, &
                        CONV_RFLU2RFLU_STEADY2STEADY         = 32, &
                        CONV_RFLU2RFLU_NONDISS2DISS          = 40, &
                        CONV_RFLU2RFLU_UNSTEADY2STEADY       = 41, &
                        CONV_RFLU2RFLU_UNSTEADY2UNSTEADY     = 42, &
                        CONV_RFLU2TETM_GRID                  = 50, &
                        CONV_RFLU2STL_GRID                   = 60
! BBR - end                       
  INTEGER :: convOption,nCells2Copy,readIter,writeIter,errorFlag,iReg
  REAL(RFREAL) :: readTime,writeTime
  TYPE(t_region), POINTER :: pRegion
  TYPE(t_global), POINTER :: global
  TYPE(t_level), POINTER :: levels(:)

! ******************************************************************************

  RCSIdentString = '$RCSfile: rfluconv.F90,v $ $Revision: 1.2 $'

! ******************************************************************************
! Start, initialize global data
! ******************************************************************************

  ALLOCATE(global,STAT=errorFlag)
  IF ( errorFlag /= ERR_NONE ) THEN 
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR - Pointer allocation failed.'
    STOP
  END IF ! errorFlag 

  casename = caseString(1:LEN(caseString))
  stamp    = stampString(1:LEN(stampString))
  
  CALL RFLU_InitGlobal(casename,verbLevel,CRAZY_VALUE_INT,global)

  CALL RegisterFunction(global,'main',__FILE__)

! ******************************************************************************
! Print header
! ******************************************************************************

  CALL RFLU_WriteVersionString(global)
  CALL RFLU_PrintHeader(global)

! ******************************************************************************
! Get conversion option
! ******************************************************************************

  WRITE(STDOUT,'(A,1X,A)')    SOLVER_NAME,'Options:'

  WRITE(STDOUT,'(A,3X,A)')    SOLVER_NAME, &
             'Convert Rocflu files from ASCII to binary little endian format:'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_ASC2BINL_GRID, &
                              ': Grid file only'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_ASC2BINL_GRIDFLOW, &
                              ': Grid and flow files'   
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_ASC2BINL_GRIDFLOWPATCH, &
                              ': Grid, flow, and patch coefficient files'
  
  WRITE(STDOUT,'(A,3X,A)')    SOLVER_NAME, &
             'Convert Rocflu files from ASCII to binary big endian format:'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_ASC2BINB_GRID, &
                              ': Grid file only'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_ASC2BINB_GRIDFLOW, &
                              ': Grid and flow files'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_ASC2BINB_GRIDFLOWPATCH, &
                              ': Grid, flow, and patch coefficient files'
 
  WRITE(STDOUT,'(A,3X,A)')    SOLVER_NAME, &
              'Convert Rocflu files from binary little endian to ASCII format:'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_BINL2ASC_GRID, &
                              ': Grid file only'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_BINL2ASC_GRIDFLOW, &
                              ': Grid and flow files' 
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_BINL2ASC_GRIDFLOWPATCH, &
                              ': Grid, flow, and patch coefficient files'

  WRITE(STDOUT,'(A,3X,A)')    SOLVER_NAME, &
              'Convert Rocflu files from binary big endian to ASCII format:'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_BINB2ASC_GRID, &
                              ': Grid file only'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_BINB2ASC_GRIDFLOW, &
                              ': Grid and flow files'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_BINB2ASC_GRIDFLOWPATCH, &
                              ': Grid, flow, and patch coefficient files'

  WRITE(STDOUT,'(A,3X,A)')    SOLVER_NAME, &
  'Convert Rocflu files from binary little endian to binary big endian format:'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_BINL2BINB_GRID, &
                              ': Grid file only'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_BINL2BINB_GRIDFLOW, &
                              ': Grid and flow files'
 WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_BINL2BINB_GRIDFLOWPATCH, &
                              ': Grid, flow, and patch coefficient files'

  WRITE(STDOUT,'(A,3X,A)')    SOLVER_NAME, &
  'Convert Rocflu files from binary big endian to binary little endian format:'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_BINB2BINL_GRID, &
                              ': Grid file only'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_BINB2BINL_GRIDFLOW, &
                              ': Grid and flow files'
 WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_BINB2BINL_GRIDFLOWPATCH, &
                              ': Grid, flow, and patch coefficient files'

  WRITE(STDOUT,'(A,3X,A)')    SOLVER_NAME, &
                              'Convert Rocflu files for different solvers:'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_DISS2NONDISS, &
                              ': Dissipative to non-dissipative' 
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_NONDISS2DISS, &
                              ': Non-dissipative to dissipative' 
  WRITE(STDOUT,'(A,3X,A)')    SOLVER_NAME, &                                  
                              'Convert Rocflu grid file to surface grid:'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2TETM_GRID, &
                              ': Tetmesh/YAMS surface grid (msh2 format)' 
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2STL_GRID, &
                              ': STL surface grid'          
  WRITE(STDOUT,'(A,3X,A)')    SOLVER_NAME, &
                              'Convert Rocflu between steady/unsteady solvers:'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_STEADY2UNSTEADY, &
                              ': Steady rocflu data to unsteady rocflu data' 
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_STEADY2STEADY, &
                              ': Steady rocflu data to steady rocflu data' 
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_UNSTEADY2STEADY, &
                              ': Unsteady rocflu data to steady rocflu data' 
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_UNSTEADY2UNSTEADY, &
                              ': Unsteady rocflu data to unsteady rocflu data' 
  WRITE(STDOUT,'(A,3X,A)')    SOLVER_NAME, &
                              'Convert Rocflu between different meshes:'
  WRITE(STDOUT,'(A,5X,I2,A)') SOLVER_NAME,CONV_RFLU2RFLU_GRID2GRID, & 
                              ': Rocflu data on larger mesh to smaller mesh' 

  WRITE(STDOUT,'(A,1X,A)')    SOLVER_NAME,'Enter selection:'

  READ(STDIN,*) convOption

! ******************************************************************************
! Read mapping file and impose serial mapping
! ******************************************************************************

  CALL RFLU_ReadRegionMappingFile(global,MAPFILE_READMODE_PEEK,global%myProcId)
  CALL RFLU_SetRegionMappingSerial(global)
  CALL RFLU_CreateRegionMapping(global,MAPTYPE_REG)
  CALL RFLU_ImposeRegionMappingSerial(global)
  
! ******************************************************************************
! Prepare data structure
! ******************************************************************************

  CALL RFLU_BuildDataStruct(global,levels) 
  CALL RFLU_ApplyRegionMapping(global,levels)
  CALL RFLU_DestroyRegionMapping(global,MAPTYPE_REG)

! ******************************************************************************
! Read input file
! ******************************************************************************

  CALL RFLU_GetUserInput(levels(1)%regions,.TRUE.) 

  IF ( global%flowType == FLOW_STEADY ) THEN
    READ(stamp,*) global%currentIter
    readIter = global%currentIter
  ELSE
    READ(stamp,*) global%currentTime
    readTime = global%currentTime
  END IF ! global%flowType  

! ******************************************************************************
! Read dimensions and allocate memory
! ****************************************************************************** 

  IF ( global%nRegions == 1 ) THEN ! single region
    pRegion => levels(1)%regions(0)
    
    CALL RFLU_ReadDimensions(pRegion)      
    CALL RFLU_CreateGrid(pRegion)       
 
    IF ( pRegion%grid%nPatches > 0 ) THEN      
      CALL RFLU_ReadBcInputFileWrapper(pRegion)    
    END IF ! pRegion%grid%nPatches     

    CALL RFLU_AllocMemSolWrapper(pRegion)
    CALL RFLU_SetVarInfoWrapper(pRegion)    
    CALL RFLU_CreatePatchCoeffs(pRegion)
  ELSE ! multiple regions
    DO iReg = 1,global%nRegionsLocal
      pRegion => levels(1)%regions(iReg)

      CALL RFLU_ReadDimensions(pRegion)         
      CALL RFLU_CreateGrid(pRegion)

      IF ( pRegion%grid%nPatches > 0 ) THEN      
        CALL RFLU_ReadBCInputFileWrapper(pRegion)    
      END IF ! pRegion%grid%nPatches
      
      CALL RFLU_AllocMemSolWrapper(pRegion)
      CALL RFLU_SetVarInfoWrapper(pRegion)             
      CALL RFLU_CreatePatchCoeffs(pRegion)
    END DO ! iReg  
  END IF ! global%nRegions

! ******************************************************************************
! Convert
! ******************************************************************************

  SELECT CASE ( convOption )

! ==============================================================================
!   Rocflu files from ASCII to binary little endian - NOTE that the gridFormat 
!   and solutFormat variables are overridden.
! ==============================================================================  
  
! ------------------------------------------------------------------------------
!   Grid file only
! ------------------------------------------------------------------------------
  
    CASE (CONV_RFLU2RFLU_ASC2BINL_GRID)        
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)          
      
        pRegion%global%gridFormat = FORMAT_ASCII           
        CALL RFLU_ReadGridWrapper(pRegion)  
      
        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
        END IF ! global%verbLevel 

        pRegion%global%gridFormat = FORMAT_BINARY_L
        CALL RFLU_WriteGridWrapper(pRegion) 
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg) 
         
          pRegion%global%gridFormat = FORMAT_ASCII              
          CALL RFLU_ReadGridWrapper(pRegion)  
      
          IF ( global%verbLevel > VERBOSE_NONE ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
          END IF ! global%verbLevel 

          pRegion%global%gridFormat = FORMAT_BINARY_L
          CALL RFLU_WriteGridWrapper(pRegion)                
        END DO ! iReg       
      END IF ! global%nRegions 

! ------------------------------------------------------------------------------
!   Grid and solution files
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_ASC2BINL_GRIDFLOW)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)          

        pRegion%global%gridFormat  = FORMAT_ASCII
        pRegion%global%solutFormat = FORMAT_ASCII   
        CALL RFLU_ReadGridWrapper(pRegion)  
        CALL RFLU_ReadFlowWrapper(pRegion)      

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_ReadAuxVarsWrapper(pRegion)
        END IF ! solverType

        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)          
        END IF ! global%verbLevel 

        pRegion%global%gridFormat  = FORMAT_BINARY_L      
        pRegion%global%solutFormat = FORMAT_BINARY_L      

        CALL RFLU_WriteGridWrapper(pRegion) 
        CALL RFLU_WriteFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! solverType 
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg) 

          pRegion%global%gridFormat  = FORMAT_ASCII
          pRegion%global%solutFormat = FORMAT_ASCII   
          CALL RFLU_ReadGridWrapper(pRegion)  
          CALL RFLU_ReadFlowWrapper(pRegion)      

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_ReadAuxVarsWrapper(pRegion)
          END IF ! solverType

          IF ( global%verbLevel > VERBOSE_NONE ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)          
          END IF ! global%verbLevel 

          pRegion%global%gridFormat  = FORMAT_BINARY_L      
          pRegion%global%solutFormat = FORMAT_BINARY_L      

          CALL RFLU_WriteGridWrapper(pRegion) 
          CALL RFLU_WriteFlowWrapper(pRegion)                 

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_WriteAuxVarsWrapper(pRegion)
          END IF ! solverType 
        END DO ! iReg       
      END IF ! global%nRegions    
      
! ------------------------------------------------------------------------------
!   Grid, solution and patch coefficient files
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_ASC2BINL_GRIDFLOWPATCH)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)          

        pRegion%global%gridFormat  = FORMAT_ASCII
        pRegion%global%solutFormat = FORMAT_ASCII   
        CALL RFLU_ReadGridWrapper(pRegion)  
        CALL RFLU_ReadFlowWrapper(pRegion)
          
        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_ReadAuxVarsWrapper(pRegion)
        END IF ! solverType

        CALL RFLU_ReadPatchCoeffsWrapper(pRegion)      

        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)          
        END IF ! global%verbLevel 

        pRegion%global%gridFormat  = FORMAT_BINARY_L      
        pRegion%global%solutFormat = FORMAT_BINARY_L      

        CALL RFLU_WriteGridWrapper(pRegion) 
        CALL RFLU_WriteFlowWrapper(pRegion) 

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! solverType 

        CALL RFLU_WritePatchCoeffsWrapper(pRegion) 
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg) 

          pRegion%global%gridFormat  = FORMAT_ASCII
          pRegion%global%solutFormat = FORMAT_ASCII   
          CALL RFLU_ReadGridWrapper(pRegion)  
          CALL RFLU_ReadFlowWrapper(pRegion)      

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_ReadAuxVarsWrapper(pRegion)
          END IF ! solverType

          CALL RFLU_ReadPatchCoeffsWrapper(pRegion)      

          IF ( global%verbLevel > VERBOSE_NONE ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)          
          END IF ! global%verbLevel 

          pRegion%global%gridFormat  = FORMAT_BINARY_L      
          pRegion%global%solutFormat = FORMAT_BINARY_L      

          CALL RFLU_WriteGridWrapper(pRegion) 
          CALL RFLU_WriteFlowWrapper(pRegion)                 

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_WriteAuxVarsWrapper(pRegion)
          END IF ! solverType 

          CALL RFLU_WritePatchCoeffsWrapper(pRegion) 
        END DO ! iReg       
      END IF ! global%nRegions    

! ==============================================================================
!   Rocflu files from ASCII to binary big endian - NOTE that the gridFormat 
!   and solutFormat variables are overridden.
! ==============================================================================  

! ------------------------------------------------------------------------------
!   Grid file only
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_ASC2BINB_GRID)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)

        pRegion%global%gridFormat = FORMAT_ASCII
        CALL RFLU_ReadGridWrapper(pRegion)

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          CALL RFLU_PrintGridInfo(pRegion)
        END IF ! global%verbLevel 

        pRegion%global%gridFormat = FORMAT_BINARY_B
        CALL RFLU_WriteGridWrapper(pRegion)
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg)

          pRegion%global%gridFormat = FORMAT_ASCII
          CALL RFLU_ReadGridWrapper(pRegion)

          IF ( global%verbLevel > VERBOSE_NONE ) THEN
            CALL RFLU_PrintGridInfo(pRegion)
          END IF ! global%verbLevel 

          pRegion%global%gridFormat = FORMAT_BINARY_B
          CALL RFLU_WriteGridWrapper(pRegion)
        END DO ! iReg       
      END IF ! global%nRegions 

! ------------------------------------------------------------------------------
!   Grid and solution files
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_ASC2BINB_GRIDFLOW)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)

        pRegion%global%gridFormat  = FORMAT_ASCII
        pRegion%global%solutFormat = FORMAT_ASCII
        CALL RFLU_ReadGridWrapper(pRegion)
        CALL RFLU_ReadFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_ReadAuxVarsWrapper(pRegion)
        END IF ! solverType

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)
        END IF ! global%verbLevel 

        pRegion%global%gridFormat  = FORMAT_BINARY_B
        pRegion%global%solutFormat = FORMAT_BINARY_B

        CALL RFLU_WriteGridWrapper(pRegion)
        CALL RFLU_WriteFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! solverType 
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg)

          pRegion%global%gridFormat  = FORMAT_ASCII
          pRegion%global%solutFormat = FORMAT_ASCII
          CALL RFLU_ReadGridWrapper(pRegion)
          CALL RFLU_ReadFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_ReadAuxVarsWrapper(pRegion)
          END IF ! solverType

          IF ( global%verbLevel > VERBOSE_NONE ) THEN
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)
          END IF ! global%verbLevel 

          pRegion%global%gridFormat  = FORMAT_BINARY_B
          pRegion%global%solutFormat = FORMAT_BINARY_B

          CALL RFLU_WriteGridWrapper(pRegion)
          CALL RFLU_WriteFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_WriteAuxVarsWrapper(pRegion)
          END IF ! solverType 
        END DO ! iReg       
      END IF ! global%nRegions 

! ------------------------------------------------------------------------------
!   Grid, solution and patch coefficient files
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_ASC2BINB_GRIDFLOWPATCH)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)

        pRegion%global%gridFormat  = FORMAT_ASCII
        pRegion%global%solutFormat = FORMAT_ASCII
        CALL RFLU_ReadGridWrapper(pRegion)
        CALL RFLU_ReadFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_ReadAuxVarsWrapper(pRegion)
        END IF ! solverType

        CALL RFLU_ReadPatchCoeffsWrapper(pRegion)

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)
        END IF ! global%verbLevel 

        pRegion%global%gridFormat  = FORMAT_BINARY_B
        pRegion%global%solutFormat = FORMAT_BINARY_B

        CALL RFLU_WriteGridWrapper(pRegion)
        CALL RFLU_WriteFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! solverType 

        CALL RFLU_WritePatchCoeffsWrapper(pRegion)
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg)

          pRegion%global%gridFormat  = FORMAT_ASCII
          pRegion%global%solutFormat = FORMAT_ASCII
          CALL RFLU_ReadGridWrapper(pRegion)
          CALL RFLU_ReadFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_ReadAuxVarsWrapper(pRegion)
          END IF ! solverType

          CALL RFLU_ReadPatchCoeffsWrapper(pRegion)

          IF ( global%verbLevel > VERBOSE_NONE ) THEN
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)
          END IF ! global%verbLevel 

          pRegion%global%gridFormat  = FORMAT_BINARY_B
          pRegion%global%solutFormat = FORMAT_BINARY_B

          CALL RFLU_WriteGridWrapper(pRegion)
          CALL RFLU_WriteFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_WriteAuxVarsWrapper(pRegion)
          END IF ! solverType 

          CALL RFLU_WritePatchCoeffsWrapper(pRegion)
        END DO ! iReg       
      END IF ! global%nRegions 

! ==============================================================================
!   Rocflu files from binary little endian to ASCII - NOTE that the gridFormat 
!   and solutFormat variables are overridden.
! ==============================================================================        

! ------------------------------------------------------------------------------
!   Grid file only
! ------------------------------------------------------------------------------
      
    CASE (CONV_RFLU2RFLU_BINL2ASC_GRID)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)          
      
        pRegion%global%gridFormat = FORMAT_BINARY_L          
        CALL RFLU_ReadGridWrapper(pRegion)  
      
        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
        END IF ! global%verbLevel 

        pRegion%global%gridFormat = FORMAT_ASCII
        CALL RFLU_WriteGridWrapper(pRegion) 
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg) 
         
          pRegion%global%gridFormat = FORMAT_BINARY_L             
          CALL RFLU_ReadGridWrapper(pRegion)  
      
          IF ( global%verbLevel > VERBOSE_NONE ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
          END IF ! global%verbLevel 

          pRegion%global%gridFormat = FORMAT_ASCII
          CALL RFLU_WriteGridWrapper(pRegion)                
        END DO ! iReg       
      END IF ! global%nRegions 
      
! ------------------------------------------------------------------------------
!   Grid and solution files
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_BINL2ASC_GRIDFLOW) 
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)          

        pRegion%global%gridFormat  = FORMAT_BINARY_L
        pRegion%global%solutFormat = FORMAT_BINARY_L  
        CALL RFLU_ReadGridWrapper(pRegion)  
        CALL RFLU_ReadFlowWrapper(pRegion)      

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_ReadAuxVarsWrapper(pRegion)
        END IF ! solverType

        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)          
        END IF ! global%verbLevel 

        pRegion%global%gridFormat  = FORMAT_ASCII      
        pRegion%global%solutFormat = FORMAT_ASCII      

        CALL RFLU_WriteGridWrapper(pRegion) 
        CALL RFLU_WriteFlowWrapper(pRegion) 

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! solverType 
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg) 

          pRegion%global%gridFormat  = FORMAT_BINARY_L
          pRegion%global%solutFormat = FORMAT_BINARY_L  
          CALL RFLU_ReadGridWrapper(pRegion)  
          CALL RFLU_ReadFlowWrapper(pRegion)      

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_ReadAuxVarsWrapper(pRegion)
          END IF ! solverType

          IF ( global%verbLevel > VERBOSE_NONE ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)          
          END IF ! global%verbLevel 

          pRegion%global%gridFormat  = FORMAT_ASCII     
          pRegion%global%solutFormat = FORMAT_ASCII     

          CALL RFLU_WriteGridWrapper(pRegion) 
          CALL RFLU_WriteFlowWrapper(pRegion)                 

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_WriteAuxVarsWrapper(pRegion)
          END IF ! solverType 
        END DO ! iReg       
      END IF ! global%nRegions    

! ------------------------------------------------------------------------------
!   Grid, solution and patch coefficient files
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_BINL2ASC_GRIDFLOWPATCH) 
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)          

        pRegion%global%gridFormat  = FORMAT_BINARY_L
        pRegion%global%solutFormat = FORMAT_BINARY_L  
        CALL RFLU_ReadGridWrapper(pRegion)  
        CALL RFLU_ReadFlowWrapper(pRegion)      

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_ReadAuxVarsWrapper(pRegion)
        END IF ! solverType

        CALL RFLU_ReadPatchCoeffsWrapper(pRegion)      

        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)          
        END IF ! global%verbLevel 

        pRegion%global%gridFormat  = FORMAT_ASCII      
        pRegion%global%solutFormat = FORMAT_ASCII      

        CALL RFLU_WriteGridWrapper(pRegion) 
        CALL RFLU_WriteFlowWrapper(pRegion) 

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! solverType 

        CALL RFLU_WritePatchCoeffsWrapper(pRegion) 
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg) 

          pRegion%global%gridFormat  = FORMAT_BINARY_L
          pRegion%global%solutFormat = FORMAT_BINARY_L  
          CALL RFLU_ReadGridWrapper(pRegion)  
          CALL RFLU_ReadFlowWrapper(pRegion)      

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_ReadAuxVarsWrapper(pRegion)
          END IF ! solverType

          CALL RFLU_ReadPatchCoeffsWrapper(pRegion)      

          IF ( global%verbLevel > VERBOSE_NONE ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)          
          END IF ! global%verbLevel 

          pRegion%global%gridFormat  = FORMAT_ASCII     
          pRegion%global%solutFormat = FORMAT_ASCII     

          CALL RFLU_WriteGridWrapper(pRegion) 
          CALL RFLU_WriteFlowWrapper(pRegion)                 

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_WriteAuxVarsWrapper(pRegion)
          END IF ! solverType 

          CALL RFLU_WritePatchCoeffsWrapper(pRegion) 
        END DO ! iReg       
      END IF ! global%nRegions    

! ==============================================================================
!   Rocflu files from binary big endian to ASCII - NOTE that the gridFormat 
!   and solutFormat variables are overridden.
! ==============================================================================        

! ------------------------------------------------------------------------------
!   Grid file only
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_BINB2ASC_GRID)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)

        pRegion%global%gridFormat = FORMAT_BINARY_B
        CALL RFLU_ReadGridWrapper(pRegion)

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          CALL RFLU_PrintGridInfo(pRegion)
        END IF ! global%verbLevel 

        pRegion%global%gridFormat = FORMAT_ASCII
        CALL RFLU_WriteGridWrapper(pRegion)
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg)

          pRegion%global%gridFormat = FORMAT_BINARY_B
          CALL RFLU_ReadGridWrapper(pRegion)

          IF ( global%verbLevel > VERBOSE_NONE ) THEN
            CALL RFLU_PrintGridInfo(pRegion)
          END IF ! global%verbLevel 

          pRegion%global%gridFormat = FORMAT_ASCII
          CALL RFLU_WriteGridWrapper(pRegion)
        END DO ! iReg       
      END IF ! global%nRegions 

! ------------------------------------------------------------------------------
!   Grid and solution files
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_BINB2ASC_GRIDFLOW)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)

        pRegion%global%gridFormat  = FORMAT_BINARY_B
        pRegion%global%solutFormat = FORMAT_BINARY_B
        CALL RFLU_ReadGridWrapper(pRegion)
        CALL RFLU_ReadFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_ReadAuxVarsWrapper(pRegion)
        END IF ! solverType

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)
        END IF ! global%verbLevel 

        pRegion%global%gridFormat  = FORMAT_ASCII
        pRegion%global%solutFormat = FORMAT_ASCII

        CALL RFLU_WriteGridWrapper(pRegion)
        CALL RFLU_WriteFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! solverType 
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg)

          pRegion%global%gridFormat  = FORMAT_BINARY_B
          pRegion%global%solutFormat = FORMAT_BINARY_B
          CALL RFLU_ReadGridWrapper(pRegion)
          CALL RFLU_ReadFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_ReadAuxVarsWrapper(pRegion)
          END IF ! solverType

          IF ( global%verbLevel > VERBOSE_NONE ) THEN
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)
          END IF ! global%verbLevel 

          pRegion%global%gridFormat  = FORMAT_ASCII
          pRegion%global%solutFormat = FORMAT_ASCII

          CALL RFLU_WriteGridWrapper(pRegion)
          CALL RFLU_WriteFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_WriteAuxVarsWrapper(pRegion)
          END IF ! solverType 
        END DO ! iReg       
      END IF ! global%nRegions 

! ------------------------------------------------------------------------------
!   Grid, solution and patch coefficient files
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_BINB2ASC_GRIDFLOWPATCH)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)

        pRegion%global%gridFormat  = FORMAT_BINARY_B
        pRegion%global%solutFormat = FORMAT_BINARY_B
        CALL RFLU_ReadGridWrapper(pRegion)
        CALL RFLU_ReadFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_ReadAuxVarsWrapper(pRegion)
        END IF ! solverType

        CALL RFLU_ReadPatchCoeffsWrapper(pRegion)

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)
        END IF ! global%verbLevel 

        pRegion%global%gridFormat  = FORMAT_ASCII
        pRegion%global%solutFormat = FORMAT_ASCII

        CALL RFLU_WriteGridWrapper(pRegion)
        CALL RFLU_WriteFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! solverType 

        CALL RFLU_WritePatchCoeffsWrapper(pRegion)
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg)

          pRegion%global%gridFormat  = FORMAT_BINARY_B
          pRegion%global%solutFormat = FORMAT_BINARY_B
          CALL RFLU_ReadGridWrapper(pRegion)
          CALL RFLU_ReadFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_ReadAuxVarsWrapper(pRegion)
          END IF ! solverType

          CALL RFLU_ReadPatchCoeffsWrapper(pRegion)

          IF ( global%verbLevel > VERBOSE_NONE ) THEN
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)
          END IF ! global%verbLevel 

          pRegion%global%gridFormat  = FORMAT_ASCII
          pRegion%global%solutFormat = FORMAT_ASCII

          CALL RFLU_WriteGridWrapper(pRegion)
          CALL RFLU_WriteFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_WriteAuxVarsWrapper(pRegion)
          END IF ! solverType 

          CALL RFLU_WritePatchCoeffsWrapper(pRegion)
        END DO ! iReg       
      END IF ! global%nRegions

! ==============================================================================
!   Rocflu files from binary little endian to big endian - NOTE that the 
!   gridFormat and solutFormat variables are overridden.
! ==============================================================================        

! ------------------------------------------------------------------------------
!   Grid file only
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_BINL2BINB_GRID)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)

        pRegion%global%gridFormat = FORMAT_BINARY_L
        CALL RFLU_ReadGridWrapper(pRegion)

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          CALL RFLU_PrintGridInfo(pRegion)
        END IF ! global%verbLevel 

        pRegion%global%gridFormat = FORMAT_BINARY_B
        CALL RFLU_WriteGridWrapper(pRegion)
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg)

          pRegion%global%gridFormat = FORMAT_BINARY_L
          CALL RFLU_ReadGridWrapper(pRegion)

          IF ( global%verbLevel > VERBOSE_NONE ) THEN
            CALL RFLU_PrintGridInfo(pRegion)
          END IF ! global%verbLevel 

          pRegion%global%gridFormat = FORMAT_BINARY_B
          CALL RFLU_WriteGridWrapper(pRegion)
        END DO ! iReg       
      END IF ! global%nRegions 

! ------------------------------------------------------------------------------
!   Grid and solution files
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_BINL2BINB_GRIDFLOW)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)

        pRegion%global%gridFormat  = FORMAT_BINARY_L
        pRegion%global%solutFormat = FORMAT_BINARY_L
        CALL RFLU_ReadGridWrapper(pRegion)
        CALL RFLU_ReadFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_ReadAuxVarsWrapper(pRegion)
        END IF ! solverType

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)
        END IF ! global%verbLevel 

        pRegion%global%gridFormat  = FORMAT_BINARY_B
        pRegion%global%solutFormat = FORMAT_BINARY_B

        CALL RFLU_WriteGridWrapper(pRegion)
        CALL RFLU_WriteFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! solverType 
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg)

          pRegion%global%gridFormat  = FORMAT_BINARY_L
          pRegion%global%solutFormat = FORMAT_BINARY_L
          CALL RFLU_ReadGridWrapper(pRegion)
          CALL RFLU_ReadFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_ReadAuxVarsWrapper(pRegion)
          END IF ! solverType

          IF ( global%verbLevel > VERBOSE_NONE ) THEN
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)
          END IF ! global%verbLevel 

          pRegion%global%gridFormat  = FORMAT_BINARY_B
          pRegion%global%solutFormat = FORMAT_BINARY_B

          CALL RFLU_WriteGridWrapper(pRegion)
          CALL RFLU_WriteFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_WriteAuxVarsWrapper(pRegion)
          END IF ! solverType 
        END DO ! iReg       
      END IF ! global%nRegions 

! ------------------------------------------------------------------------------
!   Grid, solution and patch coefficient files
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_BINL2BINB_GRIDFLOWPATCH)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)

        pRegion%global%gridFormat  = FORMAT_BINARY_L
        pRegion%global%solutFormat = FORMAT_BINARY_L
        CALL RFLU_ReadGridWrapper(pRegion)
        CALL RFLU_ReadFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_ReadAuxVarsWrapper(pRegion)
        END IF ! solverType

        CALL RFLU_ReadPatchCoeffsWrapper(pRegion)

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)
        END IF ! global%verbLevel 

        pRegion%global%gridFormat  = FORMAT_BINARY_B
        pRegion%global%solutFormat = FORMAT_BINARY_B

        CALL RFLU_WriteGridWrapper(pRegion)
        CALL RFLU_WriteFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! solverType 

        CALL RFLU_WritePatchCoeffsWrapper(pRegion)
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg)

          pRegion%global%gridFormat  = FORMAT_BINARY_L
          pRegion%global%solutFormat = FORMAT_BINARY_L
          CALL RFLU_ReadGridWrapper(pRegion)
          CALL RFLU_ReadFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_ReadAuxVarsWrapper(pRegion)
          END IF ! solverType

          CALL RFLU_ReadPatchCoeffsWrapper(pRegion)

          IF ( global%verbLevel > VERBOSE_NONE ) THEN
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)
          END IF ! global%verbLevel 

          pRegion%global%gridFormat  = FORMAT_BINARY_B
          pRegion%global%solutFormat = FORMAT_BINARY_B

          CALL RFLU_WriteGridWrapper(pRegion)
          CALL RFLU_WriteFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_WriteAuxVarsWrapper(pRegion)
          END IF ! solverType 

          CALL RFLU_WritePatchCoeffsWrapper(pRegion)
        END DO ! iReg       
      END IF ! global%nRegions

! ==============================================================================
!   Rocflu files from binary big endian to little endian - NOTE that the 
!   gridFormat and solutFormat variables are overridden.
! ==============================================================================        

! ------------------------------------------------------------------------------
!   Grid file only
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_BINB2BINL_GRID)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)

        pRegion%global%gridFormat = FORMAT_BINARY_B
        CALL RFLU_ReadGridWrapper(pRegion)

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          CALL RFLU_PrintGridInfo(pRegion)
        END IF ! global%verbLevel 

        pRegion%global%gridFormat = FORMAT_BINARY_L
        CALL RFLU_WriteGridWrapper(pRegion)
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg)

          pRegion%global%gridFormat = FORMAT_BINARY_B
          CALL RFLU_ReadGridWrapper(pRegion)

          IF ( global%verbLevel > VERBOSE_NONE ) THEN
            CALL RFLU_PrintGridInfo(pRegion)
          END IF ! global%verbLevel 

          pRegion%global%gridFormat = FORMAT_BINARY_L
          CALL RFLU_WriteGridWrapper(pRegion)
        END DO ! iReg       
      END IF ! global%nRegions 

! ------------------------------------------------------------------------------
!   Grid and solution files
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_BINB2BINL_GRIDFLOW)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)

        pRegion%global%gridFormat  = FORMAT_BINARY_B
        pRegion%global%solutFormat = FORMAT_BINARY_B
        CALL RFLU_ReadGridWrapper(pRegion)
        CALL RFLU_ReadFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_ReadAuxVarsWrapper(pRegion)
        END IF ! solverType

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)
        END IF ! global%verbLevel 

        pRegion%global%gridFormat  = FORMAT_BINARY_L
        pRegion%global%solutFormat = FORMAT_BINARY_L

        CALL RFLU_WriteGridWrapper(pRegion)
        CALL RFLU_WriteFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! solverType 
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg)

          pRegion%global%gridFormat  = FORMAT_BINARY_B
          pRegion%global%solutFormat = FORMAT_BINARY_B
          CALL RFLU_ReadGridWrapper(pRegion)
          CALL RFLU_ReadFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_ReadAuxVarsWrapper(pRegion)
          END IF ! solverType

          IF ( global%verbLevel > VERBOSE_NONE ) THEN
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)
          END IF ! global%verbLevel 

          pRegion%global%gridFormat  = FORMAT_BINARY_L
          pRegion%global%solutFormat = FORMAT_BINARY_L

          CALL RFLU_WriteGridWrapper(pRegion)
          CALL RFLU_WriteFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_WriteAuxVarsWrapper(pRegion)
          END IF ! solverType 
        END DO ! iReg       
      END IF ! global%nRegions

! ------------------------------------------------------------------------------
!   Grid, solution and patch coefficient files
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_BINB2BINL_GRIDFLOWPATCH)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)

        pRegion%global%gridFormat  = FORMAT_BINARY_B
        pRegion%global%solutFormat = FORMAT_BINARY_B
        CALL RFLU_ReadGridWrapper(pRegion)
        CALL RFLU_ReadFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_ReadAuxVarsWrapper(pRegion)
        END IF ! solverType

        CALL RFLU_ReadPatchCoeffsWrapper(pRegion)

        IF ( global%verbLevel > VERBOSE_NONE ) THEN
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)
        END IF ! global%verbLevel 

        pRegion%global%gridFormat  = FORMAT_BINARY_L
        pRegion%global%solutFormat = FORMAT_BINARY_L

        CALL RFLU_WriteGridWrapper(pRegion)
        CALL RFLU_WriteFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! solverType 

        CALL RFLU_WritePatchCoeffsWrapper(pRegion)
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg)

          pRegion%global%gridFormat  = FORMAT_BINARY_B
          pRegion%global%solutFormat = FORMAT_BINARY_B
          CALL RFLU_ReadGridWrapper(pRegion)
          CALL RFLU_ReadFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_ReadAuxVarsWrapper(pRegion)
          END IF ! solverType

          CALL RFLU_ReadPatchCoeffsWrapper(pRegion)

          IF ( global%verbLevel > VERBOSE_NONE ) THEN
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)
          END IF ! global%verbLevel 

          pRegion%global%gridFormat  = FORMAT_BINARY_L
          pRegion%global%solutFormat = FORMAT_BINARY_L

          CALL RFLU_WriteGridWrapper(pRegion)
          CALL RFLU_WriteFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_WriteAuxVarsWrapper(pRegion)
          END IF ! solverType 

          CALL RFLU_WritePatchCoeffsWrapper(pRegion)
        END DO ! iReg       
      END IF ! global%nRegions

! ==============================================================================
!   Convert Rocflu data from one solver type to another 
! ==============================================================================

! ------------------------------------------------------------------------------
!   Dissipative solver to non-dissipative solver
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_DISS2NONDISS)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)          

        CALL RFLU_AllocateMemorySolDv(pRegion)
        CALL RFLU_AllocateMemorySolGv(pRegion)

        pRegion%mixtInput%nCvOld2 = 2
        CALL RFLU_AllocateMemoryAuxVars(pRegion)

        CALL RFLU_ReadGridWrapper(pRegion)  
        CALL RFLU_ReadFlowWrapper(pRegion)      

        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)          
        END IF ! global%verbLevel 

        CALL RFLU_SetGasVars(pRegion,1,pRegion%grid%nCellsTot)
        CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCellsTot)

        CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)

        CALL RFLU_InitAuxVars(pRegion)

        CALL RFLU_HM_ConvCvD2ND(pRegion)

        CALL RFLU_WriteGridWrapper(pRegion) 
        CALL RFLU_WriteFlowWrapper(pRegion)
        CALL RFLU_WriteAuxVarsWrapper(pRegion)

        CALL RFLU_DeallocateMemorySolDv(pRegion)
        CALL RFLU_DeallocateMemorySolGv(pRegion)
        CALL RFLU_DeallocateMemoryAuxVars(pRegion)
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg) 

          CALL RFLU_AllocateMemorySolDv(pRegion)
          CALL RFLU_AllocateMemorySolGv(pRegion)

          pRegion%mixtInput%nCvOld2 = 2
          CALL RFLU_AllocateMemoryAuxVars(pRegion)

          CALL RFLU_ReadGridWrapper(pRegion)  
          CALL RFLU_ReadFlowWrapper(pRegion)      

          IF ( global%verbLevel > VERBOSE_NONE ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)          
          END IF ! global%verbLevel 

          CALL RFLU_SetGasVars(pRegion,1,pRegion%grid%nCellsTot)
          CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCellsTot)

          CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)

          CALL RFLU_InitAuxVars(pRegion)

          CALL RFLU_HM_ConvCvD2ND(pRegion)

          CALL RFLU_WriteGridWrapper(pRegion) 
          CALL RFLU_WriteFlowWrapper(pRegion)                 
          CALL RFLU_WriteAuxVarsWrapper(pRegion)

          CALL RFLU_DeallocateMemorySolDv(pRegion)
          CALL RFLU_DeallocateMemorySolGv(pRegion)
          CALL RFLU_DeallocateMemoryAuxVars(pRegion)
        END DO ! iReg       
      END IF ! global%nRegions    
      
! ------------------------------------------------------------------------------
!   Non-dissipative solver to dissipative solver
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_NONDISS2DISS)
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)          

        CALL RFLU_ReadGridWrapper(pRegion)  
        CALL RFLU_ReadFlowWrapper(pRegion)      

        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)          
        END IF ! global%verbLevel 

        CALL RFLU_HM_ConvCvND2D(pRegion)

        CALL RFLU_ConvertCvPrim2Cons(pRegion,CV_MIXT_STATE_CONS)

        CALL RFLU_WriteGridWrapper(pRegion) 
        CALL RFLU_WriteFlowWrapper(pRegion)
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg) 

          CALL RFLU_ReadGridWrapper(pRegion)  
          CALL RFLU_ReadFlowWrapper(pRegion)      

          IF ( global%verbLevel > VERBOSE_NONE ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)          
          END IF ! global%verbLevel 

          CALL RFLU_HM_ConvCvND2D(pRegion)

          CALL RFLU_ConvertCvPrim2Cons(pRegion,CV_MIXT_STATE_CONS)

          CALL RFLU_WriteGridWrapper(pRegion) 
          CALL RFLU_WriteFlowWrapper(pRegion)                 
        END DO ! iReg       
      END IF ! global%nRegions    
      
! ==============================================================================
!   Convert Rocflu grid to surface grid
! ==============================================================================

! ------------------------------------------------------------------------------
!   Tetmesh/YAMS .msh2 format
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2TETM_GRID) 
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)

        CALL RFLU_ReadGridWrapper(pRegion) 

        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)         
        END IF ! global%verbLevel         

        CALL RFLU_CreateCellMapping(pRegion)
        CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
        CALL RFLU_BuildGlob2LocCellMapping(pRegion)    

        CALL RFLU_CreateBVertexLists(pRegion)
        CALL RFLU_BuildBVertexLists(pRegion)

        CALL RFLU_CreateFaceList(pRegion)                               
        CALL RFLU_BuildFaceList(pRegion)                                         
        CALL RFLU_RenumberBFaceLists(pRegion)                  

        CALL RFLU_ConvROCFLU2TETMESH(pRegion)       
        CALL RFLU_WriteSurfGridTETMESH(pRegion)                

        CALL RFLU_DestroyFaceList(pRegion)                    
        CALL RFLU_DestroyBVertexLists(pRegion)  
        CALL RFLU_DestroyCellMapping(pRegion)                            
      ELSE ! multiple regions, must not happen
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)       
      END IF ! global%nRegions  

! ------------------------------------------------------------------------------
!   STL format
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2STL_GRID) 
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)

        CALL RFLU_ReadGridWrapper(pRegion) 

        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)         
        END IF ! global%verbLevel         

        CALL RFLU_CreateCellMapping(pRegion)
        CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
        CALL RFLU_BuildGlob2LocCellMapping(pRegion)    

        CALL RFLU_CreateBVertexLists(pRegion)
        CALL RFLU_BuildBVertexLists(pRegion)

        CALL RFLU_CreateFaceList(pRegion)                               
        CALL RFLU_BuildFaceList(pRegion)                                         
        CALL RFLU_RenumberBFaceLists(pRegion)                  

        CALL RFLU_CreateGeometry(pRegion)
        CALL RFLU_BuildGeometry(pRegion)                           

        CALL RFLU_STL_WriteSurfGridASCII(pRegion)                

        CALL RFLU_DestroyGeometry(pRegion)
        CALL RFLU_DestroyFaceList(pRegion)
        CALL RFLU_DestroyBVertexLists(pRegion)                      
        CALL RFLU_DestroyCellMapping(pRegion)                            
      ELSE ! multiple regions, must not happen
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)       
      END IF ! global%nRegions  

! ==============================================================================
!   Convert Rocflu steady/unsteady data type to steady/unsteady data type 
! ==============================================================================

! ------------------------------------------------------------------------------
!   Steady rocflu solver to unsteady rocflu solver
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_STEADY2UNSTEADY) 
      WRITE(STDOUT,'(A,1X,A)')    SOLVER_NAME,'Enter time:'

      READ(STDIN,*) writeTime

      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)          

        pRegion%global%flowType = FLOW_STEADY
        pRegion%global%currentIter = readIter

        CALL RFLU_ReadGridWrapper(pRegion)  
        CALL RFLU_ReadFlowWrapper(pRegion)      
        CALL RFLU_ReadPatchCoeffsWrapper(pRegion)      

        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)          
        END IF ! global%verbLevel 

        pRegion%global%flowType = FLOW_UNSTEADY
        pRegion%global%currentTime = writeTime

        CALL RFLU_WriteFlowWrapper(pRegion) 
        CALL RFLU_WritePatchCoeffsWrapper(pRegion) 
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg) 

          pRegion%global%flowType = FLOW_STEADY
          pRegion%global%currentIter = readIter

          CALL RFLU_ReadGridWrapper(pRegion)  
          CALL RFLU_ReadFlowWrapper(pRegion)      
          CALL RFLU_ReadPatchCoeffsWrapper(pRegion)      

          IF ( global%verbLevel > VERBOSE_NONE ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)          
          END IF ! global%verbLevel 

          pRegion%global%flowType = FLOW_UNSTEADY
          pRegion%global%currentTime = writeTime

          CALL RFLU_WriteFlowWrapper(pRegion)                 
          CALL RFLU_WritePatchCoeffsWrapper(pRegion) 
        END DO ! iReg       
      END IF ! global%nRegions    

! ------------------------------------------------------------------------------
!   Unsteady rocflu solver to steady rocflu solver
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_UNSTEADY2STEADY) 
      WRITE(STDOUT,'(A,1X,A)')    SOLVER_NAME,'Enter iteration:'

      READ(STDIN,*) writeIter

      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)          

        pRegion%global%flowType = FLOW_UNSTEADY
        pRegion%global%currentTime = readTime

        CALL RFLU_ReadGridWrapper(pRegion)  
        CALL RFLU_ReadFlowWrapper(pRegion)      
        CALL RFLU_ReadPatchCoeffsWrapper(pRegion)      

        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)          
        END IF ! global%verbLevel 

        pRegion%global%flowType = FLOW_STEADY
        pRegion%global%currentIter = writeIter

        CALL RFLU_WriteFlowWrapper(pRegion) 
        CALL RFLU_WritePatchCoeffsWrapper(pRegion) 
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg) 

          pRegion%global%flowType = FLOW_UNSTEADY
          pRegion%global%currentTime = readTime

          CALL RFLU_ReadGridWrapper(pRegion)  
          CALL RFLU_ReadFlowWrapper(pRegion)      
          CALL RFLU_ReadPatchCoeffsWrapper(pRegion)      

          IF ( global%verbLevel > VERBOSE_NONE ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)          
          END IF ! global%verbLevel 

          pRegion%global%flowType = FLOW_STEADY
          pRegion%global%currentIter = writeIter

          CALL RFLU_WriteFlowWrapper(pRegion)                 
          CALL RFLU_WritePatchCoeffsWrapper(pRegion) 
        END DO ! iReg       
      END IF ! global%nRegions    

! ------------------------------------------------------------------------------
!   Steady rocflu solver to steady rocflu solver
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_STEADY2STEADY) 
      WRITE(STDOUT,'(A,1X,A)')    SOLVER_NAME,'Enter iteration:'

      READ(STDIN,*) writeIter

      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)          

        pRegion%global%flowType = FLOW_STEADY
        pRegion%global%currentIter = readIter

        CALL RFLU_ReadGridWrapper(pRegion)  
        CALL RFLU_ReadFlowWrapper(pRegion)      
        CALL RFLU_ReadPatchCoeffsWrapper(pRegion)      

        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)          
        END IF ! global%verbLevel 

        pRegion%global%currentIter = writeIter

        CALL RFLU_WriteFlowWrapper(pRegion) 
        CALL RFLU_WritePatchCoeffsWrapper(pRegion)      
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg) 

          pRegion%global%flowType = FLOW_STEADY
          pRegion%global%currentIter = readIter

          CALL RFLU_ReadGridWrapper(pRegion)  
          CALL RFLU_ReadFlowWrapper(pRegion)      
          CALL RFLU_ReadPatchCoeffsWrapper(pRegion)      

          IF ( global%verbLevel > VERBOSE_NONE ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)          
          END IF ! global%verbLevel 

          pRegion%global%currentIter = writeIter

          CALL RFLU_WriteFlowWrapper(pRegion)                 
          CALL RFLU_WritePatchCoeffsWrapper(pRegion)      
        END DO ! iReg       
      END IF ! global%nRegions    

! ------------------------------------------------------------------------------
!   Unsteady rocflu solver to unsteady rocflu solver
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_UNSTEADY2UNSTEADY) 
      WRITE(STDOUT,'(A,1X,A)')    SOLVER_NAME,'Enter time:'

      READ(STDIN,*) writeTime

      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)          

        pRegion%global%flowType = FLOW_UNSTEADY
        pRegion%global%currentTime = readTime

        CALL RFLU_ReadGridWrapper(pRegion)  
        CALL RFLU_ReadFlowWrapper(pRegion)      
        CALL RFLU_ReadPatchCoeffsWrapper(pRegion)      

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_ReadAuxVarsWrapper(pRegion)
        END IF ! global%solverType

        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)          
        END IF ! global%verbLevel 

        pRegion%global%currentTime = writeTime

        CALL RFLU_WriteFlowWrapper(pRegion) 
        CALL RFLU_WritePatchCoeffsWrapper(pRegion)      

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! global%solverType
      ELSE ! multiple regions
        DO iReg = 1,global%nRegionsLocal
          pRegion => levels(1)%regions(iReg) 

          pRegion%global%flowType = FLOW_UNSTEADY
          pRegion%global%currentTime = readTime

          CALL RFLU_ReadGridWrapper(pRegion)  
          CALL RFLU_ReadFlowWrapper(pRegion)      
          CALL RFLU_ReadPatchCoeffsWrapper(pRegion)      

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_ReadAuxVarsWrapper(pRegion)
          END IF ! global%solverType

          IF ( global%verbLevel > VERBOSE_NONE ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
            CALL RFLU_PrintFlowInfo(pRegion)          
          END IF ! global%verbLevel 

          pRegion%global%currentTime = writeTime

          CALL RFLU_WriteFlowWrapper(pRegion)                 
          CALL RFLU_WritePatchCoeffsWrapper(pRegion)      

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_WriteAuxVarsWrapper(pRegion)
          END IF ! global%solverType
        END DO ! iReg       
      END IF ! global%nRegions    

! ==============================================================================
!   Convert Rocflu data from one mesh to another 
! ==============================================================================

! ------------------------------------------------------------------------------
!   Rocflu data from larger mesh to smaller mesh 
! ------------------------------------------------------------------------------

    CASE (CONV_RFLU2RFLU_GRID2GRID)
      WRITE(STDOUT,'(A,1X,A)')    SOLVER_NAME,'Enter nCells:'

      READ(STDIN,*) nCells2Copy

      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)          

        CALL RFLU_ReadGridWrapper(pRegion)  
        CALL RFLU_ReadFlowWrapper(pRegion)      
        CALL RFLU_ReadPatchCoeffsWrapper(pRegion)      

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_ReadAuxVarsWrapper(pRegion)
        END IF ! global%solverType

        IF ( global%verbLevel > VERBOSE_NONE ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
          CALL RFLU_PrintFlowInfo(pRegion)          
        END IF ! global%verbLevel 

        pRegion%grid%nCellsTot = nCells2Copy

        CALL RFLU_WriteFlowWrapper(pRegion) 
        CALL RFLU_WritePatchCoeffsWrapper(pRegion)      

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! global%solverType
      ELSE ! multiple regions
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! global%nRegions    

! ==============================================================================
!   Default
! ==============================================================================        

    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! convOption

! ******************************************************************************
! Deallocate memory for solution and grid
! ******************************************************************************

  IF ( global%nRegions == 1 ) THEN ! single region
    pRegion => levels(1)%regions(0)
    
    CALL RFLU_DestroyPatchCoeffs(pRegion)
    CALL RFLU_DeallocMemSolWrapper(pRegion) 
    CALL RFLU_DestroyGrid(pRegion)
  ELSE ! multiple regions
    DO iReg = 1,global%nRegionsLocal
      pRegion => levels(1)%regions(iReg)

      CALL RFLU_DestroyPatchCoeffs(pRegion)
      CALL RFLU_DeallocMemSolWrapper(pRegion) 
      CALL RFLU_DestroyGrid(pRegion)
    END DO ! iReg  
  END IF ! global%nRegions

! *****************************************************************************
! Print info about warnings
! *****************************************************************************

  CALL RFLU_PrintWarnInfo(global)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE rfluconv

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: rfluconv.F90,v $
! Revision 1.2  2015/07/23 23:11:19  brollin
! 1) The pressure coefficient of the  collision model has been changed back to its original form
! 2) New options in the format of the grid and solutions have been added. Now the user can choose the endianness, and convert from one to the over in rfluconv.
! 3) The solutions are now stored in folders named by timestamp or iteration number.
! 4) The address enty in the hashtable has been changed to an integer(8) for cases when the grid becomes very large.
! 5) RFLU_WritePM can now compute PM2 on the fly for the Macroscale problem
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.7  2010/03/15 00:29:25  mparmar
! Added CONV_RFLU2RFLU_GRID2GRID, added support for auxilliary data in UNSTEADY2UNSTEADY
!
! Revision 1.6  2009/07/08 19:12:22  mparmar
! Added aditional options
!
! Revision 1.5  2008/12/06 08:43:54  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:06  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/11/28 23:05:35  mparmar
! Reading/Writing auxilliary variables for SOLV_IMPLICIT_HM
!
! Revision 1.2  2007/06/18 18:13:08  mparmar
! Added option to convert patch coefficient files
!
! Revision 1.1  2007/04/09 18:54:50  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/02/06 23:55:54  haselbac
! Added comm argument to RFLU_InitGlobal
!
! Revision 1.3  2005/07/07 03:51:09  haselbac
! Added STL conversion, some clean-up
!
! Revision 1.2  2005/05/03 03:09:09  haselbac
! Converted to C++ reading of command-line, fixed bug in formats
!
! Revision 1.1  2005/04/18 14:57:56  haselbac
! Initial revision
!
! ******************************************************************************

