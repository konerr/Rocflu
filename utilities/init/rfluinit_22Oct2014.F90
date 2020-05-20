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
! Purpose: Driver routine for rfluinit. 
!
! Description: None.
!
! Input: 
!   caseString  String with casename
!   verbLevel   Verbosity level
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: rfluinit_22Oct2014.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE rfluinit(caseString,verbLevel)

  USE ModError
  USE ModDataTypes
  USE ModParameters  
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModMPI
#ifdef PLAG
  USE ModPartLag, ONLY: t_plag
#endif
  
  USE RFLU_ModAllocateMemory
  USE RFLU_ModAxisymmetry, ONLY: RFLU_AXI_ScaleGeometry
  USE RFLU_ModBoundLists
  USE RFLU_ModBoundXvUtils 
  USE RFLU_ModNSCBC, ONLY: RFLU_NSCBC_DecideHaveNSCBC
  USE RFLU_ModCellMapping
  USE RFLU_ModDeallocateMemory
  USE RFLU_ModDimensions
  USE RFLU_ModFaceList
  USE RFLU_ModGeometry
  USE RFLU_ModHouMahesh, ONLY: RFLU_HM_ConvCvD2ND
  USE RFLU_ModGridSpeedUtils
  USE RFLU_ModMovingFrame
  USE RFLU_ModReadBcInputFile
  USE RFLU_ModReadWriteAuxVars
  USE RFLU_ModReadWriteBcDataFile
  USE RFLU_ModReadWriteFlow
  USE RFLU_ModReadWriteGrid
  USE RFLU_ModReadWriteGridSpeeds
  USE RFLU_ModRegionMapping
  USE RFLU_ModRenumberings

#ifdef GENX
  USE RFLU_ModGENXAdmin
  USE RFLU_ModGENXIO
#endif

  USE ModInterfaces, ONLY: RFLU_AllocMemSolWrapper, & 
                           RFLU_BuildDataStruct, &   
                           RFLU_ComputeGridSpacingCyldet, &
                           RFLU_CreateGrid, & 
                           RFLU_DeallocMemSolWrapper, & 
                           RFLU_DestroyGrid, &
                           RFLU_GetUserInput, &
                           RFLU_InitAuxVars, &
                           RFLU_InitGlobal, &
                           RFLU_InitBcDataHardCode, & 
                           RFLU_InitFlowHardCodeLimWrapper, &  
                           RFLU_InitFlowHardCodeWrapper, & 
                           RFLU_InitFlowScratchWrapper, &
                           RFLU_InitFlowSerialWrapper, & 
                           RFLU_PrintFlowInfoWrapper, & 
                           RFLU_PrintGridInfo, & 
                           RFLU_PrintHeader, & 
                           RFLU_PrintWarnInfo, &
                           RFLU_RandomInit, &
                           RFLU_ReadRestartInfo, & 
                           RFLU_SetDependentVars, &
                           RFLU_SetModuleType, &  
                           RFLU_SetRestartTimeFlag, &
                           RFLU_SetVarInfoWrapper, & 
                           RFLU_WriteVersionString, &
                           ScaleRotateVector
                           
#ifdef GENX
  USE ModInterfaces, ONLY: RFLU_SetGasVars
#endif                  

#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY: PLAG_InitPatchData, &
                                     PLAG_RFLU_AllocMemSol, & 
                                     PLAG_RFLU_AllocMemSolTile, & 
                                     PLAG_RFLU_ComputeCellsContainingPcls, &
                                     PLAG_RFLU_ComputeLagPclsPerReg, &
                                     PLAG_RFLU_ComputeVolFrac, &
                                     PLAG_RFLU_DeallocMemSol, & 
                                     PLAG_RFLU_DeallocMemSolTile, &
                                     PLAG_RFLU_InitSolutionCyldet, &
                                     PLAG_RFLU_InitSolutionFile, &
                                     PLAG_RFLU_InitSolutionHardcode, &
                                     PLAG_RFLU_InitSolutionRandom, &
                                     PLAG_RFLU_InitSolutionScratch, &
                                     PLAG_RFLU_InitSolFromSerial, &
                                     PLAG_RFLU_InitSolFromSerialCopy, &
                                     PLAG_RFLU_InitSolSerialWrapper, & 
                                     PLAG_RFLU_ModifyFlowFieldVolFrac, &
                                     PLAG_RFLU_WriteSolutionASCII, &
                                     PLAG_RFLU_WriteSolutionBinary, &
                                     PLAG_RFLU_WriteUnsteadyDataASCII, &
                                     PLAG_RFLU_WriteUnsteadyDataBinary

  USE PLAG_ModDataStruct
  USE PLAG_ModDimensions,      ONLY: PLAG_RFLU_WriteDimensions, & 
                                     PLAG_SetDimensions, &
                                     PLAG_SetMaxDimensions
#endif

  IMPLICIT NONE

#ifdef GENX
  INCLUDE 'roccomf90.h'
#endif

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(*) :: caseString
  INTEGER, INTENT(IN) :: verbLevel

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: casename
#ifdef GENX
  CHARACTER(CHRLEN) :: surfWinNameInput,volWinNameInput
#endif
  INTEGER :: errorFlag,iLev,iReg,iRegLow,iRegUpp
#ifdef GENX
  INTEGER :: handleObtain
#endif
  TYPE(t_region), POINTER :: pRegion,pRegionSerial
  TYPE(t_global), POINTER :: global
  TYPE(t_level), POINTER :: levels(:)
#ifdef PLAG
  TYPE(t_plag), POINTER :: pPlag,pPlagSerial
  ! Subbu
  INTEGER :: iRegLoc,nCellsSeedPcls,nPclsSumReg
  ! Subbu
#endif

! ******************************************************************************
! Initialize global data
! ******************************************************************************  
  
  ALLOCATE(global,STAT=errorFlag)
  IF ( errorFlag /= ERR_NONE ) THEN 
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR - Pointer allocation failed.'
    STOP
  END IF ! errorFlag 

  casename = caseString(1:LEN(caseString))
  
  CALL RFLU_InitGlobal(casename,verbLevel,MPI_COMM_WORLD,global)

  CALL RegisterFunction(global,'rfluinit', & 
                        __FILE__)

  CALL RFLU_SetModuleType(global,MODULE_TYPE_INIT)

#ifdef GENX
! ******************************************************************************
! Read GENX control file, store communicator and hardcode window name
! ****************************************************************************** 

  CALL RFLU_GENX_ReadCtrlFile(global)
  CALL RFLU_GENX_StoreCommunicator(global,MPI_COMM_WORLD)
  CALL RFLU_GENX_HardCodeWindowName(global)
#endif

! ******************************************************************************
! Print header and write version string
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC ) THEN
    CALL RFLU_WriteVersionString(global)     
    IF ( global%verbLevel /= VERBOSE_NONE ) THEN
      CALL RFLU_PrintHeader(global)
    END IF ! global%verbLevel
  END IF ! global%myProcid

! ******************************************************************************
! Read mapping file, impose serial mapping, and build basic data structure
! ******************************************************************************

  CALL RFLU_ReadRegionMappingFile(global,MAPFILE_READMODE_PEEK,global%myProcId)
  CALL RFLU_SetRegionMappingSerial(global)  
  CALL RFLU_CreateRegionMapping(global,MAPTYPE_REG)
  CALL RFLU_ImposeRegionMappingSerial(global)

  CALL RFLU_BuildDataStruct(global,levels) 
  CALL RFLU_ApplyRegionMapping(global,levels)
  CALL RFLU_DestroyRegionMapping(global,MAPTYPE_REG)  

#ifdef GENX
! ******************************************************************************
! Initialize Roccom and load Rocin and Rocout
! ******************************************************************************

  CALL COM_init
! TEMPORARY
!  CALL COM_set_verbose(10)
! END TEMPORARY

  surfWinNameInput = 'RocfluInputSurf'
  volWinNameInput  = 'RocfluInputVol TEMP'

  global%winNameIn  = TRIM(global%winName)//'-IN'
  global%winNameOut = TRIM(global%winName)//'-OUT'

  CALL Rocin_load_module(TRIM(global%winNameIn))
  CALL Rocout_load_module(TRIM(global%winNameOut))
  
  handleObtain = COM_get_function_handle(TRIM(global%winNameIn)// &
                                         '.obtain_attribute')
  CALL RFLU_GENX_StoreNamesHandles(global,surfWinNameInput,volWinNameInput, & 
                                   handleObtain)
  
! ******************************************************************************
! Create windows and attributes 
! ******************************************************************************
  
  pRegionSerial => levels(1)%regions(0)

  CALL RFLU_GENX_CreateWindows(pRegionSerial) 
  CALL RFLU_GENX_CreateAttrGridSurf(pRegionSerial)
  CALL RFLU_GENX_CreateAttrFlow(pRegionSerial)    
  CALL RFLU_GENX_CreateAttrGSpeeds(pRegionSerial)
  CALL RFLU_GENX_CreateAttrInterf(pRegionSerial)
#endif

! ******************************************************************************
! Initialize random number generator. NOTE needed in order to write sensible 
! data when writing Rocpart solution files. 
! ******************************************************************************

  CALL RFLU_RandomInit(levels(1)%regions)

! ******************************************************************************
! Read input file and restart info. NOTE need restart info for GENX runs to 
! determine whether have a restart. 
! ******************************************************************************

  CALL RFLU_GetUserInput(levels(1)%regions,.TRUE.) 
  CALL RFLU_ReadRestartInfo(global)
  CALL RFLU_SetRestartTimeFlag(global)

! ******************************************************************************
! Initialize solutions
! ******************************************************************************

  IF ( global%nRegions == 1 ) THEN 
    iRegLow = 0
    iRegUpp = 0
  ELSE 
    iRegLow = 1
    iRegUpp = global%nRegions
  END IF ! global%nRegions

#ifdef PLAG
! ==============================================================================
! Temporarily disable particles so can use wrapper routines to initialize only
! Eulerian solution fields
! ==============================================================================

  global%plagUsedSave = global%plagUsed

  global%plagUsed = .FALSE.
#endif

! ==============================================================================
! Creating and init particle velocity and acceleration
! ==============================================================================

  IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
    DO iReg = iRegLow,iRegUpp
      pRegion => levels(1)%regions(iReg)

      CALL RFLU_MVF_CreatePatchVelAccel(pRegion)
      CALL RFLU_MVF_InitPatchVelAccel(pRegion)
    END DO ! iReg
  END IF ! global%mvFrameFlag

! ==============================================================================
! Initialize Eulerian solution fields
! ==============================================================================

  SELECT CASE ( global%initFlowFlag ) 
    
! ------------------------------------------------------------------------------
!   Initialize from scratch
! ------------------------------------------------------------------------------

    CASE ( INITFLOW_FROMSCRATCH )
      DO iReg = iRegLow,iRegUpp
        pRegion => levels(1)%regions(iReg)

        CALL RFLU_ReadDimensions(pRegion)                          
        CALL RFLU_CreateGrid(pRegion)

        IF ( pRegion%grid%nPatches > 0 ) THEN        
          CALL RFLU_ReadBCInputFileWrapper(pRegion)    
        END IF ! pRegion%grid%nPatches                   

        IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_ReadGridWrapper(pRegion)

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
        END IF ! RFLU_NSCBC_DecideHaveNSCBC

        CALL RFLU_AllocMemSolWrapper(pRegion)  
        CALL RFLU_SetVarInfoWrapper(pRegion) 

        CALL RFLU_InitFlowScratchWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_InitAuxVars(pRegion)
        END IF ! solverType

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_HM_ConvCvD2ND(pRegion)
        END IF ! solverType

        CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCellsTot)

! Initializing boundary array would need solution in domain to be initialized first
        IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_BXV_CreateVarsCv(pRegion)
          CALL RFLU_BXV_CreateVarsDv(pRegion)
          CALL RFLU_BXV_InitVars(pRegion)
          CALL RFLU_BXV_WriteVarsWrapper(pRegion)
          CALL RFLU_BXV_DestroyVarsCv(pRegion)
          CALL RFLU_BXV_DestroyVarsDv(pRegion)
        END IF ! RFLU_NSCBC_DecideHaveNSCBC

! Modify flow field for moving reference frame attached to moving particle
        IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
          CALL RFLU_MVF_ModifyFlowField(pRegion)
        END IF ! global%mvFrameFlag

        IF ( global%verbLevel > VERBOSE_NONE ) THEN   
          CALL RFLU_PrintFlowInfoWrapper(pRegion)    
        END IF ! global%verbLevel       

#ifdef GENX
        CALL RFLU_GENX_SetConnSize(pRegion)
        CALL RFLU_GENX_RegisterDataFlow(pRegion)               
#endif
        CALL RFLU_WriteFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! solverType

        IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_DestroyFaceList(pRegion)
          CALL RFLU_DestroyBVertexLists(pRegion)    
          CALL RFLU_DestroyCellMapping(pRegion)
          CALL RFLU_DestroyGeometry(pRegion)  
        END IF ! RFLU_NSCBC_DecideHaveNSCBC

        CALL RFLU_DeallocMemSolWrapper(pRegion)
        CALL RFLU_DestroyGrid(pRegion)
      END DO ! iReg  

! ------------------------------------------------------------------------------
!   Initialize parallel run by reading solution from serial file.
! ------------------------------------------------------------------------------
  
    CASE ( INITFLOW_FROMFILE ) 
      IF ( global%nRegions > 1 ) THEN 
        pRegionSerial => levels(1)%regions(0)

! ----- Read serial solution ---------------------------------------------------

        CALL RFLU_ReadDimensions(pRegionSerial)               
        CALL RFLU_CreateGrid(pRegionSerial)

        IF ( pRegionSerial%grid%nPatches > 0 ) THEN        
          CALL RFLU_ReadBCInputFileWrapper(pRegionSerial)    
        END IF ! pRegionSerial%grid%nPatches         

        CALL RFLU_AllocMemSolWrapper(pRegionSerial)  
        CALL RFLU_SetVarInfoWrapper(pRegionSerial)     

        CALL RFLU_ReadFlowWrapper(pRegionSerial)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_ReadAuxVarsWrapper(pRegionSerial)
        END IF ! solverType

! ----- Loop over regions and initialize ---------------------------------------

        DO iReg = 1,global%nRegions
          pRegion => levels(1)%regions(iReg)

          CALL RFLU_ReadDimensionsWrapper(pRegion)               
          CALL RFLU_CreateGrid(pRegion)

          IF ( pRegion%grid%nPatches > 0 ) THEN        
            CALL RFLU_ReadBCInputFileWrapper(pRegion)    
          END IF ! pRegion%grid%nPatches         

          CALL RFLU_AllocMemSolWrapper(pRegion)  
          CALL RFLU_SetVarInfoWrapper(pRegion)

          CALL RFLU_RNMB_CreatePC2SCMap(pRegion)    
          CALL RFLU_RNMB_CreatePV2SVMap(pRegion)
          CALL RFLU_RNMB_CreatePBF2SBFMap(pRegion)

          CALL RFLU_RNMB_ReadPxx2SxxMaps(pRegion)

          CALL RFLU_RNMB_DestroyPV2SVMap(pRegion)
          CALL RFLU_RNMB_DestroyPBF2SBFMap(pRegion)        

          CALL RFLU_InitFlowSerialWrapper(pRegion,pRegionSerial)

! Modify flow field for moving reference frame attached to moving particle
          IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
            CALL RFLU_MVF_ModifyFlowField(pRegion)
          END IF ! global%mvFrameFlag

          IF ( global%verbLevel > VERBOSE_NONE ) THEN   
            CALL RFLU_PrintFlowInfoWrapper(pRegion)    
          END IF ! global%verbLevel  

#ifdef GENX
          CALL RFLU_GENX_RegisterDataFlow(pRegion)
          CALL RFLU_SetGasVars(pRegion,1,pRegion%grid%nCellsTot)               
          CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCellsTot)
#endif
          CALL RFLU_WriteFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_WriteAuxVarsWrapper(pRegion)
          END IF ! solverType

          CALL RFLU_RNMB_DestroyPC2SCMap(pRegion)               
          CALL RFLU_DeallocMemSolWrapper(pRegion) 
          CALL RFLU_DestroyGrid(pRegion) 
        END DO ! iReg

! ----- Deallocate memory ------------------------------------------------------

        CALL RFLU_DeallocMemSolWrapper(pRegionSerial) 
        CALL RFLU_DestroyGrid(pRegionSerial)    
      END IF ! global%nRegions

! ------------------------------------------------------------------------------
!   Initialize from hardcode
! ------------------------------------------------------------------------------

    CASE ( INITFLOW_FROMHARDCODE ) 
      DO iReg = iRegLow,iRegUpp
        pRegion => levels(1)%regions(iReg)

        CALL RFLU_ReadDimensions(pRegion)               
        CALL RFLU_CreateGrid(pRegion)

        IF ( pRegion%grid%nPatches > 0 ) THEN        
          CALL RFLU_ReadBCInputFileWrapper(pRegion)    
        END IF ! pRegion%grid%nPatches       

#ifdef GENX
        CALL RFLU_GENX_ReadWindow(pRegion,GENX_WINDOW_TYPE_SURF)
        CALL RFLU_GENX_ReadWindow(pRegion,GENX_WINDOW_TYPE_VOL)
        CALL RFLU_GENX_GetDimensionsDerived(pRegion)
        CALL RFLU_GENX_CreateGridSurf(pRegion) 
        CALL RFLU_GENX_RegisterGridSurf(pRegion) 
        CALL RFLU_GENX_RegisterGridVol(pRegion)    
#endif

        CALL RFLU_ReadGridWrapper(pRegion)

        IF ( global%verbLevel > VERBOSE_LOW ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
        END IF ! global%verbLevel 

#ifdef GENX
        CALL RFLU_GENX_DestroyGridSurf(pRegion) 
#endif

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

        CALL RFLU_AllocMemSolWrapper(pRegion)  
        CALL RFLU_SetVarInfoWrapper(pRegion) 
        
        CALL RFLU_InitFlowHardCodeWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_HM_ConvCvD2ND(pRegion)
        END IF ! solverType

        IF ( RFLU_DecideReadWriteBcDataFile(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_InitBcDataHardCode(pRegion)      
        END IF ! RFLU_DecideReadWriteBcDataFile      

        CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCellsTot)

        IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_BXV_CreateVarsCv(pRegion)
          CALL RFLU_BXV_CreateVarsDv(pRegion)
          CALL RFLU_BXV_InitVars(pRegion)
          CALL RFLU_BXV_WriteVarsWrapper(pRegion)
          CALL RFLU_BXV_DestroyVarsCv(pRegion)
          CALL RFLU_BXV_DestroyVarsDv(pRegion)
        END IF ! RFLU_NSCBC_DecideHaveNSCBC

! Modify flow field for moving reference frame attached to moving particle
        IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
          CALL RFLU_MVF_ModifyFlowField(pRegion)
        END IF ! global%mvFrameFlag

        IF ( global%verbLevel > VERBOSE_NONE ) THEN   
          CALL RFLU_PrintFlowInfoWrapper(pRegion)    
        END IF ! global%verbLevel       

#ifdef GENX
        CALL RFLU_GENX_RegisterDataFlow(pRegion)              
        CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCellsTot)
#endif
        CALL RFLU_WriteFlowWrapper(pRegion)       

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! solverType

        IF ( RFLU_DecideReadWriteBcDataFile(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_WriteBcDataFile(pRegion)
        END IF ! RFLU_DecideReadWriteBcDataFile

        CALL RFLU_DestroyGeometry(pRegion)  
        CALL RFLU_DestroyFaceList(pRegion)
        CALL RFLU_DestroyBVertexLists(pRegion)    
        CALL RFLU_DestroyCellMapping(pRegion)
        CALL RFLU_DeallocMemSolWrapper(pRegion)
        CALL RFLU_DestroyGrid(pRegion)                                  
      END DO ! iReg   

! ------------------------------------------------------------------------------
!   Initialize parallel run from combo using serial solution
! ------------------------------------------------------------------------------
  
    CASE ( INITFLOW_FROMCOMBO_SERIAL ) 
      IF ( global%nRegions > 1 ) THEN 
        pRegionSerial => levels(1)%regions(0)

! ----- Read serial solution ---------------------------------------------------

        CALL RFLU_ReadDimensions(pRegionSerial)               
        CALL RFLU_CreateGrid(pRegionSerial)

        IF ( pRegionSerial%grid%nPatches > 0 ) THEN        
          CALL RFLU_ReadBCInputFileWrapper(pRegionSerial)    
        END IF ! pRegionSerial%grid%nPatches         

        CALL RFLU_AllocMemSolWrapper(pRegionSerial)  
        CALL RFLU_SetVarInfoWrapper(pRegionSerial)     

        CALL RFLU_ReadFlowWrapper(pRegionSerial)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_ReadAuxVarsWrapper(pRegionSerial)
        END IF ! solverType

! ----- Loop over regions and initialize ---------------------------------------

        DO iReg = 1,global%nRegions
          pRegion => levels(1)%regions(iReg)

          CALL RFLU_ReadDimensionsWrapper(pRegion)               
          CALL RFLU_CreateGrid(pRegion)

          IF ( pRegion%grid%nPatches > 0 ) THEN        
            CALL RFLU_ReadBCInputFileWrapper(pRegion)    
          END IF ! pRegion%grid%nPatches         

          CALL RFLU_ReadGridWrapper(pRegion)

          IF ( global%verbLevel > VERBOSE_LOW ) THEN 
            CALL RFLU_PrintGridInfo(pRegion)
          END IF ! global%verbLevel 

          CALL RFLU_AllocMemSolWrapper(pRegion)  
          CALL RFLU_SetVarInfoWrapper(pRegion)

          CALL RFLU_RNMB_CreatePC2SCMap(pRegion)    
          CALL RFLU_RNMB_CreatePV2SVMap(pRegion)
          CALL RFLU_RNMB_CreatePBF2SBFMap(pRegion)

          CALL RFLU_RNMB_ReadPxx2SxxMaps(pRegion)

          CALL RFLU_RNMB_DestroyPV2SVMap(pRegion)
          CALL RFLU_RNMB_DestroyPBF2SBFMap(pRegion)        

          CALL RFLU_InitFlowSerialWrapper(pRegion,pRegionSerial)
          
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
                  
          CALL RFLU_InitFlowHardCodeLimWrapper(pRegion)

! Modify flow field for moving reference frame attached to moving particle
          IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
            CALL RFLU_MVF_ModifyFlowField(pRegion)
          END IF ! global%mvFrameFlag

          CALL RFLU_DestroyGeometry(pRegion)  
          CALL RFLU_DestroyFaceList(pRegion)
          CALL RFLU_DestroyBVertexLists(pRegion)    
          CALL RFLU_DestroyCellMapping(pRegion)

          IF ( global%verbLevel > VERBOSE_NONE ) THEN   
            CALL RFLU_PrintFlowInfoWrapper(pRegion)    
          END IF ! global%verbLevel  

#ifdef GENX
          CALL RFLU_GENX_RegisterDataFlow(pRegion)               
          CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCellsTot)
#endif
          CALL RFLU_WriteFlowWrapper(pRegion)

          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            CALL RFLU_WriteAuxVarsWrapper(pRegion)
          END IF ! solverType

          CALL RFLU_RNMB_DestroyPC2SCMap(pRegion)               
          CALL RFLU_DeallocMemSolWrapper(pRegion) 
          CALL RFLU_DestroyGrid(pRegion) 
        END DO ! iReg

! ----- Deallocate memory ------------------------------------------------------

        CALL RFLU_DeallocMemSolWrapper(pRegionSerial) 
        CALL RFLU_DestroyGrid(pRegionSerial)    
      END IF ! global%nRegions

! ------------------------------------------------------------------------------
!   Initialize parallel run from combo using parallel solution
! ------------------------------------------------------------------------------
  
    CASE ( INITFLOW_FROMCOMBO_PARALLEL ) 

! --- Loop over regions and initialize -----------------------------------------

      DO iReg = 1,global%nRegions
        pRegion => levels(1)%regions(iReg)

        CALL RFLU_ReadDimensionsWrapper(pRegion)               
        CALL RFLU_CreateGrid(pRegion)

        IF ( pRegion%grid%nPatches > 0 ) THEN        
          CALL RFLU_ReadBCInputFileWrapper(pRegion)    
        END IF ! pRegion%grid%nPatches         

        CALL RFLU_ReadGridWrapper(pRegion)

        IF ( global%verbLevel > VERBOSE_LOW ) THEN 
          CALL RFLU_PrintGridInfo(pRegion)
        END IF ! global%verbLevel 

        CALL RFLU_AllocMemSolWrapper(pRegion)  
        CALL RFLU_SetVarInfoWrapper(pRegion)

        CALL RFLU_ReadFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_ReadAuxVarsWrapper(pRegion)
        END IF ! solverType

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

        CALL RFLU_InitFlowHardCodeLimWrapper(pRegion)

! Modify flow field for moving reference frame attached to moving particle
        IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
          CALL RFLU_MVF_ModifyFlowField(pRegion)
        END IF ! global%mvFrameFlag

        CALL RFLU_DestroyGeometry(pRegion)  
        CALL RFLU_DestroyFaceList(pRegion)
        CALL RFLU_DestroyBVertexLists(pRegion)    
        CALL RFLU_DestroyCellMapping(pRegion)

        IF ( global%verbLevel > VERBOSE_NONE ) THEN   
          CALL RFLU_PrintFlowInfoWrapper(pRegion)    
        END IF ! global%verbLevel  

#ifdef GENX
        CALL RFLU_GENX_RegisterDataFlow(pRegion)               
        CALL RFLU_SetDependentVars(pRegion,1,pRegion%grid%nCellsTot)
#endif
        CALL RFLU_WriteFlowWrapper(pRegion)
      
        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! solverType

        CALL RFLU_DeallocMemSolWrapper(pRegion) 
        CALL RFLU_DestroyGrid(pRegion) 
      END DO ! iReg    

! ------------------------------------------------------------------------------
!   Default
! ------------------------------------------------------------------------------
  
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%initFlowFlag  

! ==============================================================================
! Writing and destroying particle velocity and acceleration
! ==============================================================================

  IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
    pRegion => levels(1)%regions(iRegLow)

    IF ( pRegion%global%myProcid == MASTERPROC ) THEN
      CALL RFLU_MVF_WritePatchVelAccel(pRegion)
    END IF ! pRegion%global%myProcid

    DO iReg = iRegLow,iRegUpp
      pRegion => levels(1)%regions(iReg)

      CALL RFLU_MVF_DestroyPatchVelAccel(pRegion)
    END DO ! iReg
  END IF ! global%mvFrameFlag

#ifdef PLAG
! ==============================================================================
! Re-enable particles
! ==============================================================================

  global%plagUsed = global%plagUsedSave

! ==============================================================================
! Initialize particle solution 
! ==============================================================================

  IF ( global%plagUsed .EQV. .TRUE. ) THEN

! ==============================================================================
! Initialize if only one region
! ==============================================================================

    IF ( global%nRegions == 1 ) THEN 
      pRegion => levels(1)%regions(0)
      pPlag   => pRegion%plag
      nCellsSeedPcls = 0

      CALL RFLU_ReadDimensions(pRegion)
      CALL RFLU_CreateGrid(pRegion)

      IF ( pRegion%grid%nPatches > 0 ) THEN
        CALL RFLU_ReadBCInputFileWrapper(pRegion)
      END IF ! pRegion%grid%nPatches

      CALL RFLU_ReadGridWrapper(pRegion)

      IF ( global%verbLevel > VERBOSE_LOW ) THEN
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

! ------------------------------------------------------------------------------
!  Compute how many cells in each region can be seeded with particles
!  Sum them over all regions and thereby evaluate how many particles need to be
!  added in each cell (global%nPclsTarget)
! ------------------------------------------------------------------------------

      CALL PLAG_RFLU_ComputeCellsContainingPcls(pRegion)
      global%nPclsTarget = FLOAT(pRegion%plagInput%nPclsIni)/&
                           FLOAT(pRegion%grid%nCellsSeedPcls)

      ! Compute nPcls in each region 
      CALL PLAG_RFLU_ComputeLagPclsPerReg(pRegion)

      ! Allocating memory for particle field
      CALL PLAG_SetDimensions(pRegion,pPlag%nPcls)
      CALL PLAG_SetMaxDimensions(pRegion)
      CALL PLAG_RFLU_AllocMemSol(pRegion,pPlag)
      CALL PLAG_RFLU_AllocMemSolTile(pRegion)

      CALL RFLU_ComputeGridSpacingCyldet(pRegion)

      ! Initialize particle field - CCMT case
      nPclsSumReg = 0
      CALL PLAG_RFLU_InitSolutionCyldet(pRegion,nPclsSumReg)

      CALL PLAG_InitPatchData(pRegion)

      CALL RFLU_AllocMemSolWrapper(pRegion)
      CALL RFLU_SetVarInfoWrapper(pRegion)
      pRegion%global%plagUsed = .FALSE.
      CALL RFLU_ReadFlowWrapper(pRegion)
      pRegion%global%plagUsed = .TRUE.

      IF ( pRegion%mixtInput%axiFlag .EQV. .TRUE. ) THEN
        CALL RFLU_AXI_ScaleGeometry(pRegion)
      END IF ! pRegion%mixtInput%axiFlag

      CALL PLAG_RFLU_ComputeVolFrac(pRegion)
      CALL PLAG_RFLU_ModifyFlowFieldVolFrac(pRegion)

      CALL PLAG_SetMaxDimensions(pRegion)
      CALL PLAG_RFLU_WriteDimensions(pRegion)

      CALL RFLU_WriteFlowWrapper(pRegion)

! ------------------------------------------------------------------------------
!  Deallocate memory for grid, geometry and solution (Eulerian & Lagrangian) 
! ------------------------------------------------------------------------------

      CALL PLAG_RFLU_DeallocMemSolTile(pRegion)
      CALL PLAG_RFLU_DeallocMemSol(pRegion,pPlag)

      CALL RFLU_DeallocMemSolWrapper(pRegion)

      CALL RFLU_DestroyGeometry(pRegion)
      CALL RFLU_DestroyFaceList(pRegion)
      CALL RFLU_DestroyBVertexLists(pRegion)
      CALL RFLU_DestroyCellMapping(pRegion)

      CALL RFLU_DestroyGrid(pRegion)

    END IF ! global%nRegions=1

! ==============================================================================
! Initialize if mulitple regions
! ==============================================================================

    IF ( global%nRegions > 1 ) THEN
      ! Need to understand the purpose of the following functions & where they
      ! are being utilized
      ! Temporarily turn them off

      !CALL PLAG_DSTR_CreatePclListCSR(pRegionSerial)
      !CALL PLAG_DSTR_BuildCell2PclList(pRegionSerial)

! ------------------------------------------------------------------------------
!  Allocate memory and read dimensions, build grid, cell mapping, vertex list, 
!  face list and geometry
! ------------------------------------------------------------------------------

      nCellsSeedPcls = 0
      DO iReg = 1,global%nRegions 
        pRegion => levels(1)%regions(iReg)
        pPlag   => pRegion%plag

        CALL RFLU_ReadDimensions(pRegion)               
        CALL RFLU_CreateGrid(pRegion)       

        IF ( pRegion%grid%nPatches > 0 ) THEN        
          CALL RFLU_ReadBCInputFileWrapper(pRegion)    
        END IF ! pRegion%grid%nPatches

        CALL RFLU_ReadGridWrapper(pRegion)

        IF ( global%verbLevel > VERBOSE_LOW ) THEN 
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

! ------------------------------------------------------------------------------
!  Compute how many cells in each region can be seeded with particles
!  Sum them over all regions and thereby evaluate how many particles need to be
!  added in each cell (global%nPclsTarget)
! ------------------------------------------------------------------------------

        CALL PLAG_RFLU_ComputeCellsContainingPcls(pRegion)
        nCellsSeedPcls = nCellsSeedPcls + pRegion%grid%nCellsSeedPcls

! ------------------------------------------------------------------------------
!  Deallocate memory for grid and geometry  
! ------------------------------------------------------------------------------

        CALL RFLU_DestroyGeometry(pRegion)
        CALL RFLU_DestroyFaceList(pRegion)
        CALL RFLU_DestroyBVertexLists(pRegion)
        CALL RFLU_DestroyCellMapping(pRegion)

        CALL RFLU_DestroyGrid(pRegion)

      END DO

      global%nPclsTarget = FLOAT(pRegion%plagInput%nPclsIni)/FLOAT(nCellsSeedPcls)
      WRITE(*,*) 'nPclsTarget,nCellsSeedPcls=',global%nPclsTarget,nCellsSeedPcls
      IF (global%nPclsTarget .LT. 0.0_RFREAL) THEN
        CALL ErrorStop(global,ERR_PLAG_NPCLSTARGET_INVALID,__LINE__)
      END IF

! ------------------------------------------------------------------------------
!  Allocate memory for particle solution, initialize pcl Lagrangian field,
!  compute volume fraction and write modified Eulerian flow-field, particle
!  field including (.plag.sol, .pdim files etc.)
! ------------------------------------------------------------------------------

      DO iReg = 1,global%nRegions 
        pRegion => levels(1)%regions(iReg)
        pPlag   => pRegion%plag

        CALL RFLU_ReadDimensions(pRegion)               
        CALL RFLU_CreateGrid(pRegion)       

        IF ( pRegion%grid%nPatches > 0 ) THEN        
          CALL RFLU_ReadBCInputFileWrapper(pRegion)    
        END IF ! pRegion%grid%nPatches

        CALL RFLU_ReadGridWrapper(pRegion)

        IF ( global%verbLevel > VERBOSE_LOW ) THEN 
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
 
        ! Compute nPcls in each region 
        CALL PLAG_RFLU_ComputeLagPclsPerReg(pRegion)      
        
        ! Allocating memory for particle field
        CALL PLAG_SetDimensions(pRegion,pPlag%nPcls)
        CALL PLAG_SetMaxDimensions(pRegion)
        CALL PLAG_RFLU_AllocMemSol(pRegion,pPlag)
        CALL PLAG_RFLU_AllocMemSolTile(pRegion)

        CALL RFLU_ComputeGridSpacingCyldet(pRegion)

        ! Initialize particle field - CCMT case
        IF (pRegion%grid%initPclPresent .EQV. .TRUE.) THEN
          nPclsSumReg = 0
          IF (pRegion%iRegionGlobal > 1) THEN
            DO iRegLoc = 1,pRegion%iRegionGlobal - 1
              nPclsSumReg = nPclsSumReg + levels(1)%regions(iRegLoc)%plag%nPcls
            END DO 
          END IF
          CALL PLAG_RFLU_InitSolutionCyldet(pRegion,nPclsSumReg)
        END IF

        CALL PLAG_InitPatchData(pRegion)
  
        CALL RFLU_AllocMemSolWrapper(pRegion)  
        CALL RFLU_SetVarInfoWrapper(pRegion)
        pRegion%global%plagUsed = .FALSE.
        CALL RFLU_ReadFlowWrapper(pRegion)
        pRegion%global%plagUsed = .TRUE.

        IF ( pRegion%mixtInput%axiFlag .EQV. .TRUE. ) THEN
          CALL RFLU_AXI_ScaleGeometry(pRegion)
        END IF ! pRegion%mixtInput%axiFlag

        CALL PLAG_RFLU_ComputeVolFrac(pRegion)
        CALL PLAG_RFLU_ModifyFlowFieldVolFrac(pRegion)

        CALL PLAG_SetMaxDimensions(pRegion)
        CALL PLAG_RFLU_WriteDimensions(pRegion)

        CALL RFLU_WriteFlowWrapper(pRegion)

! ------------------------------------------------------------------------------
!  Deallocate memory for grid, geometry and solution (Eulerian & Lagrangian) 
! ------------------------------------------------------------------------------

        CALL PLAG_RFLU_DeallocMemSolTile(pRegion)
        CALL PLAG_RFLU_DeallocMemSol(pRegion,pPlag)

        CALL RFLU_DeallocMemSolWrapper(pRegion) 

        CALL RFLU_DestroyGeometry(pRegion)  
        CALL RFLU_DestroyFaceList(pRegion)
        CALL RFLU_DestroyBVertexLists(pRegion)    
        CALL RFLU_DestroyCellMapping(pRegion)

        CALL RFLU_DestroyGrid(pRegion) 
      END DO ! iReg

! ------------------------------------------------------------------------------
!  Subbu 
!   1. The following section is set up temporarily to run rflupost (merged) as 
!      it requires .pdim_00000_0.00000E+00.
!   2. Can be discarded once the paraview O/P files are fully functional. 
!   3. To run rflupost (merged) requires cyldet.*_00000 files generated from
!      rflupart - Here we assume these files still exist until we have paraview.
! ------------------------------------------------------------------------------

      IF (global%nRegions > 1) THEN
        nPclsSumReg = 0
        DO iReg = 1,global%nRegions
          pRegion => levels(1)%regions(iReg)
          pPlag   => pRegion%plag
          nPclsSumReg = nPclsSumReg + pPlag%nPcls
        END DO

        pRegionSerial => levels(1)%regions(0)
        pPlagSerial => pRegionSerial%plag
        pPlagSerial%nPcls = nPclsSumReg
        pPlagSerial%nPclsMax = NINT(1.20_RFREAL*nPclsSumReg)
        pRegionSerial%plagInput%nCont = 1
        pPlagSerial%nextIdNumber = nPclsSumReg
        CALL PLAG_RFLU_WriteDimensions(pRegionSerial)
      END IF

! ------------------------------------------------------------------------------
 
      !CALL PLAG_DSTR_DestroyPclListCSR(pRegionSerial)
      !CALL PLAG_DSTR_DestroyCell2PclList(pRegionSerial)

    END IF ! global%nRegions

  END IF ! global%plagUsed
#endif  

! ******************************************************************************
! Write grid speed files. NOTE separated from above because need number of 
! faces, which is not always known above.
! ******************************************************************************

  DO iReg = iRegLow,iRegUpp
    pRegion => levels(1)%regions(iReg)

    IF ( RFLU_DecideNeedGridSpeeds(pRegion) .EQV. .TRUE. ) THEN 
      CALL RFLU_ReadDimensions(pRegion)                          
      CALL RFLU_CreateGrid(pRegion)

      IF ( pRegion%grid%nPatches > 0 ) THEN        
        CALL RFLU_ReadBCInputFileWrapper(pRegion)    
      END IF ! pRegion%grid%nPatches                   

#ifdef GENX
      CALL RFLU_GENX_ReadWindow(pRegion,GENX_WINDOW_TYPE_SURF)
      CALL RFLU_GENX_ReadWindow(pRegion,GENX_WINDOW_TYPE_VOL)
      CALL RFLU_GENX_GetDimensionsDerived(pRegion)
      CALL RFLU_GENX_CreateGridSurf(pRegion) 
      CALL RFLU_GENX_RegisterGridSurf(pRegion) 
      CALL RFLU_GENX_RegisterGridVol(pRegion)    
#endif

      CALL RFLU_ReadGridWrapper(pRegion)

#ifdef GENX
      CALL RFLU_GENX_DestroyGridSurf(pRegion) 
#endif

      CALL RFLU_CreateCellMapping(pRegion)
      CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
      CALL RFLU_BuildGlob2LocCellMapping(pRegion)         

      CALL RFLU_CreateBVertexLists(pRegion)
      CALL RFLU_BuildBVertexLists(pRegion)

      CALL RFLU_CreateFaceList(pRegion)
      CALL RFLU_BuildFaceList(pRegion)
      CALL RFLU_RenumberBFaceLists(pRegion)

      CALL RFLU_AllocateMemoryGSpeeds(pRegion)
#ifdef GENX
      CALL RFLU_GENX_RegisterDataGSpeeds(pRegion) 
#endif
      CALL RFLU_WriteGridSpeedsWrapper(pRegion)

      CALL RFLU_DeallocateMemoryGSpeeds(pRegion) 

      CALL RFLU_DestroyFaceList(pRegion)
      CALL RFLU_DestroyBVertexLists(pRegion)    
      CALL RFLU_DestroyCellMapping(pRegion) 
      CALL RFLU_DestroyGrid(pRegion)  

#ifdef GENX
      CALL COM_delete_window(TRIM(global%volWinNameInput))
      CALL COM_delete_window(TRIM(global%surfWinNameInput))
#endif
    END IF ! RFLU_DecideNeedGridSpeeds 
  END DO ! iReg  

! ******************************************************************************
! Print info about warnings
! ******************************************************************************
 
  CALL RFLU_PrintWarnInfo(global)
                              
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE rfluinit

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: rfluinit_22Oct2014.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.2  2014/07/21 16:44:23  subbu
! Added capability to initialize lagrangian particles from ASCII file
! via subroutine PLAG_RFLU_InitSolutionFile
!
! Revision 1.1.1.1  2014/05/05 21:47:47  tmish
! Initial checkin for rocflu macro.
!
! Revision 1.7  2010/03/15 00:31:10  mparmar
! Added support for auxiliary data in all init options
!
! Revision 1.6  2008/12/06 08:43:55  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:09  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2008/05/29 01:35:31  mparmar
! Added initialization capabiity to implement step change in particle velocity
!
! Revision 1.3  2007/11/28 23:05:42  mparmar
! Added initialization of aux vars and non-dimensionalization of cv
!
! Revision 1.2  2007/06/18 18:13:43  mparmar
! Added initialization of moving reference frame data
!
! Revision 1.1  2007/04/09 18:55:41  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.19  2007/03/27 01:31:45  haselbac
! Remove superfluous USE PLAG_ModParameters statement (bad check-in)
!
! Revision 1.18  2007/03/27 00:46:14  haselbac
! Adapted to changes in RFLU_SetDimensions call
!
! Revision 1.17  2007/03/27 00:23:23  haselbac
! PLAG init completely revamped to speed up 1d cases substantially
!
! Revision 1.16  2007/03/20 17:35:16  fnajjar
! Modified USE call to streamline with new module PLAG_ModDimensions
!
! Revision 1.15  2007/03/15 22:00:58  haselbac
! Adapted to changes in PLAG init for serial runs
!
! Revision 1.14  2006/08/19 15:41:13  mparmar
! Added calls to create, init, write, and destroy patch arrays
!
! Revision 1.13  2006/05/05 18:23:47  haselbac
! Changed PLAG init so do not need serial region anymore
!
! Revision 1.12  2006/02/06 23:55:55  haselbac
! Added comm argument to RFLU_InitGlobal
!
! Revision 1.11  2005/11/10 02:44:45  haselbac
! Added support for variable properties
!
! Revision 1.10  2005/09/23 19:00:45  haselbac
! Bug fix: When init PLAG, did not know about bc
!
! Revision 1.9  2005/09/13 21:37:30  haselbac
! Added new init option
!
! Revision 1.8  2005/05/18 22:23:59  fnajjar
! Added capability of init particles
!
! Revision 1.7  2005/05/05 18:38:39  haselbac
! Removed MPI calls after bug in Rocin/out fixed
!
! Revision 1.6  2005/05/04 03:37:50  haselbac
! Commented out COM_set_verbose call
!
! Revision 1.5  2005/05/04 03:35:58  haselbac
! Added init and finalize MPI when running within GENX
!
! Revision 1.4  2005/05/03 03:10:06  haselbac
! Converted to C++ reading of command-line
!
! Revision 1.3  2005/04/22 15:20:24  haselbac
! Fixed bug in combo init: grid and geom was missing; added grid info calls
!
! Revision 1.2  2005/04/18 20:33:27  haselbac
! Removed USE RFLU_ModCommLists
!
! Revision 1.1  2005/04/15 15:08:15  haselbac
! Initial revision
!
! ******************************************************************************

