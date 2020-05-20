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
! Purpose: Driver routine for rfluextr.
!
! Description: None.
!
! Input: 
!   caseString	String with casename
!   stampString	String with iteration or time stamp
!   verbLevel	Verbosity level
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: rfluextr.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE rfluextr(caseString,stampString,verbLevel)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global 
  USE ModError
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModGrid, ONLY: t_grid
  USE ModParameters
  USE ModMPI
  
  USE RFLU_ModBoundLists
  USE RFLU_ModCellMapping
  USE RFLU_ModDimensions
  USE RFLU_ModFaceList
  USE RFLU_ModGeometry
  USE RFLU_ModGeometryTools
  USE RFLU_ModReadBcInputFile
  USE RFLU_ModReadWriteGrid
  USE RFLU_ModRegionMapping
  USE RFLU_ModRenumberings
  USE RFLU_ModInCellTest
  
  USE ModInterfaces, ONLY: RFLU_AllocateMemoryXSect, &
                           RFLU_BuildDataStruct, &
                           RFLU_CreateGrid, &
                           RFLU_DeallocateMemoryXSect, &
                           RFLU_DestroyGrid, &
                           RFLU_ExtractLineDataQuad2D, & 
                           RFLU_GetUserInput, &
                           RFLU_InitGlobal, &
                           RFLU_PrintGridInfo, &
                           RFLU_PrintHeader, &
                           RFLU_PrintWarnInfo, &
                           RFLU_ReadInputFile, & 
                           RFLU_WriteLineData, & 
                           RFLU_WritePostInfo, & 
                           RFLU_WriteVersionString
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(*) :: caseString,stampString
  INTEGER, INTENT(IN) :: verbLevel
  
! ==============================================================================
! Local variables
! ==============================================================================

  LOGICAL :: fileExists,iBBFlag
  CHARACTER(CHRLEN) :: casename,RCSIdentString,stamp
  INTEGER :: errorFlag,icg,icgStart,iReg,iRegStart
  REAL(RFREAL) :: xBeg,xMax,xMin,yBeg,yMax,yMin,zBeg,zMax,zMin
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_level), DIMENSION(:), POINTER :: levels  
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion,pRegionSerial

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: rfluextr.F90,v $ $Revision: 1.1.1.1 $'

! ******************************************************************************
! Initialize global data
! ******************************************************************************
  
  ALLOCATE(global,STAT=errorFlag)
  IF ( errorFlag /= ERR_NONE ) THEN 
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR - Pointer allocation failed.'
    STOP
  END IF ! errorFlag   

  casename = caseString(1:LEN(caseString))
  stamp    = stampString(1:LEN(stampString))
  
  CALL RFLU_InitGlobal(casename,verbLevel,CRAZY_VALUE_INT,global)

  CALL RegisterFunction(global,'rfluextr', &
                        __FILE__)

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

  CALL RFLU_GetUserInput(levels(1)%regions)  

  IF ( global%flowType == FLOW_STEADY ) THEN
    READ(stamp,*) global%currentIter
  ELSE
    READ(stamp,*) global%currentTime
  END IF ! global%flowType  

! ****************************************************************************** 
! Get initial and final positions
! ******************************************************************************

  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Enter starting position (x,y,z):'
  READ(STDIN,*) global%extrXCoordBeg,global%extrYCoordBeg,global%extrZCoordBeg

  xBeg = global%extrXCoordBeg 
  yBeg = global%extrYCoordBeg 
  zBeg = global%extrZCoordBeg 

  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Enter ending position (x,y,z):'
  READ(STDIN,*) global%extrXCoordEnd,global%extrYCoordEnd,global%extrZCoordEnd

! ****************************************************************************** 
! Allocate memory and read data
! ******************************************************************************

  pRegionSerial => levels(1)%regions(0) 

! TEMPORARY
  pRegionSerial%grid%nXSectMax = 2000
! END TEMPORARY

  CALL RFLU_AllocateMemoryXSect(pRegionSerial)

  CALL RFLU_ReadDimensions(pRegionSerial)

  IF ( pRegionSerial%grid%nCells /= pRegionSerial%grid%nHexs ) THEN 
    CALL ErrorStop(global,ERR_CELL_TYPE_INVALID,__LINE__)     
  END IF ! pRegionSerial%grid%nCells

  IF ( global%nRegions > 1 ) THEN 
    CALL RFLU_RNMB_CreateSC2RMap(pRegionSerial)
    CALL RFLU_RNMB_ReadSC2RMap(pRegionSerial)
  END IF ! global%nRegions

! ****************************************************************************** 
! Allow for multiple tries of finding initial region bcos by looking at bounding
! box, initial region is ambiguous: Starting position may lie in multiple 
! bounding boxes. At present, ambiguity can only be resolved by checking cells.
! ******************************************************************************
 
  iReg = 0
  
  outerLoop: DO 
    iRegStart = CRAZY_VALUE_INT

! ==============================================================================
!   Get initial region 
! ==============================================================================
 
    regionLoop: DO 
      iReg = iReg + 1

      pRegion => levels(1)%regions(iReg)
      pGrid   => pRegion%grid
   
      CALL RFLU_ReadDimensions(pRegion)
      CALL RFLU_CreateGrid(pRegion)

      IF ( pRegion%grid%nPatches > 0 ) THEN
        CALL RFLU_ReadBCInputFileWrapper(pRegion)
      END IF ! pRegion%grid%nPatches

      CALL RFLU_ReadGridWrapper(pRegion)

      xMin = MINVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
      xMax = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
      yMin = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
      yMax = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
      zMin = MINVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))
      zMax = MAXVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))

      IF ( RFLU_TestInBoundBox(global,xBeg,yBeg,zBeg,xMin,xMax, & 
                               yMin,yMax,zMin,zMax) .EQV. .TRUE. ) THEN 
        iRegStart = iReg  
      END IF ! RFLU_TestInBoundBox

      CALL RFLU_DestroyGrid(pRegion)

      IF ( iRegStart /= CRAZY_VALUE_INT ) THEN
        EXIT regionLoop
      END IF ! iRegStart
    END DO regionLoop

    IF ( iRegStart == CRAZY_VALUE_INT ) THEN 
! TEMPORARY
      WRITE(*,*) 'ERROR! Did not find iRegStart!'

      EXIT outerLoop
! END TEMPORARY
    ELSE 
      WRITE(STDOUT,'(A,1X,A,1X,I5.5)') SOLVER_NAME, & 
        'Region bounding box containing starting position:',iRegStart
    END IF ! iRegStart

! ==============================================================================
!   Get initial cell
! ==============================================================================
 
    pRegion => levels(1)%regions(iRegStart)
    pGrid   => pRegion%grid

    CALL RFLU_ReadDimensions(pRegion)
    CALL RFLU_CreateGrid(pRegion)

    IF ( pRegion%grid%nPatches > 0 ) THEN
      CALL RFLU_ReadBCInputFileWrapper(pRegion)
    END IF ! pRegion%grid%nPatches

    CALL RFLU_ReadGridWrapper(pRegion)

    CALL RFLU_CreateCellMapping(pRegion)
    CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
    CALL RFLU_BuildGlob2LocCellMapping(pRegion)
    
    CALL RFLU_CreateBVertexLists(pRegion)
    CALL RFLU_BuildBVertexLists(pRegion)
    
    CALL RFLU_CreateFaceList(pRegion)
    CALL RFLU_BuildFaceList(pRegion)
    CALL RFLU_RenumberBFaceLists(pRegion)

    CALL RFLU_CreateCell2FaceList(pRegion)
    CALL RFLU_BuildCell2FaceList(pRegion)
  
    CALL RFLU_CreateGeometry(pRegion)
    CALL RFLU_BuildGeometry(pRegion)

    icgStart = CRAZY_VALUE_INT

    cellLoop: DO icg = 1,pGrid%nCells
      IF ( RFLU_ICT_TestInCell(pRegion,xBeg,yBeg,zBeg,icg) .EQV. .TRUE. ) THEN 
        icgStart = icg
  
        EXIT cellLoop
      END IF ! RFLU_ICT_TestInCell
    END DO cellLoop
  
    IF ( icgStart == CRAZY_VALUE_INT ) THEN 
! TEMPORARY
      WRITE(*,*) 'ERROR! Did not find icgStart!'
! END TEMPORARY
    ELSE 
      WRITE(STDOUT,'(A,1X,A,1X,I5)') SOLVER_NAME, & 
        'Cell containing starting position:',icgStart
  
      EXIT outerLoop
    END IF ! icgStart
  
    CALL RFLU_DestroyGeometry(pRegion)
    CALL RFLU_DestroyCell2FaceList(pRegion)
    CALL RFLU_DestroyFaceList(pRegion)
    CALL RFLU_DestroyBVertexLists(pRegion)
    CALL RFLU_DestroyCellMapping(pRegion)
    CALL RFLU_DestroyGrid(pRegion)

    IF ( iReg == global%nRegions ) THEN ! NOTE here icgStart always CRAZY
      exit outerLoop 
    END IF ! iReg
  END DO outerLoop

! ****************************************************************************** 
! Extract data and deallocate memory. NOTE no extra destruction necessary as it 
! will be carried out in RFLU_ExtractLineDataQuad2D.
! ******************************************************************************

  IF ( (iRegStart /= CRAZY_VALUE_INT) .AND. (icgStart /= CRAZY_VALUE_INT) ) THEN 
    CALL RFLU_ExtractLineDataQuad2D(levels(1)%regions,iRegStart,icgStart)
    CALL RFLU_WriteLineData(levels(1)%regions)
  END IF ! iRegStart

! ****************************************************************************** 
! Deallocate memory
! ******************************************************************************

  IF ( global%nRegions > 1 ) THEN 
    CALL RFLU_RNMB_DestroySC2RMap(pRegionSerial)
  END IF ! globabl%nRegions

  CALL RFLU_DeallocateMemoryXSect(pRegionSerial)
  
! *****************************************************************************
! Print info about warnings
! *****************************************************************************

  CALL RFLU_PrintWarnInfo(global)
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE rfluextr

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: rfluextr.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.5  2008/12/06 08:43:55  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:07  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2008/01/08 19:12:52  haselbac
! Improved logic to avoid bbox problem for non-square regions
!
! Revision 1.2  2007/12/05 13:26:49  haselbac
! Some cosmetic changes
!
! Revision 1.1  2007/11/27 13:17:27  haselbac
! Initial revision
!
! ******************************************************************************

