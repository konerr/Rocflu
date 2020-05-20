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
! Purpose: Driver routine for rfluclone. 
!
! Description: None.
!
! Input: 
!   caseString	String with casename
!   verbLevel	Verbosity level
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: rfluclone.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE rfluclone(caseString,verbLevel)

  USE ModError
  USE ModDataTypes
  USE ModParameters  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModMPI
  
  USE RFLU_ModAllocateMemory
  USE RFLU_ModCellMapping
  USE RFLU_ModCommLists
  USE RFLU_ModDeallocateMemory
  USE RFLU_ModDimensions
  USE RFLU_ModReadBcInputFile
  USE RFLU_ModReadWriteGrid
  USE RFLU_ModRegionMapping


  USE ModInterfaces, ONLY: RFLU_BuildDataStruct, &
                           RFLU_CheckClonability, & 
                           RFLU_CloneCommLists, & 
                           RFLU_CloneGrid, & 
                           RFLU_CreateGrid, &
                           RFLU_DestroyGrid, &
                           RFLU_GetUserInput, &
                           RFLU_InitGlobal, & 
                           RFLU_PrintGridInfo, & 
                           RFLU_PrintHeader, & 
                           RFLU_PrintWarnInfo, &                         
                           RFLU_ReadRestartInfo, & 
                           RFLU_SetRestartTimeFlag, & 
                           RFLU_WriteVersionString
                           
  IMPLICIT NONE

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
  INTEGER :: errorFlag,iLev,iReg
  TYPE(t_region), POINTER :: pRegion,pRegionSerial
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_level), POINTER :: levels(:)
  
! ******************************************************************************
! Start, initialize global data
! ******************************************************************************  

  ALLOCATE(global,STAT=errorFlag)
  IF ( errorFlag /= ERR_NONE ) THEN 
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR - Pointer allocation failed.'
    STOP
  END IF ! errorFlag 

  casename = caseString(1:LEN(caseString))

  CALL RFLU_InitGlobal(casename,verbLevel,MPI_COMM_WORLD,global)

  CALL RegisterFunction(global,'rflupart', & 
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
! Read mapping file, impose serial mapping, and build basic data structure
! ******************************************************************************

  CALL RFLU_ReadRegionMappingFile(global,MAPFILE_READMODE_PEEK,global%myProcId)
  CALL RFLU_SetRegionMappingSerial(global)  
  CALL RFLU_CreateRegionMapping(global,MAPTYPE_REG)
  CALL RFLU_ImposeRegionMappingSerial(global)

  CALL RFLU_BuildDataStruct(global,levels) 
  CALL RFLU_ApplyRegionMapping(global,levels)
  CALL RFLU_DestroyRegionMapping(global,MAPTYPE_REG)  

! ******************************************************************************
! Read input file and restart info. NOTE need restart info for GENX runs to 
! determine whether have a restart. 
! ******************************************************************************

  CALL RFLU_GetUserInput(levels(1)%regions,.TRUE.) 
  CALL RFLU_ReadRestartInfo(global)
  CALL RFLU_SetRestartTimeFlag(global)

! ******************************************************************************
! Read dimensions, bc file, and grid
! ******************************************************************************

  pRegionSerial => levels(1)%regions(0) ! NOTE must set otherwise get core dump   

  CALL RFLU_ReadDimensionsWrapper(pRegionSerial)    
  CALL RFLU_CreateGrid(pRegionSerial)
  CALL RFLU_ReadBCInputFileWrapper(pRegionSerial)     
  CALL RFLU_ReadGridWrapper(pRegionSerial)
  
  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    CALL RFLU_PrintGridInfo(pRegionSerial)
  END IF ! global%verbLevel  

  CALL RFLU_CreateCellMapping(pRegionSerial)
  CALL RFLU_ReadLoc2GlobCellMapping(pRegionSerial)

! ******************************************************************************
! Create and read comm lists, set proc ids for borders
! ******************************************************************************

  CALL RFLU_COMM_CreateBorders(pRegionSerial,CREATE_BORDERS_MODE_DIM_KNOWN)  
  CALL RFLU_COMM_CreateCommLists(pRegionSerial)
  CALL RFLU_COMM_ReadCommLists(pRegionSerial)

! ******************************************************************************
! Check clonability
! ******************************************************************************

  CALL RFLU_CheckClonability(pRegionSerial)

! ******************************************************************************
! Clone grid, mapping, comm lists, and write out 
! ******************************************************************************
  
  DO iReg = 1,global%nRegionsLocal
    pRegion => levels(1)%regions(iReg)
    
    CALL RFLU_CloneGrid(pRegionSerial,pRegion)
    CALL RFLU_WriteGridWrapper(pRegion)
    CALL RFLU_WriteLoc2GlobCellMapping(pRegion)
    
    CALL RFLU_CloneCommLists(pRegionSerial,pRegion)
    CALL RFLU_WriteDimensionsWrapper(pRegion,WRITE_DIMENS_MODE_FORCE)    
    CALL RFLU_COMM_WriteCommLists(pRegion) 
    
    CALL RFLU_COMM_DestroyCommLists(pRegion)
    CALL RFLU_COMM_DestroyBorders(pRegion)
    CALL RFLU_DestroyCellMapping(pRegion)
    CALL RFLU_DestroyGrid(pRegion)   
  END DO ! iReg

! ******************************************************************************
! Deallocate memory 
! ******************************************************************************

  CALL RFLU_COMM_DestroyCommLists(pRegionSerial)
  CALL RFLU_COMM_DestroyBorders(pRegionSerial)  
  CALL RFLU_DestroyCellMapping(pRegionSerial)
  CALL RFLU_DestroyGrid(pRegionSerial)

! ******************************************************************************
! Print info about warnings
! ******************************************************************************
  
  CALL RFLU_PrintWarnInfo(global)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE rfluclone

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: rfluclone.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:54  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:06  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/08/07 17:13:59  haselbac
! Initial revision
!
! ******************************************************************************

