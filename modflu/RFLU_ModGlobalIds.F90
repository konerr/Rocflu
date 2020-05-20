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
! Purpose: Collection of routines to obtain global ids of virtual cells. 
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModGlobalIds.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2006-2007 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModGlobalIds

  USE ModParameters
  USE ModDataTypes  
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModBorder, ONLY: t_border
  
  USE ModMPI

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModGlobalIds.F90,v $'

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: RFLU_GID_ComputenCellsOffset, &
            RFLU_GID_ComputeGlobalIds, &
            RFLU_GID_CreatenCellsOffset, &
            RFLU_GID_CreateGlobalIds, &
            RFLU_GID_DestroynCellsOffset, &
            RFLU_GID_DestroyGlobalIds

! ==============================================================================
! Private functions
! ==============================================================================
  
  PRIVATE :: RFLU_GID_InitnCellsOffset

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
 

 




! ******************************************************************************
!
! Purpose: Compute pRegion%nCellsOffset variable. 
!
! Description: None.
!
! Input:
!   pRegion            Pointer to region data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GID_ComputenCellsOffset(regions)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), DIMENSION(:), POINTER :: regions
 
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,iProc,iReg,iReg1,nProcs,nRegions,offset
    INTEGER, DIMENSION(:), ALLOCATABLE :: globalnCellsProc,globalnCellsOffset, &
                                          localnCellsProc,localnCellsOffset
    TYPE(t_region), POINTER :: pRegion
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => regions(1)%global
    
    CALL RegisterFunction(global,'RFLU_GID_ComputenCellsOffset',__FILE__)

    nProcs   = global%nProcAlloc
    nRegions = global%nRegions

! ******************************************************************************
!   Serial region case 
! ******************************************************************************

    IF ( nRegions == 1 ) THEN
      regions(1)%nCellsOffset(0) = 0
      regions(1)%nCellsOffset(1) = 0
    ELSE

! ******************************************************************************
!     Multiple regions case, Allocate temporary memory
! ******************************************************************************

      ALLOCATE(globalnCellsProc(nProcs),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'globalnCellsProc')
      END IF ! global%error 
 
      ALLOCATE(localnCellsProc(nProcs),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'localnCellsProc')
      END IF ! global%error 
 
      ALLOCATE(globalnCellsOffset(nRegions),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'globalnCellsOffset')
      END IF ! global%error 
 
      ALLOCATE(localnCellsOffset(nRegions),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'localnCellsOffset')
      END IF ! global%error 
 
! ==============================================================================
!     Initialize local and global nCellsProc arrays 
! ==============================================================================

      DO iProc = 1,nProcs
        globalnCellsProc(iProc) = 0

        localnCellsProc(iProc)  = 0
      END DO ! iPatch

! ==============================================================================
!     Set local nCellsProc 
! ==============================================================================

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        localnCellsProc(global%myProcid+1) = localnCellsProc(global%myProcid+1) &
                                           + pRegion%grid%nCells
      END DO ! iReg

! ==============================================================================
!     Compute global nCellsProc
! ==============================================================================

      CALL MPI_AllReduce(localnCellsProc,globalnCellsProc,nProcs, &
                         MPI_INTEGER,MPI_SUM,global%mpiComm,errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
      END IF ! global%error

! ==============================================================================
!     Initialize local and global nCellsOffset arrays
! ==============================================================================

      DO iReg1 = 1,nRegions
        globalnCellsOffset(iReg1) = 0

        localnCellsOffset(iReg1)  = 0
      END DO ! iPatch

! ==============================================================================
!     Set local nCellsOffset 
! ==============================================================================

      offset = 0
      IF ( global%myProcid /= 0 ) THEN
        DO iProc=1,global%myProcid
          offset = offset + globalnCellsProc(iProc)
        END DO ! iProc
      END IF ! global%myProcid

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)
  
        localnCellsOffset(pRegion%iRegionGlobal) = offset

        offset = offset + pRegion%grid%nCells
      END DO ! iReg

! ==============================================================================
!     Compute global nCellsOffset
! ==============================================================================

      CALL MPI_AllReduce(localnCellsOffset,globalnCellsOffset,nRegions, &
                         MPI_INTEGER,MPI_SUM,global%mpiComm,errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
      END IF ! global%error

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        DO iReg1 = 1,nRegions
          pRegion%nCellsOffset(iReg1) = globalnCellsOffset(iReg1)
        END DO ! iReg1
      END DO ! iReg

! ==============================================================================
!     Deallocate temporary memory
! ==============================================================================

      DEALLOCATE(globalnCellsProc,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'globalnCellsProc')
      END IF ! global%error 
 
      DEALLOCATE(localnCellsProc,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'localnCellsProc')
      END IF ! global%error 
 
      DEALLOCATE(globalnCellsOffset,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'globalnCellsOffset')
      END IF ! global%error 
 
      DEALLOCATE(localnCellsOffset,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'localnCellsOffset')
      END IF ! global%error 

    END IF ! nRegions

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_GID_ComputenCellsOffset









! ******************************************************************************
!
! Purpose: Compute global Ids for virtual cells. 
!
! Description: None.
!
! Input:
!   pRegion            Pointer to region data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GID_ComputeGlobalIds(regions)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), DIMENSION(:), POINTER :: regions
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,iBorder,iBorder2,icg,iReg,iReg2 
    TYPE(t_global), POINTER :: global  
    TYPE(t_border), POINTER :: pBorder,pBorder2
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_region), POINTER :: pRegion,pRegion2

! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => regions(1)%global
 
    CALL RegisterFunction(global,'RFLU_GID_ComputeGlobalIds',__FILE__)

! ******************************************************************************
!   Compute global ids of virtual cells 
!   Loop over borders
! ******************************************************************************

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid     

      DO iBorder = 1,pGrid%nBorders
        pBorder => pGrid%borders(iBorder)

! ==============================================================================
!       Copy data if on same process
! ==============================================================================

        IF ( pBorder%iProc == global%myProcid ) THEN
          iReg2    = pBorder%iRegionLocal
          iBorder2 = pBorder%iBorder

          pRegion2 => regions(iReg2)   
          pBorder2 => pRegion2%grid%borders(iBorder2)

! ------------------------------------------------------------------------------
!         Copy global id of virtual cells
! ------------------------------------------------------------------------------

          DO icg=1,pBorder2%nCellsSend
            pBorder%icgRecvGlobalIds(icg) = pBorder2%icgSend(icg) &
                                  + pRegion%nCellsOffset(pRegion2%iRegionGlobal)
          END DO ! icg

          DO icg=1,pBorder2%nCellsSend
            pGrid%virt2GlobalIds(pBorder%icgRecv(icg)-pGrid%nCells) &
                                          = pBorder%icgRecvGlobalIds(icg) 
          END DO ! icg
        ELSE
! ==============================================================================
!       Read data from .com file if on different process
! ==============================================================================

          iReg2    = pBorder%iRegionGlobal
          iBorder2 = pBorder%iBorder

! ------------------------------------------------------------------------------
!         Read global id of virtual cells
! ------------------------------------------------------------------------------

          CALL RFLU_GID_ReadGlobalIds(pRegion,iReg2,iBorder2, &
                                      pBorder%nCellsRecv, &
                                      pBorder%icgRecvGlobalIds)

          DO icg=1,pBorder%nCellsRecv
            pBorder%icgRecvGlobalIds(icg) = pBorder%icgRecvGlobalIds(icg) &
                                          + pRegion%nCellsOffset(iReg2)
          END DO ! icg

          DO icg=1,pBorder%nCellsRecv
            pGrid%virt2GlobalIds(pBorder%icgRecv(icg)-pGrid%nCells) &
                                          = pBorder%icgRecvGlobalIds(icg) 
          END DO ! icg

        END IF ! pBorder
      END DO ! iBorder
    END DO ! iReg
 
! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_GID_ComputeGlobalIds







! ******************************************************************************
!
! Purpose: Allocate memory for pRegion%nCellsOffset variable. 
!
! Description: None.
!
! Input:
!   pRegion            Pointer to region data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GID_CreatenCellsOffset(pRegion)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'RFLU_GID_CreatenCellsOffset',__FILE__)

! ==============================================================================
!   Allocate memory 
! ==============================================================================

    ALLOCATE(pRegion%nCellsOffset(0:global%nRegions),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%nCellsOffset')
    END IF ! global%error 

! ==============================================================================
!   Initialize memory
! ==============================================================================

    CALL RFLU_GID_InitnCellsOffset(pRegion)
 
! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_GID_CreatenCellsOffset







! ******************************************************************************
!
! Purpose: Create array to store global ids. 
!
! Description: None.
!
! Input:
!   pRegion            Pointer to region data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GID_CreateGlobalIds(pRegion)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,iBorder
    TYPE(t_global), POINTER :: global  
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'RFLU_GID_CreateGlobalIds',__FILE__)

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    IF ( pGrid%nCellsTot > pGrid%nCells ) THEN
      ALLOCATE(pGrid%virt2GlobalIds(pGrid%nCellsTot-pGrid%nCells), &
               STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%virt2GlobalIds')
      END IF ! global%error
    ELSE
      ALLOCATE(pGrid%virt2GlobalIds(1),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%virt2GlobalIds')
      END IF ! global%error
    END IF ! pGrid%nCellsTot

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      ALLOCATE(pBorder%icgRecvGlobalIds(pBorder%nCellsRecv),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%icgRecvGlobalIds')
      END IF ! global%error
    END DO ! iBorder
 
! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_GID_CreateGlobalIds







! ******************************************************************************
!
! Purpose: Destroy pRegion%nCellsOffset variable. 
!
! Description: None.
!
! Input:
!   pRegion            Pointer to region data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GID_DestroynCellsOffset(pRegion)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'RFLU_GID_DestroynCellsOffset',__FILE__)

! ==============================================================================
!   Deallocate memory
! ==============================================================================

    DEALLOCATE(pRegion%nCellsOffset,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%nCellsOffset')
    END IF ! global%error 
 
! ==============================================================================
!   Nullify memory
! ==============================================================================

    NULLIFY(pRegion%nCellsOffset)

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_GID_DestroynCellsOffset









! ******************************************************************************
!
! Purpose: Destroy array to store global ids. 
!
! Description: None.
!
! Input:
!   pRegion            Pointer to region data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GID_DestroyGlobalIds(pRegion)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,iBorder
    TYPE(t_global), POINTER :: global  
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'RFLU_GID_DestroyGlobalIds',__FILE__)

    pGrid => pRegion%grid

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      DEALLOCATE(pBorder%icgRecvGlobalIds,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__, &
                       'pBorder%icgRecvGlobalIds')
      END IF ! global%error
    END DO ! iBorder
 
! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_GID_DestroyGlobalIds







! ******************************************************************************
!
! Purpose: Initialize pRegion%nCellsOffset variable. 
!
! Description: None.
!
! Input:
!   pRegion            Pointer to region data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GID_InitnCellsOffset(pRegion)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: iReg
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'RFLU_GID_InitnCellsOffset',__FILE__)

! ==============================================================================
!   Initialize cell offsets
! ==============================================================================

    DO iReg=0,global%nRegions
      pRegion%nCellsOffset(iReg) = 0
    END DO ! iReg
 
! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_GID_InitnCellsOffset









! ******************************************************************************
!
! Purpose: Read global ids
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GID_ReadGlobalIds(pRegion,iRegionRead,iBorderRead,nRecv, &
                                    icgRecvGlobalIds)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iBorderRead,iRegionread,nRecv
    INTEGER, DIMENSION(:), INTENT(INOUT) :: icgRecvGlobalIds 
    TYPE(t_region), POINTER :: pRegion  
   
! ==============================================================================
!   Local variables
! ==============================================================================

    LOGICAL :: readRecvGlobalIds
    INTEGER :: errorFlag,dummy,iBorder,iBorderDummy,icl,iFile,iRegionGlobal, &
               ivl,loopCounter,nBorders,nCellsRecv,nCellsSend,nVertRecv, &
               nVertSend,nVertShared
    CHARACTER(CHRLEN) :: iFileName,sectionString
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_GID_ReadGlobalIds',__FILE__)

    readRecvGlobalIds = .FALSE.

! ******************************************************************************
!   Open file
! ******************************************************************************

    iFile = IF_COMM_LISTS

    CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.com', & 
                            iRegionRead,iFileName) 

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD", & 
         IOSTAT=errorFlag)   
    global%error = errorFlag        
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   Header and general information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%myProcid

    READ(iFile,'(A)') sectionString 
    IF ( TRIM(sectionString) /= '# ROCFLU communication lists file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM
    
! ******************************************************************************
!   Dimensions
! ******************************************************************************
  
    READ(iFile,'(A)') sectionString 
    IF ( TRIM(sectionString) /= '# Dimensions' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END IF ! TRIM
    
    READ(iFile,'(I16)') nBorders

    IF ( iBorderRead > nBorders ) THEN
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__,'Border not available')    
    END IF ! iBorderRead

! ******************************************************************************
!   Rest of file
! ******************************************************************************

    loopCounter = 0

! ==============================================================================  
!   Set up infinite loop
! ==============================================================================  

    DO ! set up infinite loop
      loopCounter = loopCounter + 1
    
      READ(iFile,'(A)') sectionString
   
      SELECT CASE ( TRIM(sectionString) ) 

! ------------------------------------------------------------------------------
!       Information
! ------------------------------------------------------------------------------

        CASE ( '# Information' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Information...'
          END IF ! global%myProcid 
              
          DO iBorder = 1,nBorders          
            READ(iFile,'(2(I16))') iRegionGlobal,iBorderDummy  
          END DO ! iBorder 
    
! ------------------------------------------------------------------------------
!       Cells
! ------------------------------------------------------------------------------

        CASE ( '# Cells' ) 
              
          DO iBorder = 1,nBorders          
    
            READ(iFile,'(2(I16))') nCellsSend,nCellsRecv

            IF ( iBorder == iBorderRead ) THEN
              IF ( nCellsSend /= nRecv ) THEN 
                CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
              END IF ! nCellsSend
              
              READ(iFile,'(10(I16))') (icgRecvGlobalIds(icl),icl=1,nRecv)
              readRecvGlobalIds = .TRUE.

              READ(iFile,'(10(I16))') (dummy,icl=1,nCellsRecv)
            ELSE

              READ(iFile,'(10(I16))') (dummy,icl=1,nCellsSend)
              READ(iFile,'(10(I16))') (dummy,icl=1,nCellsRecv)
            END IF ! iBorder 
          END DO ! iBorder   

! ------------------------------------------------------------------------------
!       Vertices
! ------------------------------------------------------------------------------

        CASE ( '# Vertices' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Vertices...'
          END IF ! global%myProcid 
              
          DO iBorder = 1,nBorders          
    
            READ(iFile,'(3(I16))') nVertSend,nVertRecv,nVertShared
                        
            READ(iFile,'(10(I16))') (dummy,icl=1,nVertSend)
            READ(iFile,'(10(I16))') (dummy,icl=1,nVertRecv)
            READ(iFile,'(10(I16))') (dummy,icl=1,nVertShared)
          END DO ! iBorder   
   
! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------ 
      
        CASE ( '# End' ) 

           IF ( readRecvGlobalIds .EQV. .FALSE. ) THEN
             CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,'could not '// &
                                                               'global ids')
           END IF ! readRecvGlobalIds

           EXIT
      
! ------------------------------------------------------------------------------
!       Invalid section string
! ------------------------------------------------------------------------------ 
      
        CASE DEFAULT
          CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
      END SELECT ! TRIM
  
! ==============================================================================
!     Guard against infinite loop - might be unnecessary because of read errors?
! ==============================================================================  

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter  
    END DO ! <empty>   

! ******************************************************************************
!   Close file
! ******************************************************************************

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_GID_ReadGlobalIds






! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLU_ModGlobalIds

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModGlobalIds.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:36  brollin
! New Stable version
!
! Revision 1.5  2009/07/08 19:11:55  mparmar
! Fixed bug for serial computation with one region
!
! Revision 1.4  2008/12/06 08:43:41  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:16:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/11/29 02:03:02  mparmar
! Removed ^M at the end of lines
!
! Revision 1.1  2007/11/28 23:04:49  mparmar
! Initial revision
!
!
! ******************************************************************************

