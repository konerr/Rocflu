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
! Purpose: Suite of routines for MPI interaction.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModMPI.F90,v 1.2 2015/07/23 23:11:18 brollin Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModMPI

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt_input
  USE ModBorder, ONLY: t_border,t_border_data
  USE ModMPI

  IMPLICIT NONE

  PUBLIC :: RFLU_MPI_ClearRequestWrapper, &
            RFLU_MPI_CopyWrapper, & 
            RFLU_MPI_CreateBufferIPclSend, &        
            RFLU_MPI_CreateBuffersWrapper, &    
            RFLU_MPI_DestroyBufferIPclSend, &           
            RFLU_MPI_DestroyBuffersWrapper, &        
            RFLU_MPI_HM_CopyWrapper, & 
            RFLU_MPI_HM_ISendWrapper, & 
            RFLU_MPI_HM_RecvWrapper, & 
            RFLU_MPI_ISendWrapper, & 
            RFLU_MPI_RecreateBufferIPclSend, &
            RFLU_MPI_RecvWrapper, & 
            RFLU_MPI_PLAG_CopyWrapper, & 
            RFLU_MPI_PLAG_ISendWrapper, & 
            RFLU_MPI_PLAG_RecvWrapper, & 
            RFLU_MPI_SetTagsWrapper

  PRIVATE :: RFLU_MPI_ClearRequest, &
             RFLU_MPI_CopyCellData, &
             RFLU_MPI_CopyCellData1d, &
             RFLU_MPI_CopyCellData3d, &
             RFLU_MPI_CreateBuffers, &
             RFLU_MPI_DestroyBuffers, &
             RFLU_MPI_HM_CopyCellData, &
             RFLU_MPI_HM_ISendCellData, &
             RFLU_MPI_HM_RecvCellData, &
             RFLU_MPI_ISendCellData, &
             RFLU_MPI_ISendCellData1d, &
             RFLU_MPI_ISendCellData3d, &
             RFLU_MPI_PLAG_CopyCellData, &
             RFLU_MPI_PLAG_ISendCellData, &
             RFLU_MPI_PLAG_RecvCellData, &
             RFLU_MPI_RecvCellData, &
             RFLU_MPI_RecvCellData1d, &
             RFLU_MPI_RecvCellData3d, &
             RFLU_MPI_SetTag
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModMPI.F90,v $ $Revision: 1.2 $' 
                      
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS






! ******************************************************************************
!
! Purpose: Clearing send requests.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   request     Request (to be cleared)
!
! Output: 
!   request     Request (cleared)
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_ClearRequest(global,request)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(INOUT) :: request
    TYPE(t_global), POINTER :: global  
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag
    INTEGER :: status(MPI_STATUS_SIZE)          
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_ClearRequest',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    CALL MPI_Wait(request,status,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_ClearRequest









! ******************************************************************************
!
! Purpose: Wrapper for clearing send requests.
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

  SUBROUTINE RFLU_MPI_ClearRequestWrapper(pRegion)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iBorder
    TYPE(t_border), POINTER :: pBorder 
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid             
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_ClearRequestWrapper',__FILE__)

#ifdef ROCPROF 
    CALL FPROFILER_BEGINS("RFLU::ClearRequestWrapper")
#endif

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid
    
! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)
    
! ==============================================================================
!     Clear request if not on same process
! ==============================================================================    
    
      IF ( pBorder%iProc /= global%myProcid ) THEN 
        IF ( pBorder%nCellsSend > 0 ) THEN 

! ------------------------------------------------------------------------------
!         Mixture
! ------------------------------------------------------------------------------      

          CALL RFLU_MPI_ClearRequest(global,pBorder%mixt%sendRequest)

! ------------------------------------------------------------------------------
!         Physical modules
! ------------------------------------------------------------------------------      

#ifdef SPEC
          IF ( global%specUsed .EQV. .TRUE. ) THEN 
            CALL RFLU_MPI_ClearRequest(global,pBorder%spec%sendRequest)   
          END IF ! global%specUsed
#endif  
        END IF ! pBorder%nCellsSend
      END IF ! pBorder
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF 
    CALL FPROFILER_ENDS("RFLU::ClearRequestWrapper")
#endif

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_ClearRequestWrapper






! ******************************************************************************
!
! Purpose: Copy cell data.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pBorder     Pointer to border
!   cellData    Array with cell data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_CopyCellData(global,pBorder,pBorder2,cellData,cellData2)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cellData
    REAL(RFREAL), DIMENSION(:,:), INTENT(OUT) :: cellData2    
    TYPE(t_border), POINTER :: pBorder,pBorder2
    TYPE(t_global), POINTER :: global         
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icg2,icl,iVar,nVars        
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_CopyCellData',__FILE__)

! ******************************************************************************
!   Set variables
! ******************************************************************************

    nVars = SIZE(cellData,1)
    
    IF ( nVars /= SIZE(cellData2,1) ) THEN 
      CALL ErrorStop(global,ERR_DATADIM_MISMATCH,__LINE__)
    END IF ! nVars 
    
! ******************************************************************************
!   Copy data 
! ******************************************************************************

    DO icl = 1,pBorder%nCellsSend
      icg  = pBorder%icgSend(icl)
      icg2 = pBorder2%icgRecv(icl)

      DO iVar = 1,nVars
        cellData2(iVar,icg2) = cellData(iVar,icg) 
      END DO ! iVar      
    END DO ! icl
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_CopyCellData








! ******************************************************************************
!
! Purpose: Copy single dimensional  cell data.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pBorder     Pointer to border
!   cellData    Array with cell data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_CopyCellData1d(global,pBorder,pBorder2,cellData,cellData2)

    IMPLICIT NONE
 
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
 
! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), DIMENSION(:), INTENT(IN) :: cellData
    REAL(RFREAL), DIMENSION(:), INTENT(OUT) :: cellData2
    TYPE(t_border), POINTER :: pBorder,pBorder2
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icg2,icl,iVar
   
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_CopyCellData1d',__FILE__)

! ******************************************************************************
!   Copy data
! ******************************************************************************

    DO icl = 1,pBorder%nCellsSend
      icg  = pBorder%icgSend(icl)
      icg2 = pBorder2%icgRecv(icl)

      cellData2(icg2) = cellData(icg)
    END DO ! icl

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MPI_CopyCellData1d








! ******************************************************************************
!
! Purpose: Copy cell data.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pBorder     Pointer to border
!   cellData    Array with cell data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_CopyCellData3d(global,pBorder,pBorder2,cellData,cellData2)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), DIMENSION(:,:,:) :: cellData
    REAL(RFREAL), DIMENSION(:,:,:) :: cellData2    
    TYPE(t_border), POINTER :: pBorder,pBorder2
    TYPE(t_global), POINTER :: global         
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icg2,icl,iVar
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_CopyCellData3d',__FILE__)

! ******************************************************************************
!   Copy data 
! ******************************************************************************

    DO icl = 1,pBorder%nCellsSend
      icg  = pBorder%icgSend(icl)
      icg2 = pBorder2%icgRecv(icl)

! TEMPORARY: manoj: Need to find a permanent solution
!      DO iVar = CV_MIXT_XVEL,CV_MIXT_ZVEL
      DO iVar = 1,3 
        cellData2(XCOORD,iVar,icg2) = cellData(XCOORD,iVar,icg) 
        cellData2(YCOORD,iVar,icg2) = cellData(YCOORD,iVar,icg) 
        cellData2(ZCOORD,iVar,icg2) = cellData(ZCOORD,iVar,icg) 
      END DO ! iVar      
    END DO ! icl
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_CopyCellData3d








! ******************************************************************************
!
! Purpose: Wrapper for copying data.
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

  SUBROUTINE RFLU_MPI_CopyWrapper(regions)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), DIMENSION(:), POINTER :: regions
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iBorder,iBorder2,iReg,iReg2
    TYPE(t_border), POINTER :: pBorder,pBorder2 
    TYPE(t_global), POINTER :: global 
    TYPE(t_grid), POINTER :: pGrid     
    TYPE(t_region), POINTER :: pRegion,pRegion2
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => regions(0)%global

    CALL RegisterFunction(global,'RFLU_MPI_CopyWrapper',__FILE__)

#ifdef ROCPROF 
    CALL FPROFILER_BEGINS("RFLU::CopyWrapper")
#endif

! ******************************************************************************
!   Set pointers
! ******************************************************************************
    
! ******************************************************************************
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
!         Check dimensions
! ------------------------------------------------------------------------------

          IF ( pBorder%nCellsSend /= pBorder2%nCellsRecv ) THEN 
            CALL ErrorStop(global,ERR_BUFFERDIM_MISMATCH,__LINE__)
          END IF ! pBorder
          
! ------------------------------------------------------------------------------
!         Mixture
! ------------------------------------------------------------------------------      
      
          CALL RFLU_MPI_CopyCellData(global,pBorder,pBorder2, & 
                                     pRegion%mixt%cv,pRegion2%mixt%cv)      
                                     
! ------------------------------------------------------------------------------
!         Physical modules
! ------------------------------------------------------------------------------      

#ifdef SPEC
          IF ( global%specUsed .EQV. .TRUE. ) THEN 
            CALL RFLU_MPI_CopyCellData(global,pBorder,pBorder2, & 
                                       pRegion%spec%cv,pRegion2%spec%cv)     
          END IF ! global%specUsed
#endif           
                          
        END IF ! pBorder
      END DO ! iBorder
    END DO ! iReg

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF 
    CALL FPROFILER_ENDS("RFLU::CopyWrapper")
#endif

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_CopyWrapper








! ******************************************************************************
!
! Purpose: Create iPclSend buffers.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pBorder     Pointer to border (Optional)
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_CreateBufferIPclSend(pRegion,pBorder)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_border), POINTER, OPTIONAL :: pBorder
    TYPE(t_region), POINTER :: pRegion 
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iBorder,nVars
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid     
    TYPE(t_border), POINTER :: pBorder2
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_CreateBufferIPclSend',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Loop over borders when pBorder is not present
!    else allocate only for select pBorder
! ******************************************************************************

#ifdef PLAG
    IF ( PRESENT(pBorder) .EQV. .TRUE. ) THEN
      nVars = SIZE(pBorder%iPclSend,1)

      ALLOCATE(pBorder%iPclSend(nVars,pBorder%nPclsSendMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%iPclSend')
      END IF ! global%error

    ELSE
      DO iBorder = 1,pGrid%nBorders
        pBorder2 => pGrid%borders(iBorder)

        nVars = 2
        pBorder2%nPclsSendMax = 1000

        ALLOCATE(pBorder2%iPclSend(nVars,pBorder2%nPclsSendMax),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%iPclSend')
        END IF ! global%error
      END DO ! iBorder
    END IF ! PRESENT(pBorder)
#endif

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
      
  END SUBROUTINE RFLU_MPI_CreateBufferIPclSend








! ******************************************************************************
!
! Purpose: Create buffers.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pBorder     Pointer to border
!   borderData  Border data
!   nVars       Number of variables in borderData
!   flag        Optional flag, if present 1d buffer is created
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_CreateBuffers(global,pBorder,borderData,nVars,flag)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: nVars
    INTEGER, INTENT(IN), OPTIONAL :: flag
    TYPE(t_global), POINTER :: global
    TYPE(t_border), POINTER :: pBorder 
    TYPE(t_border_data) :: borderData
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag       
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_CreateBuffers',__FILE__)

! ******************************************************************************
!   Allocate memory
! ******************************************************************************
       
    ALLOCATE(borderData%sendBuff(nVars,pBorder%nCellsSend),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'borderData%sendBuff')
    END IF ! global%error

    ALLOCATE(borderData%recvBuff(nVars,pBorder%nCellsRecv),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'borderData%recvBuff')
    END IF ! global%error   

    IF ( PRESENT(flag) .EQV. .TRUE.) THEN
!     Allocating buffer for delP communication
      IF ( flag > 0 ) THEN
        ALLOCATE(borderData%sendBuff1d(pBorder%nCellsSend),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'borderData%sendBuff1d')
        END IF ! global%error

        ALLOCATE(borderData%recvBuff1d(pBorder%nCellsRecv),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'borderData%recvBuff1d')
        END IF ! global%error
      END IF ! flag

!     Allocating buffer for velocity gradient communication
      IF ( flag > 1 ) THEN
        ALLOCATE(borderData%sendBuff3d(3,3,pBorder%nCellsSend),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'borderData%sendBuff3d')
        END IF ! global%error

        ALLOCATE(borderData%recvBuff3d(3,3,pBorder%nCellsRecv),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'borderData%recvBuff3d')
        END IF ! global%error
      END IF ! flag
    END IF ! PRESENT(flag)

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
      
  END SUBROUTINE RFLU_MPI_CreateBuffers







! ******************************************************************************
!
! Purpose: Wrapper for creating buffers.
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

  SUBROUTINE RFLU_MPI_CreateBuffersWrapper(pRegion)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,flag,iBorder,nVars
    TYPE(t_border), POINTER :: pBorder        
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid 
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_CreateBuffersWrapper',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating buffers...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

! ==============================================================================
!     Create buffers if not on same process
! ==============================================================================    
      
      IF ( pBorder%iProc /= global%myProcid ) THEN
      
! ------------------------------------------------------------------------------
!       Mixture
! ------------------------------------------------------------------------------      

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          nVars = 1

          IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN
            flag = 2
          ELSE
            flag = 1
          END IF ! pRegion%mixtInput%flowModel

          CALL RFLU_MPI_CreateBuffers(global,pBorder,pBorder%mixt,nVars,flag)
        ELSE
          nVars = pRegion%mixtInput%nCv
#ifdef PLAG
! ------- Create 1d buffer for communicating particle volume fraction ----------
          IF ( global%plagUsed .EQV. .TRUE. ) THEN 
            flag = 1
            CALL RFLU_MPI_CreateBuffers(global,pBorder,pBorder%mixt,nVars,flag)
          ELSE
            CALL RFLU_MPI_CreateBuffers(global,pBorder,pBorder%mixt,nVars)
          END IF ! global%plagUsed        
#else
          CALL RFLU_MPI_CreateBuffers(global,pBorder,pBorder%mixt,nVars)
#endif                           
        END IF ! global%solverType

! ------------------------------------------------------------------------------
!       Physical modules
! ------------------------------------------------------------------------------      

#ifdef SPEC
        IF ( global%specUsed .EQV. .TRUE. ) THEN 
          CALL RFLU_MPI_CreateBuffers(global,pBorder,pBorder%spec, & 
                                      pRegion%specInput%nSpecies)   

        END IF ! global%specUsed        
#endif                           
      END IF ! pBorder%iProc
    END DO ! iBorder
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating buffers done.'
    END IF ! global%verbLevel
      
  END SUBROUTINE RFLU_MPI_CreateBuffersWrapper








! ******************************************************************************
!
! Purpose: Destroy iPclSend buffers.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pBorder     Pointer to border (optional)
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_DestroyBufferIPclSend(pRegion,pBorder)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_border), POINTER, OPTIONAL :: pBorder
    TYPE(t_region), POINTER :: pRegion 
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iBorder      
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid     
    TYPE(t_border), POINTER :: pBorder2
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_DestroyBufferIPclSend',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Loop over borders when pBorder is not present
!    else deallocate only for select pBorder
! ******************************************************************************

#ifdef PLAG
    IF ( PRESENT(pBorder) .EQV. .TRUE. ) THEN
      DEALLOCATE(pBorder%iPclSend,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pBorder%iPclSend')
      END IF ! global%error

    ELSE
      DO iBorder = 1,pGrid%nBorders
        pBorder2 => pGrid%borders(iBorder)

        DEALLOCATE(pBorder2%iPclSend,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pBorder%iPclSend')
        END IF ! global%error
      END DO ! iBorder
    END IF ! PRESENT(pBorder)
#endif

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
      
  END SUBROUTINE RFLU_MPI_DestroyBufferIPclSend








! ******************************************************************************
!
! Purpose: Destroy buffers.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pBorder     Pointer to border
!   borderData  Border data
!   flag        Optional flag, if equal to 1 then 1d buffers are destroyed
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_DestroyBuffers(global,pBorder,borderData,flag)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN), OPTIONAL :: flag
    TYPE(t_global), POINTER :: global
    TYPE(t_border), POINTER :: pBorder 
    TYPE(t_border_data) :: borderData
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag       
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_DestroyBuffers',__FILE__)

! ******************************************************************************
!   Allocate memory
! ******************************************************************************
       
    DEALLOCATE(borderData%sendBuff,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'borderData%sendBuff')
    END IF ! global%error

    DEALLOCATE(borderData%recvBuff,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'borderData%recvBuff')
    END IF ! global%error   

    IF ( PRESENT(flag) .EQV. .TRUE. ) THEN
      IF ( flag > 0 ) THEN
        DEALLOCATE(borderData%sendBuff1d,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'borderData%sendBuff1d')
        END IF ! global%error

        DEALLOCATE(borderData%recvBuff1d,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'borderData%recvBuff1d')
        END IF ! global%error
      END IF ! flag  

      IF ( flag > 1 ) THEN
        DEALLOCATE(borderData%sendBuff3d,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'borderData%sendBuff3d')
        END IF ! global%error

        DEALLOCATE(borderData%recvBuff3d,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'borderData%recvBuff3d')
        END IF ! global%error
      END IF ! flag  
    END IF ! PRESENT(flag)

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
      
  END SUBROUTINE RFLU_MPI_DestroyBuffers







! ******************************************************************************
!
! Purpose: Wrapper for destroying buffers.
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

  SUBROUTINE RFLU_MPI_DestroyBuffersWrapper(pRegion)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,flag,iBorder   
    TYPE(t_border), POINTER :: pBorder        
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid 
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_DestroyBuffersWrapper',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying buffers...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

! ==============================================================================
!     Destroy buffers if not on same process
! ==============================================================================    
      
      IF ( pBorder%iProc /= global%myProcid ) THEN
      
! ------------------------------------------------------------------------------
!       Mixture
! ------------------------------------------------------------------------------      

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN
            flag = 2
          ELSE
            flag = 1
          END IF ! pRegion%mixtInput%flowModel

          CALL RFLU_MPI_DestroyBuffers(global,pBorder,pBorder%mixt,flag)
        ELSE
#ifdef PLAG
! ------- Destroy 1d buffer for communicating particle volume fraction ---------
          IF ( global%plagUsed .EQV. .TRUE. ) THEN 
            flag = 1
            CALL RFLU_MPI_DestroyBuffers(global,pBorder,pBorder%mixt,flag)
          ELSE
            CALL RFLU_MPI_DestroyBuffers(global,pBorder,pBorder%mixt)
          END IF ! global%plagUsed        
#else
          CALL RFLU_MPI_DestroyBuffers(global,pBorder,pBorder%mixt)
#endif                           
        END IF ! global%solverType

! ------------------------------------------------------------------------------
!       Physical modules
! ------------------------------------------------------------------------------      

#ifdef SPEC
        IF ( global%specUsed .EQV. .TRUE. ) THEN 
          CALL RFLU_MPI_DestroyBuffers(global,pBorder,pBorder%spec)   
        END IF ! global%specUsed        
#endif       
      END IF ! pBorder%iProc
    END DO ! iBorder
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying buffers done.'
    END IF ! global%verbLevel
      
  END SUBROUTINE RFLU_MPI_DestroyBuffersWrapper








! ******************************************************************************
!
! Purpose: Copy cell data needed for HM solver.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pBorder     Pointer to border
!   iVar        variable to be copied
!   cellData    Array with cell data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_HM_CopyCellData(global,pBorder,pBorder2,iVar,cellData, &
                                      cellData2)

    IMPLICIT NONE
 
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
 
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iVar
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cellData
    REAL(RFREAL), DIMENSION(:,:), INTENT(OUT) :: cellData2
    TYPE(t_border), POINTER :: pBorder,pBorder2
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icg2,icl,nVars
   
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_HM_CopyCellData',__FILE__)

! ******************************************************************************
!   Set variables
! ******************************************************************************

    nVars = SIZE(cellData,1)

    IF ( nVars /= SIZE(cellData2,1) ) THEN
      CALL ErrorStop(global,ERR_DATADIM_MISMATCH,__LINE__)
    END IF ! nVars

! ******************************************************************************
!   Copy data
! ******************************************************************************

    DO icl = 1,pBorder%nCellsSend
      icg  = pBorder%icgSend(icl)
      icg2 = pBorder2%icgRecv(icl)

      cellData2(iVar,icg2) = cellData(iVar,icg)
    END DO ! icl

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MPI_HM_CopyCellData








! ******************************************************************************
!
! Purpose: Send cell data related to Hou-Mahesh solver.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   pBorder             Pointer to border
!   iVar                Variable to be sent
!   cellDataBuff        Buffer array
!   cellData            Data array
!   tag                 Tag
!
! Output: 
!   request             Request
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_HM_ISendCellData(global,pBorder,iVar,cellDataBuff, &
                                       cellData,tag,request)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iVar,tag
    INTEGER, INTENT(OUT) :: request
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cellData
    REAL(RFREAL), DIMENSION(:), INTENT(OUT) :: cellDataBuff
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global         
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icl,nVars        
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_HM_ISendCellData',__FILE__)
    
! ******************************************************************************
!   Pack data into buffer
! ******************************************************************************

    DO icl = 1,pBorder%nCellsSend
      icg = pBorder%icgSend(icl)

      cellDataBuff(icl) = cellData(iVar,icg)
    END DO ! icl

! ******************************************************************************
!   Send data  
! ******************************************************************************

    IF ( pBorder%nCellsSend > 0 ) THEN 
      CALL MPI_ISend(cellDataBuff,pBorder%nCellsSend,MPI_RFREAL, & 
                     pBorder%iProc,tag,global%mpiComm,request,errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
      END IF ! global%error                                      
    END IF ! pBorder%nCellsSend

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_HM_ISendCellData








! ******************************************************************************
!
! Purpose: Send cell data.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   pBorder             Pointer to border
!   cellDataBuff        Buffer array
!   cellData            Data array
!   tag                 Tag
!
! Output: 
!   request             Request
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_ISendCellData(global,pBorder,cellDataBuff,cellData,tag, & 
                                    request)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: tag
    INTEGER, INTENT(OUT) :: request
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cellData
    REAL(RFREAL), DIMENSION(:,:), INTENT(OUT) :: cellDataBuff
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global         
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icl,iVar,nVars        
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_ISendCellData',__FILE__)
    
! ******************************************************************************
!   Set variables
! ******************************************************************************

    nVars = SIZE(cellData,1)

! ******************************************************************************
!   Pack data into buffer
! ******************************************************************************

    DO icl = 1,pBorder%nCellsSend
      icg = pBorder%icgSend(icl)

      DO iVar = 1,nVars
        cellDataBuff(iVar,icl) = cellData(iVar,icg)
      END DO ! iVar      
    END DO ! icl

! ******************************************************************************
!   Send data  
! ******************************************************************************

    IF ( pBorder%nCellsSend > 0 ) THEN 
      CALL MPI_ISend(cellDataBuff,pBorder%nCellsSend*nVars,MPI_RFREAL, & 
                     pBorder%iProc,tag,global%mpiComm,request,errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
      END IF ! global%error                                      
    END IF ! pBorder%nCellsSend

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_ISendCellData








! ******************************************************************************
!
! Purpose: Send one dimensional cell data.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   pBorder             Pointer to border
!   cellDataBuff        Buffer array
!   cellData            Data array
!   tag                 Tag
!
! Output: 
!   request             Request
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_ISendCellData1d(global,pBorder,cellDataBuff,cellData, &
                                      tag,request)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: tag
    INTEGER, INTENT(OUT) :: request
    REAL(RFREAL), DIMENSION(:), INTENT(IN) :: cellData
    REAL(RFREAL), DIMENSION(:), INTENT(OUT) :: cellDataBuff
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global         
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icl,iVar
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_ISendCellData1d',__FILE__)
    
! ******************************************************************************
!   Pack data into buffer
! ******************************************************************************

    DO icl = 1,pBorder%nCellsSend
      icg = pBorder%icgSend(icl)

      cellDataBuff(icl) = cellData(icg)
    END DO ! icl

! ******************************************************************************
!   Send data  
! ******************************************************************************

    IF ( pBorder%nCellsSend > 0 ) THEN 
      CALL MPI_ISend(cellDataBuff,pBorder%nCellsSend,MPI_RFREAL, & 
                     pBorder%iProc,tag,global%mpiComm,request,errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
      END IF ! global%error                                      
    END IF ! pBorder%nCellsSend

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_ISendCellData1d







! ******************************************************************************
!
! Purpose: Send cell data.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   pBorder             Pointer to border
!   cellDataBuff        Buffer array
!   cellData            Data array
!   tag                 Tag
!
! Output: 
!   request             Request
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_ISendCellData3d(global,pBorder,cellDataBuff,cellData,tag, & 
                                    request)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: tag
    INTEGER, INTENT(OUT) :: request
    REAL(RFREAL), DIMENSION(:,:,:), INTENT(IN) :: cellData
    REAL(RFREAL), DIMENSION(:,:,:), INTENT(OUT) :: cellDataBuff
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global         
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icl,iVar,iVar1     
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_ISendCellData3d',__FILE__)
    
! ******************************************************************************
!   Pack data into buffer
! ******************************************************************************

    DO icl = 1,pBorder%nCellsSend
      icg = pBorder%icgSend(icl)

      iVar1 = 1
! TEMPORARY: manoj: Need to find a permanent solution
!      DO iVar = CV_MIXT_XVEL,CV_MIXT_ZVEL
      DO iVar = 1,3 
        cellDataBuff(XCOORD,iVar1,icl) = cellData(XCOORD,iVar,icg)
        cellDataBuff(YCOORD,iVar1,icl) = cellData(YCOORD,iVar,icg)
        cellDataBuff(ZCOORD,iVar1,icl) = cellData(ZCOORD,iVar,icg)

        iVar1 = iVar1 + 1
      END DO ! iVar      
    END DO ! icl

! ******************************************************************************
!   Send data  
! ******************************************************************************

    IF ( pBorder%nCellsSend > 0 ) THEN 
      CALL MPI_ISend(cellDataBuff,pBorder%nCellsSend*9,MPI_RFREAL, & 
                     pBorder%iProc,tag,global%mpiComm,request,errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
      END IF ! global%error                                      
    END IF ! pBorder%nCellsSend

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_ISendCellData3d








! ******************************************************************************
!
! Purpose: Recreate iPclSend buffers.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pBorder     Pointer to border (Optional)
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_RecreateBufferIPclSend(pRegion,pBorder)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_border), POINTER, OPTIONAL :: pBorder
    TYPE(t_region), POINTER :: pRegion 
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iPcl,iVar,nPclsSendMax,nPclsSendMaxOld,nVars
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iPclSendTemp     
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid     
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_RecreateBufferIPclSend',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

#ifdef PLAG
! ******************************************************************************
!   Set variables 
! ******************************************************************************

    nVars           = SIZE(pBorder%iPclSend,1)
    nPclsSendMaxOld = SIZE(pBorder%iPclSend,2)
     
    pBorder%nPclsSendMax = & 
      NINT(1.20_RFREAL*REAL(pBorder%nPclsSend,KIND=RFREAL))
      ! Subbu - Debug
      !NINT(10.20_RFREAL*REAL(pBorder%nPclsSend,KIND=RFREAL))
      ! Subbu - End Debug
! ******************************************************************************
!   Allocate temporary array
! ******************************************************************************

    ALLOCATE(iPclSendTemp(nVars,pBorder%nPclsSendMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'iPclSendTemp')
    END IF ! global%error

    
! ******************************************************************************
!   Copy array
! ******************************************************************************
   
    DO iPcl = 1,nPclsSendMaxOld
      DO iVar = 1,nVars
         iPclSendTemp(iVar,iPcl) = pBorder%iPclSend(iVar,iPcl)
      END DO ! iVar
    END DO ! iPcl 

! ******************************************************************************
!   Deallocate iPclSend array
! ******************************************************************************      

    CALL RFLU_MPI_DestroyBufferIPclSend(pRegion,pBorder)

! ******************************************************************************
!   Reallocate iPclSend array
! ******************************************************************************      

    CALL RFLU_MPI_CreateBufferIPclSend(pRegion,pBorder)   

! ******************************************************************************
!   Copy array
! ******************************************************************************
    
    DO iPcl = 1,pBorder%nPclsSend
      DO iVar = 1,nVars
        pBorder%iPclSend(iVar,iPcl) = iPclSendTemp(iVar,iPcl)
      END DO ! iVar
    END DO ! iPcl     

! ******************************************************************************
!   Deallocate temporary array
! ******************************************************************************

    DEALLOCATE(iPclSendTemp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'iPclSendTemp')
    END IF ! global%error
#endif

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
      
  END SUBROUTINE RFLU_MPI_RecreateBufferIPclSend






! ******************************************************************************
!
! Purpose: Wrapper for copying data needed for HM solver.
!
! Description: None.
!           
! Input:    
!   pRegion     Pointer to region
!   iVar        Variable to be copied
!           
! Output: None.
!           
! Notes: None.
!           
! ******************************************************************************
  
  SUBROUTINE RFLU_MPI_HM_CopyWrapper(regions,iVar)
  
    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
  
! ==============================================================================
!   Arguments
! ==============================================================================
  
    INTEGER, INTENT(IN) :: iVar
    TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
!   Local variables
! ==============================================================================
  
    INTEGER :: errorFlag,iBorder,iBorder2,iReg,iReg2
    TYPE(t_border), POINTER :: pBorder,pBorder2
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_region), POINTER :: pRegion,pRegion2

! ******************************************************************************
!   Start
! ******************************************************************************

    global => regions(0)%global

    CALL RegisterFunction(global,'RFLU_MPI_HM_CopyWrapper',__FILE__)

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::HM_CopyWrapper")
#endif

! ******************************************************************************
!   Set pointers
! ******************************************************************************

! ******************************************************************************
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
!         Check dimensions
! ------------------------------------------------------------------------------

          IF ( pBorder%nCellsSend /= pBorder2%nCellsRecv ) THEN
            CALL ErrorStop(global,ERR_BUFFERDIM_MISMATCH,__LINE__)
          END IF ! pBorder

! ------------------------------------------------------------------------------
!         Mixture
! ------------------------------------------------------------------------------

          SELECT CASE (iVar)
            CASE (1,2,3,4,5)
              CALL RFLU_MPI_HM_CopyCellData(global,pBorder,pBorder2,iVar, &
                                            pRegion%mixt%cv,pRegion2%mixt%cv)
            CASE (6)
              CALL RFLU_MPI_CopyCellData1d(global,pBorder,pBorder2, &
                                           pRegion%mixt%delP,pRegion2%mixt%delP)
            CASE (7)
              CALL RFLU_MPI_HM_CopyCellData(global,pBorder,pBorder2,1, &
                                            pRegion%mixt%cvOld, &
                                            pRegion2%mixt%cvOld)
            CASE (8)
              CALL RFLU_MPI_HM_CopyCellData(global,pBorder,pBorder2,5, &
                                            pRegion%mixt%cvOld, &
                                            pRegion2%mixt%cvOld)
            CASE (9)
              CALL RFLU_MPI_CopyCellData3d(global,pBorder,pBorder2, &
                                            pRegion%mixt%gradCell, &
                                            pRegion2%mixt%gradCell)
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! iVar

        END IF ! pBorder
      END DO ! iBorder
    END DO ! iReg

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::HM_CopyWrapper")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MPI_HM_CopyWrapper








! ******************************************************************************
!
! Purpose: Wrapper for sending data for HM solver.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   iVar        variable to be sent
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_HM_ISendWrapper(pRegion,iVar)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iVar
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iBorder
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_HM_ISendWrapper',__FILE__)

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::HM_ISendWrapper")
#endif

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

! ==============================================================================
!     Send data if not on same process
! ==============================================================================

      IF ( pBorder%iProc /= global%myProcid ) THEN

! ------------------------------------------------------------------------------
!       Mixture
! ------------------------------------------------------------------------------

        SELECT CASE (iVar)
          CASE (1,2,3,4,5)
            CALL RFLU_MPI_HM_ISendCellData(global,pBorder,iVar, &
                                           pBorder%mixt%sendBuff1d, &
                                           pRegion%mixt%cv, &
                                           pBorder%mixt%tag, &
                                           pBorder%mixt%sendRequest)
          CASE (6)
            CALL RFLU_MPI_ISendCellData1d(global,pBorder, &
                                          pBorder%mixt%sendBuff1d, &
                                          pRegion%mixt%delP, &
                                          pBorder%mixt%tag, &
                                          pBorder%mixt%sendRequest)
          CASE (7)
            CALL RFLU_MPI_HM_ISendCellData(global,pBorder,1, &
                                           pBorder%mixt%sendBuff1d, &
                                           pRegion%mixt%cvOld, &
                                           pBorder%mixt%tag, &
                                           pBorder%mixt%sendRequest)
          CASE (8)
            CALL RFLU_MPI_HM_ISendCellData(global,pBorder,5, &
                                           pBorder%mixt%sendBuff1d, &
                                           pRegion%mixt%cvOld, &
                                           pBorder%mixt%tag, &
                                           pBorder%mixt%sendRequest)
          CASE (9)
            CALL RFLU_MPI_ISendCellData3d(global,pBorder, &
                                          pBorder%mixt%sendBuff3d, &
                                          pRegion%mixt%gradCell, &
                                          pBorder%mixt%tag, &
                                          pBorder%mixt%sendRequest)
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! iVar

      END IF ! pBorder
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::HM_ISendWrapper")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MPI_HM_ISendWrapper








! ******************************************************************************
!
! Purpose: Wrapper for receiving data for HM solver.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   iVar        variable to be sent
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_HM_RecvWrapper(pRegion,iVar)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iVar
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iBorder
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_HM_RecvWrapper',__FILE__)

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::HM_RecvWrapper")
#endif

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

! ==============================================================================
!     Receive data if not on same process
! ==============================================================================

      IF ( pBorder%iProc /= global%myProcid ) THEN

! ------------------------------------------------------------------------------
!       Mixture
! ------------------------------------------------------------------------------

        SELECT CASE (iVar)
          CASE (1,2,3,4,5)
            CALL RFLU_MPI_HM_RecvCellData(global,pBorder,iVar, &
                                          pBorder%mixt%recvBuff1d, &
                                          pRegion%mixt%cv, &
                                          pBorder%mixt%tag)
          CASE (6)
            CALL RFLU_MPI_RecvCellData1d(global,pBorder, &
                                         pBorder%mixt%recvBuff1d, &
                                         pRegion%mixt%delP, &
                                         pBorder%mixt%tag)
          CASE (7)
            CALL RFLU_MPI_HM_RecvCellData(global,pBorder,1, &
                                          pBorder%mixt%recvBuff1d, &
                                          pRegion%mixt%cvOld, &
                                          pBorder%mixt%tag)
          CASE (8)
            CALL RFLU_MPI_HM_RecvCellData(global,pBorder,5, &
                                          pBorder%mixt%recvBuff1d, &
                                          pRegion%mixt%cvOld, &
                                          pBorder%mixt%tag)
          CASE (9)
            CALL RFLU_MPI_RecvCellData3d(global,pBorder, &
                                         pBorder%mixt%recvBuff3d, &
                                         pRegion%mixt%gradCell, &
                                         pBorder%mixt%tag)
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! iVar

      END IF ! pBorder
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::HM_RecvWrapper")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MPI_HM_RecvWrapper








! ******************************************************************************
!
! Purpose: Receive cell data related to Hou-Mahesh solver.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   pBorder             Pointer to border
!   iVar                variable to be recieved
!   cellDataBuff        Buffer array
!   cellData            Data array
!   tag                 Tag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_HM_RecvCellData(global,pBorder,iVar,cellDataBuff, &
                                      cellData,tag)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iVar,tag
    REAL(RFREAL), DIMENSION(:), INTENT(IN) :: cellDataBuff
    REAL(RFREAL), DIMENSION(:,:), INTENT(OUT) :: cellData
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global         
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icl,nVars
    INTEGER :: status(MPI_STATUS_SIZE)        
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_HM_RecvCellData',__FILE__)

! ******************************************************************************
!   Set variables
! ******************************************************************************

    nVars = SIZE(cellData,1)
    
! ******************************************************************************
!   Recv data 
! ******************************************************************************

    IF ( pBorder%nCellsRecv > 0 ) THEN 
      CALL MPI_Recv(cellDataBuff,pBorder%nCellsRecv,MPI_RFREAL, & 
                    pBorder%iProc,tag,global%mpiComm,status,errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
      END IF ! global%error  
    END IF ! pBorder%nCellsRecv

! ******************************************************************************
!   Unpack data from buffer
! ******************************************************************************

    DO icl = 1,pBorder%nCellsRecv
      icg = pBorder%icgRecv(icl)

      cellData(iVar,icg) = cellDataBuff(icl) 
    END DO ! icl
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_HM_RecvCellData






! ******************************************************************************
!
! Purpose: Wrapper for sending data.
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

  SUBROUTINE RFLU_MPI_ISendWrapper(pRegion)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iBorder
    TYPE(t_border), POINTER :: pBorder 
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid             
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_ISendWrapper',__FILE__)

#ifdef ROCPROF 
    CALL FPROFILER_BEGINS("RFLU::ISendWrapper")
#endif

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid
    
! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)
    
! ==============================================================================
!     Send data if not on same process
! ==============================================================================    
    
      IF ( pBorder%iProc /= global%myProcid ) THEN 
      
! ------------------------------------------------------------------------------
!       Mixture
! ------------------------------------------------------------------------------      

        CALL RFLU_MPI_ISendCellData(global,pBorder,pBorder%mixt%sendBuff, & 
                                    pRegion%mixt%cv,pBorder%mixt%tag, &
                                    pBorder%mixt%sendRequest)

! ------------------------------------------------------------------------------
!       Physical modules
! ------------------------------------------------------------------------------      

#ifdef SPEC
        IF ( global%specUsed .EQV. .TRUE. ) THEN 
          CALL RFLU_MPI_ISendCellData(global,pBorder,pBorder%spec%sendBuff, & 
                                      pRegion%spec%cv,pBorder%spec%tag, &
                                      pBorder%spec%sendRequest)        
        END IF ! global%specUsed
#endif  
      END IF ! pBorder
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF 
    CALL FPROFILER_ENDS("RFLU::ISendWrapper")
#endif

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_ISendWrapper








! ******************************************************************************
!
! Purpose: Receive cell data.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   pBorder             Pointer to border
!   cellDataBuff        Buffer array
!   cellData            Data array
!   tag                 Tag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_RecvCellData(global,pBorder,cellDataBuff,cellData,tag)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: tag
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cellDataBuff
    REAL(RFREAL), DIMENSION(:,:), INTENT(OUT) :: cellData
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global         
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icl,iVar,nVars
    INTEGER :: status(MPI_STATUS_SIZE)        
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_RecvCellData',__FILE__)

! ******************************************************************************
!   Set variables
! ******************************************************************************

    nVars = SIZE(cellData,1)
    
! ******************************************************************************
!   Recv data 
! ******************************************************************************

    IF ( pBorder%nCellsRecv > 0 ) THEN 
      CALL MPI_Recv(cellDataBuff,pBorder%nCellsRecv*nVars,MPI_RFREAL, & 
                    pBorder%iProc,tag,global%mpiComm,status,errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
      END IF ! global%error  
    END IF ! pBorder%nCellsRecv

! ******************************************************************************
!   Unpack data from buffer
! ******************************************************************************

    DO icl = 1,pBorder%nCellsRecv
      icg = pBorder%icgRecv(icl)

      DO iVar = 1,nVars
        cellData(iVar,icg) = cellDataBuff(iVar,icl) 
      END DO ! iVar      
    END DO ! icl
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_RecvCellData







! ******************************************************************************
!
! Purpose: Receive one dimensional cell data.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   pBorder             Pointer to border
!   cellDataBuff        Buffer array
!   cellData            Data array
!   tag                 Tag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_RecvCellData1d(global,pBorder,cellDataBuff,cellData,tag)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: tag
    REAL(RFREAL), DIMENSION(:), INTENT(IN) :: cellDataBuff
    REAL(RFREAL), DIMENSION(:), INTENT(OUT) :: cellData
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global         
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icl
    INTEGER :: status(MPI_STATUS_SIZE)        
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_RecvCellData1d',__FILE__)

! ******************************************************************************
!   Recv data 
! ******************************************************************************

    IF ( pBorder%nCellsRecv > 0 ) THEN 
      CALL MPI_Recv(cellDataBuff,pBorder%nCellsRecv,MPI_RFREAL, & 
                    pBorder%iProc,tag,global%mpiComm,status,errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
      END IF ! global%error  
    END IF ! pBorder%nCellsRecv

! ******************************************************************************
!   Unpack data from buffer
! ******************************************************************************

    DO icl = 1,pBorder%nCellsRecv
      icg = pBorder%icgRecv(icl)

      cellData(icg) = cellDataBuff(icl) 
    END DO ! icl
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_RecvCellData1d






! ******************************************************************************
!
! Purpose: Receive cell data.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   pBorder             Pointer to border
!   cellDataBuff        Buffer array
!   cellData            Data array
!   tag                 Tag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_RecvCellData3d(global,pBorder,cellDataBuff,cellData,tag)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: tag
    REAL(RFREAL), DIMENSION(:,:,:), INTENT(IN) :: cellDataBuff
    REAL(RFREAL), DIMENSION(:,:,:), INTENT(OUT) :: cellData
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global         
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icl,iVar,iVar1
    INTEGER :: status(MPI_STATUS_SIZE)        
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_RecvCellData3d',__FILE__)

! ******************************************************************************
!   Recv data 
! ******************************************************************************

    IF ( pBorder%nCellsRecv > 0 ) THEN 
      CALL MPI_Recv(cellDataBuff,pBorder%nCellsRecv*9,MPI_RFREAL, & 
                    pBorder%iProc,tag,global%mpiComm,status,errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
      END IF ! global%error  
    END IF ! pBorder%nCellsRecv

! ******************************************************************************
!   Unpack data from buffer
! ******************************************************************************

    DO icl = 1,pBorder%nCellsRecv
      icg = pBorder%icgRecv(icl)

      iVar1 = 1
! TEMPORARY: manoj: Need to find a permanent solution
!      DO iVar = CV_MIXT_XVEL,CV_MIXT_ZVEL
      DO iVar = 1,3 
        cellData(XCOORD,iVar,icg) = cellDataBuff(XCOORD,iVar1,icl) 
        cellData(YCOORD,iVar,icg) = cellDataBuff(YCOORD,iVar1,icl) 
        cellData(ZCOORD,iVar,icg) = cellDataBuff(ZCOORD,iVar1,icl) 

        iVar1 = iVar1 + 1
      END DO ! iVar      
    END DO ! icl
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_RecvCellData3d







! ******************************************************************************
!
! Purpose: Wrapper for receiving data.
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

  SUBROUTINE RFLU_MPI_RecvWrapper(pRegion)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iBorder
    TYPE(t_border), POINTER :: pBorder 
    TYPE(t_global), POINTER :: global 
    TYPE(t_grid), POINTER :: pGrid     
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_RecvWrapper',__FILE__)

#ifdef ROCPROF 
    CALL FPROFILER_BEGINS("RFLU::RecvWrapper")
#endif

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid
    
! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

! ==============================================================================
!     Receive data if not on same process
! ==============================================================================    

      IF ( pBorder%iProc /= global%myProcid ) THEN 

! ------------------------------------------------------------------------------
!       Mixture
! ------------------------------------------------------------------------------      

        CALL RFLU_MPI_RecvCellData(global,pBorder,pBorder%mixt%recvBuff, & 
                                   pRegion%mixt%cv,pBorder%mixt%tag)      

! ------------------------------------------------------------------------------
!       Physical modules
! ------------------------------------------------------------------------------      

#ifdef SPEC
        IF ( global%specUsed .EQV. .TRUE. ) THEN 
          CALL RFLU_MPI_RecvCellData(global,pBorder,pBorder%spec%recvBuff, & 
                                     pRegion%spec%cv,pBorder%spec%tag)           
        END IF ! global%specUsed
#endif         
      END IF ! pBorder
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF 
    CALL FPROFILER_ENDS("RFLU::RecvWrapper")
#endif

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_RecvWrapper





! ******************************************************************************
!
! Purpose: Set tag.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   iReg1       Index of first region
!   iReg2       Index of second region
!   iMsg        Index of message
!   tagMax      Maximum allowed value of tag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  INTEGER FUNCTION RFLU_MPI_SetTag(global,iReg1,iReg2,iMsg,tagMax)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iMsg,iReg1,iReg2,tagMax
    TYPE(t_global), POINTER :: global 
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: iRegMax,iRegMin  
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_SetTag',__FILE__)

! ******************************************************************************
!   Set tag
! ******************************************************************************

    iRegMax = MAX(iReg1,iReg2)
    iRegMin = MIN(iReg1,iReg2)

    !RFLU_MPI_SetTag = iRegMin + (iRegMax-1)*global%nRegions &
    !                + (iRegMax-1)*(iMsg-1)*global%nRegions*global%nRegions 

    ! BBR - begin - fixing MPI_SetTag formula & Error test
    RFLU_MPI_SetTag = iRegMin*(2*(global%nRegions-1) - iRegMin + 1)/2 & 
                      + (iRegMax - iRegMin) &
                      + (iMsg-1)*global%nRegions*(global%nRegions-1)/2

!    IF ( RFLU_MPI_SetTag > tagMax ) THEN 
    IF ( RFLU_MPI_SetTag > tagMax .OR. RFLU_MPI_SetTag < 0 ) THEN
      CALL ErrorStop(global,ERR_MPI_TAGMAX,__LINE__)
    END IF ! RFLU_MPI_SetTag

    ! BBR - end
 
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END FUNCTION RFLU_MPI_SetTag










! ******************************************************************************
!
! Purpose: Wrapper for setting tags.
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

  SUBROUTINE RFLU_MPI_SetTagsWrapper(pRegion)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion
       
! ==============================================================================
!   Local variables
! ==============================================================================

    LOGICAL :: dummyLogical
    INTEGER :: errorFlag,iBorder,iMsg,tagMax
    TYPE(t_border), POINTER :: pBorder 
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid             
  
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_SetTagsWrapper',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Get maximum allowed value of tag. NOTE must always use MPI_COMM_WORLD - 
!   cannot use global%mpiComm because may be split off in GENx computations
!   and get zero tagMax as a result.
! ******************************************************************************

    !BBR - Begin - Swap CALL name becasue of deprecated MPI function
    ! MPI_Attr_get --> deprecated but, working on all machine, in particular
    ! VULCAN
    CALL MPI_Attr_get(MPI_COMM_WORLD,MPI_TAG_UB,tagMax,dummyLogical,errorFlag)
    ! MPI_MPI_Comm_get_attr --> up to date call, but does not work with MPI
    ! version on VULCAN
    !CALL MPI_Comm_get_attr(MPI_COMM_WORLD,MPI_TAG_UB,tagMax &
    !                      ,dummyLogical,errorFlag)

    !BBR - End
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
    END IF ! global%error

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)
    
! ==============================================================================
!     Mixture
! ==============================================================================
            
      iMsg = 1

      pBorder%mixt%tag = RFLU_MPI_SetTag(global,pRegion%iRegionGlobal, & 
                                         pBorder%iRegionGlobal,iMsg,tagMax)

! ==============================================================================
!     Physical modules
! ==============================================================================

#ifdef SPEC
      iMsg = iMsg + 1

      pBorder%spec%tag = RFLU_MPI_SetTag(global,pRegion%iRegionGlobal, & 
                                         pBorder%iRegionGlobal,iMsg,tagMax)            
#endif        

#ifdef PLAG
      iMsg = iMsg + 1

      pBorder%plag%tagCount = RFLU_MPI_SetTag(global,pRegion%iRegionGlobal, & 
                                              pBorder%iRegionGlobal,iMsg,tagMax)            
      iMsg = iMsg + 1

      pBorder%plag%tagInt = RFLU_MPI_SetTag(global,pRegion%iRegionGlobal, & 
                                            pBorder%iRegionGlobal,iMsg,tagMax)            
      iMsg = iMsg + 1

      pBorder%plag%tag = RFLU_MPI_SetTag(global,pRegion%iRegionGlobal, & 
                                         pBorder%iRegionGlobal,iMsg,tagMax)            
#endif        
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_SetTagsWrapper







! ******************************************************************************
!
! Purpose: Copy particle volume fraction data.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pBorder     Pointer to border
!   cellData    Array with cell data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_PLAG_CopyCellData(global,pBorder,pBorder2,cellData, &
                                      cellData2)

    IMPLICIT NONE
 
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
 
! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cellData
    REAL(RFREAL), DIMENSION(:,:), INTENT(OUT) :: cellData2
    TYPE(t_border), POINTER :: pBorder,pBorder2
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icg2,icl,nVars
   
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_PLAG_CopyCellData',__FILE__)

! ******************************************************************************
!   Set variables
! ******************************************************************************

    nVars = SIZE(cellData,1)

    IF ( nVars /= SIZE(cellData2,1) ) THEN
      CALL ErrorStop(global,ERR_DATADIM_MISMATCH,__LINE__)
    END IF ! nVars

! ******************************************************************************
!   Copy data
! ******************************************************************************

    DO icl = 1,pBorder%nCellsSend
      icg  = pBorder%icgSend(icl)
      icg2 = pBorder2%icgRecv(icl)

      cellData2(1,icg2) = cellData(1,icg)
    END DO ! icl

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MPI_PLAG_CopyCellData







! ******************************************************************************
!
! Purpose: Wrapper for copying particle volume fraction data.
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
  
  SUBROUTINE RFLU_MPI_PLAG_CopyWrapper(regions)
  
    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
  
! ==============================================================================
!   Arguments
! ==============================================================================
  
    TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
!   Local variables
! ==============================================================================
  
    INTEGER :: errorFlag,iBorder,iBorder2,iReg,iReg2
    TYPE(t_border), POINTER :: pBorder,pBorder2
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_region), POINTER :: pRegion,pRegion2

! ******************************************************************************
!   Start
! ******************************************************************************

    global => regions(0)%global

    CALL RegisterFunction(global,'RFLU_MPI_PLAG_CopyWrapper',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

! ******************************************************************************
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
!         Check dimensions
! ------------------------------------------------------------------------------

          IF ( pBorder%nCellsSend /= pBorder2%nCellsRecv ) THEN
            CALL ErrorStop(global,ERR_BUFFERDIM_MISMATCH,__LINE__)
          END IF ! pBorder

! ------------------------------------------------------------------------------
!         Lagrangian particle volume fraction
! ------------------------------------------------------------------------------

          CALL RFLU_MPI_PLAG_CopyCellData(global,pBorder,pBorder2, &
                                          pRegion%plag%vFracE, &
                                          pRegion2%plag%vFracE)
        END IF ! pBorder
      END DO ! iBorder
    END DO ! iReg

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MPI_PLAG_CopyWrapper







! ******************************************************************************
!
! Purpose: Send cell data related to PLAG.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   pBorder             Pointer to border
!   cellDataBuff        Buffer array
!   cellData            Data array
!   tag                 Tag
!
! Output: 
!   request             Request
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_PLAG_ISendCellData(global,pBorder,cellDataBuff,cellData, &
                                         tag,request)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: tag
    INTEGER, INTENT(OUT) :: request
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cellData
    REAL(RFREAL), DIMENSION(:), INTENT(OUT) :: cellDataBuff
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global         
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icl,nVars        
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_PLAG_ISendCellData',__FILE__)
    
! ******************************************************************************
!   Pack data into buffer
! ******************************************************************************

    DO icl = 1,pBorder%nCellsSend
      icg = pBorder%icgSend(icl)

      cellDataBuff(icl) = cellData(1,icg)
    END DO ! icl

! ******************************************************************************
!   Send data  
! ******************************************************************************

    IF ( pBorder%nCellsSend > 0 ) THEN 
      CALL MPI_ISend(cellDataBuff,pBorder%nCellsSend,MPI_RFREAL, & 
                     pBorder%iProc,tag,global%mpiComm,request,errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
      END IF ! global%error                                      
    END IF ! pBorder%nCellsSend

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_PLAG_ISendCellData



  



! ******************************************************************************
!
! Purpose: Wrapper for sending data for PLAG.
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

  SUBROUTINE RFLU_MPI_PLAG_ISendWrapper(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iBorder
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_PLAG_ISendWrapper',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

! ==============================================================================
!     Send data if not on same process
! ==============================================================================

      IF ( pBorder%iProc /= global%myProcid ) THEN

! ------------------------------------------------------------------------------
!       Mixture
! ------------------------------------------------------------------------------

        CALL RFLU_MPI_PLAG_ISendCellData(global,pBorder, &
                                         pBorder%mixt%sendBuff1d, &
                                         pRegion%plag%vFracE, &
                                         pBorder%mixt%tag, &
                                         pBorder%mixt%sendRequest)
      END IF ! pBorder
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MPI_PLAG_ISendWrapper








! ******************************************************************************
!
! Purpose: Receive cell data related to PLAG.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   pBorder             Pointer to border
!   cellDataBuff        Buffer array
!   cellData            Data array
!   tag                 Tag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_PLAG_RecvCellData(global,pBorder,cellDataBuff, &
                                      cellData,tag)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: tag
    REAL(RFREAL), DIMENSION(:), INTENT(IN) :: cellDataBuff
    REAL(RFREAL), DIMENSION(:,:), INTENT(OUT) :: cellData
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global         
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icl,nVars
    INTEGER :: status(MPI_STATUS_SIZE)        
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_PLAG_RecvCellData',__FILE__)

! ******************************************************************************
!   Set variables
! ******************************************************************************

    nVars = SIZE(cellData,1)
    
! ******************************************************************************
!   Recv data 
! ******************************************************************************

    IF ( pBorder%nCellsRecv > 0 ) THEN 
      CALL MPI_Recv(cellDataBuff,pBorder%nCellsRecv,MPI_RFREAL, & 
                    pBorder%iProc,tag,global%mpiComm,status,errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
      END IF ! global%error  
    END IF ! pBorder%nCellsRecv

! ******************************************************************************
!   Unpack data from buffer
! ******************************************************************************

    DO icl = 1,pBorder%nCellsRecv
      icg = pBorder%icgRecv(icl)

      cellData(1,icg) = cellDataBuff(icl) 
    END DO ! icl
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_PLAG_RecvCellData






  

! ******************************************************************************
!
! Purpose: Wrapper for receiving data for PLAG.
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

  SUBROUTINE RFLU_MPI_PLAG_RecvWrapper(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iBorder
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_PLAG_RecvWrapper',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

! ==============================================================================
!     Receive data if not on same process
! ==============================================================================

      IF ( pBorder%iProc /= global%myProcid ) THEN

! ------------------------------------------------------------------------------
!       Mixture
! ------------------------------------------------------------------------------

        CALL RFLU_MPI_PLAG_RecvCellData(global,pBorder, &
                                        pBorder%mixt%recvBuff1d, &
                                        pRegion%plag%vFracE, &
                                        pBorder%mixt%tag)
      END IF ! pBorder
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MPI_PLAG_RecvWrapper








! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModMPI


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModMPI.F90,v $
! Revision 1.2  2015/07/23 23:11:18  brollin
! 1) The pressure coefficient of the  collision model has been changed back to its original form
! 2) New options in the format of the grid and solutions have been added. Now the user can choose the endianness, and convert from one to the over in rfluconv.
! 3) The solutions are now stored in folders named by timestamp or iteration number.
! 4) The address enty in the hashtable has been changed to an integer(8) for cases when the grid becomes very large.
! 5) RFLU_WritePM can now compute PM2 on the fly for the Macroscale problem
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:36  brollin
! New Stable version
!
! Revision 1.6  2009/09/21 18:06:07  mparmar
! Bug fix in reading gradient data for virtual cells
!
! Revision 1.5  2008/12/06 08:43:41  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:16:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/12/05 13:22:12  haselbac
! Added PRIVATE to RCSIdentString, ifort compiler on vonkarman complained
!
! Revision 1.2  2007/11/28 23:05:24  mparmar
! Added routines for communication in SOLV_IMPLICIT_HM
!
! Revision 1.1  2007/04/09 18:49:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:41  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.14  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.13  2006/02/09 03:37:44  haselbac
! Bug fix: Must always use MPI_COMM_WORLD to get max tag
!
! Revision 1.12  2005/12/14 21:50:32  haselbac
! Cosmetics
!
! Revision 1.11  2005/12/14 21:20:28  fnajjar
! Added subroutine and made changes for dynamic allocation of iPclsSend
!
! Revision 1.10  2005/12/13 23:30:53  haselbac
! Cosmetics
!
! Revision 1.9  2005/12/13 23:06:56  fnajjar
! Added defs of tag for PLAG
!
! Revision 1.8  2005/12/08 03:01:01  haselbac
! Major bug fix: spec tag was not set
!
! Revision 1.7  2005/12/03 19:48:23  haselbac
! Bug fix: Only clear request if have indeed sent message, cosmetics
!
! Revision 1.6  2005/09/19 18:40:37  haselbac
! Added IFs for border sizes before send and recv
!
! Revision 1.5  2005/07/08 15:01:29  haselbac
! Added profiling calls
!
! Revision 1.4  2005/05/26 22:01:01  haselbac
! Fixed two serious bugs: status is array in recv and wait
!
! Revision 1.3  2005/05/18 22:12:04  fnajjar
! ACH: Added routines to create and destroy iPclSend buffers
!
! Revision 1.2  2005/04/29 00:05:47  haselbac
! Added MPI_Wait routines
!
! Revision 1.1  2005/04/15 15:06:43  haselbac
! Initial revision
!
! ******************************************************************************

