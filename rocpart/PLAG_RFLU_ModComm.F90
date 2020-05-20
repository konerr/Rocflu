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
! Purpose: Suite of routines for communication between regions with RFLU solver.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: PLAG_RFLU_ModComm.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE PLAG_RFLU_ModComm

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt_input
  USE ModBorder, ONLY: t_border,t_border_data
  USE ModPartLag, ONLY: t_plag,t_plag_input
  USE ModMPI
  
  USE PLAG_ModParameters
  USE ModInterfaces, ONLY: RFLU_DecidePrint
  USE PLAG_ModInterfaces, ONLY: PLAG_UpdateDataStruct
  USE PLAG_ModReallocateMemory
  USE PLAG_RFLU_ModFindCells

! TEMPORARY
  USE PLAG_ModInterfaces, ONLY: PLAG_CalcDerivedVariables
! END TEMPORARY

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: PLAG_RFLU_CommDriver, &
            PLAG_RFLU_InitSendCounters
   
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: PLAG_RFLU_ModComm.F90,v $ $Revision: 1.1.1.1 $' 
                      
  INTEGER, PARAMETER, PRIVATE :: REQUEST_TYPE_COUNTER = 1, &
                                 REQUEST_TYPE_DATA    = 2

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

  SUBROUTINE PLAG_RFLU_ClearRequest(global,request)

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

    CALL RegisterFunction(global,'PLAG_RFLU_ClearRequest',__FILE__)

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
  
  END SUBROUTINE PLAG_RFLU_ClearRequest









! ******************************************************************************
!
! Purpose: Wrapper for clearing send requests.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   iReqFlag    Request flag
!
! Output: None.
!
! Notes: 
!   1. When clearing the data request, the sending side initiates
!      communication; hence, we need to check if nPclsSend is not null.
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLU_ClearRequestWrapper(pRegion,iReqFlag)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iReqFlag
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

    CALL RegisterFunction(global,'PLAG_RFLU_ClearRequestWrapper',__FILE__)

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
        SELECT CASE (iReqFlag)
          CASE(REQUEST_TYPE_COUNTER)
            CALL PLAG_RFLU_ClearRequest(global,pBorder%plag%sendRequestCount)  
          
          CASE(REQUEST_TYPE_DATA)          
            IF ( pBorder%nPclsSend == 0 ) CYCLE
            CALL PLAG_RFLU_ClearRequest(global,pBorder%plag%sendRequest)
            CALL PLAG_RFLU_ClearRequest(global,pBorder%plag%sendRequestInt)  
        END SELECT ! iReqFlag
      END IF ! pBorder

    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE PLAG_RFLU_ClearRequestWrapper









! ******************************************************************************
!
! Purpose: Driver routine for communication.
!
! Description: None.
!
! Input:
!   regions     Pointer to all regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLU_CommDriver(regions)

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

    LOGICAL :: doPrint
    INTEGER :: errorFlag,iBorder,iPclBeg,iPclEnd,iReg,istage,loopCounter, &
               nPclsRecvPerReg
    REAL(RFREAL), PARAMETER :: PLAG_REALLOC_FACT = 1.20_RFREAL
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global 
    TYPE(t_grid), POINTER :: pGrid     
    TYPE(t_region), POINTER :: pRegion
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => regions(0)%global

    CALL RegisterFunction(global,'PLAG_RFLU_CommDriver',__FILE__)

! ******************************************************************************
!   Initialize variables
! ******************************************************************************
    
#ifdef ROCPROF
  CALL FPROFILER_BEGINS("PLAG_RFLU::PLAG_RFLU_CommDriver")
#endif

    loopCounter = 0

! ******************************************************************************
!   Loop till all particles have been communicated
! ******************************************************************************

    commLoop: DO 

! ==============================================================================
!     Find if any region need to communicate particles 
!      and exit commLoop if null
! ==============================================================================    

#ifdef ROCPROF
      CALL FPROFILER_BEGINS("PLAG_RFLU::TotalnPclsComm")
#endif

      CALL PLAG_RFLU_TotalnPclsComm(regions)

!BBR  -- making the code less chatty
!      IF ( global%myProcid == MASTERPROC   .AND. &
!           global%verbLevel > VERBOSE_LOW  .AND. &
!           global%nPclsCommTot /= 0              ) THEN
!        WRITE(STDOUT,'(A,3X,A,1X,1PE12.5,2(1X,I10))') &
!              SOLVER_NAME,'Current Loop Counter & nPclsComm:',&
!              global%currentTime,loopCounter,global%nPclsCommTot
!      END IF ! global%myProcid

#ifdef ROCPROF
      CALL FPROFILER_ENDS("PLAG_RFLU::TotalnPclsComm")
#endif

      IF ( global%nPclsCommTot == 0 ) THEN
        EXIT commLoop
      END IF ! nPclsSendTot

! ------------------------------------------------------------------------------
!     Initialize receive counters
! ------------------------------------------------------------------------------
 
#ifdef ROCPROF
      CALL FPROFILER_BEGINS("PLAG_RFLU::InitRecvCount")
#endif

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)
        CALL PLAG_RFLU_InitRecvCounters(pRegion)
      END DO ! iReg

#ifdef ROCPROF
      CALL FPROFILER_ENDS("PLAG_RFLU::InitRecvCount")
#endif

! ------------------------------------------------------------------------------
!     Send counters
! ------------------------------------------------------------------------------
 
#ifdef ROCPROF
      CALL FPROFILER_BEGINS("PLAG_RFLU::ISendCount")
#endif

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)
        CALL PLAG_RFLU_ISendCounters(pRegion)
      END DO ! iReg

#ifdef ROCPROF
      CALL FPROFILER_ENDS("PLAG_RFLU::ISendCount")
#endif

! ------------------------------------------------------------------------------
!     Copy counters on same processors
! ------------------------------------------------------------------------------
 
#ifdef ROCPROF
      CALL FPROFILER_BEGINS("PLAG_RFLU::CopyCount")
#endif

      CALL PLAG_RFLU_CopyCounters(regions)

#ifdef ROCPROF
      CALL FPROFILER_ENDS("PLAG_RFLU::CopyCount")
#endif

! ------------------------------------------------------------------------------
!     Create and fill send data buffers 
! ------------------------------------------------------------------------------

#ifdef ROCPROF
      CALL FPROFILER_BEGINS("PLAG_RFLU::CreateBuffSend")
#endif

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)
        CALL PLAG_RFLU_CreateBuffersSend(pRegion)
        CALL PLAG_RFLU_LoadBuffersSend(pRegion)
      END DO ! iReg

#ifdef ROCPROF
      CALL FPROFILER_ENDS("PLAG_RFLU::CreateBuffSend")
#endif

! ------------------------------------------------------------------------------
!     Receive counters 
!      Note: separate from clear requests due MPI issues on uP
! ------------------------------------------------------------------------------

#ifdef ROCPROF
      CALL FPROFILER_BEGINS("PLAG_RFLU::RecvCount")
#endif

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)
        CALL PLAG_RFLU_RecvCounters(pRegion)
      END DO ! iReg

#ifdef ROCPROF
      CALL FPROFILER_ENDS("PLAG_RFLU::RecvCount")
#endif

! ------------------------------------------------------------------------------
!     Clear requests and create receive data buffers 
! ------------------------------------------------------------------------------

#ifdef ROCPROF
      CALL FPROFILER_BEGINS("PLAG_RFLU::CreateBuffRecv")
#endif

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)
        CALL PLAG_RFLU_ClearRequestWrapper(pRegion,REQUEST_TYPE_COUNTER)
        CALL PLAG_RFLU_CreateBuffersRecv(pRegion)
      END DO ! iReg

#ifdef ROCPROF
      CALL FPROFILER_ENDS("PLAG_RFLU::CreateBuffRecv")
#endif

! ------------------------------------------------------------------------------
!    Send communication buffers
! ------------------------------------------------------------------------------

#ifdef ROCPROF
      CALL FPROFILER_BEGINS("PLAG_RFLU::ISendData")
#endif

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)
        CALL PLAG_RFLU_ISendData(pRegion)
      END DO ! iReg

#ifdef ROCPROF
      CALL FPROFILER_ENDS("PLAG_RFLU::ISendData")
#endif

! ------------------------------------------------------------------------------
!     Communicate data buffers for borders on same processor
! ------------------------------------------------------------------------------
 
#ifdef ROCPROF
      CALL FPROFILER_BEGINS("PLAG_RFLU::CopyData")
#endif

      CALL PLAG_RFLU_CopyData(regions)

#ifdef ROCPROF
      CALL FPROFILER_ENDS("PLAG_RFLU::CopyData")
#endif

! ------------------------------------------------------------------------------
!     Update datastructure after buffers have been sent
! ------------------------------------------------------------------------------

#ifdef ROCPROF
      CALL FPROFILER_BEGINS("PLAG_RFLU::UpdateDataPar")
#endif

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        IF ( pRegion%plag%nPcls > 0 ) THEN
          CALL PLAG_UpdateDataStruct(pRegion)
        END IF ! nPcls
      END DO ! iReg

#ifdef ROCPROF
      CALL FPROFILER_ENDS("PLAG_RFLU::UpdateDataPar")
#endif

! ------------------------------------------------------------------------------
!     Receive data buffers
! ------------------------------------------------------------------------------

#ifdef ROCPROF
      CALL FPROFILER_BEGINS("PLAG_RFLU::RecvData")
#endif

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)
        CALL PLAG_RFLU_RecvData(pRegion)
      END DO ! iReg

#ifdef ROCPROF
      CALL FPROFILER_ENDS("PLAG_RFLU::RecvData")
#endif

! ------------------------------------------------------------------------------
!     Clear requests for data buffers
! ------------------------------------------------------------------------------

#ifdef ROCPROF
      CALL FPROFILER_BEGINS("PLAG_RFLU::ClearReqsData")
#endif

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)
        CALL PLAG_RFLU_ClearRequestWrapper(pRegion,REQUEST_TYPE_DATA)
      END DO ! iReg

#ifdef ROCPROF
      CALL FPROFILER_ENDS("PLAG_RFLU::ClearReqsData")
#endif

! ------------------------------------------------------------------------------
!     Unload receive buffers, deallocate buffers, and reallocate memory
! ------------------------------------------------------------------------------

#ifdef ROCPROF
      CALL FPROFILER_BEGINS("PLAG_RFLU::UnloadBuffRecv")
#endif

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        pRegion%plag%nPclsPrev = pRegion%plag%nPcls

        ! Subbu - Reallocate memory before copying the nPclsRecv info
        ! into array ranging nPcls+1 < index < nPclsMax
        ! Reallocate memory if nPclsRecv + nPcls > nPclsMax in a region
        nPclsRecvPerReg = 0
        DO iBorder = 1,pRegion%grid%nBorders
          pBorder => pRegion%grid%borders(iBorder)
          nPclsRecvPerReg = nPclsRecvPerReg + pBorder%nPclsRecv
        END DO   
        IF (nPclsRecvPerReg .GT. (pRegion%plag%nPclsMax - pRegion%plag%nPcls)) THEN
          pRegion%plag%nPclsMax = NINT(PLAG_REALLOC_FACT* &
                                       REAL((pRegion%plag%nPcls+nPclsRecvPerReg),KIND=RFREAL))
          CALL PLAG_ReallocMemWrappernPclsPrev(pRegion) 
        END IF
        ! Subbu - End reallocate memory for nPclsRecv + nPcls > nPclsMax
        CALL PLAG_RFLU_UnloadBuffersRecv(pRegion)
        CALL PLAG_RFLU_DestroyBuffersSend(pRegion)
        CALL PLAG_RFLU_DestroyBuffersRecv(pRegion)
        CALL PLAG_ReallocMemWrapper(pRegion) 
      END DO ! iReg

#ifdef ROCPROF
      CALL FPROFILER_ENDS("PLAG_RFLU::UnloadBuffRecv")
#endif

! ------------------------------------------------------------------------------
!     Initialize send counters
! ------------------------------------------------------------------------------

#ifdef ROCPROF
      CALL FPROFILER_BEGINS("PLAG_RFLU::InitSendCountPar")
#endif

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)
        CALL PLAG_RFLU_InitSendCounters(pRegion)
      END DO ! iReg

#ifdef ROCPROF
      CALL FPROFILER_ENDS("PLAG_RFLU::InitSendCountPar")
#endif

! ------------------------------------------------------------------------------
!     Continue tracking particles if remaining trajectory distance is not zero
!       Tracking performed when particles are added to datastructure, i.e.
!       iPclEnd greater or equal to iPclBeg.
! ------------------------------------------------------------------------------

#ifdef ROCPROF
      CALL FPROFILER_BEGINS("PLAG_RFLU::FindCellsPar")
#endif

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)
        pGrid   => pRegion%grid      
 
        iPclBeg = pRegion%plag%nPclsPrev +1
        iPclEnd = pRegion%plag%nPcls

        IF ( iPclEnd < iPclBeg ) CYCLE

! TEMPORARY
        CALL PLAG_CalcDerivedVariables(pRegion)
! END TEMPORARY

        SELECT CASE ( pRegion%plagInput%findPclMethod )
          CASE ( FIND_PCL_METHOD_TRAJ_FAST ) 
             CALL PLAG_RFLU_FindCellsTrajFast(pRegion,iPclBeg,iPclEnd)      
          CASE ( FIND_PCL_METHOD_TRAJ_SAFE ) 
            CALL PLAG_RFLU_FindCellsTrajSafe(pRegion,iPclBeg,iPclEnd)  
          ! Subbu - Add hardcode pcl tracking call routine
          CASE ( FIND_PCL_METHOD_HARDCODE ) 
            CALL PLAG_RFLU_FindCellsHardCode(pRegion,iPclBeg,iPclEnd)  
          ! Subbu - End add hardcode pcl tracking call routine
         CASE DEFAULT 
           CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pRegion%plagInput%findPclMethod
      END DO ! iReg

#ifdef ROCPROF
      CALL FPROFILER_ENDS("PLAG_RFLU::FindCellsPar")
#endif

! ==============================================================================
!     Update loopCounter 
! ============================================================================== 

      loopCounter = loopCounter + 1

! ==============================================================================
!     Guard against infinite loop 
! ============================================================================== 

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter        
    END DO commLoop

! ******************************************************************************
!   Update datastructure after communication is completed
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("PLAG_RFLU::UpdateDataPar2")
#endif

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)

      IF ( pRegion%plag%nPcls > 0 ) THEN
        CALL PLAG_UpdateDataStruct(pRegion)
      END IF ! nPcls
    END DO ! iReg

#ifdef ROCPROF
    CALL FPROFILER_ENDS("PLAG_RFLU::UpdateDataPar2")
#endif

! ******************************************************************************
! Print number of particles in each region
! ******************************************************************************

!  DO iReg = 1,global%nRegionsLocal
!     pRegion => regions(iReg)
!     istage = pRegion%irkStep 
!
!!     IF ( istage == global%nrkSteps ) THEN 
!       doPrint = RFLU_DecidePrint(global)
!       IF ( (global%verbLevel > VERBOSE_NONE) .AND. (doPrint .EQV. .TRUE.)  ) & 
!         WRITE(STDOUT,'(A,I4,I8)') 'iRegGlobal: nPcls = ',pRegion%iRegionGlobal,pRegion%plag%nPcls
!!     END IF ! istage
!  END DO ! iReg 
  
#ifdef ROCPROF
    CALL FPROFILER_ENDS("PLAG_RFLU::PLAG_RFLU_CommDriver")
#endif

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE PLAG_RFLU_CommDriver









! ******************************************************************************
!
! Purpose: Copy counters.
!
! Description: None.
!
! Input:
!   regions     Pointer to all regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLU_CopyCounters(regions)

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

    CALL RegisterFunction(global,'PLAG_RFLU_CopyCounters',__FILE__)
 
! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid      

      DO iBorder = 1,pGrid%nBorders
        pBorder => pGrid%borders(iBorder)

        IF ( pBorder%iProc == global%myProcid ) THEN
          iReg2    = pBorder%iRegionLocal
          iBorder2 = pBorder%iBorder

          pRegion2 => regions(iReg2)     
          pBorder2 => pRegion2%grid%borders(iBorder2) 

          pBorder2%nPclsRecv = pBorder%nPclsSend

        END IF ! pBorder
      END DO ! iBorder
    END DO ! iReg

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE PLAG_RFLU_CopyCounters









! ******************************************************************************
!
! Purpose: Copy data buffers for regions on same processor.
!
! Description: None.
!
! Input:
!   regions    Pointer to all regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLU_CopyData(regions)

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

    CALL RegisterFunction(global,'PLAG_RFLU_CopyData',__FILE__)

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

          IF ( pBorder%nPclsSend /= pBorder2%nPclsRecv ) THEN 
            CALL ErrorStop(global,ERR_BUFFERDIM_MISMATCH,__LINE__)
          END IF ! pBorder%nPclsSend

          IF ( pBorder%nPclsSend == 0 ) CYCLE

! ------------------------------------------------------------------------------
!         Real data
! ------------------------------------------------------------------------------      

          CALL PLAG_RFLU_CopyDataReal(global,pBorder%plag%sendBuff, & 
                                      pBorder2%plag%recvBuff)      
                                     
! ------------------------------------------------------------------------------
!         Integer data
! ------------------------------------------------------------------------------      

          CALL PLAG_RFLU_CopyDataInt(global,pBorder%plag%sendBuffInt, & 
                                     pBorder2%plag%recvBuffInt)      
                 
        END IF ! pBorder
      END DO ! iBorder
    END DO ! iReg

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE PLAG_RFLU_CopyData








! ******************************************************************************
!
! Purpose: Copy particle data for integer variables.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pclData     Integer array with particle data
!
! Output:
!   pclData2    Integer array with particle data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLU_CopyDataInt(global,pclData,pclData2)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, DIMENSION(:,:), INTENT(IN) :: pclData
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: pclData2    
    TYPE(t_global), POINTER :: global         
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iPcl,iVar,nPclsSend,nVars        
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'PLAG_RFLU_CopyDataInt',__FILE__)

! ******************************************************************************
!   Set variables
! ******************************************************************************

    nVars = SIZE(pclData,1)
    
    IF ( nVars /= SIZE(pclData2,1) ) THEN 
      CALL ErrorStop(global,ERR_DATADIM_MISMATCH,__LINE__)
    END IF ! nVars 

    nPclsSend = SIZE(pclData,2)
    
    IF ( nPclsSend /= SIZE(pclData2,2) ) THEN 
      CALL ErrorStop(global,ERR_DATADIM_MISMATCH,__LINE__)
    END IF ! nVars 
        
! ******************************************************************************
!   Copy data 
! ******************************************************************************

    DO iPcl = 1,nPclsSend
      DO iVar = 1,nVars
        pclData2(iVar,iPcl) = pclData(iVar,iPcl) 
      END DO ! iVar      
    END DO ! iPcl
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE PLAG_RFLU_CopyDataInt









! ******************************************************************************
!
! Purpose: Copy particle data for real variables.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pclData     Real array with particle data
!
! Output: 
!   pclData2    Real array with particle data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLU_CopyDataReal(global,pclData,pclData2)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: pclData
    REAL(RFREAL), DIMENSION(:,:), INTENT(OUT) :: pclData2    
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global         
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iPcl,iVar,nPclsSend,nVars        
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'PLAG_RFLU_CopyDataReal',__FILE__)

! ******************************************************************************
!   Set variables
! ******************************************************************************

    nVars = SIZE(pclData,1)
    
    IF ( nVars /= SIZE(pclData2,1) ) THEN 
      CALL ErrorStop(global,ERR_DATADIM_MISMATCH,__LINE__)
    END IF ! nVars 

    nPclsSend = SIZE(pclData,2)
    
    IF ( nPclsSend /= SIZE(pclData2,2) ) THEN 
      CALL ErrorStop(global,ERR_DATADIM_MISMATCH,__LINE__)
    END IF ! nVars 
        
! ******************************************************************************
!   Copy data 
! ******************************************************************************

    DO iPcl = 1,nPclsSend
      DO iVar = 1,nVars
        pclData2(iVar,iPcl) = pclData(iVar,iPcl) 
      END DO ! iVar      
    END DO ! iPcl
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE PLAG_RFLU_CopyDataReal











! ******************************************************************************
!
! Purpose: Create receive buffers.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. aiv, arv cv, cvOld and rhsSum arrays are communicated.
!   2. cvOld array needs to be communicated for proper trajectory tracking.
!   3. TEMPORARY: Manoj: Adding dudtMixt, dudtPlag
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLU_CreateBuffersRecv(pRegion)

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

    INTEGER :: errorFlag,iBorder,nUnsteadyData,nVals,nVarsInt,nVarsReal
    TYPE(t_border), POINTER :: pBorder        
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid 
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_RFLU_CreateBuffersRecv',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Set dimensions
!    Add for the integer variable an additional value 
!    to account for cell mapping
! ******************************************************************************

    nUnsteadyData = pRegion%plagInput%nUnsteadyData

    nVarsReal = 3*pRegion%plag%nCv  +pRegion%plag%nArv  +2*3*nUnsteadyData
    nVarsInt  = pRegion%plag%nAiv +1    

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      nVals = pBorder%nPclsRecv

      IF ( pBorder%nPclsRecv == 0 ) CYCLE

! ==============================================================================
!     Allocate memory
! ==============================================================================    
      
      ALLOCATE(pBorder%plag%recvBuff(nVarsReal,nVals),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%plag%recvBuff')
      END IF ! global%error

      ALLOCATE(pBorder%plag%recvBuffInt(nVarsInt,nVals),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%plag%recvBuffInt')
      END IF ! global%error

    END DO ! iBorder
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_RFLU_CreateBuffersRecv









! ******************************************************************************
!
! Purpose: Create send buffers.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. aiv, arv cv, cvOld and rhsSum arrays are communicated.
!   2. cvOld array needs to be communicated for proper trajectory tracking.
!   3. TEMPORARY: Manoj: Adding dudtMixt, dudtPlag
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLU_CreateBuffersSend(pRegion)

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

    INTEGER :: errorFlag,iBorder,nUnsteadyData,nVals,nVarsInt,nVarsReal
    TYPE(t_border), POINTER :: pBorder        
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid 
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_RFLU_CreateBuffersSend',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Set dimensions
!    Note: Add for the integer variable an additional value 
!    to account for cell mapping
! ******************************************************************************

    nUnsteadyData = pRegion%plagInput%nUnsteadyData

    nVarsReal = 3*pRegion%plag%nCv  +pRegion%plag%nArv   +2*3*nUnsteadyData 
    nVarsInt  = pRegion%plag%nAiv +1    

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      nVals = pBorder%nPclsSend

      IF ( pBorder%nPclsSend == 0 ) CYCLE

! ==============================================================================
!     Allocate memory
! ==============================================================================    

      ALLOCATE(pBorder%plag%sendBuff(nVarsReal,nVals),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%plag%sendBuff')
      END IF ! global%error
      
      ALLOCATE(pBorder%plag%sendBuffInt(nVarsInt,nVals),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%plag%sendBuffInt')
      END IF ! global%error

    END DO ! iBorder
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
      
  END SUBROUTINE PLAG_RFLU_CreateBuffersSend









! ******************************************************************************
!
! Purpose: Destroy receive buffers.
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

  SUBROUTINE PLAG_RFLU_DestroyBuffersRecv(pRegion)

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

    CALL RegisterFunction(global,'PLAG_RFLU_DestroyBuffersRecv',__FILE__)

!    IF ( global%myProcid == MASTERPROC .AND. &
!         global%verbLevel > VERBOSE_NONE ) THEN
!      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying receive buffers...'
!      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
!                                       pRegion%iRegionGlobal
!    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      IF ( pBorder%nPclsRecv == 0 ) CYCLE

! ==============================================================================
!     Deallocate memory
! ==============================================================================    
      
      DEALLOCATE(pBorder%plag%recvBuff,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pBorder%plag%recvBuff')
      END IF ! global%error
      
      DEALLOCATE(pBorder%plag%recvBuffInt,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pBorder%plag%recvBuffInt')
      END IF ! global%error

    END DO ! iBorder
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

!    IF ( global%myProcid == MASTERPROC .AND. &
!         global%verbLevel > VERBOSE_NONE ) THEN
!      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying receive buffers done.'
!    END IF ! global%verbLevel
      
  END SUBROUTINE PLAG_RFLU_DestroyBuffersRecv









! ******************************************************************************
!
! Purpose: Destroy send buffers.
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

  SUBROUTINE PLAG_RFLU_DestroyBuffersSend(pRegion)

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

    CALL RegisterFunction(global,'PLAG_RFLU_DestroyBuffersSend',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      IF ( pBorder%nPclsSend == 0 ) CYCLE

! ==============================================================================
!     Deallocate memory
! ==============================================================================    

      DEALLOCATE(pBorder%plag%sendBuff,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pBorder%plag%sendBuff')
      END IF ! global%error
      
      DEALLOCATE(pBorder%plag%sendBuffInt,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pBorder%plag%sendBuffInt')
      END IF ! global%error

    END DO ! iBorder
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
   
  END SUBROUTINE PLAG_RFLU_DestroyBuffersSend









! ******************************************************************************
!
! Purpose: Initialize receive counter.
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

  SUBROUTINE PLAG_RFLU_InitRecvCounters(pRegion)

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
    TYPE(t_grid),   POINTER :: pGrid     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_RFLU_InitRecvCounters',__FILE__)
 
! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid   => pRegion%grid      
 
! ******************************************************************************
!   Initialize counter
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)
      pBorder%nPclsRecv = 0
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE PLAG_RFLU_InitRecvCounters










! ******************************************************************************
!
! Purpose: Initialize send counter.
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

  SUBROUTINE PLAG_RFLU_InitSendCounters(pRegion)

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
    TYPE(t_grid),   POINTER :: pGrid     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_RFLU_InitSendCounters',__FILE__)
 
! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid   => pRegion%grid      
 
! ******************************************************************************
!   Initialize counter
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)
      pBorder%nPclsSend = 0
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE PLAG_RFLU_InitSendCounters









! ******************************************************************************
!
! Purpose: Send counters.
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

  SUBROUTINE PLAG_RFLU_ISendCounters(pRegion)

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

    INTEGER :: errorFlag,iBorder,nVars,tag        

    TYPE(t_border), POINTER :: pBorder 
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_RFLU_ISendCounters',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

    nVars = 1

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)
    
! ==============================================================================
!     Send counter if not on same processor
! ==============================================================================    
    
      IF ( pBorder%iProc /= global%myProcid ) THEN     
        tag = pBorder%plag%tagCount
        
        CALL MPI_ISend(pBorder%nPclsSend,nVars,MPI_INTEGER,pBorder%iProc,tag, & 
                       global%mpiComm,pBorder%plag%sendRequestCount,errorFlag)
        
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
        END IF ! global%error                                      
      END IF ! pBorder
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE PLAG_RFLU_ISendCounters









! ******************************************************************************
!
! Purpose: Send data buffers.
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

  SUBROUTINE PLAG_RFLU_ISendData(pRegion)

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

    INTEGER :: errorFlag,iBorder,nBuffsInt,nBuffsReal,nVals,nVarsInt, &
               nVarsReal,tagInt,tagReal        

    TYPE(t_border), POINTER :: pBorder 
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_RFLU_ISendData',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      IF ( pBorder%nPclsSend == 0 ) CYCLE
  
      nVarsReal = SIZE(pBorder%plag%sendBuff,1)
      nVarsInt  = SIZE(pBorder%plag%sendBuffInt,1) 
      nVals     = SIZE(pBorder%plag%sendBuff,2)

      nBuffsReal = nVarsReal *nVals
      nBuffsInt  = nVarsInt  *nVals 

      IF ( nVals /= pBorder%nPclsSend ) THEN 
        CALL ErrorStop(global,ERR_DATADIM_MISMATCH,__LINE__)
      END IF ! nVarsSend

! ==============================================================================
!     Send data if not on same processor and none null size
! ==============================================================================    
    
      IF ( pBorder%iProc /= global%myProcid ) THEN     
      
! ------------------------------------------------------------------------------
!       Integer data buffers
! ------------------------------------------------------------------------------    

        tagInt  = pBorder%plag%tagInt

        CALL MPI_ISend(pBorder%plag%sendBuffInt,nBuffsInt,MPI_INTEGER, &
                       pBorder%iProc,tagInt,global%mpiComm,             &
                       pBorder%plag%sendRequestInt,errorFlag            )        
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
        END IF ! global%error                                      

! ------------------------------------------------------------------------------
!       Real data buffers
! ------------------------------------------------------------------------------    

        tagReal = pBorder%plag%tag
       
        CALL MPI_ISend(pBorder%plag%sendBuff,nBuffsReal,MPI_RFREAL, &
                       pBorder%iProc,tagReal,global%mpiComm,        &
                       pBorder%plag%sendRequest,errorFlag           )        
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
        END IF ! global%error 

      END IF ! pBorder
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE PLAG_RFLU_ISendData










! ******************************************************************************
!
! Purpose: Load particle data into communication buffers 
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

  SUBROUTINE PLAG_RFLU_LoadBuffersSend(pRegion)

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

    INTEGER :: errorFlag,iBorder,icg,iLoc,iPcl,iPcl2,iVar,iVarBuff,&
               nAiv,nArv,nCv,nUnsteadyData
 
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global         
    TYPE(t_grid),   POINTER :: pGrid
    TYPE(t_plag),   POINTER :: pPlag
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_RFLU_LoadBuffersSend',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid
    pPlag => pRegion%plag

! ******************************************************************************
!   Set dimensions
!    Add for the integer variable an additional value 
!    to account for cell mapping
! ******************************************************************************

    nAiv = pRegion%plag%nAiv
    nArv = pRegion%plag%nArv 
    nCv  = pRegion%plag%nCv

    nUnsteadyData = pRegion%plagInput%nUnsteadyData

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      IF ( pBorder%nPclsSend == 0 ) CYCLE

! ==============================================================================
!     Load send buffers from particle datastructure
! ==============================================================================

      DO iPcl = 1,pBorder%nPclsSend

! ------------------------------------------------------------------------------
!       Get particle and cell mappings
! ------------------------------------------------------------------------------

        iPcl2 = pBorder%iPclSend(1,iPcl)

        IF ( pPlag%aiv(AIV_PLAG_STATUS,iPcl2) /= PLAG_STATUS_COMM ) THEN
         WRITE(*,*) ' PLAG_RFLU_LoadBuffersSend: PLAG_STATUS Mismatch'
         WRITE(*,*) '  iPcl2     = ',iPcl2
         WRITE(*,*) '  aivStatus = ' , pPlag%aiv(AIV_PLAG_STATUS,iPcl2)
!TO DO   CALL ErrorStop(global,ERR_STATUS_MISMATCH,__LINE__)
        ENDIF ! aiv
        
        icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcl2)

        IF ( icg < pGrid%nCells+1 ) THEN
          WRITE(*,*) ' PLAG_RFLU_LoadBuffersSend: Cell Bound Mismatch' 
          WRITE(*,*) '  icg      = ',icg
          WRITE(*,*) '  nCells+1 = ' ,  pGrid%nCells+1
          STOP  
!TO DO   CALL ErrorStop(global,ERR_CELLBOUND_MISMATCH,__LINE__)
        ENDIF ! icg
        
! ------------------------------------------------------------------------------
!       Load array location
! ------------------------------------------------------------------------------
        
        iLoc = pBorder%iPclSend(2,iPcl)

        IF ( iLoc > pBorder%nCellsRecv ) THEN
          WRITE(*,*) ' PLAG_RFLU_LoadBuffersSend: Cell Bound Mismatch on Recv Side' 
          WRITE(*,*) '  iLoc       = ',iLoc
          WRITE(*,*) '  nCellsRecv = ' ,  pBorder%nCellsRecv 
          STOP
!TO DO   CALL ErrorStop(global,ERR_CELLBOUND_MISMATCH,__LINE__)
        ENDIF ! iLoc
        
! ------------------------------------------------------------------------------
!       Load communication buffers for real data
! ------------------------------------------------------------------------------

        iVarBuff = 0
        
        DO iVar = 1,nCv
          iVarBuff = iVarBuff+1
          pBorder%plag%sendBuff(iVarBuff,iPcl) = pPlag%cv(iVar,iPcl2)
        END DO ! iVar

        DO iVar = 1,nCv
          iVarBuff = iVarBuff+1
          pBorder%plag%sendBuff(iVarBuff,iPcl) = pPlag%cvOld(iVar,iPcl2)
        END DO ! iVar   

        DO iVar = 1,nCv
          iVarBuff = iVarBuff+1
          pBorder%plag%sendBuff(iVarBuff,iPcl) = pPlag%rhsSum(iVar,iPcl2)
        END DO ! iVar 
        
        DO iVar = 1,nArv
          iVarBuff = iVarBuff+1
          pBorder%plag%sendBuff(iVarBuff,iPcl) = pPlag%arv(iVar,iPcl2)
        END DO ! iVar

        DO iVar = 1,nUnsteadyData
          iVarBuff = iVarBuff+1
          pBorder%plag%sendBuff(iVarBuff,iPcl) = pPlag%dudtMixt(XCOORD,iVar,iPcl2)
        END DO ! iVar

        DO iVar = 1,nUnsteadyData
          iVarBuff = iVarBuff+1
          pBorder%plag%sendBuff(iVarBuff,iPcl) = pPlag%dudtMixt(YCOORD,iVar,iPcl2)
        END DO ! iVar

        DO iVar = 1,nUnsteadyData
          iVarBuff = iVarBuff+1
          pBorder%plag%sendBuff(iVarBuff,iPcl) = pPlag%dudtMixt(ZCOORD,iVar,iPcl2)
        END DO ! iVar

        DO iVar = 1,nUnsteadyData
          iVarBuff = iVarBuff+1
          pBorder%plag%sendBuff(iVarBuff,iPcl) = pPlag%dudtPlag(XCOORD,iVar,iPcl2)
        END DO ! iVar

        DO iVar = 1,nUnsteadyData
          iVarBuff = iVarBuff+1
          pBorder%plag%sendBuff(iVarBuff,iPcl) = pPlag%dudtPlag(YCOORD,iVar,iPcl2)
        END DO ! iVar

        DO iVar = 1,nUnsteadyData
          iVarBuff = iVarBuff+1
          pBorder%plag%sendBuff(iVarBuff,iPcl) = pPlag%dudtPlag(ZCOORD,iVar,iPcl2)
        END DO ! iVar

! ------------------------------------------------------------------------------
!       Load communication buffers for integer data
! ------------------------------------------------------------------------------

        DO iVar = 1,nAiv
          pBorder%plag%sendBuffInt(iVar,iPcl)  = pPlag%aiv(iVar,iPcl2)
        END DO ! iVar
        
        pBorder%plag%sendBuffInt(nAiv+1,iPcl)  = iLoc

! ------------------------------------------------------------------------------
!       Overwrite status arrays
! ------------------------------------------------------------------------------
        
        pBorder%plag%sendBuffInt(AIV_PLAG_STATUS,iPcl ) = PLAG_STATUS_KEEP
                       pPlag%aiv(AIV_PLAG_STATUS,iPcl2) = PLAG_STATUS_DELETE     

      END DO ! iPcl
     END DO ! iBorder        

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE PLAG_RFLU_LoadBuffersSend











! ******************************************************************************
!
! Purpose: Receive counters.
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

  SUBROUTINE PLAG_RFLU_RecvCounters(pRegion)

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

    INTEGER :: errorFlag,iBorder,nVars,tag
    INTEGER :: status(MPI_STATUS_SIZE)

    TYPE(t_border), POINTER :: pBorder 
    TYPE(t_global), POINTER :: global 
    TYPE(t_grid), POINTER :: pGrid     
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_RFLU_RecvCounters',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid
    
    nVars = 1

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

! ==============================================================================
!     Receive data if not on same process
! ==============================================================================    

      IF ( pBorder%iProc /= global%myProcid ) THEN 
        tag = pBorder%plag%tagCount

        CALL MPI_Recv(pBorder%nPclsRecv,nVars,MPI_INTEGER,pBorder%iProc,tag, & 
                      global%mpiComm,status,errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
        END IF ! global%error 
      
      END IF ! pBorder
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE PLAG_RFLU_RecvCounters






! ******************************************************************************
!
! Purpose: Receive data buffers.
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

  SUBROUTINE PLAG_RFLU_RecvData(pRegion)

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

    INTEGER :: errorFlag,iBorder,nBuffsInt,nBuffsReal,nVals,nVarsInt, &
               nVarsReal,tagInt,tagReal
    INTEGER :: status(MPI_STATUS_SIZE),statusInt(MPI_STATUS_SIZE)

    TYPE(t_border), POINTER :: pBorder 
    TYPE(t_global), POINTER :: global 
    TYPE(t_grid), POINTER :: pGrid     
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_RFLU_RecvData',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      IF ( pBorder%nPclsRecv == 0 ) CYCLE
  
      nVarsReal = SIZE(pBorder%plag%recvBuff,1)
      nVarsInt  = SIZE(pBorder%plag%recvBuffInt,1)
      nVals     = SIZE(pBorder%plag%recvBuff,2)
      
      nBuffsReal = nVarsReal *nVals
      nBuffsInt  = nVarsInt  *nVals

      IF ( nVals /= pBorder%nPclsRecv ) THEN 
        CALL ErrorStop(global,ERR_DATADIM_MISMATCH,__LINE__)
      END IF ! nVarsSend

! ==============================================================================
!     Receive data if not on same process
! ==============================================================================    

      IF ( pBorder%iProc /= global%myProcid ) THEN 
      
! ------------------------------------------------------------------------------
!       Integer data buffers
! ------------------------------------------------------------------------------    

        tagInt = pBorder%plag%tagInt

        CALL MPI_Recv(pBorder%plag%recvBuffInt,nBuffsInt,MPI_INTEGER,      &
                      pBorder%iProc,tagInt,global%mpiComm,statusInt,errorFlag )
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
        END IF ! global%error 
      
! ------------------------------------------------------------------------------
!       Real data buffers
! ------------------------------------------------------------------------------    

        tagReal = pBorder%plag%tag

        CALL MPI_Recv(pBorder%plag%recvBuff,nBuffsReal,MPI_RFREAL,          &
                      pBorder%iProc,tagReal,global%mpiComm,status,errorFlag )
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
        END IF ! global%error       

      END IF ! pBorder
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE PLAG_RFLU_RecvData







! ******************************************************************************
!
! Purpose: Calculate the total number of particles to be communicated for all 
!   regions on all processors.
!
! Description: First determine the total number of particles for all regions on 
!   the same processor. 
!
! Input: 
!   regions           Data of all regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLU_TotalnPclsComm(regions)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: regions(:)
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iBorder,iReg,nPclsCommGlobal,nPclsCommLocal
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global 
    TYPE(t_grid), POINTER :: pGrid     
    TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
!   Start
! ******************************************************************************

    global => regions(0)%global

    CALL RegisterFunction(global,'PLAG_RFLU_TotalnPclsComm',__FILE__)

! ******************************************************************************
!   Initialize variables
! ******************************************************************************

    nPclsCommLocal = 0

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pRegion%global%nPclsCommTot = 0
    END DO ! iReg

! ******************************************************************************
!   Compute sum of communicated particles for all regions on the same processor
! ******************************************************************************

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid 
      
      DO iBorder = 1,pGrid%nBorders
        pBorder => pGrid%borders(iBorder)
        
        nPclsCommLocal = nPclsCommLocal + pBorder%nPclsSend
      END DO ! iBorder 
    END DO ! iReg

! ******************************************************************************
!   Determine global sum of particles communicated over all processors and store
! ******************************************************************************

    CALL MPI_Allreduce(nPclsCommLocal,nPclsCommGlobal,1,MPI_INTEGER,MPI_SUM, & 
                       global%mpiComm,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
    END IF ! global%errorFlag
    
    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)

      pRegion%global%nPclsCommTot = nPclsCommGlobal 
      global%nPclsCommTot         = nPclsCommGlobal 
    END DO ! iReg

! ******************************************************************************
!   Print information
! ******************************************************************************

! TEMPORARY
!    IF ( global%myProcid == MASTERPROC   .AND. & 
!         global%verbLevel > VERBOSE_LOW  .AND. &
!         global%nPclsCommTot /= 0              ) THEN
!      WRITE(STDOUT,'(A,1X,A)') &
!        SOLVER_NAME,'Printing total communicated particles information...'
!      WRITE(STDOUT,'(A,3X,A,1X,I8)') &
!         SOLVER_NAME,'Total communicated particles:',global%nPclsCommTot
!      WRITE(STDOUT,'(A,1X,A)') &
!        SOLVER_NAME,'Printing total communicated particles information done.'
!     END IF ! global%myProcid
! END TEMPORARY

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_RFLU_TotalnPclsComm







! ******************************************************************************
!
! Purpose: Unload communication buffers into particle datastructure.
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

  SUBROUTINE PLAG_RFLU_UnloadBuffersRecv(pRegion)

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

    INTEGER :: errorFlag,iBorder,icg,iLoc,iPcl,iPcl2,iVar,iVarBuff,nAiv,nArv,&
               nCv,nPclsSend,nUnsteadyData
  
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global         
    TYPE(t_grid),   POINTER :: pGrid
    TYPE(t_plag),   POINTER :: pPlag
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_RFLU_UnloadBuffersRecv',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid
    pPlag => pRegion%plag

! ******************************************************************************
!   Set dimensions
!    Add for the integer variable an additional value 
!    to account for cell mapping
! ******************************************************************************

    nAiv = pRegion%plag%nAiv
    nArv = pRegion%plag%nArv 
    nCv  = pRegion%plag%nCv

    nUnsteadyData = pRegion%plagInput%nUnsteadyData
    
    iPcl2 = pPlag%nPcls

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      IF ( pBorder%nPclsRecv == 0 ) CYCLE

! ==============================================================================
!     Append particle datastructure with receive buffers
! ==============================================================================

      DO iPcl = 1,pBorder%nPclsRecv
        iPcl2 = iPcl2 +1

! ------------------------------------------------------------------------------
!       Load array location and get mapping
! ------------------------------------------------------------------------------
        
        iLoc = pBorder%plag%recvBuffInt(nAiv+1,iPcl)
        
        IF ( iLoc > pBorder%nCellsSend ) THEN
          WRITE(*,*) ' PLAG_RFLU_UnloadBuffersRecv: Cell Bound Mismatch on Send Side' 
          WRITE(*,*) '  iLoc       = ',iLoc
          WRITE(*,*) '  nCellsSend = ' ,  pBorder%nCellsSend
          STOP
!TO DO   CALL ErrorStop(global,ERR_CELLBOUND_MISMATCH,__LINE__)
        ENDIF ! iLoc

        icg  = pBorder%icgSend(iLoc)

! ------------------------------------------------------------------------------
!       Load data structure from communication buffers for real data
! ------------------------------------------------------------------------------

        iVarBuff = 0
       
        ! Subbu - Debug
        IF (iPcl2 .GT. SIZE(pPlag%cv,2)) THEN 
          WRITE(*,*) '********* pPlag%cv is out of bounds **********'
          CALL ErrorStop(global,ERR_EXCEEDS_DECL_MEM,__LINE__)
        END IF
        ! Subbu - End Debug

        DO iVar = 1,nCv
          iVarBuff = iVarBuff+1
          pPlag%cv(iVar,iPcl2)     = pBorder%plag%recvBuff(iVarBuff,iPcl)
        END DO ! iVar
        
        DO iVar = 1,nCv
          iVarBuff = iVarBuff+1
          pPlag%cvOld(iVar,iPcl2)  = pBorder%plag%recvBuff(iVarBuff,iPcl)
        END DO ! iVar
        
        DO iVar = 1,nCv
          iVarBuff = iVarBuff+1
          pPlag%rhsSum(iVar,iPcl2) = pBorder%plag%recvBuff(iVarBuff,iPcl)
        END DO ! iVar 
        
        DO iVar = 1,nArv
          iVarBuff = iVarBuff+1
          pPlag%arv(iVar,iPcl2)    = pBorder%plag%recvBuff(iVarBuff,iPcl)
        END DO ! iVar

        DO iVar = 1,nUnsteadyData
          iVarBuff = iVarBuff+1
          pPlag%dudtMixt(XCOORD,iVar,iPcl2) = pBorder%plag%recvBuff(iVarBuff,iPcl)
        END DO ! iVar

        DO iVar = 1,nUnsteadyData
          iVarBuff = iVarBuff+1
          pPlag%dudtMixt(YCOORD,iVar,iPcl2) = pBorder%plag%recvBuff(iVarBuff,iPcl)
        END DO ! iVar

        DO iVar = 1,nUnsteadyData
          iVarBuff = iVarBuff+1
          pPlag%dudtMixt(ZCOORD,iVar,iPcl2) = pBorder%plag%recvBuff(iVarBuff,iPcl)
        END DO ! iVar

        DO iVar = 1,nUnsteadyData
          iVarBuff = iVarBuff+1
          pPlag%dudtPlag(XCOORD,iVar,iPcl2) = pBorder%plag%recvBuff(iVarBuff,iPcl)
        END DO ! iVar

        DO iVar = 1,nUnsteadyData
          iVarBuff = iVarBuff+1
          pPlag%dudtPlag(YCOORD,iVar,iPcl2) = pBorder%plag%recvBuff(iVarBuff,iPcl)
        END DO ! iVar

        DO iVar = 1,nUnsteadyData
          iVarBuff = iVarBuff+1
          pPlag%dudtPlag(ZCOORD,iVar,iPcl2) = pBorder%plag%recvBuff(iVarBuff,iPcl)
        END DO ! iVar

! ------------------------------------------------------------------------------
!       Load data structure from communication buffers for integer data
!        aivOld is needed for proper trajectory tracking
! ------------------------------------------------------------------------------

        DO iVar = 1,nAiv
          pPlag%aiv(iVar,iPcl2)    = pBorder%plag%recvBuffInt(iVar,iPcl)
        END DO ! iVar
        
        pPlag%aiv(AIV_PLAG_ICELLS,iPcl2)    = icg
        pPlag%aivOld(AIV_PLAG_ICELLS,iPcl2) = icg
      END DO ! iPcl
     END DO ! iBorder   

! ******************************************************************************
!   Update particle size
! ******************************************************************************      
     pRegion%plag%nPcls = iPcl2

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE PLAG_RFLU_UnloadBuffersRecv








! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE PLAG_RFLU_ModComm


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_ModComm.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:52  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/06/01 13:16:57  haselbac
! Major bug fix: Region 0 should not have been included in loops
!
! Revision 1.1  2007/04/09 18:50:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:34  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.9  2006/08/18 21:10:34  fnajjar
! Added PROF calls, enabled serial periodicity, cosmetics
!
! Revision 1.8  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.7  2005/12/30 16:29:02  fnajjar
! Cosmetic cleanups, streamlined nPclsTotComm routine and added call to CalcDv 
! for proper datastruc defs
!
! Revision 1.6  2005/12/19 16:46:31  fnajjar
! Bug fix for proper calling of PLAG_UpdateDataStruc
!
! Revision 1.5  2005/12/13 23:10:47  fnajjar
! Major bugs fixes for clear requests, tags and computing communicated particles, 
! and removed debug statements
!
! Revision 1.4  2005/11/10 22:24:27  fnajjar
! Removed DEBUG STOP
!
! Revision 1.3  2005/07/18 20:48:54  fnajjar
! Added MPI routines and started initial testing
!
! Revision 1.2  2005/05/27 00:57:43  haselbac
! Fixed bug in clearing requests
!
! Revision 1.1  2005/05/18 22:27:45  fnajjar
! Initial revision
!
! ******************************************************************************

