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
! Purpose: Integrate the governing equations in time; move/regenerate the grid.
!
! Description: None.
!
! Input: 
!   dTimeSystem         Total solution time (unsteady flow)
!   dIterSystem         Total number of iterations (steady flow)
!   regions             Data for all grid regions
!
! Output: 
!   regions             Data for all grid regions
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_TimeStepping.F90,v 1.8 2016/02/08 17:00:56 fred Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_TimeStepping(dTimeSystem,dIterSystem,regions)

  USE ModDataTypes
  USE ModError  
  USE ModParameters
  USE ModGlobal, ONLY: t_global  
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input    
  USE ModMPI
  
  USE RFLU_ModDimensions
  USE RFLU_ModForcesMoments
  USE RFLU_ModGeometry
  USE RFLU_ModHouMahesh, ONLY: RFLU_HM_PredCorrMP
  USE RFLU_ModMovingFrame, ONLY: RFLU_MVF_WritePatchVelAccel
  USE RFLU_ModPatchCoeffs
  USE RFLU_ModProbes
  USE RFLU_ModReadWriteAuxVars
  USE RFLU_ModReadWriteFlow
  USE RFLU_ModReadWriteGrid
  USE RFLU_ModReadWriteGridSpeeds
  USE RFLU_ModThrustSpecImpulse
  USE RFLU_ModWeights
  USE RFLU_ModBoundXvUtils, ONLY: RFLU_BXV_WriteVarsWrapper
  

  USE RFLU_ModTimeZoom, ONLY: RFLU_UnZoomGridSpeeds, &
                              RFLU_ZoomGridSpeeds

#ifdef PLAG
  USE PLAG_ModDimensions, ONLY: PLAG_CalcNPclsGlobal, &
                                PLAG_PrintNPclsGlobal, &
                                PLAG_RFLU_WriteDimensions
  USE PLAG_ModSurfStats, ONLY: PLAG_DecideHaveSurfStats, & 
                               PLAG_WriteSurfStatsWrapper
#endif  
  
  USE ModInterfaces, ONLY: DummyRungeKuttaMP, &
                           IntegrateSourceTermsMP, &
                           RFLU_ComputeGridSpeeds, &
                           RFLU_ComputeIntegralValues, & 
                           RFLU_CheckGridSpeeds, &
                           RFLU_DecidePrint, &
                           RFLU_DecideNeedBGradFace, & 
                           RFLU_DecideWrite, &
                           !begin BBR 
                           RFLU_DecideSmallWrite, &
                           !end BBR 
                           RFLU_ExplicitMultiStage, &
                           RFLU_GetDeformationWrapper, & 
                           RFLU_MinimumTimeStep, &
                           RFLU_MoveGridWrapper, &
                           RFLU_PrintChangeInfo, &  
                           RFLU_PrintGridInfo, &
                           RFLU_PrintFlowInfoWrapper, &
                           RFLU_PrintWriteConvergence, &
                           RFLU_ResidualNorm, &
                           RFLU_TimeStepInviscid, &
                           RFLU_TimeStepViscous, &
                           !begin BBR
                           RFLU_WriteInteg, &
                           RFLU_WritePM, &
                           !end BBR
                           RFLU_WriteRestartInfo, &
                           RFLU_WriteStatsFileOLES, &
                           RungeKuttaMP, &
                           WriteProbe, &
                           WriteTotalMass

#ifdef GENX
  USE ModInterfaces, ONLY: RFLU_PutBoundaryValues
#endif
#ifdef STATS
  USE ModInterfaces, ONLY: RFLU_WriteStat
  USE ModStatsRoutines, ONLY: GetStatistics
#endif

  IMPLICIT NONE

#ifdef GENX
#include "roccomf90.h"
#endif

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER :: dIterSystem
  REAL(RFREAL) :: dTimeSystem
  TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iPatch,iReg,stophandle,dtlimcount
  LOGICAL :: doPrint,doProbe,doFomWrite,doWrite,finished,ftermNew,moveGrid, &
             residFterm,stopFileExists
  TYPE(t_global), POINTER :: global
  TYPE(t_mixt_input), POINTER :: pMixtInput   
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion,pRegionSerial

! BBR - begin
  LOGICAL :: doSmallWrite,first,donesafewrite
  INTEGER :: t1,t2,clock_rate,clock_max,safew_l,safew
  INTEGER :: timeroutineat,errorFlag
  REAL(RFREAL) :: tc1,tc2,iter,cumulreal,cumulcpu,elapsedtime
  CHARACTER(CHRLEN) :: NameSolDir
! BBR - end

!CRN - begin
  REAL(RFREAL) :: timerStart, timerEnd
!CRN - end

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_TimeStepping.F90,v $ $Revision: 1.8 $'

! ******************************************************************************
! Set pointers and variables
! ****************************************************************************** 
  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_TimeStepping',__FILE__)

  finished = .FALSE. ! run not finished yet

  ! BBR - begin
  donesafewrite = .FALSE.
  safew_l = 0; safew = 0;
  ! BBR - end

#ifdef GENX
  global%timeSinceRestart = 0.0_RFREAL
#endif

! ==============================================================================
! No multigrid here
! ==============================================================================

  ftermNew   = .FALSE. ! no new forcing term
  residFterm = .FALSE. ! do not add forcing term to residual

! ==============================================================================
! Determine whether have moving grids
! ==============================================================================
  dtlimcount = 0
  moveGrid = .FALSE.

  DO iReg = 1,global%nRegionsLocal
    IF ( regions(iReg)%mixtInput%moveGrid .EQV. .TRUE. ) THEN 
      moveGrid = .TRUE.
    END IF ! regions
  END DO ! iReg

! ******************************************************************************
! Loop over iterations/time steps
! ****************************************************************************** 

! BBR - begin
#ifdef PLAG
  IF ( global%plagUsed .EQV. .TRUE. ) THEN
     CALL PLAG_CalcNPclsGlobal(regions)
  END IF ! global%plagUsed
#endif

  IF(global%casename .EQ. "cyldet")THEN
    !CALL RFLU_WriteRProf(regions)
    CALL RFLU_WritePM(regions)
  END IF 

  DO iReg=1,global%nRegionsLocal
     pRegion => regions(iReg)

     CALL RFLU_WriteFlowWrapper(pRegion)

  END DO

  iter = 0.0_RFREAL
  cumulreal = 0.0_RFREAL
  cumulcpu = 0.0_RFREAL

  ! Set timeSince* variable in case of a safe restart (out of sync)
  ! WARNING : printTime and Probe time get out of sync relative to 
  ! run restarted from
  IF ( global%currentTime .NE. 0.0_RFREAL ) THEN
    global%timeSinceFomWrite = global%currentTime - global%previousTime 
    global%timeSinceWrite = global%currentTime - global%previousTime
  END IF

#ifdef FOLDER
  ! create a folder to store probes data
  IF ( global%nProbes > 0 ) THEN
    IF ( global%myProcid == 0 ) THEN
      CALL system ( 'mkdir Probes')
    END IF
  END IF
#endif 
!rahul folder debug

! BBR - end

!CRN - Begin
!Set initial time before starting loop
 CALL MPI_Barrier(global%mpiComm,errorFlag)
 IF(global%myProcid == 0) THEN
    timerStart=MPI_Wtime()
 END IF
!CRN - end

  DO

! BBR - begin
   
  !CRN - Begin
  IF(global%myProcid == 0) THEN
    timerEnd=MPI_Wtime()
    elapsedtime = timerEnd - timerStart
  END IF
  !CRN - end


  !CALL SYSTEM_CLOCK ( timeroutineat, clock_rate, clock_max )
  !elapsedtime = REAL( timeroutineat - global%timeruninit ) / REAL(clock_rate)

  IF ( global%TimeCycle ) THEN
    IF ( global%myProcid .EQ. 1 ) THEN
      IF ( iter .EQ. 0.0_RFREAL ) THEN
        OPEN(IF_TIME,file='TimePerCycle.tim',status='unknown')
      END IF
      CALL SYSTEM_CLOCK ( t1, clock_rate, clock_max )
      CALL CPU_TIME ( tc1 )
      ! BBR - testing safe restart timing works properly
      !WRITE(999,*)global%wallTime,elapsedtime,global%wallTime-elapsedtime
    END IF
  END IF
! BBR - end

! ==============================================================================
!   Check for non-zero system time or iteration step. This is needed because 
!   of changed restart mechanism. If restarting code when final time was 
!   reached, can get zero time step which gives indeterminate or infinite 
!   grid speeds for moving grid computations. 
! ==============================================================================

    IF ( global%flowType == FLOW_UNSTEADY ) THEN 
      IF ( dTimeSystem == 0.0_RFREAL ) THEN
        global%warnCounter = global%warnCounter + 1       
      
        IF ( global%myProcid == MASTERPROC .AND. & 
             global%verbLevel > VERBOSE_NONE ) THEN 
          WRITE(STDOUT,'(A,2(1X,A))') SOLVER_NAME,'*** WARNING *** ', & 
                'Nothing to be done. Returning to calling procedure.'
        END IF ! global
        
        EXIT
      END IF ! dTimeSystem
    ELSE
      IF ( dIterSystem == 0 ) THEN
        global%warnCounter = global%warnCounter + 1       
      
        IF ( global%myProcid == MASTERPROC .AND. & 
             global%verbLevel > VERBOSE_NONE ) THEN 
          WRITE(STDOUT,'(A,2(1X,A))') SOLVER_NAME,'*** WARNING *** ', & 
                'Nothing to be done. Returning to calling procedure.'
        END IF ! global
        
        EXIT
      END IF ! dIterSystem
    END IF ! global%flowType

! ==============================================================================
!   Update iteration counter for steady flow. NOTE does not have to be done
!   here, but is more consistent, otherwise can get output for iteration 0 
!   from RFLU_DecidePrint, and again from iteration 1 because of how MOD 
!   function works. If update iteration counters here, only get output once.
!   Need iteration counters for case of non-dissipative implicit solver also.
! ==============================================================================

    IF ( global%flowType == FLOW_STEADY .OR. &
         global%solverType == SOLV_IMPLICIT_HM ) THEN
      global%currentIter      = global%currentIter      + 1
      global%iterSinceRestart = global%iterSinceRestart + 1
    END IF ! global%flowType

! ==============================================================================
!   Compute time step. For unsteady flow, compute minimum time step and make 
!   sure do not run over maximum time.
!   For non-dissipative solver dtMin is fixed.
! ==============================================================================

    IF ( global%solverType /= SOLV_IMPLICIT_HM ) THEN
      DO iReg = 1,global%nRegionsLocal  
        pRegion => regions(iReg)
        
        IF ( pRegion%mixtInput%flowModel == FLOW_EULER ) THEN      
          CALL RFLU_TimeStepInviscid(pRegion)
        ELSE 
          CALL RFLU_TimeStepViscous(pRegion)
        END IF ! pRegion
      END DO ! iReg
    
#ifdef PLAG
      IF ( global%plagUsed .AND. (global%flowType == FLOW_UNSTEADY) ) THEN
        IF ( regions(1)%plagInput%flagStability ) THEN
          CALL DummyRungeKuttaMP(regions)
        END IF ! regions(1)%plagInput%flagStability
      END IF ! global%plagUsed
#endif

      IF ( global%flowType == FLOW_UNSTEADY ) THEN         
         CALL RFLU_MinimumTimeStep(regions)
#ifdef GENX
        IF ( global%dtMin < global%dtImposed .AND. &
             global%dtMin < global%dtMinLimit ) THEN
          dtlimcount = dtlimcount + 1
          IF ( dtlimcount > 2) THEN
            dtlimcount = 0
            stophandle   = COM_get_function_handle('Rocman.interrupt')
            IF(stophandle <= 0) THEN
              WRITE(*,*) 'Could not get Rocman.stop function handle.'
            ELSE
              CALL COM_call_function(stophandle,2,3,'Rocflu dt < limit, &
                                     requesting remesh')
            ENDIF
          ENDIF
        ENDIF
#endif
        IF ( global%timeSinceRestart + global%dtMin > dTimeSystem ) THEN 
           global%dtMin = dTimeSystem - global%timeSinceRestart
           finished     = .TRUE.                     
        END IF ! global%timeSinceRestart
      END IF ! global%flowType            
    ELSE
      global%dtMin = global%dtImposed

      IF ( global%flowType == FLOW_UNSTEADY ) THEN         
        IF ( global%timeSinceRestart + global%dtMin > dTimeSystem ) THEN 
           finished     = .TRUE.                     
        END IF ! global%timeSinceRestart
      END IF ! global%flowType            
    END IF ! solverType
     
! ==============================================================================
!   Move or generate new grid
! ==============================================================================

    IF ( global%flowType == FLOW_UNSTEADY ) THEN    
      IF ( moveGrid .EQV. .TRUE. ) THEN     
        CALL RFLU_GetDeformationWrapper(regions)        
        CALL RFLU_MoveGridWrapper(regions)
            
        DO iReg=1,global%nRegionsLocal
          pRegion => regions(iReg)
                   
          CALL RFLU_BuildGeometry(pRegion)          
          CALL RFLU_ComputeGridSpeeds(pRegion)          

          IF(global%zoomFactor > 1.0_RFREAL) THEN
             CALL RFLU_UnZoomGridSpeeds(pRegion)
          ENDIF

          IF ( global%checkLevel == CHECK_HIGH ) THEN           
            CALL RFLU_CheckGridSpeeds(pRegion)               
          END IF ! global%checkLevel
        END DO ! iReg
      END IF ! moveGrid     
    END IF ! global%flowType  

! ==============================================================================
!   Recompute weights
! ==============================================================================
#ifndef PROPONLY

    IF ( global%flowType == FLOW_UNSTEADY ) THEN    
      IF ( moveGrid .EQV. .TRUE. ) THEN     
        DO iReg=1,global%nRegionsLocal
          pRegion => regions(iReg)
          pMixtInput => pRegion%mixtInput

          IF ( pMixtInput%spaceOrder > 1 ) THEN 
            CALL RFLU_ComputeWtsC2CWrapper(pRegion,pMixtInput%spaceOrder-1)
          END IF ! pMixtInput%spaceOrder      

          IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN
            CALL RFLU_ComputeWtsF2CWrapper(pRegion,pMixtInput%spaceOrder-1)
          END IF ! pMixtInput%flowModel           
           
          DO iPatch = 1,pRegion%grid%nPatches
            pPatch => pRegion%patches(iPatch)

            IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
              CALL RFLU_ComputeWtsBF2CWrapper(pRegion,pPatch,pPatch%spaceOrder)
            END IF ! RFLU_DecideNeedBGradFace
          END DO ! iPatch
        END DO ! iReg
      END IF ! moveGrid     
    END IF ! global%flowType 

! ==============================================================================
!   Relocate probes
! ==============================================================================
       
    IF ( global%flowType == FLOW_UNSTEADY ) THEN    
      IF ( moveGrid .EQV. .TRUE. ) THEN     
        IF ( global%nProbes > 0 ) THEN 
          DO iReg = 1,global%nRegionsLocal
            pRegion => regions(iReg)
            CALL RFLU_FindProbeCells(pRegion)
          END DO ! iReg

          CALL RFLU_PrintProbeInfo(global)
        END IF ! global         
      END IF ! moveGrid     
    END IF ! global%flowType             
            
! ==============================================================================
!   Compute new solution
! ==============================================================================

    !Fred - Restting JWL Logical flags for storing initial iterative guess for
    !each cell
    IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
    pRegion%mixt%set_p_JWL(:) = .FALSE.
    pRegion%mixt%set_e_JWL(:) = .FALSE.
    END IF !End reset

    global%forceX  = 0.0_RFREAL
    global%forceY  = 0.0_RFREAL
    global%forceZ  = 0.0_RFREAL
    global%massIn  = 0.0_RFREAL
    global%massOut = 0.0_RFREAL  

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
        CALL RFLU_HM_PredCorrMP(regions)
      ELSE 
        CALL RungeKuttaMP(regions)      
! TEMPORARY - At present, do not use operator-split integration of chemistry
!             source terms based on modifications by Luca. This means that 
!             call source term residual function directly in SourceTermsMP.F90
!      CALL IntegrateSourceTermsMP(regions)  
! END TEMPORARY
      END IF ! global%solverType
    ELSE
      CALL RFLU_ExplicitMultiStage(regions)
    END IF ! global%flowType
 
#ifdef STATS
! ==============================================================================
!   Get statistics
! ==============================================================================

    CALL GetStatistics(regions)
#endif

! ==============================================================================
!   Reset global%timeSince* variables. NOTE this must be done right before
!   the update of the various time variables, otherwise calling the routines
!   RFLU_Decide* from within rungeKuttaMP or explicitMultiStage will not work
!   correctly.
! ==============================================================================

    IF ( RFLU_DecidePrint(global) .EQV. .TRUE. ) THEN 
      IF ( global%flowType == FLOW_UNSTEADY ) THEN    
        IF ( global%iterSinceRestart > 1 ) THEN 
          global%timeSincePrint = 0.0_RFREAL
        END IF ! global%iterSinceRestart
      END IF ! global%flowType      
    END IF ! RFLU_DecidePrint

    IF ( RFLU_DecideWrite(global) .EQV. .TRUE. ) THEN 
      IF ( global%flowType == FLOW_UNSTEADY ) THEN    
        global%timeSinceWrite = 0.0_RFREAL
      END IF ! global%flowType      
    END IF ! RFLU_DecideWrite

    IF ( RFLU_DecideWriteFom(global) .EQV. .TRUE. ) THEN
      IF ( global%flowType == FLOW_UNSTEADY ) THEN
        global%timeSinceFomWrite = 0.0_RFREAL
      END IF ! global%flowType
    END IF ! RFLU_DecideWriteFom

    IF ( RFLU_DecideWriteProbes(global) .EQV. .TRUE. ) THEN 
      IF ( global%flowType == FLOW_UNSTEADY ) THEN    
        IF ( global%iterSinceRestart > 1 ) THEN 
          global%timeSinceProbe = 0.0_RFREAL
        END IF ! global%iterSinceRestart 
      END IF ! global%flowType      
    END IF ! RFLU_DecideWriteProbes   

    !begin BBR
    IF ( RFLU_DecideSmallWrite(global) .EQV. .TRUE. ) THEN
      IF ( global%flowType == FLOW_UNSTEADY ) THEN
        global%halftimeSinceWrite = 0.0_RFREAL
      END IF ! global%flowType      
    END IF ! RFLU_DecideSmallWrite
    !end BBR
! ==============================================================================
!   Update times for unsteady flow
! ==============================================================================

! ENDIF PROPONLY
#endif

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      global%currentTime      = global%currentTime      + global%dtMin      
      global%timeSinceRestart = global%timeSinceRestart + global%dtMin
      
      global%timeSincePrint = global%timeSincePrint + global%dtMin
      global%timeSinceFomWrite = global%timeSinceFomWrite + global%dtMin
      global%timeSinceWrite = global%timeSinceWrite + global%dtMin
      global%timeSinceProbe = global%timeSinceProbe + global%dtMin      
      !begin BBR
      global%halftimeSinceWrite = global%halftimeSinceWrite + global%dtMin
      !end BBR     
 
      global%iterSinceRestart = global%iterSinceRestart + 1
    END IF ! global%flowType

! ==============================================================================
!   Decide whether to print/write convergence, data, and probes
! ==============================================================================

    doFomWrite = RFLU_DecideWriteFom(global)
    doPrint = RFLU_DecidePrint(global)
    doWrite = RFLU_DecideWrite(global)
    doProbe = RFLU_DecideWriteProbes(global)
    ! begin BBR
    doSmallWrite = RFLU_DecideSmallWrite(global)
    ! end BBR 

! BBR - begin
! ==============================================================================
!  Create time/iteration-based directory for the next dataset
! ==============================================================================

#ifdef FOLDER   
 IF( ((global%currentIter .GT. 0) .OR. (global%currentTime .GT. 0.0_RFREAL)) &
    .AND. ((doFomWrite .EQV. .TRUE.) .OR. (doWrite .EQV. .TRUE.) ))THEN! .OR. &
    !(doSmallWrite .EQV. .TRUE.)) )THEN
    
    IF ( global%flowType == FLOW_UNSTEADY ) THEN
    WRITE(NameSolDir,'(1PE11.5)') global%currentTime
    ELSEIF ( global%flowType == FLOW_STEADY ) THEN
    WRITE(NameSolDir,'(i10)') global%currentIter
    END IF
    IF(global%myProcid==MASTERPROC) THEN
      CALL system( 'mkdir -p SOL_'//NameSolDir )
      CALL system( 'mkdir -p PARAVIEW_'//NameSolDir )
    END IF   

 END IF
#endif
!rahul folder debug   

! BBR - end

! BBR - begin
! ==============================================================================
!   Check whether it is close enough to walltime to force a safe restart dump
! ==============================================================================

! CRN - Edited MPI call
!Have the master processor check the walltime and then pass a variable to all
!other processors. The variable is safew_l (that is lowercase L). A 0 means to
!keep going. A one means to stop and print the final solution before exiting.

IF(global%myProcid == 0) THEN

   IF ( ( (global%wallTime - elapsedtime) .LE. global%safeWriteTime  )  &
     .AND. ( donesafewrite .EQV. .FALSE. ) ) THEN
       safew_l = 1

   END IF

END IF

!Broadcast the data to all processors from the master processor
CALL MPI_Bcast(safew_l, 1, MPI_INT, 0,global%mpiComm,global%mpierr)

!All processors check the broadcasted variable that they recieved
IF(safew_l == 1) THEN

  doWrite = .TRUE.
  donesafewrite = .TRUE.
  global%SafeWrite = safew_l

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    WRITE(NameSolDir,'(E11.5)') global%currentTime
  ELSEIF ( global%flowType == FLOW_STEADY ) THEN
    WRITE(NameSolDir,'(i10)') global%currentIter
  END IF
#ifdef FOLDER
  IF(global%myProcid==MASTERPROC) THEN
      CALL system( 'mkdir -p SOL_'//NameSolDir )
      CALL system( 'mkdir -p PARAVIEW_'//NameSolDir )
  END IF
#endif
!rahul folder debug
END IF

! END CRN


! BBR - end

! ==============================================================================
!   Check for stop file
! ==============================================================================

    INQUIRE(FILE="STOP",EXIST=stopFileExists)
    IF ( stopFileExists .EQV. .TRUE. ) THEN 
      IF ( global%myProcid  == MASTERPROC .AND. &
           global%verbLevel /= VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A)')      SOLVER_NAME
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Stop file detected!'
        WRITE(STDOUT,'(A)')      SOLVER_NAME        
      END IF ! global

      finished = .TRUE.
    END IF ! stopFileExists
  
! ==============================================================================
!   Check for end of time stepping
! ==============================================================================

    IF ( global%flowType == FLOW_UNSTEADY ) THEN    
      IF ( global%timeSinceRestart >= dTimeSystem ) THEN 
        finished = .TRUE.
      END IF ! global%timeSinceRestart
    ELSE
      IF ( (doPrint .EQV. .TRUE.) .OR. (global%iterSinceRestart >= dIterSystem) ) THEN 
        CALL RFLU_ResidualNorm(regions)    
        IF ( (global%iterSinceRestart >= dIterSystem) .OR. &
             (global%residual/global%resInit <= global%resTol) ) THEN 
           finished = .TRUE.
        END IF ! global%iterSinceRestart      
      END IF ! doPrint
    END IF ! global%flowType

! ==============================================================================
!   Write convergence (file & screen) and total mass (file) history
! ==============================================================================

    IF ( (doPrint .EQV. .TRUE.) .OR. (finished .EQV. .TRUE.) ) THEN
      CALL RFLU_PrintWriteConvergence(global)
    
      !BBR - begin
      !IF (global%casename .EQ. "cyldet") THEN 
       !CALL RFLU_WritePM(regions)
       !CALL RFLU_WriteInteg(regions)
      !END IF ! casename
      !BBR - end   
   
#ifndef GENX      
      IF ( moveGrid .EQV. .TRUE. ) THEN
        CALL RFLU_ComputeIntegralValues(regions) 
        CALL WriteTotalMass(regions)
      END IF ! moveGrid
#endif    
      
#ifndef PROPONLY
      DO iReg = 1,global%nRegionsLocal
        IF ( regions(iReg)%mixtInput%spaceDiscr == DISCR_OPT_LES ) THEN 
          CALL RFLU_WriteStatsFileOLES(global)
        END IF ! mixtInput
      END DO ! iReg
#endif
    END IF ! doPrint

#ifndef PROPONLY
! ==============================================================================
!   Compute forces and mass flow
! ==============================================================================

    IF ( global%forceFlag .EQV. .TRUE. ) THEN 
      IF ( (doFomWrite .EQV. .TRUE.) .OR. (doWrite .EQV. .TRUE.) &
           .OR. (finished .EQV. .TRUE.) ) THEN 
        DO iReg = 1,global%nRegionsLocal
          pRegion => regions(iReg)

          CALL RFLU_ComputeLocalForcesMoments(pRegion)     
        END DO ! iReg

        CALL RFLU_ComputeGlobalForcesMoments(regions)

        DO iReg = 1,global%nRegionsLocal
          pRegion => regions(iReg)

          CALL RFLU_TSI_ComputeGlobalThrustSI(pRegion)     
        END DO ! iReg

        pRegion => regions(1)

        IF ( pRegion%global%myProcid == MASTERPROC ) THEN 
          CALL RFLU_PrintGlobalForcesMoments(pRegion)
          CALL RFLU_WriteGlobalForcesMoments(pRegion)
          CALL RFLU_TSI_PrintGlobalVals(pRegion)
          CALL RFLU_TSI_WriteGlobalVals(pRegion)
        END IF ! pRegion
      END IF ! doFomWrite
    END IF ! global%forceFlag

! ==============================================================================
!   Write Particle velocity and acceleration
! ==============================================================================

    IF ( (doFomWrite .EQV. .TRUE.) .OR. (doWrite .EQV. .TRUE.) &
         .OR. (finished .EQV. .TRUE.) ) THEN
      pRegion => regions(1)
      global => pRegion%global
 
      IF ( pRegion%global%myProcid == MASTERPROC ) THEN
        IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
          CALL RFLU_MVF_WritePatchVelAccel(pRegion)
        END IF ! global%mvFrameFlag
      END IF ! pRegion
    END IF ! doFomWrite

! ==============================================================================
!   Store probe data
! ==============================================================================

    IF ( global%nProbes > 0 ) THEN
      IF ( doProbe .EQV. .TRUE. ) THEN 
        DO iReg = 1,global%nRegionsLocal
          CALL WriteProbe(regions,iReg)
        END DO ! iReg
      END IF ! doProbe
    END IF ! global

#ifndef GENX
! ==============================================================================
!   Store flow field (and grid if moving). Write restart info file after 
!   flow (and grid) files so that should those not be written due to reaching
!   the time limit, the restart file will not contain the time level of the 
!   incomplete flow (and grid) files.
! ==============================================================================

    !BBR - begin - writing light output to track sanity of run 
    ! to be replaced by some Catalyst function later
   
    IF (global%casename .EQ. "cyldet" ) THEN
      IF ( (doSmallWrite .EQV. .TRUE.) .AND. (finished .EQV. .FALSE.) &
        .AND. (doWrite .EQV. .FALSE.) ) THEN

        !CALL RFLU_WriteSlice(pRegion)
        !CALL RFLU_WriteRProf(regions)
        CALL RFLU_WritePM(regions)
        !CALL RFLU_WriteInteg(regions)
 
      END IF ! doSmallWrite
     END IF !casename
    !BBR - end

    IF( (doWrite .EQV. .TRUE.) .AND. (finished .EQV. .FALSE.) )THEN

    !BBR - Begin - write Prediction Metrics (PM)
    ! to be replaced by some Catalyst function later

    IF (global%casename .EQ. "cyldet" ) THEN
      CALL RFLU_WritePM(regions)
    END IF

    !BBR - End 

#ifdef PLAG
      IF ( global%plagUsed .EQV. .TRUE. ) THEN
        CALL PLAG_CalcNPclsGlobal(regions)

        IF ( global%myProcid == MASTERPROC ) THEN 
          pRegionSerial => regions(0)

          CALL PLAG_RFLU_WriteDimensions(pRegionSerial)
        END IF ! global%myProcid
      END IF ! global%plagUsed
#endif

      DO iReg=1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_WriteDimensionsWrapper(pRegion,WRITE_DIMENS_MODE_MAYBE)
                
        IF ( moveGrid .EQV. .TRUE. ) THEN        
          CALL RFLU_WriteGridWrapper(pRegion)
          CALL RFLU_WriteGridSpeedsWrapper(pRegion)
          
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_NONE ) THEN           
            CALL RFLU_PrintGridInfo(pRegion)
          END IF ! global%myProcid
        END IF ! moveGrid
                       
        CALL RFLU_WriteFlowWrapper(pRegion)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          CALL RFLU_WriteAuxVarsWrapper(pRegion)
        END IF ! solverType

        CALL RFLU_BXV_WriteVarsWrapper(pRegion)

        IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN 
          CALL RFLU_WritePatchCoeffsWrapper(pRegion)
        END IF ! global%patchCoeffFlag
        
#ifdef PLAG
        IF ( global%plagUsed .EQV. .TRUE. ) THEN 
          IF ( PLAG_DecideHaveSurfStats(pRegion) .EQV. .TRUE. ) THEN
            CALL PLAG_WriteSurfStatsWrapper(pRegion)
          END IF ! PLAG_DecideHaveSurfStats
        END IF ! global%plagUsed
#endif  
        
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_NONE ) THEN          
          CALL RFLU_PrintFlowInfoWrapper(pRegion)

          IF ( global%solverType /= SOLV_IMPLICIT_HM ) THEN          
            IF ( global%verbLevel > VERBOSE_LOW ) THEN 
              CALL RFLU_PrintChangeInfo(pRegion)
            END IF ! global%verbLevel
          END IF ! solverType

#ifdef PLAG
          IF ( global%plagUsed .EQV. .TRUE. ) THEN 
            CALL PLAG_PrintNPclsGlobal(pRegion)
          END IF ! global%plagUsed
#endif
        END IF ! global%myProcid        
      END DO ! iReg
    
      !BBR - begin 
      IF (donesafewrite .EQV. .FALSE. ) THEN ! prevent writing twice in *.rin
                                             ! RFLU_EndFlowSolver will do it 
        CALL RFLU_WriteRestartInfo(global)
      END IF      
      !BBR - end

    END IF ! doWrite

#ifdef STATS
! ==============================================================================
!   Output statistics
! ==============================================================================

    IF ( (doWrite .EQV. .TRUE.) .AND. (finished .EQV. .FALSE.) .AND. &
         (global%doStat == ACTIVE) ) THEN
      IF (global%myProcid==MASTERPROC .AND. &
          global%verbLevel/=VERBOSE_NONE) THEN
        WRITE(STDOUT,'(A)') SOLVER_NAME,'Saving statistics ...'
      ENDIF
      
      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)
        CALL RFLU_WriteStat(pRegion)
      END DO ! iReg
    END IF ! iReg
#endif
#endif
#endif

! ==============================================================================
!   If run finished, update GENX buffers and exit
! ==============================================================================

    IF ( finished .EQV. .TRUE. ) THEN 
#ifdef GENX
       DO iReg = 1,global%nRegionsLocal
         CALL RFLU_PutBoundaryValues(regions(iReg))
       END DO ! iReg
       
       global%timeStamp = global%currentTime       
#endif    
      EXIT
    END IF ! finished

! BBR - Begin - timing

  IF ( global%TimeCycle ) THEN
    IF ( global%myProcid .EQ. 1)THEN
      iter = iter + 1.0_RFREAL

      CALL SYSTEM_CLOCK ( t2, clock_rate, clock_max )
     ! write ( *, * ) 'Elapsed real time = ', real ( t2 - t1 ) / real ( clock_rate )

      CALL CPU_TIME ( tc2 )
     ! write ( *, * ) 'Elapsed CPU time = ', tc2 - tc1

      cumulreal = cumulreal + REAL ( t2 - t1 ) / REAL ( clock_rate )
      cumulcpu = cumulcpu + tc2 - tc1

      WRITE(IF_TIME,'(5E13.4)') iter, REAL ( t2 - t1 ) / REAL ( clock_rate ) &
              ,tc2 - tc1, cumulreal/iter, cumulcpu/iter

    END IF !global%myProcid
  END IF ! global%TimeCycle

! BBR - End - timing

! BBR - Begin - Exit main loop if saferestart done

! ==============================================================================
!   If run saferestart done
!   Total runtime close to walltime requested. Restart dump accomplished. 
!   Exiting run cleanly
! ==============================================================================

  IF ( global%SafeWrite .EQ. 1 ) THEN
    CALL MPI_Barrier(global%mpiComm,errorFlag)
    EXIT  ! exiting loop
  END IF
! BBR - End - Exit main loop

  END DO ! loop over time/iterations



  IF (global%myProcid==MASTERPROC .AND. (global%SafeWrite .EQ. 1) ) THEN
    WRITE(STDOUT,'(A)') SOLVER_NAME,'Safe Restart Done, exiting this run ...'
  END IF

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TimeStepping

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_TimeStepping.F90,v $
! Revision 1.8  2016/02/08 17:00:56  fred
! Fixing Vulcan compiling issue
!
! Revision 1.7  2016/02/06 17:22:18  fred
! Adding JWL Initial condition flags
!
! Revision 1.6  2015/12/19 00:39:46  rahul
! Added compiler flag FOLDER.
!
! Revision 1.5  2015/07/24 00:07:10  neal
! changed from using SYSTEM_CLOCK() to the MPI_Wtime() function for timing runs due to issues with SYSTEM_CLOCK() during long simulations.
!
! Revision 1.4  2015/07/23 23:11:18  brollin
! 1) The pressure coefficient of the  collision model has been changed back to its original form
! 2) New options in the format of the grid and solutions have been added. Now the user can choose the endianness, and convert from one to the over in rfluconv.
! 3) The solutions are now stored in folders named by timestamp or iteration number.
! 4) The address enty in the hashtable has been changed to an integer(8) for cases when the grid becomes very large.
! 5) RFLU_WritePM can now compute PM2 on the fly for the Macroscale problem
!
! Revision 1.3  2015/04/15 20:26:38  neal
! Changed the safe restart algorithm to fix a bug that caused code to hang.
!
! Revision 1.2  2015/02/11 16:58:26  brollin
! New fix to Safe Restart Feature. Now Safe restart communicated to all ranks
! once any rank hit the Safe restart condition.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.7  2009/08/28 18:29:48  mtcampbe
! RocfluMP integration with Rocstar and some makefile tweaks.  To build
! Rocstar with new Rocflu:
! make ROCFLU=RocfluMP
! To build Rocstar with the new RocfluND:
! make ROCFLU=RocfluMP HYPRE=/the/hypre/install/path
!
! Revision 1.6  2008/12/06 08:43:48  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2008/01/19 20:19:45  haselbac
! Added calls to PLAG_DecideHaveSurfStats
!
! Revision 1.3  2007/11/28 23:05:33  mparmar
! Modifying timestepping according to SOLV_IMPLICIT_HM
!
! Revision 1.2  2007/06/18 18:12:24  mparmar
! Writing files for moving reference frame
!
! Revision 1.1  2007/04/09 18:49:57  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:02  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.70  2007/03/31 23:54:34  haselbac
! Bug fix: Writing of PLAG dims should only be called by serial region on master proc
!
! Revision 1.69  2007/03/27 00:41:17  haselbac
! Added calls to calculate, write, and print particle dimensions
!
! Revision 1.68  2007/02/18 03:17:56  mtcampbe
! Added proponly for Rocstar proponly runs
!
! Revision 1.67  2006/10/20 21:33:07  mparmar
! Added calls to write thrust/specific impulse and again added code for NSCBC implementation
!
! Revision 1.66  2006/09/12 14:58:44  mtcampbe
! Moved include of Roccom after IMPLICIT NONE
!
! Revision 1.65  2006/09/11 15:43:11  mtcampbe
! Added Rocman interrupt call to support automatic remeshing.
!
! Revision 1.64  2006/08/19 15:40:04  mparmar
! Changed because of NSCBC implementation
!
! Revision 1.63  2006/04/07 16:04:03  haselbac
! Adapted to changes in bf2c wts computation
!
! Revision 1.62  2006/04/07 15:19:22  haselbac
! Removed tabs
!
! Revision 1.61  2006/04/07 14:55:06  haselbac
! Adapted to changes in bface stencil routine
!
! Revision 1.60  2006/03/09 14:10:31  haselbac
! Now call wrapper routine for F2C weights
!
! Revision 1.59  2006/01/06 22:16:13  haselbac
! Adapted to name changes, removed commented-out call to ExplMS routine
!
! Revision 1.58  2005/11/10 16:51:29  fnajjar
! Added plagUsed IF statement around PLAG routines
!
! Revision 1.57  2005/10/25 19:39:23  haselbac
! Added IF on forceFlag
!
! Revision 1.56  2005/10/05 14:20:35  haselbac
! Added call to recompute bface wts for moving grids
!
! Revision 1.55  2005/08/09 00:59:56  haselbac
! Enclosed writing of patch coeffs within IF (patchCoeffFlag)
!
! Revision 1.54  2005/06/06 14:23:35  haselbac
! Adapted to Lucas changes
!
! Revision 1.53  2005/05/16 20:44:55  haselbac
! Now compute time step also for steady flow, call RFLU_ExplicitMultiStage
!
! Revision 1.52  2005/04/29 12:50:00  haselbac
! Added USE RFLU_ModProbes, removed interfaces for probe routines
!
! Revision 1.51  2005/04/20 14:44:04  haselbac
! Removed CHECK_UNIFLOW code section
!
! Revision 1.50  2005/04/15 15:07:26  haselbac
! Now use RFLU_PrintWriteConvergence, fixed bug in writing forces
!
! Revision 1.49  2004/12/28 20:28:20  wasistho
! moved statistics routines into module ModStatsRoutines
!
! Revision 1.48  2004/12/21 15:05:11  fnajjar
! Included calls for PLAG surface statistics
!
! Revision 1.47  2004/11/29 17:17:29  wasistho
! use ModInterfacesStatistics
!
! Revision 1.46  2004/11/03 17:05:53  haselbac
! Removed HACK_PERIODIC ifdef
!
! Revision 1.45  2004/10/19 19:29:36  haselbac
! Added writing of grid speeds, changed time step interfaces, cosmetics
!
! Revision 1.44  2004/07/06 15:14:53  haselbac
! Adapted to changes in libflu and modflu, cosmetics
!
! Revision 1.43  2004/06/16 20:01:15  haselbac
! Added forces and moments stuff, cosmetics
!
! Revision 1.42  2004/04/14 02:10:07  haselbac
! Removed call to DescaleGridSpeeds
!
! Revision 1.41  2004/04/01 21:30:35  haselbac
! Added IntegrateSourceTermsMP, cosmetic changes
!
! Revision 1.40  2004/03/17 04:28:25  haselbac
! Adapted call to RFLU_WriteDimensionsWrapper
!
! Revision 1.39  2004/03/11 16:34:32  fnajjar
! ACH: Changed call to RFLU_WriteDimWrapper bcos of Rocpart
!
! Revision 1.38  2004/01/31 03:59:03  haselbac
! Improved updating of iteration counters
!
! Revision 1.37  2004/01/29 22:59:24  haselbac
! Added setting of timeSince* vars for improved useability of RFLU_Decide* 
! funcs
!
! Revision 1.36  2003/12/04 03:30:09  haselbac
! Adapted recomputation of weights
!
! Revision 1.35  2003/11/25 21:04:46  haselbac
! Added call to RFLU_PrintFlowInfoWrapper, cosmetic changes
!
! Revision 1.34  2003/11/04 01:35:27  haselbac
! Bug fix: added init for timeSinceRestart for GENX
!
! Revision 1.33  2003/10/29 21:39:57  haselbac
! Bug fix and clean-up: Writing of data to screen and files
!
! Revision 1.32  2003/10/15 02:44:21  haselbac
! Removed unnecessary initialization of dtMin
!
! Revision 1.31  2003/10/03 20:47:21  haselbac
! Changed from RungeKutta to RungeKuttaMP
!
! Revision 1.30  2003/08/13 20:28:49  haselbac
! Fixed bug with writing probe data within GENx
!
! Revision 1.29  2003/07/22 02:11:40  haselbac
! Added global%warnCounter
!
! Revision 1.28  2003/07/09 22:37:48  haselbac
! Added check to avoid wrong restart, find probes on moving grids
!
! Revision 1.27  2003/06/20 22:36:23  haselbac
! Added call to RFLU_WriteRestartInfo, fixed bug in steady conv check
!
! Revision 1.26  2003/04/24 15:43:35  haselbac
! Adapted interface to RFLU_PutBoundaryValues
!
! Revision 1.25  2003/04/07 14:27:44  haselbac
! Added interval writing of probe info
!
! Revision 1.24  2003/03/31 16:18:17  haselbac
! Replaced MoveGrid call by call to wrapper
!
! Revision 1.23  2003/03/15 18:59:56  haselbac
! Added calls to DescaleGridSpeeds and GetDeformationWrapper
!
! Revision 1.22  2003/02/25 21:47:39  haselbac
! Added call to DescaleGridSpeeds
!
! Revision 1.21  2003/01/28 14:52:08  haselbac
! Changes to be consistent with rewrite of RFLU_InitFlowSolver
!
! Revision 1.20  2002/12/20 23:21:33  haselbac
! Fixed output bug: no output for verbosity=0
!
! Revision 1.19  2002/11/15 21:26:49  haselbac
! Added RFLU_ComputeIntegralValues
!
! Revision 1.18  2002/11/08 21:36:25  haselbac
! Fixed bug in grid-speed comp, added RFLU_CheckGridSpeeds and WriteTotalMass
!
! Revision 1.17  2002/11/02 02:04:49  wasistho
! Added TURB statistics
!
! Revision 1.16  2002/10/27 19:19:31  haselbac
! Logical modifications of last part, cosmetic redesign
!
! Revision 1.15  2002/10/16 21:18:50  haselbac
! Added call to RFLU_NewGrid, some rearrangement of calls
!
! Revision 1.14  2002/10/12 15:00:44  haselbac
! Added statements for time check for GENX runs
!
! Revision 1.13  2002/10/05 19:35:19  haselbac
! Call wrapper routines for output, write probe data
!
! Revision 1.12  2002/09/09 15:51:56  haselbac
! global and mixtInput under regions, adapated calls to RFLU_CheckCalcGrad, 
! OLES routines removed, added viscous calls
!
! Revision 1.11  2002/07/25 14:23:40  haselbac
! Added OLES calls, gradient check, and HACK_PERIODIC segment
!
! Revision 1.10  2002/06/18 00:26:20  wasistho
! Added prefix SOLVER NAME to satistics STDOutput
!
! Revision 1.9  2002/06/17 13:34:12  haselbac
! Prefixed SOLVER_NAME to all screen output
!
! Revision 1.8  2002/06/14 22:26:08  wasistho
! update statistics
!
! Revision 1.7  2002/06/14 21:54:35  wasistho
! Added time avg statistics
!
! Revision 1.6  2002/06/14 20:19:46  haselbac
! Deleted ModLocal, changed local%nRegions to global%nRegionsLocal
!
! Revision 1.5  2002/06/10 21:35:50  haselbac
! Changed order of convergence and file writing, changed flag to CHECK_UNIFLOW
!
! Revision 1.4  2002/06/05 18:59:12  haselbac
! Miscellaneous small functionality changes
!
! Revision 1.3  2002/05/28 14:02:04  haselbac
! Activated unsteady routines
!
! Revision 1.2  2002/05/04 17:13:15  haselbac
! Many changes to enable flow solution
!
! Revision 1.1  2002/04/11 18:57:06  haselbac
! Initial revision, commented out RFLO stuff
!
! ******************************************************************************

