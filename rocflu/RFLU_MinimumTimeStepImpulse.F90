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
! Purpose: Calculate the minimum time step for all regions on all processors.
!
! Description: First determine smallest time step and region in which it 
!   occurs. If this smallest time step is smaller than that imposed by the 
!   user and smaller than the user-specified lower limit, then print out 
!   also in which cell the smallest time step occurs. This can be useful 
!   to analyze problems, particularly with moving grids. 
!
! Input: 
!   regions           Data of all regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_MinimumTimeStepImpulse.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_MinimumTimeStepImpulse(regions)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  USE ModMPI

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,iCellDtMin,iReg,iRegDtMin
  INTEGER, DIMENSION(:), ALLOCATABLE :: globalValsInt,localValsInt  
  REAL(RFREAL) :: dtMinImpulse
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: globalValsReal,localValsReal
  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_MinimumTimeStepImpulse.F90,v $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_MinimumTimeStepImpulse',__FILE__)

! ******************************************************************************
! Compute actual time step for each region. NOTE take into account CFL number 
! and user-specified minimum time step
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => regions(iReg)

    IF ( pRegion%mixtInput%frozenFlag .EQV. .FALSE. ) THEN 
      pRegion%dtMinImpulse = pRegion%plagInput%cfl*pRegion%dtMinImpulse       
      pRegion%dtMinImpulse = MIN(global%dtImposed,pRegion%dtMinImpulse)    
    ELSE 
      pRegion%dtMinImpulse = global%dtImposed   
    END IF ! pRegion%mixtInput%frozenFlag
  END DO ! iReg

! ******************************************************************************
! Determine minimum time step over all processors
! ******************************************************************************

! ==============================================================================
! Allocate temporary memory
! ==============================================================================

  ALLOCATE(localValsReal(0:global%nRegions),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'localValsReal')
  END IF ! global%error

  ALLOCATE(globalValsReal(0:global%nRegions),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'globalValsReal')
  END IF ! global%error

! ==============================================================================
! Peform reduction operation. NOTE need to include region index 0 to make sure
! this works properly for serial runs.
! ==============================================================================

  DO iReg = 0,global%nRegions
    localValsReal(iReg) = HUGE(1.0_RFREAL)
  END DO ! iReg

  DO iReg = 1,global%nRegionsLocal
    pRegion => regions(iReg)

    localValsReal(pRegion%iRegionGlobal) = pRegion%dtMinImpulse
  END DO ! iReg    

  CALL MPI_Allreduce(localValsReal(0:global%nRegions), &
                     globalValsReal(0:global%nRegions),global%nRegions+1, &
                     MPI_RFREAL,MPI_MIN,global%mpiComm,errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
  END IF ! global%errorFlag
  
! ==============================================================================
! Find smallest time step and store as global time step
! ==============================================================================

  dtMinImpulse = HUGE(1.0_RFREAL) 

  DO iReg = 0,global%nRegions
    IF ( globalValsReal(iReg) < dtMinImpulse ) THEN
      dtMinImpulse = globalValsReal(iReg)
      iRegDtMin = iReg
    END IF ! globalValsReal 
  END DO ! iReg

  DO iReg = 1,global%nRegionsLocal
    regions(iReg)%global%dtMinImpulse = dtMinImpulse
  END DO ! iReg

! ==============================================================================
! Deallocate temporary memory
! ==============================================================================

  DEALLOCATE(localValsReal,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'localValsReal')
  END IF ! global%error 

  DEALLOCATE(globalValsReal,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'globalValsReal')
  END IF ! global%error  

#ifndef GENX
! ******************************************************************************
! Store minimum time step as system time step for coupled simulations
! ******************************************************************************

  global%dTimeSystem = global%dtMinImpulse
#endif

! ******************************************************************************
! If smallest time step is smaller than user-specified time step and smaller 
! than user-specified limit, determine in which cell this occurs. NOTE the 
! region in which the smallest time step occurs has already been determined 
! above.
! ******************************************************************************

  IF ( global%dtMinImpulse < global%dtImposed .AND. & 
       global%dtMinImpulse < global%dtMinLimit ) THEN 

! ==============================================================================
!   Allocate temporary memory
! ==============================================================================

    ALLOCATE(localValsInt(0:global%nRegions),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'localValsInt')
    END IF ! global%error

    ALLOCATE(globalValsInt(0:global%nRegions),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'globalValsInt')
    END IF ! global%error
    
! ==============================================================================
!   Peform reduction operations. NOTE need to include region index 0 to make 
!   sure this works properly for serial runs.
! ==============================================================================
  
    DO iReg = 0,global%nRegions
      localValsInt(iReg) = -HUGE(1)
    END DO ! iReg  
    
    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)      

      localValsInt(pRegion%iRegionGlobal) = pRegion%dtMinImpulseLoc
    END DO ! iReg      

    CALL MPI_Allreduce(localValsInt(0:global%nRegions), &
                       globalValsInt(0:global%nRegions),global%nRegions+1, &
                       MPI_INTEGER,MPI_MAX,global%mpiComm,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
    END IF ! global%errorFlag

    IF ( global%myProcid == MASTERPROC ) THEN 
      iCellDtMin = globalValsInt(iRegDtMin)
    END IF ! global%myProcid
    
! ==============================================================================
!   Deallocate temporary memory
! ==============================================================================

    DEALLOCATE(localValsInt,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'localValsInt')
    END IF ! global%error 
  
    DEALLOCATE(globalValsInt,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'globalValsInt')
    END IF ! global%error  

! ==============================================================================
!   Print information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing time step information...'
      WRITE(STDOUT,'(A,3X,A,1X,E16.9)') SOLVER_NAME,'Smallest time step:', & 
                                        global%dtMinImpulse
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Location of smallest time step:'
      WRITE(STDOUT,'(A,5X,A,4X,I6)') SOLVER_NAME,'Region:',iRegDtMin
      WRITE(STDOUT,'(A,5X,A,1X,I9)') SOLVER_NAME,'Cell:  ',iCellDtMin
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Printing time step information done.'
    END IF ! global%myProcid
  END IF ! global%dtMinImpulse

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_MinimumTimeStepImpulse

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_MinimumTimeStepImpulse.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
! ******************************************************************************
