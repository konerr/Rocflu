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
! Purpose: Collection of routines related to Hypre Solver implementation. 
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModHypre.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2006-2007 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModHypre

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModMPI

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModHypre.F90,v $'

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: RFLU_HYPRE_AssembleMatrixVector, &
            RFLU_HYPRE_CreateObjects, &
            RFLU_HYPRE_DestroyObjects, &
            RFLU_HYPRE_InitializeMatrixVector, &
            RFLU_HYPRE_SolveMatrix

! ==============================================================================
! Private functions
! ==============================================================================
  

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  





! ******************************************************************************
!
! Purpose: Assemble Hypre matrix and vector. 
!
! Description: None.
!
! Input:
!   regions            Pointer to region data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_HYPRE_AssembleMatrixVector(regions)

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
    
    INTEGER :: iErrHypre 
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => regions(1)%global
 
    CALL RegisterFunction(global,'RFLU_HYPRE_AssembleMatrixVector',__FILE__)

! ==============================================================================
!   Assemble matrix and vector
! ==============================================================================

#ifdef HYPRE
    CALL HYPRE_IJMatrixAssemble(regions(1)%AHypre,iErrHypre)
    CALL HYPRE_IJVectorAssemble(regions(1)%RhsHypre,iErrHypre)
#endif   
  
! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HYPRE_AssembleMatrixVector











! ******************************************************************************
!
! Purpose: Create Hypre matrix, vector and solver. 
!
! Description: None.
!
! Input:
!   regions            Pointer to region data
!
! Notes: 
!   1. All the local regions on current processor contribute to one block
!      of rows.
!   2. Matrix and vector structure (non-zeros locations) would remain same
!      during whole computation. So this subroutine is called only once 
!      for all time stepping cycles.
!   3. Solver can be choosen differently in each time step but for now it is
!      kept same for all time steps.
! ******************************************************************************

  SUBROUTINE RFLU_HYPRE_CreateObjects(regions)

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
    
    INTEGER :: c,c1,c2,colId,curPos,diagCount,dummy,errorFlag,globalCellno,i, &
               icg,ifg,iLowerHypre,iReg,iReg2,iErrHypre,iRow,iUpperHypre,j, &
               maxIter,nEntries,procOffset,sizeHypre
    INTEGER, ALLOCATABLE, DIMENSION(:) :: cols,diagSizesHypre,growsHypre, &
                                          istrt,ncols,offdSizesHypre,ncolsHypre
    REAL(RFREAL) :: amgTol
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_region), POINTER :: pRegion
    TYPE(t_global), POINTER :: global  
     
! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => regions(1)%global
 
    CALL RegisterFunction(global,'RFLU_HYPRE_CreateObjects',__FILE__)

! ******************************************************************************
!   Find number of rows on currect processor
! ******************************************************************************

    sizeHypre = 0

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid

      sizeHypre = sizeHypre + pGrid%nCells
    END DO ! iReg

! ==============================================================================
!   Initialize HYPRE variables
!   Note: 
!     1. For serial runs, iLower = 1 and iUpper = nCells
!     2. For parallel runs, iLower(myProcId) = iUpper(myProcId-1)+1  
! ==============================================================================

    ALLOCATE(ncolsHypre(sizeHypre),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'ncolsHypre')
    END IF ! global%error

    ALLOCATE(diagSizesHypre(sizeHypre),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'diagSizesHypre')
    END IF ! global%error
    
    ALLOCATE(offdSizesHypre(sizeHypre),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'offdSizesHypre')
    END IF ! global%error
    
    procOffset  = regions(1)%nCellsOffset(regions(1)%iRegionGlobal)
    iLowerHypre = 1 + procOffset
    iUpperHypre = iLowerHypre + sizeHypre - 1

! ==============================================================================
!   Build CSR format, Set diagonal and off-diagonal sizes 
! ==============================================================================

    DO icg = 1,sizeHypre
      ncolsHypre(icg) = 1
      diagSizesHypre(icg) = 1
      offdSizesHypre(icg) = 0
    END DO ! icg

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid

      iReg2 = pRegion%iRegionGlobal

! --- Loop over actual faces ---------------------------------------------------

      DO ifg = 1,pGrid%nFaces-pGrid%nFacesAV
        c1 = pGrid%f2c(1,ifg)
        c2 = pGrid%f2c(2,ifg)

        ncolsHypre(c1+pRegion%nCellsOffset(iReg2)-procOffset) &
         = ncolsHypre(c1+pRegion%nCellsOffset(iReg2)-procOffset) + 1
        ncolsHypre(c2+pRegion%nCellsOffset(iReg2)-procOffset) &
         = ncolsHypre(c2+pRegion%nCellsOffset(iReg2)-procOffset) + 1

        diagSizesHypre(c1+pRegion%nCellsOffset(iReg2)-procOffset) &
         = diagSizesHypre(c1+pRegion%nCellsOffset(iReg2)-procOffset) + 1
        diagSizesHypre(c2+pRegion%nCellsOffset(iReg2)-procOffset) &
         = diagSizesHypre(c2+pRegion%nCellsOffset(iReg2)-procOffset) + 1
      END DO ! ifg

! --- Loop over actual to virtual faces ----------------------------------------

      DO ifg = pGrid%nFaces-pGrid%nFacesAV+1,pGrid%nFaces
        c1 = pGrid%f2c(1,ifg)
        c2 = pGrid%f2c(2,ifg)

        IF ( c1 > pGrid%nCells ) THEN ! c1 is virtual
          dummy = c1
          c1    = c2
          c2    = dummy
        END IF ! c1

! ----- c1 is actual cell and c2 is virtual cell -------------------------------

        ncolsHypre(c1+pRegion%nCellsOffset(iReg2)-procOffset) &
        = ncolsHypre(c1+pRegion%nCellsOffset(iReg2)-procOffset) + 1

! ----- Determine if c2 lies on same processor or different processor ----------

        globalCellno = pGrid%virt2GlobalIds(c2-pGrid%nCells)

        IF ( (globalCellno .LE. procOffset) .OR. &
             (globalCellno .GT. procOffset + sizeHypre) ) THEN
          offdSizesHypre(c1+pRegion%nCellsOffset(iReg2)-procOffset) &
          = offdSizesHypre(c1+pRegion%nCellsOffset(iReg2)-procOffset) + 1
        ELSE
          diagSizesHypre(c1+pRegion%nCellsOffset(iReg2)-procOffset) &
          = diagSizesHypre(c1+pRegion%nCellsOffset(iReg2)-procOffset) + 1
        END IF ! globalCellno
      END DO ! ifg
    END DO ! iReg

! ==============================================================================
!   Initialize Hypre Matrices
! ==============================================================================

#ifdef HYPRE
    CALL HYPRE_IJMatrixCreate(MPI_COMM_WORLD,iLowerHypre,iUpperHypre,&
                              iLowerHypre,iUpperHypre,regions(1)%AHypre, &          
                              iErrHypre)                              
    IF ( iErrHypre /= 0 ) THEN 
      WRITE(*,*) 'HYPRE Error in IJMatrixCreate'
      CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                     'RFLU_HYPRE_CreateObjects')
    END IF ! iErrHypre

    CALL HYPRE_IJMatrixSetObjectType(regions(1)%AHypre,HYPRE_PARCSR,iErrHypre)
    IF ( iErrHypre /= 0 ) THEN 
      WRITE(*,*) 'HYPRE Error in IJMatrixSetObjectType'
      CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                     'RFLU_HYPRE_CreateObjects')
    END IF ! iErrHypre

    CALL HYPRE_IJMatrixSetRowSizes(regions(1)%AHypre,nColsHypre,iErrHypre)
    IF ( iErrHypre /= 0 ) THEN 
      WRITE(*,*) 'HYPRE Error in IJMatrixSetRowSizes'
      CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                     'RFLU_HYPRE_CreateObjects')
    END IF ! iErrHypre
    
    CALL HYPRE_IJMatrixSetDiagOffdSizes(regions(1)%AHypre,diagSizesHypre, &
                                        offdSizesHypre,iErrHypre)
    IF ( iErrHypre /= 0 ) THEN 
      WRITE(*,*) 'HYPRE Error in IJMatrixSetDiagOffdSizes'
      CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                     'RFLU_HYPRE_CreateObjects')
    END IF ! iErrHypre
    
    CALL HYPRE_IJMatrixInitialize(regions(1)%AHypre,iErrHypre)
    IF ( iErrHypre /= 0 ) THEN 
      WRITE(*,*) 'HYPRE Error in IJMatrixInitialize'
      CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                     'RFLU_HYPRE_CreateObjects')
    END IF ! iErrHypre
      
    CALL HYPRE_IJMatrixGetObject(regions(1)%AHypre,regions(1)%parAHypre, &
                                 iErrHypre)
    IF ( iErrHypre /= 0 ) THEN 
      WRITE(*,*) 'HYPRE Error in IJMatrixGetObject'
      CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                     'RFLU_HYPRE_CreateObjects')
    END IF ! iErrHypre
#endif    

! ==============================================================================
!   Initialize Hypre RHS
! ==============================================================================

#ifdef HYPRE
    CALL HYPRE_IJVectorCreate(MPI_COMM_WORLD,iLowerHypre,iUpperHypre, &
                              regions(1)%RhsHypre,iErrHypre)
    IF ( iErrHypre /= 0 ) THEN
      WRITE(*,*) 'HYPRE Error in IJVectorCreate'
      CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                     'RFLU_HYPRE_CreateObjects')
    END IF ! iErrHypre
  
    CALL HYPRE_IJVectorSetObjectType(regions(1)%RhsHypre,HYPRE_PARCSR,iErrHypre)
    IF ( iErrHypre /= 0 ) THEN
      WRITE(*,*) 'HYPRE Error in IJVectorSetObjectType'
      CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                     'RFLU_HYPRE_CreateObjects')
    END IF ! iErrHypre

    CALL HYPRE_IJVectorInitialize(regions(1)%RhsHypre,iErrHypre)
    IF ( iErrHypre /= 0 ) THEN
      WRITE(*,*) 'HYPRE Error in IJVectorInitialize'
      CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                     'RFLU_HYPRE_CreateObjects')
    END IF ! iErrHypre
  
    CALL HYPRE_IJVectorGetObject(regions(1)%RhsHypre,regions(1)%parRhsHypre, &
                                 iErrHypre)
    IF ( iErrHypre /= 0 ) THEN
      WRITE(*,*) 'HYPRE Error in IJVectorGetObject'
      CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                     'RFLU_HYPRE_CreateObjects')
    END IF ! iErrHypre
#endif  
  
! ==============================================================================
!   Initialize Hypre Solution
! ==============================================================================
 
#ifdef HYPRE 
    CALL HYPRE_IJVectorCreate(MPI_COMM_WORLD,iLowerHypre,iUpperHypre, &
                              regions(1)%SolHypre,iErrHypre)
    IF ( iErrHypre /= 0 ) THEN
      WRITE(*,*) 'HYPRE Error in IJVectorCreate-SolHypre'
      CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                     'RFLU_HYPRE_CreateObjects')
    END IF ! iErrHypre
    
    CALL HYPRE_IJVectorSetObjectType(regions(1)%SolHypre,HYPRE_PARCSR,iErrHypre)
    IF ( iErrHypre /= 0 ) THEN
      WRITE(*,*) 'HYPRE Error in IJVectorSetObjectType-SolHypre'
      CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                     'RFLU_HYPRE_CreateObjects')
    END IF ! iErrHypre

    CALL HYPRE_IJVectorInitialize(regions(1)%SolHypre,iErrHypre)
    IF ( iErrHypre /= 0 ) THEN
      WRITE(*,*) 'HYPRE Error in IJVectorInitialize-SolHypre'
      CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                     'RFLU_HYPRE_CreateObjects')
    END IF ! iErrHypre
  
    CALL HYPRE_IJVectorGetObject(regions(1)%SolHypre,regions(1)%parSolHypre, &
                                 iErrHypre) 
    IF ( iErrHypre /= 0 ) THEN
      WRITE(*,*) 'HYPRE Error in IJVectorGetObject-SolHypre'
      CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                     'RFLU_HYPRE_CreateObjects')
    END IF ! iErrHypre
#endif

! ==============================================================================
!   Initialize Hypre Solver
! ==============================================================================

#ifdef HYPRE
    SELECT CASE( global%hypreSolver )
      CASE ( HYPRE_SOLV_BAMG )
        CALL HYPRE_BoomerAMGCreate(regions(1)%hypreSolver,iErrHypre) 

        amgTol = global%tolHypre
        CALL HYPRE_BoomerAMGSetTol(regions(1)%hypreSolver,amgTol,iErrHypre)

        maxIter = global%maxIterHypre
        CALL HYPRE_BoomerAMGSetMaxIter(regions(1)%hypreSolver,maxIter,iErrHypre)
      CASE ( HYPRE_SOLV_GMRES )
        CALL HYPRE_ParCSRGMRESCreate(global%mpiComm,regions(1)%hypreSolver,iErrHypre) 

        amgTol = global%tolHypre
        CALL HYPRE_ParCSRGMRESSetTol(regions(1)%hypreSolver,amgTol,iErrHypre)

        maxIter = global%maxIterHypre
        CALL HYPRE_ParCSRGMRESSetMaxIter(regions(1)%hypreSolver,maxIter,iErrHypre)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_HYPRE_SOLVER_INVALID,__LINE__)
    END SELECT ! pPatch%bcType
#endif

! ==============================================================================
!   Following are various hypre routines which can be used in future to set
!   many other parameters which are left default for now.
! ==============================================================================

!    CALL HYPRE_BoomerAMGSetCoarsenType(hypreSolver,(hybrid*coarsenType), &
!                                       iErrHypre)
!    CALL HYPRE_BoomerAMGSetMeasureType(hypreSolver,measureType,iErrHypre)
!    CALL HYPRE_BoomerAMGSetStrongThreshold(hypreSolver,strongThresh,iErrHypre)
!    CALL HYPRE_BoomerAMGSetTruncFactor(hypreSolver,truncFact,iErrHypre)
!    CALL HYPRE_BoomerAMGSetCycleType(hypreSolver,cycType, iErrHypre)
!    CALL HYPRE_BoomerAMGSetNumGridSweeps(hypreSolver, nGridSweeps,iErrHypre)
!    CALL HYPRE_BoomerAMGSetGridRelaxType(hypreSolver, gridRelaxType,iErrHypre)
!    CALL HYPRE_BoomerAMGSetRelaxWeight(hypreSolver, relaxWt,iErrHypre)
!    CALL HYPRE_BoomerAMGSetOmega(hypreSolver, omega,iErrHypre)
!    CALL HYPRE_BoomerAMGSetSmoothType(hypreSolver, smoothType,iErrHypre)
!    CALL HYPRE_BoomerAMGSetSmoothNumLevels(hypreSolver, smoothNLev,iErrHypre)
!    CALL HYPRE_BoomerAMGSetSmoothNumSweeps(hypreSolver, smoothNSweep,iErrHypre)
!    CALL HYPRE_BoomerAMGSetGridRelaxPoints(hypreSolver, grdRelaxPts,iErrHypre)
!    CALL HYPRE_BoomerAMGSetMaxLevels(hypreSolver, maxLev,iErrHypre)
!    CALL HYPRE_BoomerAMGSetMaxRowSum(hypreSolver, maxRowSum,iErrHypre)
!    CALL HYPRE_BoomerAMGSetNumFunctions(hypreSolver, numFunc,iErrHypre)
!    CALL HYPRE_BoomerAMGSetVariant(hypreSolver, variant,iErrHypre)
!    CALL HYPRE_BoomerAMGSetOverlap(hypreSolver, overlap,iErrHypre)
!    CALL HYPRE_BoomerAMGSetDomainType(hypreSolver, domainType,iErrHypre)

! ==============================================================================
!   Copy Hypre objects to other local regions
! ==============================================================================

    IF ( global%nRegionsLocal > 1 ) THEN
      DO iReg = 2,global%nRegionsLocal
        pRegion => regions(iReg)

        pRegion%AHypre      = regions(1)%AHypre
        pRegion%parAHypre   = regions(1)%parAHypre
        pRegion%RhsHypre    = regions(1)%RhsHypre
        pRegion%parRhsHypre = regions(1)%parRhsHypre
        pRegion%SolHypre    = regions(1)%SolHypre
        pRegion%parSolHypre = regions(1)%parSolHypre
      END DO ! iReg
    END IF ! global%nRegionLocal

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HYPRE_CreateObjects









! ******************************************************************************
!
! Purpose: Destroy Hypre matrix,vector and solver objects. 
!
! Description: None.
!
! Input:
!   regions            Pointer to region data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_HYPRE_DestroyObjects(regions)

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
    
    INTEGER :: c1,c2,iErrHypre,ifg,iRow
    REAL(RFREAL) :: zeroValue 
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => regions(1)%global
 
    CALL RegisterFunction(global,'RFLU_HYPRE_DestroyObjects',__FILE__)

! ==============================================================================
!   Destroy Hypre matrix,vectors and solver.
! ==============================================================================

#ifdef HYPRE
    CALL HYPRE_IJMatrixDestroy(regions(1)%AHypre,iErrHypre)
    IF ( iErrHypre /= 0 ) THEN 
      WRITE(*,*) 'HYPRE Error in IJMatrixDestroy-AHypre'
      CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                     'RFLU_HYPRE_DestroyObjects')
    END IF ! iErrHypre
    
    CALL HYPRE_IJVectorDestroy(regions(1)%RhsHypre,iErrHypre)
    IF ( iErrHypre /= 0 ) THEN
      WRITE(*,*) 'HYPRE Error in IJVectorDestroy-RhsHypre'
      CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                     'RFLU_HYPRE_DestroyObjects')
    END IF ! iErrHypre
  
    CALL HYPRE_IJVectorDestroy(regions(1)%SolHypre,iErrHypre)
    IF ( iErrHypre /= 0 ) THEN
      WRITE(*,*) 'HYPRE Error in IJVectorDestroy-SolHypre'
      CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                     'RFLU_HYPRE_DestroyObjects')
    END IF ! iErrHypre
  
    SELECT CASE( global%hypreSolver )
      CASE ( HYPRE_SOLV_BAMG )
        CALL HYPRE_BoomerAMGDestroy(regions(1)%hypreSolver,iErrHypre) 
        IF ( iErrHypre /= 0 ) THEN
          WRITE(*,*) 'HYPRE Error in HYPRE_BoomerAMGDestroy-hypreSolver'
          CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                         'RFLU_HYPRE_DestroyObjects')
        END IF ! iErrHypre
      CASE ( HYPRE_SOLV_GMRES )
      CALL HYPRE_ParCSRGMRESDestroy(regions(1)%hypreSolver,iErrHypre) 
        IF ( iErrHypre /= 0 ) THEN
          WRITE(*,*) 'HYPRE Error in HYPRE_ParCSRGMRESDestroy-hypreSolver'
          CALL ErrorStop(global,ERR_HYPRE_OUTPUT,__LINE__, &
                       'RFLU_HYPRE_DestroyObjects')
        END IF ! iErrHypre
      CASE DEFAULT
        CALL ErrorStop(global,ERR_HYPRE_SOLVER_INVALID,__LINE__)
    END SELECT ! pPatch%bcType
#endif  
  
! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HYPRE_DestroyObjects









! ******************************************************************************
!
! Purpose: Initialize Hypre matrix and vector. 
!
! Description: None.
!
! Input:
!   regions            Pointer to region data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_HYPRE_InitializeMatrixVector(regions)

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
    
    INTEGER :: c1,c2,dummy,iColGlobal,iErrHypre,ifg,iReg,iReg1,iRow, &
               iRowGlobal,procOffset,sizeHypre
    REAL(RFREAL) :: zeroValue 
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => regions(1)%global
 
    CALL RegisterFunction(global,'RFLU_HYPRE_InitializeMatrixVector',__FILE__)

    zeroValue = 0.0_RFREAL
     
! *****************************************************************************
!   Find number of rows on currect processor
! *****************************************************************************

    sizeHypre = 0

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid

      sizeHypre = sizeHypre + pGrid%nCells
    END DO ! iReg

    procOffset  = regions(1)%nCellsOffset(regions(1)%iRegionGlobal)

! ==============================================================================
!   Initialization of Hypre matrix and vectors already done
!   Set matrix and vector to Zero before assembling them
! ==============================================================================

    DO iRow = 1,sizeHypre
      iRowGlobal = iRow + procOffset
      iColGlobal = iRowGlobal
      
#ifdef HYPRE      
      CALL HYPRE_IJMatrixSetValues(pRegion%AHypre,1,1,iRowGlobal,iColGlobal, &
                                   zeroValue,iErrHypre)

      CALL HYPRE_IJVectorSetValues(pRegion%RhsHypre,1,iRowGlobal,zeroValue, &
                                   iErrHypre)
#endif                                                                      
    END DO ! iRow

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid

      iReg1 = pRegion%iRegionGlobal

! --- Loop over actual faces ---------------------------------------------------
      DO ifg = 1,pGrid%nFaces-pGrid%nFacesAV
        c1 = pGrid%f2c(1,ifg)
        c2 = pGrid%f2c(2,ifg)

        iRowGlobal = c1 + pRegion%nCellsOffset(iReg1)
        iColGlobal = c2 + pRegion%nCellsOffset(iReg1)
        
#ifdef HYPRE        
        CALL HYPRE_IJMatrixSetValues(pRegion%AHypre,1,1,iRowGlobal,iColGlobal, &
                                     zeroValue,iErrHypre)
#endif

        iRowGlobal = c2 + pRegion%nCellsOffset(iReg1)
        iColGlobal = c1 + pRegion%nCellsOffset(iReg1) 
        
#ifdef HYPRE        
        CALL HYPRE_IJMatrixSetValues(pRegion%AHypre,1,1,iRowGlobal,iColGlobal, &
                                     zeroValue,iErrHypre)
#endif                                    
      END DO ! ifg

! --- Loop over actual to virtual faces ----------------------------------------
      DO ifg = pGrid%nFaces-pGrid%nFacesAV+1,pGrid%nFaces
        c1 = pGrid%f2c(1,ifg)
        c2 = pGrid%f2c(2,ifg)

        IF ( c1 .GT. pGrid%nCells ) THEN ! c1 is virtual
          dummy = c1
          c1    = c2
          c2    = dummy
        END IF ! c1

        iRowGlobal = c1 + pRegion%nCellsOffset(iReg1)
        iColGlobal = pGrid%virt2GlobalIds(c2 - pGrid%nCells)
        
#ifdef HYPRE        
        CALL HYPRE_IJMatrixSetValues(pRegion%AHypre,1,1,iRowGlobal,iColGlobal, &
                                     zeroValue,iErrHypre)
#endif                                     
      END DO ! ifg
    END DO ! iReg

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HYPRE_InitializeMatrixVector






! ******************************************************************************
!
! Purpose: Solve Hypre matrix. 
!
! Description: None.
!
! Input:
!   regions            Pointer to region data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_HYPRE_SolveMatrix(regions)

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
    
    INTEGER :: iErrHypre,nIterHypre
    REAL(RFREAL) :: resNormHypre 
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => regions(1)%global
 
    CALL RegisterFunction(global,'RFLU_HYPRE_SolveMatrix',__FILE__)

! ==============================================================================
!   Solve system of equation using HYPRE 
! ==============================================================================

#ifdef HYPRE
    SELECT CASE( global%hypreSolver )
      CASE ( HYPRE_SOLV_BAMG )
        CALL HYPRE_BoomerAMGSetup(regions(1)%hypreSolver,regions(1)%parAHypre, &
                                  regions(1)%parRhsHypre,regions(1)%parSolHypre, &
                                  iErrHypre)

        CALL HYPRE_BoomerAMGSolve(regions(1)%hypreSolver,regions(1)%parAHypre, &
                                  regions(1)%parRhsHypre,regions(1)%parSolHypre, &
                                  iErrHypre)
  
        CALL HYPRE_BoomerAMGGetNumIterations(regions(1)%hypreSolver,nIterHypre, &
                                             iErrHypre)
      CASE ( HYPRE_SOLV_GMRES )
        CALL HYPRE_ParCSRGMRESSetup(regions(1)%hypreSolver,regions(1)%parAHypre, &
                                    regions(1)%parRhsHypre,regions(1)%parSolHypre, &
                                    iErrHypre)

        CALL HYPRE_ParCSRGMRESSolve(regions(1)%hypreSolver,regions(1)%parAHypre, &
                                    regions(1)%parRhsHypre,regions(1)%parSolHypre, &
                                    iErrHypre)

!        CALL HYPRE_ParCSRGMRESGetNumIterations(regions(1)%hypreSolver,nIterHypre, &
!                                               iErrHypre)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_HYPRE_SOLVER_INVALID,__LINE__)
    END SELECT ! pPatch%bcType
#endif

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HYPRE_SolveMatrix





END MODULE RFLU_ModHypre

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModHypre.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:36  brollin
! New Stable version
!
! Revision 1.6  2010/03/14 23:52:32  mparmar
! Added option for Hypre GMRES solver
!
! Revision 1.5  2008/12/06 08:43:41  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:16:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/12/04 20:02:22  haselbac
! Enclosed HYPRE calls within ifdef HYPRE sections
!
! Revision 1.2  2007/11/29 02:03:16  mparmar
! Removed ^M at the end of lines
!
! Revision 1.1  2007/11/28 23:04:45  mparmar
! Initial revision
!
! ******************************************************************************

