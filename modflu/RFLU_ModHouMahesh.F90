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
! Purpose: Collection of routines for non-dissipative implicit scheme by
!          Hou and Mahesh. 
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModHouMahesh.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2006-2007 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModHouMahesh

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModError
  USE ModMPI

  USE RFLU_ModDifferentiationCells, ONLY: RFLU_ComputeGradCellsGGScalar, &
                                          RFLU_ComputeGradCellsGGVector, &
                                          RFLU_AXI_ComputeGradCellsGGScalar, &
                                          RFLU_AXI_ComputeGradCellsGGVector
  USE RFLU_ModHypre
  USE RFLU_ModMovingFrame, ONLY: RFLU_MVF_ComputeAcceleration, &
                                 RFLU_MVF_SetVelocity, &
                                 RFLU_MVF_UpdateBC

  USE RFLU_ModRelatedPatches, ONLY: RFLU_RELP_TransformWrapper  
  
  USE ModInterfaces, ONLY: ReflectVector, &
                           RFLU_SetVarsWrapper, &
                           MixtPerf_HM_T_DGMP

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModHouMahesh.F90,v $'

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: RFLU_HM_ComputeGlobalDelP, &
            RFLU_HM_ConvCvD2ND, &
            RFLU_HM_ConvCvND2D, &
            RFLU_HM_ConvDvD2ND, &
            RFLU_HM_ConvDvND2D, &
            RFLU_HM_ConvGridCoordD2ND, &
            RFLU_HM_ConvTvD2ND, &
            RFLU_HM_ConvTvND2D, &
            RFLU_HM_PredictFaceNormalVelocity, &  
            RFLU_HM_PredCorrMP, &  
            RFLU_HM_TestSolver

! ==============================================================================
! Private functions
! ==============================================================================
  
  PRIVATE :: RFLU_HM_CorrectFaceNormalVelocity, &
             RFLU_HM_CorrectVelocity, &
             RFLU_HM_InitPredCorrMP, &
             RFLU_HM_SolveContinuityEq, &
             RFLU_HM_SolveEnergyEq, &
             RFLU_HM_SolveMomentumEq, &
             RFLU_HM_SolveMomentumEqWrapper, &
             RFLU_HM_UpdatePressure

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  




! ******************************************************************************
!
! Purpose: Compute local and global maximum delP.
!
! Description: None.
!
! Input:
!   regions             Region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_HM_ComputeGlobalDelP(regions,maxDelP)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), DIMENSION(:), POINTER :: regions
    REAL(RFREAL), INTENT(OUT) :: maxDelP 
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,i,ipos,iReg,iRegionGlobal,nVals
    REAL(RFREAL) :: maxValue,minValue,value
    REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: globalVals,localVals
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_region), POINTER :: pRegion
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => regions(1)%global
    
    CALL RegisterFunction(global,'RFLU_HM_ComputeGlobalDelP',__FILE__)

    nvals = global%nRegions

! ******************************************************************************
!   Check for serial run 
! ******************************************************************************

    IF ( nVals == 1 ) THEN
      ipos = 0
      maxDelP = 0.0_RFREAL
      DO i=1,regions(1)%grid%nCells
        IF ( maxDelP < ABS(regions(1)%mixt%delP(i)) ) THEN
          maxDelP = ABS(regions(1)%mixt%delP(i))
          ipos  = i
        END IF ! value      
      END DO ! i
    ELSE

! ******************************************************************************
!     Compute maxDelP over all regions
!     Allocate temporary memory
! ******************************************************************************
        
      ALLOCATE(globalVals(nvals),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'globalVals')
      END IF ! global%error        
        
      ALLOCATE(localVals(nvals),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'localVals')
      END IF ! global%error

! ******************************************************************************
!     Compute global max delP 
! ******************************************************************************

      DO iReg = 1,global%nRegions
        globalVals(iReg) = 0.0_RFREAL

        localVals(iReg) = 0.0_RFREAL      
      END DO ! iPatch

! ==============================================================================
!     Set local max delP
! ==============================================================================
        
      DO iReg = 1,global%nRegionsLocal   
        pRegion => regions(iReg)
        pGrid   => pRegion%grid
        
        iRegionGlobal = pRegion%iRegionGlobal
          
        ipos = 0 
        value = 0.0_RFREAL 
        DO i=1,pGrid%nCells
          IF ( value < ABS(pRegion%mixt%delP(i)) ) THEN
            value = pRegion%mixt%delP(i)
            ipos  = i
          END IF ! value      
        END DO ! i

        localVals(iRegionGlobal) = value  
      END DO ! iReg

! ==============================================================================
!     Compute global max delp 
! ==============================================================================

      CALL MPI_AllReduce(localVals,globalVals,nvals,MPI_RFREAL,MPI_SUM, &
                         global%mpiComm,errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
      END IF ! global%error

! ==============================================================================
!     Find max delp
! ==============================================================================

      maxDelP = 0.0
      DO iReg = 1,global%nRegions
     
        IF ( maxDelP < globalVals(iReg) ) THEN
          maxDelP = globalVals(iReg)
        END IF ! maxDelP
      END DO ! iReg
    
! ******************************************************************************
!     Deallocate temporary memory
! ******************************************************************************
        
      DEALLOCATE(globalVals,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'globalVals')
      END IF ! global%error        
        
      DEALLOCATE(localVals,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'localVals')
      END IF ! global%error

    END IF ! nVals

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_ComputeGlobalDelP









! ******************************************************************************
!
! Purpose: Non-dimensionalize flow variables.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!
! Notes: None.
!
! ******************************************************************************
   
  SUBROUTINE RFLU_HM_ConvCvD2ND(pRegion)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,icg
    REAL(RFREAL) :: irRef,iru2Ref,iuRef,pRef
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pCvOld
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid
 
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'RFLU_HM_ConvCvD2ND',__FILE__)

    pGrid => pRegion%grid

    pCv    => pRegion%mixt%cv
    pCvOld => pRegion%mixt%cvOld

    irRef   = 1.0_RFREAL/global%refDensity
    pRef    = global%refPressure
    iuRef   = 1.0_RFREAL/global%refVelocity
    iru2Ref = irRef*iuRef*iuRef 

! ==============================================================================
!   Loop over cells and non-dimensionalize flow variables 
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot
      pCv(CV_MIXT_DENS,icg) = pCv(CV_MIXT_DENS,icg)*irRef
      pCv(CV_MIXT_XVEL,icg) = pCv(CV_MIXT_XVEL,icg)*iuRef
      pCv(CV_MIXT_YVEL,icg) = pCv(CV_MIXT_YVEL,icg)*iuRef
      pCv(CV_MIXT_ZVEL,icg) = pCv(CV_MIXT_ZVEL,icg)*iuRef
      pCv(CV_MIXT_PRES,icg) = (pCv(CV_MIXT_PRES,icg)-pRef)*iru2Ref

      pCvOld(CV_MIXT_DENS,icg) = pCvOld(CV_MIXT_DENS,icg)*irRef
      pCvOld(CV_MIXT_PRES,icg) = (pCvOld(CV_MIXT_PRES,icg)-pRef)*iru2Ref
    END DO ! icg

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_ConvCvD2ND










! ******************************************************************************
!
! Purpose: Dimensionalize flow variables.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!
! Notes: None.
!
! ******************************************************************************
   
  SUBROUTINE RFLU_HM_ConvCvND2D(pRegion)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,icg
    REAL(RFREAL) :: rRef,ru2Ref,pRef,uRef 
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pCvOld
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid
 
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'RFLU_HM_ConvCvND2D',__FILE__)

    pGrid  => pRegion%grid
    pCv    => pRegion%mixt%cv
    pCvOld => pRegion%mixt%cvOld

    rRef   = global%refDensity
    pRef   = global%refPressure
    uRef   = global%refVelocity
    ru2Ref = rRef*uRef*uRef

! ==============================================================================
!   Loop over cells and dimensionalize flow variables 
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot
      pCv(CV_MIXT_DENS,icg) = pCv(CV_MIXT_DENS,icg)*rRef
      pCv(CV_MIXT_XVEL,icg) = pCv(CV_MIXT_XVEL,icg)*uRef
      pCv(CV_MIXT_YVEL,icg) = pCv(CV_MIXT_YVEL,icg)*uRef
      pCv(CV_MIXT_ZVEL,icg) = pCv(CV_MIXT_ZVEL,icg)*uRef
      pCv(CV_MIXT_PRES,icg) = (pCv(CV_MIXT_PRES,icg)*ru2Ref)+pRef

      pCvOld(CV_MIXT_DENS,icg) = pCvOld(CV_MIXT_DENS,icg)*rRef
      pCvOld(CV_MIXT_PRES,icg) = (pCvOld(CV_MIXT_PRES,icg)*ru2Ref)+pRef
    END DO ! icg

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_ConvCvND2D










! ******************************************************************************
!
! Purpose: Non-dimensionalize dependent flow variables.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!
! Notes: None.
!
! ******************************************************************************
   
  SUBROUTINE RFLU_HM_ConvDvD2ND(pRegion)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,icg
    REAL(RFREAL) :: iRefT
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pDv
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid
 
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'RFLU_HM_ConvDvD2ND',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Non-dimensionalizing dv...'
    END IF ! global%verbLevel

    pGrid => pRegion%grid
    pDv   => pRegion%mixt%dv

    iRefT = 1.0_RFREAL/global%refTemperature

! ==============================================================================
!   Loop over cells and non-dimensionalize flow variables 
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot
      pDv(DV_MIXT_TEMP,icg) = iRefT*pDv(DV_MIXT_TEMP,icg) 
    END DO ! icg

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Non-dimensionalizing dv done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_ConvDvD2ND










! ******************************************************************************
!
! Purpose: Dimensionalize dependent flow variables.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!
! Notes: None.
!
! ******************************************************************************
   
  SUBROUTINE RFLU_HM_ConvDvND2D(pRegion)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,icg
    REAL(RFREAL) :: refT
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pDv
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid
 
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'RFLU_HM_ConvDvND2D',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Dimensionalizing dv...'
    END IF ! global%verbLevel

    pGrid => pRegion%grid
    pDv   => pRegion%mixt%dv

    refT = global%refTemperature

! ==============================================================================
!   Loop over cells and dimensionalize flow variables 
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot
      pDv(DV_MIXT_TEMP,icg) = refT*pDv(DV_MIXT_TEMP,icg)
    END DO ! icg

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Dimensionalizing dv done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_ConvDvND2D










! ******************************************************************************
!
! Purpose: Non-dimensionalize Grid Coordinates.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!
! Notes: None.
!
! ******************************************************************************
   
  SUBROUTINE RFLU_HM_ConvGridCoordD2ND(pRegion)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,ivg
    REAL(RFREAL) :: iRefLength
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid
 
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************

    global => pRegion%global
    pGrid  => pRegion%grid
 
    CALL RegisterFunction(global,'RFLU_HM_ConvGridCoordD2ND',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Non-dimensionalizing grid...'
    END IF ! global%verbLevel

    iRefLength = 1.0_RFREAL/global%refLength

! ==============================================================================
!   Loop over vertices and non-dimensionalize grid coordinates 
! ==============================================================================

    DO ivg = 1,pGrid%nVert
      pGrid%xyz(XCOORD,ivg) = pGrid%xyz(XCOORD,ivg)*iRefLength
      pGrid%xyz(YCOORD,ivg) = pGrid%xyz(YCOORD,ivg)*iRefLength
      pGrid%xyz(ZCOORD,ivg) = pGrid%xyz(ZCOORD,ivg)*iRefLength
    END DO ! icg

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Non-dimensionalizing flow done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_ConvGridCoordD2ND









! ******************************************************************************
!
! Purpose: Non-dimensionalize transport flow variables.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!
! Notes: None.
!
! ******************************************************************************
   
  SUBROUTINE RFLU_HM_ConvTvD2ND(pRegion)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,icg
    REAL(RFREAL) :: iRefKappa,iRefVisc
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pTv
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid
 
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'RFLU_HM_ConvTvD2ND',__FILE__)

    pGrid => pRegion%grid
    pTv   => pRegion%mixt%tv

    iRefVisc  = 1.0_RFREAL/global%refVisc
    iRefKappa = global%prLam/(global%refCp*global%refVisc)

! ==============================================================================
!   Loop over cells and non-dimensionalize flow variables 
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot
      pTv(TV_MIXT_MUEL,icg) =  iRefVisc*pTv(TV_MIXT_MUEL,icg) 
      pTv(TV_MIXT_TCOL,icg) = iRefKappa*pTv(TV_MIXT_TCOL,icg) 
    END DO ! icg

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_ConvTvD2ND










! ******************************************************************************
!
! Purpose: Dimensionalize transport flow variables.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!
! Notes: None.
!
! ******************************************************************************
   
  SUBROUTINE RFLU_HM_ConvTvND2D(pRegion)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,icg
    REAL(RFREAL) :: refKappa,refVisc
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pTv
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid
 
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'RFLU_HM_ConvTvND2D',__FILE__)

    pGrid => pRegion%grid
    pTv   => pRegion%mixt%tv

    refVisc  = global%refVisc
    refKappa = global%refCp*global%refVisc/global%prLam

! ==============================================================================
!   Loop over cells and dimensionalize flow variables 
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot
      pTv(TV_MIXT_MUEL,icg) =  refVisc*pTv(TV_MIXT_MUEL,icg) 
      pTv(TV_MIXT_TCOL,icg) = refKappa*pTv(TV_MIXT_TCOL,icg) 
    END DO ! icg

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_ConvTvND2D










! ******************************************************************************
!
! Purpose: Correct face-normal velocity with pressure-correction.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!
! Notes: None.
!
! ******************************************************************************
   
  SUBROUTINE RFLU_HM_CorrectFaceNormalVelocity(pRegion)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: c1,c2,errorFlag,ifg
    REAL(RFREAL) :: ddelPdn,deln,delT,rho_k_n1 
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pCvOld
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid
 
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************

    global => pRegion%global
    pGrid  => pRegion%grid
    pCv    => pRegion%mixt%cv
    pCvOld => pRegion%mixt%cvOld
 
    CALL RegisterFunction(global,'RFLU_HM_CorrectFaceNormalVelocity',__FILE__)

    delT = global%dtMin*global%refVelocity/global%refLength

! ==============================================================================
!   Loop over interior faces and correct face-normal velocity, 
!   Note : no need to loop over boundary faces bcoz boundary face-normal
!   velocity is zero.
! ==============================================================================

    DO ifg = 1,pGrid%nFaces
      c1 = pGrid%f2c(1,ifg)
      c2 = pGrid%f2c(2,ifg)

      deln = ( (pGrid%cofg(XCOORD,c2)-pGrid%cofg(XCOORD,c1)) &
              *(pGrid%cofg(XCOORD,c2)-pGrid%cofg(XCOORD,c1)) &
             + (pGrid%cofg(YCOORD,c2)-pGrid%cofg(YCOORD,c1)) &
              *(pGrid%cofg(YCOORD,c2)-pGrid%cofg(YCOORD,c1)) &
             + (pGrid%cofg(ZCOORD,c2)-pGrid%cofg(ZCOORD,c1)) &
              *(pGrid%cofg(ZCOORD,c2)-pGrid%cofg(ZCOORD,c1)) )**(0.5_RFREAL)

      rho_k_n1 = ( pCv(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c1) &
                 + pCv(CV_MIXT_DENS,c2)+pCvOld(CV_MIXT_DENS,c2) )/4.0_RFREAL

      ddelPdn = (pRegion%mixt%delP(c2) - pRegion%mixt%delP(c1))/deln

      pRegion%mixt%vfMixt(ifg) = pRegion%mixt%vfMixt(ifg) &
                               - delT*ddelPdn/(4.0_RFREAL*rho_k_n1)
    END DO ! ifg

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_CorrectFaceNormalVelocity










! ******************************************************************************
!
! Purpose: Correct velocity with pressure-correction.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!
! Output:
!
! Notes: None.
!
! ******************************************************************************
   
  SUBROUTINE RFLU_HM_CorrectVelocity(pRegion)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,flag,icg
    REAL(RFREAL) :: delT,rho_n1 
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pCvOld,var
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
    TYPE(t_global), POINTER :: global  
 
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************

    global => pRegion%global
    pCv    => pRegion%mixt%cv
    pCvOld => pRegion%mixt%cvOld
 
    CALL RegisterFunction(global,'RFLU_HM_CorrectVelocity',__FILE__)

    delT = global%dtMin*global%refVelocity/global%refLength

! ==============================================================================
!   Allocate temporary memory allocated for delp gradient computation 
! ==============================================================================

    ALLOCATE(grad(XCOORD:ZCOORD,1,pRegion%grid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'grad')
    END IF ! global%error

    ALLOCATE(var(1,pRegion%grid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'var')
    END IF ! global%error

! ==============================================================================
!   call gradient computing routine for delp
! ==============================================================================

    DO icg = 1,pRegion%grid%nCellsTot
      var(1,icg) = pRegion%mixt%delP(icg)
    END DO ! icg

    flag = 3
    IF ( pRegion%mixtInput%axiFlag .EQV. .TRUE. ) THEN
      CALL RFLU_AXI_ComputeGradCellsGGScalar(pRegion,1,1,1,1,var,grad,flag)
    ELSE
      CALL RFLU_ComputeGradCellsGGScalar(pRegion,1,1,1,1,var,grad,flag)
    END IF ! pRegion%mixtInput%axiFlag


! ==============================================================================
!   Loop over cells and add velocity corrections 
! ==============================================================================

    DO icg = 1,pRegion%grid%nCells
      rho_n1 = 0.5_RFREAL*(pCv(CV_MIXT_DENS,icg)+pCvOld(CV_MIXT_DENS,icg)) 

      pCv(CV_MIXT_XVEL,icg) = pCv(CV_MIXT_XVEL,icg) &
                            - delT*grad(XCOORD,1,icg)/(4.0_RFREAL*rho_n1)
      pCv(CV_MIXT_YVEL,icg) = pCv(CV_MIXT_YVEL,icg) &
                            - delT*grad(YCOORD,1,icg)/(4.0_RFREAL*rho_n1)
      pCv(CV_MIXT_ZVEL,icg) = pCv(CV_MIXT_ZVEL,icg) &
                            - delT*grad(ZCOORD,1,icg)/(4.0_RFREAL*rho_n1)
    END DO ! icg

! ==============================================================================
!   Deallocate temporary memory allocated for delp gradient computation 
! ==============================================================================

    DEALLOCATE(grad,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'grad')
    END IF ! global%error

    DEALLOCATE(var,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'var')
    END IF ! global%error

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_CorrectVelocity








! ******************************************************************************
!
! Purpose: Initialize first iteration (q=1) for predictor-corrector scheme
!          by Hou-Mahesh. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!
! Output:
!
! Notes: None.
!
! ******************************************************************************
   
  SUBROUTINE RFLU_HM_InitPredCorrMP(pRegion)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,icg,ifg,ifl,iPatch,iVar,iVarBeg,iVarEnd
    REAL(RFREAL), DIMENSION(:), POINTER :: pVfMixt,pVfMixtOld
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pCvOld,pCvOld2
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGradCell,pGradCellOld, &
                                               pGradCellOld2
    TYPE(t_global), POINTER :: global  
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_grid), POINTER :: pGrid  
 
! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'RFLU_HM_InitPredCorrMP',__FILE__)

    pGrid         => pRegion%grid
    pVfMixt       => pRegion%mixt%vfMixt
    pVfMixtOld    => pRegion%mixt%vfMixtOld
    pGradCell     => pRegion%mixt%gradCell
    pGradCellOld  => pRegion%mixt%gradCellOld
    pGradCellOld2 => pRegion%mixt%gradCellOld2
    pCv           => pRegion%mixt%cv
    pCvOld        => pRegion%mixt%cvOld
    pCvOld2       => pRegion%mixt%cvOld2

! ==============================================================================
!   Loop over cells and initialize cell variables (cv, grad)
! ==============================================================================

! ------------------------------------------------------------------------------
!  Initialize cv variables 
! ------------------------------------------------------------------------------

    DO icg = 1,pGrid%nCellsTot
      pCvOld2(CV_OLD2_MIXT_DENS,icg) = pCvOld(CV_MIXT_DENS,icg)
      pCvOld2(CV_OLD2_MIXT_PRES,icg) = pCvOld(CV_MIXT_PRES,icg)

      pCvOld(CV_MIXT_DENS,icg) = pCv(CV_MIXT_DENS,icg)
      pCvOld(CV_MIXT_XVEL,icg) = pCv(CV_MIXT_XVEL,icg)
      pCvOld(CV_MIXT_YVEL,icg) = pCv(CV_MIXT_YVEL,icg)
      pCvOld(CV_MIXT_ZVEL,icg) = pCv(CV_MIXT_ZVEL,icg)
      pCvOld(CV_MIXT_PRES,icg) = pCv(CV_MIXT_PRES,icg)
    END DO ! icg

! ------------------------------------------------------------------------------
!  Initialize gradient variables 
! ------------------------------------------------------------------------------

    iVar = GRC_MIXT_PRES
    DO icg = 1,pGrid%nCellsTot
      pGradCellOld2(XCOORD,iVar,icg) = pGradCellOld(XCOORD,iVar,icg)
      pGradCellOld2(YCOORD,iVar,icg) = pGradCellOld(YCOORD,iVar,icg)
      pGradCellOld2(ZCOORD,iVar,icg) = pGradCellOld(ZCOORD,iVar,icg)
    END DO ! icg

    IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN
      iVarBeg = GRC_MIXT_XVEL
      iVarEnd = GRC_MIXT_PRES
    ELSE
      iVarBeg = GRC_MIXT_PRES
      iVarEnd = GRC_MIXT_PRES
    END IF ! pRegion%mixtInput%flowModel

    DO iVar=iVarBeg,iVarEnd
      DO icg = 1,pGrid%nCellsTot
        pGradCellOld(XCOORD,iVar,icg) = pGradCell(XCOORD,iVar,icg)
        pGradCellOld(YCOORD,iVar,icg) = pGradCell(YCOORD,iVar,icg)
        pGradCellOld(ZCOORD,iVar,icg) = pGradCell(ZCOORD,iVar,icg)
      END DO ! icg
    END DO ! iVar

! ==============================================================================
!   Loop over faces and initialize face variables
! ==============================================================================

    DO ifg = 1,pGrid%nFaces
      pVfMixtOld(ifg) = pVfMixt(ifg) 
    END DO ! ifg

! ==============================================================================
!   Compute temperature from cv (earlier cvOld was used for Temp)
! ==============================================================================

    CALL RFLU_SetVarsWrapper(pRegion,1,pRegion%grid%nCellsTot)

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_InitPredCorrMP











! ******************************************************************************
!
! Purpose: Calculate solution at a new time level.
!
! Description: The governing equations are integrated in time using
!              predictor corrector scheme proposed by Hou and Mahesh.
!
! Input: 
!   regions        Data of all regions.
!
! Notes: 
!   1. Routines to communicate variables at region boundaries needs information
!      of which variable to communicate. Not all variables needed to be 
!      communicated at a time in this predictor-corrector solver unlike in 
!      dissipative solver. Information about which variable to communicate is 
!      passed through variable 'iVar' which can have values from 1 to 9 which 
!      has following meaning,
!      1          density from cv
!      2,3,4      x,y,z velocities from cv
!      5          pressure from cv
!      6          delp
!      7          density from cvOld
!      8          pressure from cvold
!      9          gradients
!
! ******************************************************************************

  SUBROUTINE RFLU_HM_PredCorrMP( regions )

    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE ModGlobal, ONLY     : t_global
    USE ModError
    USE ModParameters
    USE ModMPI
    USE RFLU_ModMPI

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,flag,flowModel,icg,iReg,iRegLocal,iter,iVar, &
               MAXITERATION
    REAL(RFREAL) :: maxDelP,minDelP,tolerance 
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad 
    TYPE(t_global), POINTER :: global
    TYPE(t_region), POINTER :: pRegion
 
! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => regions(1)%global

    CALL RegisterFunction( global,'RFLU_HM_PredCorrMP',__FILE__ )

! =============================================================================
!   Initialize solution for first iteration 
! =============================================================================

    DO iReg=1,global%nRegionsLocal
      pRegion => regions(iReg)

      CALL RFLU_HM_InitPredCorrMP(pRegion)
    END DO ! iRegLocal

! =============================================================================
! Compute particle accelerations and velocities 
! =============================================================================

    IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
      CALL RFLU_MVF_ComputeAcceleration(regions)

      DO iRegLocal=1,global%nRegionsLocal
        iReg = iRegLocal

        pRegion => regions(iReg)

        IF ( global%mvfAccFlag .EQV. .TRUE. ) THEN
          CALL RFLU_MVF_SetVelocity(pRegion)
        END IF ! global%mvfAccFlag

        CALL RFLU_MVF_UpdateBC(pRegion)
      END DO ! iRegLocal
    END IF ! global%mvFrameFlag

! =============================================================================
!   Test hypre solver: Useful routine to test hypre solver 
! =============================================================================
!
!    CALL RFLU_HM_TestSolver(regions)

! =============================================================================
!   Loop over iterations and regions 
! =============================================================================

    tolerance = global%tolPC 
    MAXITERATION = global%maxIterPC

    iter = 0

    DO WHILE ( iter < MAXITERATION ) 
      iter = iter + 1

! -----------------------------------------------------------------------------
!     Solve continuity equation
! -----------------------------------------------------------------------------

      CALL RFLU_HM_SolveContinuityEq(regions)

! === Update virtual cells ====================================================

      iVar = 1
      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_MPI_HM_ISendWrapper(pRegion,iVar)
      END DO ! iReg

      CALL RFLU_MPI_HM_CopyWrapper(regions,iVar)
    
      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_MPI_HM_RecvWrapper(pRegion,iVar)
      END DO ! iReg

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_MPI_ClearRequestWrapper(pRegion)
      END DO ! iReg

! -----------------------------------------------------------------------------
!     Compute pressure gradient for n+3/2,q
! -----------------------------------------------------------------------------

      DO iReg=1,global%nRegionsLocal
        pRegion => regions(iReg)

        flag = 1
        IF ( pRegion%mixtInput%axiFlag .EQV. .TRUE. ) THEN
          CALL RFLU_AXI_ComputeGradCellsGGScalar(pRegion,CV_MIXT_PRES, &
                                                 CV_MIXT_PRES,GRC_MIXT_PRES, &
                                                 GRC_MIXT_PRES,pRegion%mixt%cv, &
                                                 pRegion%mixt%gradCell,flag)
        ELSE
          CALL RFLU_ComputeGradCellsGGScalar(pRegion,CV_MIXT_PRES,CV_MIXT_PRES, &
                                             GRC_MIXT_PRES,GRC_MIXT_PRES, &
                                             pRegion%mixt%cv, &
                                             pRegion%mixt%gradCell,flag)
        END IF ! pRegion%mixtInput%axiFlag
      END DO ! iReg

! -----------------------------------------------------------------------------
!     Solve momentum equations
! -----------------------------------------------------------------------------

      CALL RFLU_HM_SolveMomentumEqWrapper(regions)

! === Update virtual cells ====================================================

      DO iVar = 2,4
        DO iReg = 1,global%nRegionsLocal
          pRegion => regions(iReg)

          CALL RFLU_MPI_HM_ISendWrapper(pRegion,iVar)
        END DO ! iReg

        CALL RFLU_MPI_HM_CopyWrapper(regions,iVar)
    
        DO iReg = 1,global%nRegionsLocal
          pRegion => regions(iReg)

          CALL RFLU_MPI_HM_RecvWrapper(pRegion,iVar)
        END DO ! iReg

        DO iReg = 1,global%nRegionsLocal
          pRegion => regions(iReg)

          CALL RFLU_MPI_ClearRequestWrapper(pRegion)
        END DO ! iReg
      END DO ! iVar

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_RELP_TransformWrapper(pRegion)
      END DO ! iReg

! -----------------------------------------------------------------------------
!     Predict face-normal velocity
! -----------------------------------------------------------------------------

      DO iReg=1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_HM_PredictFaceNormalVelocity(pRegion,pRegion%mixt%cv, &
                                               pRegion%mixt%vfMixt)
      END DO ! iReg

! -----------------------------------------------------------------------------
!     Solve energy equation
! -----------------------------------------------------------------------------

      CALL RFLU_HM_SolveEnergyEq(regions)

! === Update virtual cells ====================================================

      iVar = 6
      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_MPI_HM_ISendWrapper(pRegion,iVar)
      END DO ! iReg

      CALL RFLU_MPI_HM_CopyWrapper(regions,iVar)
    
      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_MPI_HM_RecvWrapper(pRegion,iVar)
      END DO ! iReg

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_MPI_ClearRequestWrapper(pRegion)
      END DO ! iReg

! -----------------------------------------------------------------------------
!     Check convergence, compute max delp
! -----------------------------------------------------------------------------

      CALL RFLU_HM_ComputeGlobalDelP(regions,maxDelP)

      IF ( iter >= MAXITERATION ) THEN
        IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
          WRITE(*,*) 'Convergence not reached upto q-level=',iter, &
                     ' maxDelP=',maxDelP
          CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__, &
                         'RFLU_HM_PredCorrMP')
        END IF ! global
      ELSE
        IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
          WRITE(*,*) 'q-level=',iter,' maxDelP=',maxDelP
        END IF ! global
      END IF

      IF ( maxDelP < tolerance ) THEN
        EXIT
      END IF ! maxDelP

! -----------------------------------------------------------------------------
!     Corrector step
! -----------------------------------------------------------------------------

      DO iReg=1,global%nRegionsLocal
        pRegion => regions(iReg)

! -----------------------------------------------------------------------------
!       Update pressure
! -----------------------------------------------------------------------------

        CALL RFLU_HM_UpdatePressure(pRegion)

! -----------------------------------------------------------------------------
!       Correct velocities
! -----------------------------------------------------------------------------

        CALL RFLU_HM_CorrectVelocity(pRegion)
      END DO ! iReg

! === Update virtual cells ====================================================

      DO iVar = 2,5
        DO iReg = 1,global%nRegionsLocal
          pRegion => regions(iReg)

          CALL RFLU_MPI_HM_ISendWrapper(pRegion,iVar)
        END DO ! iReg

        CALL RFLU_MPI_HM_CopyWrapper(regions,iVar)
    
        DO iReg = 1,global%nRegionsLocal
          pRegion => regions(iReg)

          CALL RFLU_MPI_HM_RecvWrapper(pRegion,iVar)
        END DO ! iReg

        DO iReg = 1,global%nRegionsLocal
          pRegion => regions(iReg)

          CALL RFLU_MPI_ClearRequestWrapper(pRegion)
        END DO ! iReg
      END DO ! iVar

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_RELP_TransformWrapper(pRegion)
      END DO ! iReg

! -----------------------------------------------------------------------------
!     Correct face-normal velocities
! -----------------------------------------------------------------------------

      DO iReg=1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_HM_CorrectFaceNormalVelocity(pRegion)
      END DO ! iReg

! -----------------------------------------------------------------------------
!     Compute velocity gradient for n+1,q+1
! -----------------------------------------------------------------------------

      DO iReg=1,global%nRegionsLocal
        pRegion => regions(iReg)

        IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN
          IF ( pRegion%mixtInput%axiFlag .EQV. .TRUE. ) THEN
            CALL RFLU_AXI_ComputeGradCellsGGVector(pRegion,CV_MIXT_XVEL, &
                                                   CV_MIXT_ZVEL,GRC_MIXT_XVEL, &
                                                   GRC_MIXT_ZVEL, &
                                                   pRegion%mixt%cv, &
                                                   pRegion%mixt%gradCell)
          ELSE
            CALL RFLU_ComputeGradCellsGGVector(pRegion,CV_MIXT_XVEL, &
                                               CV_MIXT_ZVEL,GRC_MIXT_XVEL, &
                                               GRC_MIXT_ZVEL,pRegion%mixt%cv, &
                                               pRegion%mixt%gradCell)
          END IF ! pRegion%mixtInput%axiFlag
        END IF ! pRegion%mixtInput%flowModel
      END DO ! iReg

! === Update virtual cells ====================================================

      IF ( regions(1)%mixtInput%flowModel == FLOW_NAVST ) THEN
        iVar = 9
        DO iReg = 1,global%nRegionsLocal
          pRegion => regions(iReg)

          CALL RFLU_MPI_HM_ISendWrapper(pRegion,iVar)
        END DO ! iReg

        CALL RFLU_MPI_HM_CopyWrapper(regions,iVar)
    
        DO iReg = 1,global%nRegionsLocal
          pRegion => regions(iReg)

          CALL RFLU_MPI_HM_RecvWrapper(pRegion,iVar)
        END DO ! iReg

        DO iReg = 1,global%nRegionsLocal
          pRegion => regions(iReg)

          CALL RFLU_MPI_ClearRequestWrapper(pRegion)
        END DO ! iReg
      END IF ! pRegion%mixtInput%flowModel

! =============================================================================
!   End loop over iterations
! =============================================================================

    END DO ! max(delP) .OR. iter

! ******************************************************************************
!   Finalize
! ******************************************************************************

    CALL DeregisterFunction( global )

  END SUBROUTINE RFLU_HM_PredCorrMP








! ******************************************************************************
!
! Purpose: Compute face-normal velocity from cell centered velocities.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   pCv         Pointer to cv,cvOld, face-normal is computed from velocities
!               in either cv or cvOld.
!
! Output:
!   pVfMixt     Face-normal velocity.
!
! Notes: None.
!
! ******************************************************************************
   
  SUBROUTINE RFLU_HM_PredictFaceNormalVelocity(pRegion,pCv,pVfMixt)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), DIMENSION(:), POINTER :: pVfMixt
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: c1,c2,errorFlag,ifg
    REAL(RFREAL) :: nm,nx,ny,nz
    TYPE(t_global), POINTER :: global  
 
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'RFLU_HM_PredictFaceNormalVelocity',__FILE__)

! ==============================================================================
!   Loop over interior faces and compute face-normal velocity
! ==============================================================================

    DO ifg = 1,pRegion%grid%nFaces
      c1 = pRegion%grid%f2c(1,ifg)
      c2 = pRegion%grid%f2c(2,ifg)

      nx = pRegion%grid%fn(XCOORD,ifg)
      ny = pRegion%grid%fn(YCOORD,ifg)
      nz = pRegion%grid%fn(ZCOORD,ifg)
      nm = pRegion%grid%fn(XYZMAG,ifg)

      pVfMixt(ifg) = 0.5_RFREAL* &
                   ( (pCv(CV_MIXT_XVEL,c1) + pCv(CV_MIXT_XVEL,c2))*nx &
                   + (pCv(CV_MIXT_YVEL,c1) + pCv(CV_MIXT_YVEL,c2))*ny &
                   + (pCv(CV_MIXT_ZVEL,c1) + pCv(CV_MIXT_ZVEL,c2))*nz )
    END DO ! ifg

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_PredictFaceNormalVelocity









! ******************************************************************************
!
!   Purpose: Solve continuity equation in non-dissipative scheme by Hou and 
!            Mahesh.
!
!   Description: None.
!
!   Input:
!     regions            Pointer to region data
!
!   Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_HM_SolveContinuityEq(regions)

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
    
    INTEGER :: c1,c2,distrib,dummy,errorFlag,icg,iColGlobal,iErrHypre,ifg,ifl, &
               iPatch,iRow,iRowGlobal,iReg,iReg1,nCellsOffsetRegion
    REAL(RFREAL) :: aRef,delT,direction,iDt,nm,refM2,value,vfnmi4,voliDt,volume
    REAL(RFREAL) :: aoa,aos,mf,minj,nx,ny,nz,uo,vo,wo,ro,po,sigma, &
                    SpeedOfSound,tf,VelX,VelY,VelZ
    REAL(RFREAL), DIMENSION(:), POINTER :: pVfMixt
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pCvOld
    TYPE(t_global), POINTER :: global 
    TYPE(t_patch), POINTER :: pPatch 
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => regions(1)%global
 
    CALL RegisterFunction(global,'RFLU_HM_SolveContinuityEq',__FILE__)

    delT      = global%dtMin*global%refVelocity/global%refLength
    iDt       = 1.0_RFREAL/delT

! *****************************************************************************
!   Initialize Hypre matrix and vectors
! *****************************************************************************

    CALL RFLU_HYPRE_InitializeMatrixVector(regions)

! *****************************************************************************
!   Loop over regions and assign pointers
! *****************************************************************************

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid

      pVfMixt => pRegion%mixt%vfMixt
      pCv     => pRegion%mixt%cv
      pCvOld  => pRegion%mixt%cvOld

      iReg1 = pRegion%iRegionGlobal
      nCellsOffsetRegion = pRegion%nCellsOffset(iReg1) 

! ==============================================================================
!     Loop over cells and initialize the coefficients
! ==============================================================================

      DO icg = 1,pGrid%nCells
        volume = pGrid%vol(icg)

        voliDt = volume*iDt

        value = voliDt
        iRowGlobal = icg + nCellsOffsetRegion 
        iColGlobal = iRowGlobal
        
#ifdef HYPRE        
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
#endif

        value = pCvOld(CV_MIXT_DENS,icg)*voliDt
        iRowGlobal = icg + nCellsOffsetRegion
        
#ifdef HYPRE           
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif                                       
      END DO ! icg

! ==============================================================================
!     Loop over interior faces and accumulate coefficients
! ==============================================================================

      DO ifg = 1,pGrid%nFaces-pGrid%nFacesAV
        c1 = pGrid%f2c(1,ifg)
        c2 = pGrid%f2c(2,ifg)

        nm = pGrid%fn(XYZMAG,ifg)

        vfnmi4 = pVfMixt(ifg)*nm/4.0_RFREAL

        value =  vfnmi4
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = iRowGlobal       
#ifdef HYPRE           
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
#endif
        value = -vfnmi4
        iRowGlobal = c2 + nCellsOffsetRegion
        iColGlobal = iRowGlobal        
#ifdef HYPRE        
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
#endif
        value =  vfnmi4 
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = c2 + nCellsOffsetRegion
#ifdef HYPRE           
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
#endif
        value = -vfnmi4
        iRowGlobal = c2 + nCellsOffsetRegion
        iColGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE           
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
#endif
        value = - (pCvOld(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c2))*vfnmi4
        iRowGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE           
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif
        value = + (pCvOld(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c2))*vfnmi4
        iRowGlobal = c2 + nCellsOffsetRegion
#ifdef HYPRE           
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif                                       
      END DO ! ifg


      DO ifg = pGrid%nFaces-pGrid%nFacesAV+1,pGrid%nFaces
        c1 = pGrid%f2c(1,ifg)
        c2 = pGrid%f2c(2,ifg)

        direction = 1.0_RFREAL

        IF ( c1 .GT. pGrid%nCells ) THEN ! c1 is virtual
          dummy = c1
          c1    = c2
          c2    = dummy

          direction = -1.0_RFREAL
        END IF

        nm = pGrid%fn(XYZMAG,ifg)

        vfnmi4 = direction*pVfMixt(ifg)*nm/4.0_RFREAL

        value =  vfnmi4 
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = iRowGlobal
#ifdef HYPRE   
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
#endif
        value =  vfnmi4 
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = pGrid%virt2GlobalIds(c2 - pGrid%nCells)
#ifdef HYPRE   
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
#endif
        value = - (pCvOld(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c2))*vfnmi4
        iRowGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE   
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif                                       
      END DO ! ifg

! ==============================================================================
!     Loop over boundary faces. Continuity eqn.
! ==============================================================================

      DO iPatch = 1,pRegion%grid%nPatches
        pPatch => pRegion%patches(iPatch)
        
        distrib = pPatch%mixt%distrib
 
        SELECT CASE( pPatch%bcType )
        CASE (BC_SLIPWALL,BC_NOSLIPWALL_HFLUX,BC_NOSLIPWALL_TEMP)
        CASE (BC_OUTFLOW)
          DO ifg = 1,pPatch%nBFaces
            c1 = pPatch%bf2c(ifg)

            nx = pPatch%fn(XCOORD,ifg)
            ny = pPatch%fn(YCOORD,ifg)
            nz = pPatch%fn(ZCOORD,ifg)
            nm = pPatch%fn(XYZMAG,ifg)

            vfnmi4 = (pCv(CV_MIXT_XVEL,c1)*nx+pCv(CV_MIXT_YVEL,c1)*ny &
                     +pCv(CV_MIXT_ZVEL,c1)*nz)*nm/4.0_RFREAL

            value =  vfnmi4
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = iRowGlobal       
#ifdef HYPRE           
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)
#endif
            value =  vfnmi4 
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE           
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)
#endif
            value = - (pCvOld(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c1))*vfnmi4
            iRowGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE           
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                           iErrHypre)
#endif
          END DO ! ifg

        CASE (BC_FARFIELD)
          DO ifg = 1,pPatch%nBFaces
            c1 = pPatch%bf2c(ifg)
    
            nx = pPatch%fn(XCOORD,ifg)
            ny = pPatch%fn(YCOORD,ifg)
            nz = pPatch%fn(ZCOORD,ifg)
            nm = pPatch%fn(XYZMAG,ifg)

            mf  = pPatch%mixt%vals(BCDAT_FARF_MACH,distrib*ifg)
            aoa = pPatch%mixt%vals(BCDAT_FARF_ATTACK,distrib*ifg)
            aos = pPatch%mixt%vals(BCDAT_FARF_SLIP,distrib*ifg)
            tf  = pPatch%mixt%vals(BCDAT_FARF_TEMP,distrib*ifg)

            SpeedOfSound = SQRT(tf/global%refTemperature)

            VelX = mf*SpeedOfSound*COS(aoa)*COS(aos)
            VelY = mf*SpeedOfSound*SIN(aoa)*COS(aos)
            VelZ = mf*SpeedOfSound*SIN(aos)

            vfnmi4 = (nx*VelX+ny*VelY+nz*VelZ)*nm/4.0_RFREAL 

            value =  vfnmi4 
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = iRowGlobal
#ifdef HYPRE   
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)
#endif                                       

            value = - vfnmi4*pRegion%mixt%cv(CV_MIXT_DENS,c1)
            iRowGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE   
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                           iErrHypre)
#endif                                       

            value = - (pCvOld(CV_MIXT_DENS,c1)+ &
                       pRegion%mixt%cvOld(CV_MIXT_DENS,c1))*vfnmi4
            iRowGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE   
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                           iErrHypre)
#endif                                       
          END DO ! ifg

        CASE (BC_INFLOW_VELTEMP)
          DO ifg = 1,pPatch%nBFaces
            c1 = pPatch%bf2c(ifg)

            nm = pPatch%fn(XYZMAG,ifg)

            minj = -pPatch%mixt%vals(BCDAT_INJECT_MFRATE,distrib*ifg)
            minj = minj/(global%refDensity*global%refVelocity)

            value = -minj*nm
            iRowGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE   
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif                                       
          END DO ! ifg

        CASE (BC_INJECTION)
          DO ifg = 1,pPatch%nBFaces
            c1 = pPatch%bf2c(ifg)

            nm = pPatch%fn(XYZMAG,ifg)

            minj = -pPatch%mixt%vals(BCDAT_INJECT_MFRATE,distrib*ifg)
            minj = minj/(global%refDensity*global%refVelocity)

            value = -minj*nm
            iRowGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE   
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif                                       
          END DO ! ifg

        CASE (BC_PERIODIC,BC_SYMMETRY,BC_VIRTUAL)
        CASE DEFAULT                
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pPatch%bcType
      END DO ! iPatch

! ==============================================================================
!     Compute coefficients due to absorbing boundary condition. Continuity eqn. 
! ==============================================================================

      IF ( (global%abcFlag .EQV. .TRUE.) .AND. (global%abcKind == 0) ) THEN
        DO icg = 1,pGrid%nCells
          distrib = global%abcDistrib
          ro = pRegion%mixt%cvRef(CV_MIXT_DENS,distrib*icg)
          uo = pRegion%mixt%cvRef(CV_MIXT_XVEL,distrib*icg)
          vo = pRegion%mixt%cvRef(CV_MIXT_YVEL,distrib*icg)
          wo = pRegion%mixt%cvRef(CV_MIXT_ZVEL,distrib*icg)
          po = pRegion%mixt%cvRef(CV_MIXT_PRES,distrib*icg)

          volume = pGrid%vol(icg)

          sigma = pRegion%mixt%sigma(icg)*global%refLength/global%refVelocity 

          value = sigma*volume*0.5_RFREAL
          iRowGlobal = icg + nCellsOffsetRegion 
          iColGlobal = iRowGlobal
#ifdef HYPRE        
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)
#endif

          value = -(pCvOld(CV_MIXT_DENS,icg)*0.5_RFREAL-ro)*sigma*volume
          iRowGlobal = icg + nCellsOffsetRegion
#ifdef HYPRE           
          CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                         iErrHypre)
#endif                                       
        END DO ! icg
      END IF ! global%abcFlag

    END DO ! iReg

! ==============================================================================
!   Assemble matrix and vector
! ==============================================================================

    CALL RFLU_HYPRE_AssembleMatrixVector(regions)

! ==============================================================================
!   Solve system of equation for density using HYPRE 
! ==============================================================================

    CALL RFLU_HYPRE_SolveMatrix(regions)

! ==============================================================================
!   Extract solution from HYPRE solution vector parSolHypre
! ==============================================================================

    global => regions(1)%global

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid

      pCv     => pRegion%mixt%cv

      iReg1 = pRegion%iRegionGlobal
      nCellsOffsetRegion = pRegion%nCellsOffset(iReg1) 

! --- Extract actual cell data -------------------------------------------------
      DO icg = 1,pGrid%nCells
        iRowGlobal = icg + nCellsOffsetRegion  
#ifdef HYPRE           
        CALL HYPRE_IJVectorGetValues(regions(1)%SolHypre,1,iRowGlobal,value, &
                                     iErrHypre)
#endif                                     
        pCv(CV_MIXT_DENS,icg) = value
      END DO ! icg
    END DO ! iReg

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_SolveContinuityEq










! ******************************************************************************
!
! Purpose: Solve energy equation in non-dissipative scheme by Hou and Mahesh.
!
! Description: None.
!
! Input:
!   regions            Pointer to region data
!
! Notes:
! 1. Divergence of velocity at no-slip wall is same in axisymmetric computation
!    and in 2D computation as velocity is zero at no-slip wall. The case of
!    accelerating particle is simulated by transforming governing equations in
!    in moving reference frame attached to particle, thus velocity at no-slip
!    is again zero. Hence, friction coefficients at no-slip wall are same in
!    axisymmetric computation as in 2D computations. 
!
! ******************************************************************************

  SUBROUTINE RFLU_HM_SolveEnergyEq(regions)

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
    
    INTEGER :: c1,c2,distrib,dummy,errorFlag,icg,iColGlobal,iErrHypre,ifg, &
               iPatch,iRow,iRowGlobal,iReg,iReg1,nCellsOffsetRegion,numIter
    REAL(RFREAL), PARAMETER :: ONE_THIRD = 1.0_RFREAL/3.0_RFREAL
    REAL(RFREAL), PARAMETER :: TWO_THIRD = 2.0_RFREAL/3.0_RFREAL
    REAL(RFREAL), PARAMETER :: EPS = 1.0E-15_RFREAL
    REAL(RFREAL) :: aRef,aini,aiui_nHalf,deln,delT,direction,du1dx1,du1dx2, &
                    du1dx3,du2dx1,du2dx2,du2dx3,du3dx1,du3dx2,du3dx3,iDt, &
                    iRefMu,iRefReNum,mu,mVoliDt,nx,ny,nz,nm,Pr,pRef,p_tilda, &
                    pi_k,qn,refM,refM2,rho_c1_n1,rho_c2_n1,rho_k,rho_n, &
                    rho_k_n1,rho_n1,refGamma,relRes,rRef,t1,t1D,t1N,t1T,t2, &
                    t2D,t2N,t2T,t3,t3D,t3N,t3T,Temp1,Temp1Old,Temp2,Temp2Old, &
                    term1,term2,term3,Twall,uc1,uc1Old,uc2,uc2Old,uic1_ni, &
                    uic2_ni,u1,u2,u3,uRef,uu_k_nHalf,uu_n,uu_n1,uysk,v_nHalf, &
                    value,vc1,vc1Old,vc2,vc2Old,wc1,wc1Old,wc2,wc2Old,volume, &
                    xc,yc,zc
    REAL(RFREAL) :: aoa,aos,mf,minj,tf,tinj,uo,vo,wo,ro,po,sigma,SpeedOfsound, &
                    uu_np1,sigma_c1,sigma_c2,uini_c1,uini_c2,VelX,VelY,VelZ
    REAL(RFREAL), DIMENSION(:), POINTER :: pVfMixt,pVfMixtOld
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pCvOld,pCvOld2,pDv,pGradP, &
                                             pGradPOld,pGradPOld2,pGv,pTv
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGradCell,pGradCellOld, &
                                               pGradCellOld2
    TYPE(t_region), POINTER :: pRegion
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => regions(1)%global
 
    CALL RegisterFunction(global,'RFLU_HM_SolveEnergyEq',__FILE__)

    delT = global%dtMin*global%refVelocity/global%refLength
    iDt = 1.0_RFREAL/delT

    refGamma = global%refGamma
    rRef     = global%refDensity
    pRef     = global%refPressure
    uRef     = global%refVelocity 
    aRef     = (refGamma*pRef/rRef)**0.5_RFREAL
    refM     = uRef/aRef
    refM2    = (uRef/aRef)**2.0_RFREAL
    iRefMu   = 1.0_RFREAL/global%refVisc
    iRefReNum = 1.0_RFREAL/global%refReNum

! *****************************************************************************
!   Initialize Hypre matrix and vectors
! *****************************************************************************

    CALL RFLU_HYPRE_InitializeMatrixVector(regions)

! *****************************************************************************
!   Loop over regions and assign pointers
! *****************************************************************************

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid

      pVfMixt    => pRegion%mixt%vfMixt
      pVfMixtOld => pRegion%mixt%vfMixtOld
      pCv        => pRegion%mixt%cv
      pCvOld     => pRegion%mixt%cvOld
      pCvOld2    => pRegion%mixt%cvOld2
      pGradP     => pRegion%mixt%gradCell(:,GRC_MIXT_PRES,:)
      pGradPOld  => pRegion%mixt%gradCellOld(:,GRC_MIXT_PRES,:)
      pGradPOld2 => pRegion%mixt%gradCellOld2(:,GRC_MIXT_PRES,:)
      pGradCell     => pRegion%mixt%gradCell
      pGradCellOld  => pRegion%mixt%gradCellOld
      pGradCellOld2 => pRegion%mixt%gradCellOld2

      IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN
        pDv => pRegion%mixt%dv
        pGv => pRegion%mixt%gv
        pTv => pRegion%mixt%tv
      END IF ! pRegion%mixtInput%flowModel

      iReg1 = pRegion%iRegionGlobal
      nCellsOffsetRegion = pRegion%nCellsOffset(iReg1)

! ==============================================================================
!     Loop over cells and initialize the coefficients
! ==============================================================================

      DO icg = 1,pGrid%nCells
        volume = pGrid%vol(icg)

        rho_n  = (pCvOld(CV_MIXT_DENS,icg)+pCvOld2(CV_OLD2_MIXT_DENS,icg)) &
                 /2.0_RFREAL
        rho_n1 = (pCv(CV_MIXT_DENS,icg)+pCvOld(CV_MIXT_DENS,icg))/2.0_RFREAL
        uu_n   = pCvOld(CV_MIXT_XVEL,icg)*pCvOld(CV_MIXT_XVEL,icg) &
               + pCvOld(CV_MIXT_YVEL,icg)*pCvOld(CV_MIXT_YVEL,icg) &
               + pCvOld(CV_MIXT_ZVEL,icg)*pCvOld(CV_MIXT_ZVEL,icg)
        uu_n1  = pCv(CV_MIXT_XVEL,icg)*pCv(CV_MIXT_XVEL,icg) &
               + pCv(CV_MIXT_YVEL,icg)*pCv(CV_MIXT_YVEL,icg) &
               + pCv(CV_MIXT_ZVEL,icg)*pCv(CV_MIXT_ZVEL,icg)

        mVoliDt = refM2*volume*iDt

        value = mVoliDt/(2.0_RFREAL)
        iRowGlobal = icg + nCellsOffsetRegion
        iColGlobal = iRowGlobal       
#ifdef HYPRE           
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
#endif
        value = -mVoliDt &
              *( 0.5_RFREAL*(refGamma-1.0_RFREAL)*(rho_n1*uu_n1-rho_n*uu_n) &
               + 0.5_RFREAL*(pCv(CV_MIXT_PRES,icg) &
                            -pCvOld2(CV_OLD2_MIXT_PRES,icg)))
        iRowGlobal = icg + nCellsOffsetRegion
#ifdef HYPRE           
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif                                       
      END DO ! icg

! ==============================================================================
!     Loop over interior faces and accumulate coefficients
! ==============================================================================

      DO ifg = 1,pGrid%nFaces-pGrid%nFacesAV
        c1 = pGrid%f2c(1,ifg)
        c2 = pGrid%f2c(2,ifg)

        nx = pGrid%fn(XCOORD,ifg)
        ny = pGrid%fn(YCOORD,ifg)
        nz = pGrid%fn(ZCOORD,ifg)
        nm = pGrid%fn(XYZMAG,ifg)

        deln = ( (pGrid%cofg(XCOORD,c2)-pGrid%cofg(XCOORD,c1)) &
                *(pGrid%cofg(XCOORD,c2)-pGrid%cofg(XCOORD,c1)) &
               + (pGrid%cofg(YCOORD,c2)-pGrid%cofg(YCOORD,c1)) &
                *(pGrid%cofg(YCOORD,c2)-pGrid%cofg(YCOORD,c1)) &
               + (pGrid%cofg(ZCOORD,c2)-pGrid%cofg(ZCOORD,c1)) &
                *(pGrid%cofg(ZCOORD,c2)-pGrid%cofg(ZCOORD,c1)) )**(0.5_RFREAL)

        rho_k = (pCvOld(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c2))/2.0_RFREAL

        rho_k_n1 = ( pCv(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c1) &
                   + pCv(CV_MIXT_DENS,c2)+pCvOld(CV_MIXT_DENS,c2) )/4.0_RFREAL

        v_nHalf = (pVfMixt(ifg)+pVfMixtOld(ifg))/2.0_RFREAL

        p_tilda = ( pCvOld2(CV_OLD2_MIXT_PRES,c1) &
                  + 2.0_RFREAL*pCvOld(CV_MIXT_PRES,c1) &
                  + pCv(CV_MIXT_PRES,c1) &
                  + pCvOld2(CV_OLD2_MIXT_PRES,c2) &
                  + 2.0_RFREAL*pCvOld(CV_MIXT_PRES,c2) &
                  + pCv(CV_MIXT_PRES,c2) )/8.0_RFREAL

        uu_k_nHalf = (( pCv(CV_MIXT_XVEL,c1)+pCvOld(CV_MIXT_XVEL,c1) &
                      + pCv(CV_MIXT_XVEL,c2)+pCvOld(CV_MIXT_XVEL,c2) ) &
                     *( pCv(CV_MIXT_XVEL,c1)+pCvOld(CV_MIXT_XVEL,c1) &
                      + pCv(CV_MIXT_XVEL,c2)+pCvOld(CV_MIXT_XVEL,c2) ) &
                    + ( pCv(CV_MIXT_YVEL,c1)+pCvOld(CV_MIXT_YVEL,c1) &
                      + pCv(CV_MIXT_YVEL,c2)+pCvOld(CV_MIXT_YVEL,c2) ) &
                     *( pCv(CV_MIXT_YVEL,c1)+pCvOld(CV_MIXT_YVEL,c1) &
                      + pCv(CV_MIXT_YVEL,c2)+pCvOld(CV_MIXT_YVEL,c2) ) &
                    + ( pCv(CV_MIXT_ZVEL,c1)+pCvOld(CV_MIXT_ZVEL,c1) &
                      + pCv(CV_MIXT_ZVEL,c2)+pCvOld(CV_MIXT_ZVEL,c2) ) &       
                     *( pCv(CV_MIXT_ZVEL,c1)+pCvOld(CV_MIXT_ZVEL,c1) &
                      + pCv(CV_MIXT_ZVEL,c2)+pCvOld(CV_MIXT_ZVEL,c2) )) &
                    /16.0_RFREAL       

        pi_k = refGamma*p_tilda &
             + 0.5_RFREAL*(refGamma-1.0_RFREAL)*rho_k*uu_k_nHalf

        uic1_ni = pCv(CV_MIXT_XVEL,c1)*nx &
                + pCv(CV_MIXT_YVEL,c1)*ny &
                + pCv(CV_MIXT_ZVEL,c1)*nz

        uic2_ni = - ( pCv(CV_MIXT_XVEL,c2)*nx &
                    + pCv(CV_MIXT_YVEL,c2)*ny &
                    + pCv(CV_MIXT_ZVEL,c2)*nz )

        term1 = v_nHalf*refGamma*nm/8.0_RFREAL
        term2 = delT*nm/(8.0_RFREAL*rho_k_n1*deln)
        term3 = (refGamma-1.0_RFREAL)*delT*rho_k*v_nHalf*v_nHalf*nm &
                /(8.0_RFREAL*deln*rho_k_n1)

        value = refM2*(term1 + pi_k*term2 + term3)
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = iRowGlobal
#ifdef HYPRE           
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
                                       
        value = refM2*( - term1 + pi_k*term2 + term3)
        iRowGlobal = c2 + nCellsOffsetRegion
        iColGlobal = iRowGlobal          
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value = refM2*(term1 - pi_k*term2 - term3)
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = c2 + pRegion%nCellsOffset(iReg1)         
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value = refM2*( - term1 - pi_k*term2 - term3)
        iRowGlobal = c2 + nCellsOffsetRegion
        iColGlobal = c1 + nCellsOffsetRegion       
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value = term2 
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = iRowGlobal         
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
                                       
        value = term2
        iRowGlobal = c2 + nCellsOffsetRegion
        iColGlobal = iRowGlobal    
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value = - (refGamma-1.0_RFREAL)*refM2*(uic1_ni)*nm/8.0_RFREAL
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = c2 + nCellsOffsetRegion       
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
                                       
        value = - (refGamma-1.0_RFREAL)*refM2*(uic2_ni)*nm/8.0_RFREAL
        iRowGlobal = c2 + nCellsOffsetRegion
        iColGlobal = c1 + nCellsOffsetRegion         
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value = - term2
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = c2 + nCellsOffsetRegion   
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value = - term2
        iRowGlobal = c2 + nCellsOffsetRegion
        iColGlobal = c1 + nCellsOffsetRegion
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value = - refM2*pi_k*v_nHalf*nm
        iRowGlobal = c1 + nCellsOffsetRegion 
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)

        value =  refM2*pi_k*v_nHalf*nm
        iRowGlobal = c2 + nCellsOffsetRegion
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)

        value = - v_nHalf*nm
        iRowGlobal = c1 + nCellsOffsetRegion
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)

        value =  v_nHalf*nm
        iRowGlobal = c2 + nCellsOffsetRegion
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif
      END DO ! ifg

      DO ifg = pGrid%nFaces-pGrid%nFacesAV+1,pGrid%nFaces
        c1 = pGrid%f2c(1,ifg)
        c2 = pGrid%f2c(2,ifg)

        direction = 1.0_RFREAL

        IF ( c1 .GT. pGrid%nCells ) THEN ! c1 is virtual
          dummy = c1
          c1    = c2
          c2    = dummy
 
          direction = -1.0_RFREAL
        END IF

        nx = direction*pGrid%fn(XCOORD,ifg)
        ny = direction*pGrid%fn(YCOORD,ifg)
        nz = direction*pGrid%fn(ZCOORD,ifg)
        nm = pGrid%fn(XYZMAG,ifg)

        deln = ( (pGrid%cofg(XCOORD,c2)-pGrid%cofg(XCOORD,c1)) &
                *(pGrid%cofg(XCOORD,c2)-pGrid%cofg(XCOORD,c1)) &
               + (pGrid%cofg(YCOORD,c2)-pGrid%cofg(YCOORD,c1)) &
                *(pGrid%cofg(YCOORD,c2)-pGrid%cofg(YCOORD,c1)) &
               + (pGrid%cofg(ZCOORD,c2)-pGrid%cofg(ZCOORD,c1)) &
                *(pGrid%cofg(ZCOORD,c2)-pGrid%cofg(ZCOORD,c1)) )**(0.5_RFREAL)

        rho_k = (pCvOld(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c2))/2.0_RFREAL

        rho_k_n1 = ( pCv(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c1) &
                   + pCv(CV_MIXT_DENS,c2)+pCvOld(CV_MIXT_DENS,c2) )/4.0_RFREAL

        v_nHalf = direction*(pVfMixt(ifg)+pVfMixtOld(ifg))/2.0_RFREAL

        p_tilda = ( pCvOld2(CV_OLD2_MIXT_PRES,c1) &
                  + 2.0_RFREAL*pCvOld(CV_MIXT_PRES,c1) &
                  + pCv(CV_MIXT_PRES,c1) &
                  + pCvOld2(CV_OLD2_MIXT_PRES,c2) &
                  + 2.0_RFREAL*pCvOld(CV_MIXT_PRES,c2) &
                  + pCv(CV_MIXT_PRES,c2) )/8.0_RFREAL

        uu_k_nHalf = (( pCv(CV_MIXT_XVEL,c1)+pCvOld(CV_MIXT_XVEL,c1) &
                      + pCv(CV_MIXT_XVEL,c2)+pCvOld(CV_MIXT_XVEL,c2) ) &
                     *( pCv(CV_MIXT_XVEL,c1)+pCvOld(CV_MIXT_XVEL,c1) &
                      + pCv(CV_MIXT_XVEL,c2)+pCvOld(CV_MIXT_XVEL,c2) ) &
                    + ( pCv(CV_MIXT_YVEL,c1)+pCvOld(CV_MIXT_YVEL,c1) &
                      + pCv(CV_MIXT_YVEL,c2)+pCvOld(CV_MIXT_YVEL,c2) ) &
                     *( pCv(CV_MIXT_YVEL,c1)+pCvOld(CV_MIXT_YVEL,c1) &
                      + pCv(CV_MIXT_YVEL,c2)+pCvOld(CV_MIXT_YVEL,c2) ) &
                    + ( pCv(CV_MIXT_ZVEL,c1)+pCvOld(CV_MIXT_ZVEL,c1) &
                      + pCv(CV_MIXT_ZVEL,c2)+pCvOld(CV_MIXT_ZVEL,c2) ) &       
                     *( pCv(CV_MIXT_ZVEL,c1)+pCvOld(CV_MIXT_ZVEL,c1) &
                      + pCv(CV_MIXT_ZVEL,c2)+pCvOld(CV_MIXT_ZVEL,c2) )) &
                   /16.0_RFREAL       

        pi_k = refGamma*p_tilda &
             + 0.5_RFREAL*(refGamma-1.0_RFREAL)*rho_k*uu_k_nHalf

        uic1_ni = (pCv(CV_MIXT_XVEL,c1)*nx &
                 + pCv(CV_MIXT_YVEL,c1)*ny &
                 + pCv(CV_MIXT_ZVEL,c1)*nz)

        term1 = v_nHalf*refGamma*nm/8.0_RFREAL
        term2 = delT*nm/(8.0_RFREAL*rho_k_n1*deln)
        term3 = (refGamma-1.0_RFREAL)*delT*rho_k*v_nHalf*v_nHalf*nm &
                /(8.0_RFREAL*deln*rho_k_n1)

        value = refM2*(term1 + pi_k*term2 + term3)
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = iRowGlobal
#ifdef HYPRE   
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value = refM2*(term1 - pi_k*term2 - term3)
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = pGrid%virt2GlobalIds(c2 - pGrid%nCells)
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value = term2 
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = iRowGlobal
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value = - (refGamma-1.0_RFREAL)*refM2*(uic1_ni)*nm/8.0_RFREAL
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = pGrid%virt2GlobalIds(c2 - pGrid%nCells)
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value = - term2
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = pGrid%virt2GlobalIds(c2 - pGrid%nCells)
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value = - refM2*pi_k*v_nHalf*nm
        iRowGlobal = c1 + nCellsOffsetRegion
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)

        value = - v_nHalf*nm
        iRowGlobal = c1 + nCellsOffsetRegion
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif                                       
      END DO ! ifg

! ------------------------------------------------------------------------------
!     Loop over boundary faces. Energy eqn.
! ------------------------------------------------------------------------------
                          
      DO iPatch = 1,pRegion%grid%nPatches
        pPatch => pRegion%patches(iPatch)
                          
        distrib = pPatch%mixt%distrib
 
        SELECT CASE( pPatch%bcType )
        CASE (BC_SLIPWALL,BC_NOSLIPWALL_HFLUX,BC_NOSLIPWALL_TEMP)
          DO ifg = 1,pPatch%nBFaces
            c1 = pPatch%bf2c(ifg)
                         
            nx = pPatch%fn(XCOORD,ifg)
            ny = pPatch%fn(YCOORD,ifg)
            nz = pPatch%fn(ZCOORD,ifg)
            nm = pPatch%fn(XYZMAG,ifg)
                     
            uic1_ni =  pCv(CV_MIXT_XVEL,c1)*nx + pCv(CV_MIXT_YVEL,c1)*ny &
                     + pCv(CV_MIXT_ZVEL,c1)*nz

            value = - (refGamma-1.0_RFREAL)*refM2*(uic1_ni)*nm/8.0_RFREAL
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = iRowGlobal
#ifdef HYPRE   
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)
#endif                                       

            pPatch%cp(ifg) = 2.0_RFREAL*pCv(CV_MIXT_PRES,c1)
          END DO ! ifg

        CASE (BC_OUTFLOW)
          DO ifg = 1,pPatch%nBFaces
            c1 = pPatch%bf2c(ifg)
                         
            xc = pPatch%fc(XCOORD,ifg)
            yc = pPatch%fc(YCOORD,ifg)
            zc = pPatch%fc(ZCOORD,ifg)

            nx = pPatch%fn(XCOORD,ifg)
            ny = pPatch%fn(YCOORD,ifg)
            nz = pPatch%fn(ZCOORD,ifg)
            nm = pPatch%fn(XYZMAG,ifg)
                     
            deln = 2.0_RFREAL &
                   *( (xc-pGrid%cofg(XCOORD,c1))*(xc-pGrid%cofg(XCOORD,c1)) &
                    + (yc-pGrid%cofg(YCOORD,c1))*(yc-pGrid%cofg(YCOORD,c1)) &
                    + (zc-pGrid%cofg(ZCOORD,c1))*(zc-pGrid%cofg(ZCOORD,c1)) &
                    )**(0.5_RFREAL)

            rho_k = (pCvOld(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c1))/2.0_RFREAL

            rho_k_n1 = ( pCv(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c1) &
                   + pCv(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c1) )/4.0_RFREAL

            v_nHalf = (pVfMixt(ifg)+pVfMixtOld(ifg))/2.0_RFREAL

            p_tilda = pPatch%mixt%vals(BCDAT_OUTFLOW_PRESS,distrib*ifg)
            p_tilda = (p_tilda-global%refPressure)/(global%refDensity &
                      *global%refVelocity*global%refVelocity)

            uu_k_nHalf = (( pCv(CV_MIXT_XVEL,c1)+pCvOld(CV_MIXT_XVEL,c1) &
                          + pCv(CV_MIXT_XVEL,c1)+pCvOld(CV_MIXT_XVEL,c1) ) &
                         *( pCv(CV_MIXT_XVEL,c1)+pCvOld(CV_MIXT_XVEL,c1) &
                          + pCv(CV_MIXT_XVEL,c1)+pCvOld(CV_MIXT_XVEL,c1) ) &
                        + ( pCv(CV_MIXT_YVEL,c1)+pCvOld(CV_MIXT_YVEL,c1) &
                          + pCv(CV_MIXT_YVEL,c1)+pCvOld(CV_MIXT_YVEL,c1) ) &
                         *( pCv(CV_MIXT_YVEL,c1)+pCvOld(CV_MIXT_YVEL,c1) &
                          + pCv(CV_MIXT_YVEL,c1)+pCvOld(CV_MIXT_YVEL,c1) ) &
                        + ( pCv(CV_MIXT_ZVEL,c1)+pCvOld(CV_MIXT_ZVEL,c1) &
                          + pCv(CV_MIXT_ZVEL,c1)+pCvOld(CV_MIXT_ZVEL,c1) ) &       
                         *( pCv(CV_MIXT_ZVEL,c1)+pCvOld(CV_MIXT_ZVEL,c1) &
                          + pCv(CV_MIXT_ZVEL,c1)+pCvOld(CV_MIXT_ZVEL,c1) )) &
                        /16.0_RFREAL       

            pi_k = refGamma*p_tilda &
                 + 0.5_RFREAL*(refGamma-1.0_RFREAL)*rho_k*uu_k_nHalf

            uic1_ni = pCv(CV_MIXT_XVEL,c1)*nx &
                    + pCv(CV_MIXT_YVEL,c1)*ny &
                    + pCv(CV_MIXT_ZVEL,c1)*nz

            uic2_ni = - ( pCv(CV_MIXT_XVEL,c2)*nx &
                        + pCv(CV_MIXT_YVEL,c2)*ny &
                        + pCv(CV_MIXT_ZVEL,c2)*nz )

            term1 = v_nHalf*refGamma*nm/8.0_RFREAL
            term2 = delT*nm/(8.0_RFREAL*rho_k_n1*deln)
            term3 = (refGamma-1.0_RFREAL)*delT*rho_k*v_nHalf*v_nHalf*nm &
                    /(8.0_RFREAL*deln*rho_k_n1)

            value = refM2*(pi_k*term2 + term3)
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = iRowGlobal
#ifdef HYPRE           
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)
#endif
                                       
            value = -refM2*(- pi_k*term2 - term3)
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = c1 + pRegion%nCellsOffset(iReg1)         
#ifdef HYPRE           
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)
#endif

            value = term2 
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = iRowGlobal         
#ifdef HYPRE           
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)
#endif
                                       
            value = (refGamma-1.0_RFREAL)*refM2*(uic1_ni)*nm/8.0_RFREAL
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = c1 + nCellsOffsetRegion       
#ifdef HYPRE           
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)
#endif
                                       
            value = term2
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = c1 + nCellsOffsetRegion   
#ifdef HYPRE           
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)
#endif

            value = - refM2*pi_k*v_nHalf*nm
            iRowGlobal = c1 + nCellsOffsetRegion 
#ifdef HYPRE           
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                           iErrHypre)
#endif

            value = - v_nHalf*nm
            iRowGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE           
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                           iErrHypre)
#endif

            pPatch%cp(ifg) = 2.0_RFREAL*pCv(CV_MIXT_PRES,c1)
          END DO ! ifg

        CASE (BC_FARFIELD)
          DO ifg = 1,pPatch%nBFaces
            c1 = pPatch%bf2c(ifg)
                         
            xc = pPatch%fc(XCOORD,ifg)
            yc = pPatch%fc(YCOORD,ifg)
            zc = pPatch%fc(ZCOORD,ifg)

            nx = pPatch%fn(XCOORD,ifg)
            ny = pPatch%fn(YCOORD,ifg)
            nz = pPatch%fn(ZCOORD,ifg)
            nm = pPatch%fn(XYZMAG,ifg)
                     
            deln = 2.0_RFREAL &
                   *( (xc-pGrid%cofg(XCOORD,c1))*(xc-pGrid%cofg(XCOORD,c1)) &
                    + (yc-pGrid%cofg(YCOORD,c1))*(yc-pGrid%cofg(YCOORD,c1)) &
                    + (zc-pGrid%cofg(ZCOORD,c1))*(zc-pGrid%cofg(ZCOORD,c1)) &
                    )**(0.5_RFREAL)

            rho_k = pCvOld(CV_MIXT_DENS,c1)
            
            rho_k_n1 = (pCv(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c1))/2.0_RFREAL

            mf  = pPatch%mixt%vals(BCDAT_FARF_MACH,distrib*ifg)
            aoa = pPatch%mixt%vals(BCDAT_FARF_ATTACK,distrib*ifg)
            aos = pPatch%mixt%vals(BCDAT_FARF_SLIP,distrib*ifg)
            tf  = pPatch%mixt%vals(BCDAT_FARF_TEMP,distrib*ifg)

            SpeedOfSound = SQRT(tf/global%refTemperature)

            VelX = mf*SpeedOfSound*COS(aoa)*COS(aos)
            VelY = mf*SpeedOfSound*SIN(aoa)*COS(aos)
            VelZ = mf*SpeedOfSound*SIN(aos)

            v_nHalf = (nx*VelX+ny*VelY+nz*VelZ)

            p_tilda = ( pCvOld2(CV_OLD2_MIXT_PRES,c1) &
                      + 2.0_RFREAL*pCvOld(CV_MIXT_PRES,c1) &
                      + pCv(CV_MIXT_PRES,c1))/4.0_RFREAL

            uu_k_nHalf = (( pCv(CV_MIXT_XVEL,c1)+pCvOld(CV_MIXT_XVEL,c1) &
                          + VelX + VelX ) &
                         *( pCv(CV_MIXT_XVEL,c1)+pCvOld(CV_MIXT_XVEL,c1) &
                          + VelX + VelX ) &
                        + ( pCv(CV_MIXT_YVEL,c1)+pCvOld(CV_MIXT_YVEL,c1) &
                          + VelY + VelY ) &
                         *( pCv(CV_MIXT_YVEL,c1)+pCvOld(CV_MIXT_YVEL,c1) &
                          + VelY + VelY ) &
                        + ( pCv(CV_MIXT_ZVEL,c1)+pCvOld(CV_MIXT_ZVEL,c1) &
                          + VelZ + VelZ ) &
                         *( pCv(CV_MIXT_ZVEL,c1)+pCvOld(CV_MIXT_ZVEL,c1) &
                          + VelZ + VelZ ) &
                        )/16.0_RFREAL

            pi_k = refGamma*p_tilda &
                 + 0.5_RFREAL*(refGamma-1.0_RFREAL)*rho_k*uu_k_nHalf

            uic1_ni =  pCv(CV_MIXT_XVEL,c1)*nx + pCv(CV_MIXT_YVEL,c1)*ny &
                     + pCv(CV_MIXT_ZVEL,c1)*nz

            term1 = v_nHalf*refGamma*nm/8.0_RFREAL
            term2 = delT*nm/(8.0_RFREAL*rho_k_n1*deln)
            term3 = (refGamma-1.0_RFREAL)*delT*rho_k*v_nHalf*v_nHalf*nm &
                    /(8.0_RFREAL*deln*rho_k_n1)

#ifdef HYPRE   
            value = refM2*(term1 + pi_k*term2 + term3)
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = iRowGlobal
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)

! ------------------------------------------------------------------------------
! TEMPORARY:
!  Dont know how to deal with a_beta terms as can not be put on rhs as delp_beta
!  is not known. As of now commenting out all a_beta terms that means
!  treating delp_boundary=0.
!  For Slipwall, delp_wall = delp_interior, that means a_beta terms go to
!  a_alpha
! END TEMPORARY
! ------------------------------------------------------------------------------

            value = term2
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = iRowGlobal
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)

            value = - refM2*pi_k*v_nHalf*nm
            iRowGlobal = c1 + nCellsOffsetRegion
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal, &
                                           value,iErrHypre)

            value = - v_nHalf*nm
            iRowGlobal = c1 + nCellsOffsetRegion
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal, &
                                           value,iErrHypre)
#endif

            pPatch%cp(ifg) = 2.0_RFREAL*pCv(CV_MIXT_PRES,c1)
          END DO ! ifg

        CASE (BC_INFLOW_VELTEMP)
          DO ifg = 1,pPatch%nBFaces
            c1 = pPatch%bf2c(ifg)

            nx = pPatch%fn(XCOORD,ifg)
            ny = pPatch%fn(YCOORD,ifg)
            nz = pPatch%fn(ZCOORD,ifg)
            nm = pPatch%fn(XYZMAG,ifg)

            minj = -pPatch%mixt%vals(BCDAT_INJECT_MFRATE,distrib*ifg)
            minj = minj/(global%refDensity*global%refVelocity)
            tinj = pPatch%mixt%vals(BCDAT_INJECT_TEMP,distrib*ifg)
            tinj = tinj/global%refTemperature

            v_nHalf = minj*tinj/(refGamma*refM2*pCvOld(CV_MIXT_PRES,c1)+1.0_RFREAL)

            p_tilda = ( pCvOld2(CV_OLD2_MIXT_PRES,c1) &
                      + 2.0_RFREAL*pCvOld(CV_MIXT_PRES,c1) &
                      + pCv(CV_MIXT_PRES,c1) )/4.0_RFREAL

            pi_k = refGamma*p_tilda &
                 + 0.5_RFREAL*(refGamma-1.0_RFREAL)*minj*v_nHalf

            uic1_ni = pCv(CV_MIXT_XVEL,c1)*nx + pCv(CV_MIXT_YVEL,c1)*ny &
                    + pCv(CV_MIXT_ZVEL,c1)*nz

            value = -(refGamma-1.0_RFREAL)*refM2*uic1_ni*nm/8.0_RFREAL 
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = iRowGlobal
#ifdef HYPRE   
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)
#endif                                       

            value = refGamma*refM2*v_nHalf*nm/4.0_RFREAL 
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = iRowGlobal
#ifdef HYPRE   
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)
#endif                                       

            value = -(refM2*pi_k + 1.0_RFREAL)*v_nHalf*nm
            iRowGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE   
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                           iErrHypre)
#endif                                       

            pPatch%cp(ifg) = 2.0_RFREAL*pCv(CV_MIXT_PRES,c1)
          END DO ! ifg

        CASE (BC_INJECTION)
          DO ifg = 1,pPatch%nBFaces
            c1 = pPatch%bf2c(ifg)

            nx = pPatch%fn(XCOORD,ifg)
            ny = pPatch%fn(YCOORD,ifg)
            nz = pPatch%fn(ZCOORD,ifg)
            nm = pPatch%fn(XYZMAG,ifg)

            minj = -pPatch%mixt%vals(BCDAT_INJECT_MFRATE,distrib*ifg)
            minj = minj/(global%refDensity*global%refVelocity)
            tinj = pPatch%mixt%vals(BCDAT_INJECT_TEMP,distrib*ifg)
            tinj = tinj/global%refTemperature

            v_nHalf = minj*tinj/(refGamma*refM2*pCvOld(CV_MIXT_PRES,c1)+1.0_RFREAL)

            p_tilda = ( pCvOld2(CV_OLD2_MIXT_PRES,c1) &
                      + 2.0_RFREAL*pCvOld(CV_MIXT_PRES,c1) &
                      + pCv(CV_MIXT_PRES,c1) )/4.0_RFREAL

            pi_k = refGamma*p_tilda &
                 + 0.5_RFREAL*(refGamma-1.0_RFREAL)*minj*v_nHalf

            uic1_ni = pCv(CV_MIXT_XVEL,c1)*nx + pCv(CV_MIXT_YVEL,c1)*ny &
                    + pCv(CV_MIXT_ZVEL,c1)*nz

            value = -(refGamma-1.0_RFREAL)*refM2*uic1_ni*nm/8.0_RFREAL 
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = iRowGlobal
#ifdef HYPRE   
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)
#endif                                       

            value = refGamma*refM2*v_nHalf*nm/4.0_RFREAL 
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = iRowGlobal
#ifdef HYPRE   
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)
#endif                                       

            value = -(refM2*pi_k + 1.0_RFREAL)*v_nHalf*nm
            iRowGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE   
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                           iErrHypre)
#endif                                       

            pPatch%cp(ifg) = 2.0_RFREAL*pCv(CV_MIXT_PRES,c1)
          END DO ! ifg

        CASE (BC_PERIODIC,BC_SYMMETRY,BC_VIRTUAL)
          DO ifg = 1,pPatch%nBFaces
            c1 = pPatch%bf2c(ifg)

            pPatch%cp(ifg) = 2.0_RFREAL*pCv(CV_MIXT_PRES,c1)
          END DO ! ifg
        CASE DEFAULT                
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pPatch%bcType
      END DO ! iPatch

! ==============================================================================
!     Compute coefficients due to moving reference frame source terms. EnergyEq.
! ==============================================================================

      global => regions(1)%global

      IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN

      DO icg = 1,pGrid%nCells
        volume = pGrid%vol(icg)

        aiui_nHalf = (pRegion%mvfAcc(XCOORD)*(pCv(CV_MIXT_XVEL,icg) &
                                         +pCvOld(CV_MIXT_XVEL,icg)) &
                     +pRegion%mvfAcc(YCOORD)*(pCv(CV_MIXT_YVEL,icg) &
                                         +pCvOld(CV_MIXT_YVEL,icg)) &
                     +pRegion%mvfAcc(ZCOORD)*(pCv(CV_MIXT_ZVEL,icg) &
                                         +pCvOld(CV_MIXT_ZVEL,icg)))*0.5_RFREAL

#ifdef HYPRE           
        value = -(refGamma-1.0_RFREAL)*refM2*volume*pCvOld(CV_MIXT_DENS,icg) &
                *aiui_nHalf
        iRowGlobal = icg + nCellsOffsetRegion
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif                                       
      END DO ! icg

      DO ifg = 1,pGrid%nFaces-pGrid%nFacesAV
        c1 = pGrid%f2c(1,ifg)
        c2 = pGrid%f2c(2,ifg)

        nx = pGrid%fn(XCOORD,ifg)
        ny = pGrid%fn(YCOORD,ifg)
        nz = pGrid%fn(ZCOORD,ifg)
        nm = pGrid%fn(XYZMAG,ifg)

        rho_c1_n1 = (pCv(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c1))/2.0_RFREAL
        rho_c2_n1 = (pCv(CV_MIXT_DENS,c2)+pCvOld(CV_MIXT_DENS,c2))/2.0_RFREAL

        aini = pRegion%mvfAcc(XCOORD)*nx+pRegion%mvfAcc(YCOORD)*ny &
              +pRegion%mvfAcc(ZCOORD)*nz

#ifdef HYPRE           
        value = -(refGamma-1.0_RFREAL)*refM2*delT*aini*nm &
                 *pCvOld(CV_MIXT_DENS,c1)/(16_RFREAL*rho_c1_n1)
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = c2 + nCellsOffsetRegion       
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
#endif                                       

#ifdef HYPRE           
        value = (refGamma-1.0_RFREAL)*refM2*delT*aini*nm &
                 *pCvOld(CV_MIXT_DENS,c2)/(16_RFREAL*rho_c2_n1)
        iRowGlobal = c2 + nCellsOffsetRegion
        iColGlobal = c1 + nCellsOffsetRegion       
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
#endif                                       
      END DO ! ifg

      DO ifg = pGrid%nFaces-pGrid%nFacesAV+1,pGrid%nFaces
        c1 = pGrid%f2c(1,ifg)
        c2 = pGrid%f2c(2,ifg)

        direction = 1.0_RFREAL

        IF ( c1 .GT. pGrid%nCells ) THEN ! c1 is virtual
          dummy = c1
          c1    = c2
          c2    = dummy
 
          direction = -1.0_RFREAL
        END IF

        nx = direction*pGrid%fn(XCOORD,ifg)
        ny = direction*pGrid%fn(YCOORD,ifg)
        nz = direction*pGrid%fn(ZCOORD,ifg)
        nm = pGrid%fn(XYZMAG,ifg)

        rho_c1_n1 = (pCv(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c1))/2.0_RFREAL

        aini = pRegion%mvfAcc(XCOORD)*nx+pRegion%mvfAcc(YCOORD)*ny &
              +pRegion%mvfAcc(ZCOORD)*nz

#ifdef HYPRE           
        value = -(refGamma-1.0_RFREAL)*refM2*delT*aini*nm &
                 *pCvOld(CV_MIXT_DENS,c1)/(16_RFREAL*rho_c1_n1)
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = pGrid%virt2GlobalIds(c2 - pGrid%nCells)
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
#endif                                       
      END DO ! ifg

      DO iPatch = 1,pRegion%grid%nPatches
        pPatch => pRegion%patches(iPatch)

        distrib = pPatch%mixt%distrib

      SELECT CASE( pPatch%bcType )
      CASE (BC_SLIPWALL,BC_NOSLIPWALL_HFLUX,BC_NOSLIPWALL_TEMP)
        DO ifg = 1,pPatch%nBFaces
          c1 = pPatch%bf2c(ifg)

          nx = pPatch%fn(XCOORD,ifg)
          ny = pPatch%fn(YCOORD,ifg)
          nz = pPatch%fn(ZCOORD,ifg)
          nm = pPatch%fn(XYZMAG,ifg)

          rho_c1_n1 = (pCv(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c1))/2.0_RFREAL
          
          aini = pRegion%mvfAcc(XCOORD)*nx+pRegion%mvfAcc(YCOORD)*ny &
                +pRegion%mvfAcc(ZCOORD)*nz

#ifdef HYPRE           
          value = -(refGamma-1.0_RFREAL)*refM2*delT*aini*nm &
                   *pCvOld(CV_MIXT_DENS,c1)/(16_RFREAL*rho_c1_n1)
          iRowGlobal = c1 + nCellsOffsetRegion
          iColGlobal = iRowGlobal 
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)
#endif  
        END DO ! ifg

      CASE (BC_OUTFLOW)
      CASE (BC_FARFIELD)
      CASE (BC_INFLOW_VELTEMP)
      CASE (BC_INJECTION)
      CASE (BC_PERIODIC,BC_SYMMETRY,BC_VIRTUAL)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%bcType
      END DO ! iPatch
      END IF ! global%mvFrameFlag

! ==============================================================================
!     Compute coefficients due to absorbing boundary condition. Energy eqn.
! ==============================================================================

      IF ( (global%abcFlag .EQV. .TRUE.) .AND. (global%abcKind == 0) ) THEN
        DO icg = 1,pGrid%nCells
          distrib = global%abcDistrib
          ro = pRegion%mixt%cvRef(CV_MIXT_DENS,distrib*icg)
          uo = pRegion%mixt%cvRef(CV_MIXT_XVEL,distrib*icg)
          vo = pRegion%mixt%cvRef(CV_MIXT_YVEL,distrib*icg)
          wo = pRegion%mixt%cvRef(CV_MIXT_ZVEL,distrib*icg)
          po = pRegion%mixt%cvRef(CV_MIXT_PRES,distrib*icg)

          volume = pGrid%vol(icg)

          sigma = pRegion%mixt%sigma(icg)*global%refLength/global%refVelocity 

          uu_np1 = pCv(CV_MIXT_XVEL,icg)*pCv(CV_MIXT_XVEL,icg) &
                  +pCv(CV_MIXT_YVEL,icg)*pCv(CV_MIXT_YVEL,icg) &
                  +pCv(CV_MIXT_ZVEL,icg)*pCv(CV_MIXT_ZVEL,icg)

          uu_n = pCvOld(CV_MIXT_XVEL,icg)*pCvOld(CV_MIXT_XVEL,icg) &
                +pCvOld(CV_MIXT_YVEL,icg)*pCvOld(CV_MIXT_YVEL,icg) &
                +pCvOld(CV_MIXT_ZVEL,icg)*pCvOld(CV_MIXT_ZVEL,icg)

          p_tilda = ( pCvOld2(CV_OLD2_MIXT_PRES,icg) &
                    + 2.0_RFREAL*pCvOld(CV_MIXT_PRES,icg) &
                    + pCv(CV_MIXT_PRES,icg) )/4.0_RFREAL

          value = refM2*sigma*volume*0.25_RFREAL
          iRowGlobal = icg + nCellsOffsetRegion 
          iColGlobal = iRowGlobal
#ifdef HYPRE        
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)
#endif

          value = -(pCvOld(CV_MIXT_DENS,icg)*(uu_np1+uu_n) &
                    *(refGamma-1.0_RFREAL)/4.0_RFREAL &
                    + p_tilda &
                    - ro*uo*uo*(refGamma-1.0_RFREAL)/2.0_RFREAL - po) &
                  *refM2*sigma*volume
          iRowGlobal = icg + nCellsOffsetRegion
#ifdef HYPRE           
          CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                         iErrHypre)
#endif                                       
        END DO ! icg

        DO ifg = 1,pGrid%nFaces-pGrid%nFacesAV
          c1 = pGrid%f2c(1,ifg)
          c2 = pGrid%f2c(2,ifg)

          nx = pGrid%fn(XCOORD,ifg)
          ny = pGrid%fn(YCOORD,ifg)
          nz = pGrid%fn(ZCOORD,ifg)
          nm = pGrid%fn(XYZMAG,ifg)

          sigma_c1 = pRegion%mixt%sigma(c1)*global%refLength/global%refVelocity 
          sigma_c2 = pRegion%mixt%sigma(c2)*global%refLength/global%refVelocity

          rho_c1_n1 = (pCv(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c1))/2.0_RFREAL
          rho_c2_n1 = (pCv(CV_MIXT_DENS,c2)+pCvOld(CV_MIXT_DENS,c2))/2.0_RFREAL

          uini_c1 = pCv(CV_MIXT_XVEL,c1)*nx + pCv(CV_MIXT_YVEL,c1)*ny &
                  + pCv(CV_MIXT_ZVEL,c1)*nz
          uini_c2 = pCv(CV_MIXT_XVEL,c2)*nx + pCv(CV_MIXT_YVEL,c2)*ny &
                  + pCv(CV_MIXT_ZVEL,c2)*nz

#ifdef HYPRE           
          value = -(refGamma-1.0_RFREAL)*refM2*delT*sigma_c1*nm &
                 *pCvOld(CV_MIXT_DENS,c1)*uini_c1/(16_RFREAL*rho_c1_n1)
          iRowGlobal = c1 + nCellsOffsetRegion
          iColGlobal = c2 + nCellsOffsetRegion       
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)
#endif                                       

#ifdef HYPRE           
          value = (refGamma-1.0_RFREAL)*refM2*delT*sigma_c2*nm &
                 *pCvOld(CV_MIXT_DENS,c2)*uini_c2/(16_RFREAL*rho_c2_n1)
          iRowGlobal = c2 + nCellsOffsetRegion
          iColGlobal = c1 + nCellsOffsetRegion       
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)
#endif                                       
        END DO ! ifg

        DO ifg = pGrid%nFaces-pGrid%nFacesAV+1,pGrid%nFaces
          c1 = pGrid%f2c(1,ifg)
          c2 = pGrid%f2c(2,ifg)

          direction = 1.0_RFREAL

          IF ( c1 .GT. pGrid%nCells ) THEN ! c1 is virtual
            dummy = c1
            c1    = c2
            c2    = dummy
 
            direction = -1.0_RFREAL
          END IF

          nx = direction*pGrid%fn(XCOORD,ifg)
          ny = direction*pGrid%fn(YCOORD,ifg)
          nz = direction*pGrid%fn(ZCOORD,ifg)
          nm = pGrid%fn(XYZMAG,ifg)

          sigma_c1 = pRegion%mixt%sigma(c1)*global%refLength/global%refVelocity 

          rho_c1_n1 = (pCv(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c1))/2.0_RFREAL

          uini_c1 = pCv(CV_MIXT_XVEL,c1)*nx + pCv(CV_MIXT_YVEL,c1)*ny &
                  + pCv(CV_MIXT_ZVEL,c1)*nz

#ifdef HYPRE           
          value = -(refGamma-1.0_RFREAL)*refM2*delT*sigma_c1*nm &
                 *pCvOld(CV_MIXT_DENS,c1)*uini_c1/(16_RFREAL*rho_c1_n1)
          iRowGlobal = c1 + nCellsOffsetRegion
          iColGlobal = pGrid%virt2GlobalIds(c2 - pGrid%nCells)
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)
#endif                                       
        END DO ! ifg

! ERROR: manoj: Need to loop over boundary face also.
      END IF ! global%abcFlag

! ==============================================================================
!     Compute viscous flux
!     Loop over interior faces and accumulate coefficients
! ==============================================================================

      global => regions(1)%global

      IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN
        DO ifg = 1,pGrid%nFaces-pGrid%nFacesAV
          c1 = pGrid%f2c(1,ifg)
          c2 = pGrid%f2c(2,ifg)

          nx = pGrid%fn(XCOORD,ifg)
          ny = pGrid%fn(YCOORD,ifg)
          nz = pGrid%fn(ZCOORD,ifg)
          nm = pGrid%fn(XYZMAG,ifg)

          deln = ( (pGrid%cofg(XCOORD,c2)-pGrid%cofg(XCOORD,c1)) &
                  *(pGrid%cofg(XCOORD,c2)-pGrid%cofg(XCOORD,c1)) &
                 + (pGrid%cofg(YCOORD,c2)-pGrid%cofg(YCOORD,c1)) &
                  *(pGrid%cofg(YCOORD,c2)-pGrid%cofg(YCOORD,c1)) &
                 + (pGrid%cofg(ZCOORD,c2)-pGrid%cofg(ZCOORD,c1)) &
                  *(pGrid%cofg(ZCOORD,c2)-pGrid%cofg(ZCOORD,c1)) )**(0.5_RFREAL)

          mu = 0.5_RFREAL*(pTv(TV_MIXT_MUEL,c1)+pTv(TV_MIXT_MUEL,c2)) 
          Pr = global%prLam 

          u1 = 0.25_RFREAL*(pCv(CV_MIXT_XVEL,c1)+pCv(CV_MIXT_XVEL,c2) &
                           +pCvOld(CV_MIXT_XVEL,c1)+pCvOld(CV_MIXT_XVEL,c2))

          u2 = 0.25_RFREAL*(pCv(CV_MIXT_YVEL,c1)+pCv(CV_MIXT_YVEL,c2) &
                           +pCvOld(CV_MIXT_YVEL,c1)+pCvOld(CV_MIXT_YVEL,c2))

          u3 = 0.25_RFREAL*(pCv(CV_MIXT_ZVEL,c1)+pCv(CV_MIXT_ZVEL,c2) &
                           +pCvOld(CV_MIXT_ZVEL,c1)+pCvOld(CV_MIXT_ZVEL,c2))

          du1dx1 = 0.25_RFREAL*(pGradCell(XCOORD,GRC_MIXT_XVEL,c1) + &
                                pGradCell(XCOORD,GRC_MIXT_XVEL,c2) + &
                                pGradCellOld(XCOORD,GRC_MIXT_XVEL,c1) + &
                                pGradCellOld(XCOORD,GRC_MIXT_XVEL,c2))

          du2dx1 = 0.25_RFREAL*(pGradCell(XCOORD,GRC_MIXT_YVEL,c1) + &
                                pGradCell(XCOORD,GRC_MIXT_YVEL,c2) + &
                                pGradCellOld(XCOORD,GRC_MIXT_YVEL,c1) + &
                                pGradCellOld(XCOORD,GRC_MIXT_YVEL,c2))

          du3dx1 = 0.25_RFREAL*(pGradCell(XCOORD,GRC_MIXT_ZVEL,c1) + &
                                pGradCell(XCOORD,GRC_MIXT_ZVEL,c2) + &
                                pGradCellOld(XCOORD,GRC_MIXT_ZVEL,c1) + &
                                pGradCellOld(XCOORD,GRC_MIXT_ZVEL,c2))

          du1dx2 = 0.25_RFREAL*(pGradCell(YCOORD,GRC_MIXT_XVEL,c1) + &
                                pGradCell(YCOORD,GRC_MIXT_XVEL,c2) + &
                                pGradCellOld(YCOORD,GRC_MIXT_XVEL,c1) + &
                                pGradCellOld(YCOORD,GRC_MIXT_XVEL,c2))

          du2dx2 = 0.25_RFREAL*(pGradCell(YCOORD,GRC_MIXT_YVEL,c1) + &
                                pGradCell(YCOORD,GRC_MIXT_YVEL,c2) + &
                                pGradCellOld(YCOORD,GRC_MIXT_YVEL,c1) + &
                                pGradCellOld(YCOORD,GRC_MIXT_YVEL,c2))

          du3dx2 = 0.25_RFREAL*(pGradCell(YCOORD,GRC_MIXT_ZVEL,c1) + &
                                pGradCell(YCOORD,GRC_MIXT_ZVEL,c2) + &
                                pGradCellOld(YCOORD,GRC_MIXT_ZVEL,c1) + &
                                pGradCellOld(YCOORD,GRC_MIXT_ZVEL,c2))

          du1dx3 = 0.25_RFREAL*(pGradCell(ZCOORD,GRC_MIXT_XVEL,c1) + &
                                pGradCell(ZCOORD,GRC_MIXT_XVEL,c2) + &
                                pGradCellOld(ZCOORD,GRC_MIXT_XVEL,c1) + &
                                pGradCellOld(ZCOORD,GRC_MIXT_XVEL,c2))

          du2dx3 = 0.25_RFREAL*(pGradCell(ZCOORD,GRC_MIXT_YVEL,c1) + &
                                pGradCell(ZCOORD,GRC_MIXT_YVEL,c2) + &
                                pGradCellOld(ZCOORD,GRC_MIXT_YVEL,c1) + &
                                pGradCellOld(ZCOORD,GRC_MIXT_YVEL,c2))

          du3dx3 = 0.25_RFREAL*(pGradCell(ZCOORD,GRC_MIXT_ZVEL,c1) + &
                                pGradCell(ZCOORD,GRC_MIXT_ZVEL,c2) + &
                                pGradCellOld(ZCOORD,GRC_MIXT_ZVEL,c1) + &
                                pGradCellOld(ZCOORD,GRC_MIXT_ZVEL,c2))

          t1N = 0.5_RFREAL*mu*(pCv(CV_MIXT_XVEL,c2)-pCv(CV_MIXT_XVEL,c1)+&
                    pCvOld(CV_MIXT_XVEL,c2)-pCvOld(CV_MIXT_XVEL,c1))/deln
          t1D = ONE_THIRD*mu*(du1dx1 + du2dx2 + du3dx3)*nx
          t1T = mu*((du2dx1*ny-du2dx2*nx) + (du3dx1*nz-du3dx3*nx))
          t2N = 0.5_RFREAL*mu*(pCv(CV_MIXT_YVEL,c2)-pCv(CV_MIXT_YVEL,c1)+&
                    pCvOld(CV_MIXT_YVEL,c2)-pCvOld(CV_MIXT_YVEL,c1))/deln
          t2D = ONE_THIRD*mu*(du1dx1 + du2dx2 + du3dx3)*ny
          t2T = mu*((du1dx2*nx-du1dx1*ny) + (du3dx2*nz-du3dx3*ny))
          t3N = 0.5_RFREAL*mu*(pCv(CV_MIXT_ZVEL,c2)-pCv(CV_MIXT_ZVEL,c1)+&
                    pCvOld(CV_MIXT_ZVEL,c2)-pCvOld(CV_MIXT_ZVEL,c1))/deln
          t3D = ONE_THIRD*mu*(du1dx1 + du2dx2 + du3dx3)*nz
          t3T = mu*((du1dx3*nx-du1dx1*nz) + (du2dx3*ny-du2dx2*nz))

          t1 = t1N + t1D + t1T
          t2 = t2N + t2D + t2T
          t3 = t3N + t3D + t3T

! ------------------------------------------------------------------------------
!           Viscous flux
! ------------------------------------------------------------------------------

#ifdef HYPRE          
          value = (refGamma-1.0_RFREAL)*refM2* &
                    (t1*u1+t2*u2+t3*u3)*iRefReNum*nm 
          iRowGlobal = c1 + nCellsOffsetRegion
          CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                         iErrHypre)

          value = -(refGamma-1.0_RFREAL)*refM2* &
                     (t1*u1+t2*u2+t3*u3)*iRefReNum*nm 
          iRowGlobal = c2 + nCellsOffsetRegion
          CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                         iErrHypre)
#endif                                         
                                         
! ------------------------------------------------------------------------------
!           Thermal flux
! ------------------------------------------------------------------------------

#ifdef HYPRE             
          value = 0.5_RFREAL*refGamma*refM2*mu*nm*iRefReNum &
                  /(Pr*deln*pCv(CV_MIXT_DENS,c1))
          iRowGlobal = c1 + nCellsOffsetRegion
          iColGlobal = iRowGlobal
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)

          value = 0.5_RFREAL*refGamma*refM2*mu*nm*iRefReNum &
                  /(Pr*deln*pCv(CV_MIXT_DENS,c2))
          iRowGlobal = c2 + nCellsOffsetRegion
          iColGlobal = iRowGlobal
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)

          value = -0.5_RFREAL*refGamma*refM2*mu*nm*iRefReNum &
                  /(Pr*deln*pCv(CV_MIXT_DENS,c2))
          iRowGlobal = c1 + nCellsOffsetRegion
          iColGlobal = c2 + pRegion%nCellsOffset(iReg1)
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)

          value = -0.5_RFREAL*refGamma*refM2*mu*nm*iRefReNum &
                  /(Pr*deln*pCv(CV_MIXT_DENS,c1))
          iRowGlobal = c2 + nCellsOffsetRegion
          iColGlobal = c1 + nCellsOffsetRegion
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)
#endif                                         

          Temp1    = MixtPerf_HM_T_DGMP(pCv(CV_MIXT_DENS,c1),refGamma,refM, &
                     pCv(CV_MIXT_PRES,c1))
          Temp2    = MixtPerf_HM_T_DGMP(pCv(CV_MIXT_DENS,c2),refGamma,refM, &
                     pCv(CV_MIXT_PRES,c2))
          Temp1Old = MixtPerf_HM_T_DGMP(pCvOld(CV_MIXT_DENS,c1),refGamma,refM, &
                     pCvOld(CV_MIXT_PRES,c1))
          Temp2Old = MixtPerf_HM_T_DGMP(pCvOld(CV_MIXT_DENS,c2),refGamma,refM, &
                     pCvOld(CV_MIXT_PRES,c2))

#ifdef HYPRE             
          value = (Temp2+Temp2Old-Temp1-Temp1Old)*mu*nm*iRefReNum &
                 /(Pr*deln*2.0_RFREAL) 
          iRowGlobal = c1 + nCellsOffsetRegion         
          CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                         iErrHypre)

          value = -(Temp2+Temp2Old-Temp1-Temp1Old)*mu*nm*iRefReNum &
                  /(Pr*deln*2.0_RFREAL) 
          iRowGlobal = c2 + nCellsOffsetRegion
          CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                         iErrHypre)
#endif                                         
        END DO ! ifg
      
        DO ifg = pGrid%nFaces-pGrid%nFacesAV+1,pGrid%nFaces
          c1 = pGrid%f2c(1,ifg)
          c2 = pGrid%f2c(2,ifg)

          direction = 1.0_RFREAL

          IF ( c1 .GT. pGrid%nCells ) THEN ! c1 is virtual
            dummy = c1
            c1    = c2
            c2    = dummy
   
            direction = -1.0_RFREAL
          END IF

          nx = direction*pGrid%fn(XCOORD,ifg)
          ny = direction*pGrid%fn(YCOORD,ifg)
          nz = direction*pGrid%fn(ZCOORD,ifg)
          nm = pGrid%fn(XYZMAG,ifg)

          deln = ( (pGrid%cofg(XCOORD,c2)-pGrid%cofg(XCOORD,c1)) &
                  *(pGrid%cofg(XCOORD,c2)-pGrid%cofg(XCOORD,c1)) &
                 + (pGrid%cofg(YCOORD,c2)-pGrid%cofg(YCOORD,c1)) &
                  *(pGrid%cofg(YCOORD,c2)-pGrid%cofg(YCOORD,c1)) &
                 + (pGrid%cofg(ZCOORD,c2)-pGrid%cofg(ZCOORD,c1)) &
                  *(pGrid%cofg(ZCOORD,c2)-pGrid%cofg(ZCOORD,c1)) )**(0.5_RFREAL)

          mu = 0.5_RFREAL*(pTv(TV_MIXT_MUEL,c1)+pTv(TV_MIXT_MUEL,c2)) 
          Pr = global%prLam 

          u1 = 0.25_RFREAL*(pCv(CV_MIXT_XVEL,c1)+pCv(CV_MIXT_XVEL,c2) &
                           +pCvOld(CV_MIXT_XVEL,c1)+pCvOld(CV_MIXT_XVEL,c2))

          u2 = 0.25_RFREAL*(pCv(CV_MIXT_YVEL,c1)+pCv(CV_MIXT_YVEL,c2) &
                           +pCvOld(CV_MIXT_YVEL,c1)+pCvOld(CV_MIXT_YVEL,c2))

          u3 = 0.25_RFREAL*(pCv(CV_MIXT_ZVEL,c1)+pCv(CV_MIXT_ZVEL,c2) &
                           +pCvOld(CV_MIXT_ZVEL,c1)+pCvOld(CV_MIXT_ZVEL,c2))

          du1dx1 = 0.25_RFREAL*(pGradCell(XCOORD,GRC_MIXT_XVEL,c1) + &
                                pGradCell(XCOORD,GRC_MIXT_XVEL,c2) + &
                                pGradCellOld(XCOORD,GRC_MIXT_XVEL,c1) + &
                                pGradCellOld(XCOORD,GRC_MIXT_XVEL,c2))

          du2dx1 = 0.25_RFREAL*(pGradCell(XCOORD,GRC_MIXT_YVEL,c1) + &
                                pGradCell(XCOORD,GRC_MIXT_YVEL,c2) + &
                                pGradCellOld(XCOORD,GRC_MIXT_YVEL,c1) + &
                                pGradCellOld(XCOORD,GRC_MIXT_YVEL,c2))

          du3dx1 = 0.25_RFREAL*(pGradCell(XCOORD,GRC_MIXT_ZVEL,c1) + &
                                pGradCell(XCOORD,GRC_MIXT_ZVEL,c2) + &
                                pGradCellOld(XCOORD,GRC_MIXT_ZVEL,c1) + &
                                pGradCellOld(XCOORD,GRC_MIXT_ZVEL,c2))

          du1dx2 = 0.25_RFREAL*(pGradCell(YCOORD,GRC_MIXT_XVEL,c1) + &
                                pGradCell(YCOORD,GRC_MIXT_XVEL,c2) + &
                                pGradCellOld(YCOORD,GRC_MIXT_XVEL,c1) + &
                                pGradCellOld(YCOORD,GRC_MIXT_XVEL,c2))

          du2dx2 = 0.25_RFREAL*(pGradCell(YCOORD,GRC_MIXT_YVEL,c1) + &
                                pGradCell(YCOORD,GRC_MIXT_YVEL,c2) + &
                                pGradCellOld(YCOORD,GRC_MIXT_YVEL,c1) + &
                                pGradCellOld(YCOORD,GRC_MIXT_YVEL,c2))

          du3dx2 = 0.25_RFREAL*(pGradCell(YCOORD,GRC_MIXT_ZVEL,c1) + &
                                pGradCell(YCOORD,GRC_MIXT_ZVEL,c2) + &
                                pGradCellOld(YCOORD,GRC_MIXT_ZVEL,c1) + &
                                pGradCellOld(YCOORD,GRC_MIXT_ZVEL,c2))

          du1dx3 = 0.25_RFREAL*(pGradCell(ZCOORD,GRC_MIXT_XVEL,c1) + &
                                pGradCell(ZCOORD,GRC_MIXT_XVEL,c2) + &
                                pGradCellOld(ZCOORD,GRC_MIXT_XVEL,c1) + &
                                pGradCellOld(ZCOORD,GRC_MIXT_XVEL,c2))

          du2dx3 = 0.25_RFREAL*(pGradCell(ZCOORD,GRC_MIXT_YVEL,c1) + &
                                pGradCell(ZCOORD,GRC_MIXT_YVEL,c2) + &
                                pGradCellOld(ZCOORD,GRC_MIXT_YVEL,c1) + &
                                pGradCellOld(ZCOORD,GRC_MIXT_YVEL,c2))

          du3dx3 = 0.25_RFREAL*(pGradCell(ZCOORD,GRC_MIXT_ZVEL,c1) + &
                                pGradCell(ZCOORD,GRC_MIXT_ZVEL,c2) + &
                                pGradCellOld(ZCOORD,GRC_MIXT_ZVEL,c1) + &
                                pGradCellOld(ZCOORD,GRC_MIXT_ZVEL,c2))

          t1N = 0.5_RFREAL*mu*(pCv(CV_MIXT_XVEL,c2)-pCv(CV_MIXT_XVEL,c1)+&
                    pCvOld(CV_MIXT_XVEL,c2)-pCvOld(CV_MIXT_XVEL,c1))/deln
          t1D = ONE_THIRD*mu*(du1dx1 + du2dx2 + du3dx3)*nx
          t1T = mu*((du2dx1*ny-du2dx2*nx) + (du3dx1*nz-du3dx3*nx))
          t2N = 0.5_RFREAL*mu*(pCv(CV_MIXT_YVEL,c2)-pCv(CV_MIXT_YVEL,c1)+&
                    pCvOld(CV_MIXT_YVEL,c2)-pCvOld(CV_MIXT_YVEL,c1))/deln
          t2D = ONE_THIRD*mu*(du1dx1 + du2dx2 + du3dx3)*ny
          t2T = mu*((du1dx2*nx-du1dx1*ny) + (du3dx2*nz-du3dx3*ny))
          t3N = 0.5_RFREAL*mu*(pCv(CV_MIXT_ZVEL,c2)-pCv(CV_MIXT_ZVEL,c1)+&
                    pCvOld(CV_MIXT_ZVEL,c2)-pCvOld(CV_MIXT_ZVEL,c1))/deln
          t3D = ONE_THIRD*mu*(du1dx1 + du2dx2 + du3dx3)*nz
          t3T = mu*((du1dx3*nx-du1dx1*nz) + (du2dx3*ny-du2dx2*nz))

          t1 = t1N + t1D + t1T
          t2 = t2N + t2D + t2T
          t3 = t3N + t3D + t3T

! ------------------------------------------------------------------------------
!           Viscous flux
! ------------------------------------------------------------------------------

          value = (refGamma-1.0_RFREAL)*refM2* &
                    (t1*u1+t2*u2+t3*u3)*iRefReNum*nm 
          iRowGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE             
          CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                         iErrHypre)
#endif

! ------------------------------------------------------------------------------
!           Thermal flux
! ------------------------------------------------------------------------------

#ifdef HYPRE             
          value = 0.5_RFREAL*refGamma*refM2*mu*nm*iRefReNum &
                  /(Pr*deln*pCv(CV_MIXT_DENS,c1))
          iRowGlobal = c1 + nCellsOffsetRegion
          iColGlobal = iRowGlobal
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre) 
                                
          value = -0.5_RFREAL*refGamma*refM2*mu*nm*iRefReNum &
                  /(Pr*deln*pCv(CV_MIXT_DENS,c2))
          iRowGlobal = c1 + nCellsOffsetRegion
          iColGlobal = pGrid%virt2GlobalIds(c2 - pGrid%nCells)
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)
#endif                                         

          Temp1    = MixtPerf_HM_T_DGMP(pCv(CV_MIXT_DENS,c1),refGamma,refM, &
                     pCv(CV_MIXT_PRES,c1))
          Temp2    = MixtPerf_HM_T_DGMP(pCv(CV_MIXT_DENS,c2),refGamma,refM, &
                     pCv(CV_MIXT_PRES,c2))
          Temp1Old = MixtPerf_HM_T_DGMP(pCvOld(CV_MIXT_DENS,c1),refGamma,refM, &
                     pCvOld(CV_MIXT_PRES,c1))
          Temp2Old = MixtPerf_HM_T_DGMP(pCvOld(CV_MIXT_DENS,c2),refGamma,refM, &
                     pCvOld(CV_MIXT_PRES,c2))

#ifdef HYPRE             
          value = (Temp2+Temp2Old-Temp1-Temp1Old)*mu*nm*iRefReNum &
                 /(Pr*deln*2.0_RFREAL) 
          iRowGlobal = c1 + nCellsOffsetRegion
          CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                         iErrHypre)
#endif                                         
        END DO ! ifg

! ==============================================================================
!       Loop over boundary faces. Energy eqn. Viscous.
! ==============================================================================

        DO iPatch = 1,pRegion%grid%nPatches
          pPatch => pRegion%patches(iPatch)

          distrib = pPatch%mixt%distrib
         
        SELECT CASE( pPatch%bcType )
        CASE (BC_SLIPWALL)
        CASE (BC_NOSLIPWALL_HFLUX,BC_NOSLIPWALL_TEMP)
          DO ifg = 1,pPatch%nBFaces
            c1 = pPatch%bf2c(ifg)

            xc = pPatch%fc(XCOORD,ifg)
            yc = pPatch%fc(YCOORD,ifg)
            zc = pPatch%fc(ZCOORD,ifg)

            nx = pPatch%fn(XCOORD,ifg)
            ny = pPatch%fn(YCOORD,ifg)
            nz = pPatch%fn(ZCOORD,ifg)
            nm = pPatch%fn(XYZMAG,ifg)

            deln = 2.0_RFREAL &
                   *( (xc-pGrid%cofg(XCOORD,c1))*(xc-pGrid%cofg(XCOORD,c1)) &
                    + (yc-pGrid%cofg(YCOORD,c1))*(yc-pGrid%cofg(YCOORD,c1)) &
                    + (zc-pGrid%cofg(ZCOORD,c1))*(zc-pGrid%cofg(ZCOORD,c1)) &
                    )**(0.5_RFREAL)

            mu = pTv(TV_MIXT_MUEL,c1)
            Pr = global%prLam 

            u1 = 0.0_RFREAL
            u2 = 0.0_RFREAL
            u3 = 0.0_RFREAL

            du1dx1 = 0.5_RFREAL*(pGradCell(XCOORD,GRC_MIXT_XVEL,c1) + &
                                 pGradCellOld(XCOORD,GRC_MIXT_XVEL,c1))

            du2dx1 = 0.5_RFREAL*(pGradCell(XCOORD,GRC_MIXT_YVEL,c1) + &
                                 pGradCellOld(XCOORD,GRC_MIXT_YVEL,c1))

            du3dx1 = 0.5_RFREAL*(pGradCell(XCOORD,GRC_MIXT_ZVEL,c1) + &
                                 pGradCellOld(XCOORD,GRC_MIXT_ZVEL,c1))

            du1dx2 = 0.5_RFREAL*(pGradCell(YCOORD,GRC_MIXT_XVEL,c1) + &
                                 pGradCellOld(YCOORD,GRC_MIXT_XVEL,c1))

            du2dx2 = 0.5_RFREAL*(pGradCell(YCOORD,GRC_MIXT_YVEL,c1) + &
                                 pGradCellOld(YCOORD,GRC_MIXT_YVEL,c1))

            du3dx2 = 0.5_RFREAL*(pGradCell(YCOORD,GRC_MIXT_ZVEL,c1) + &
                                 pGradCellOld(YCOORD,GRC_MIXT_ZVEL,c1))

            du1dx3 = 0.5_RFREAL*(pGradCell(ZCOORD,GRC_MIXT_XVEL,c1) + &
                                 pGradCellOld(ZCOORD,GRC_MIXT_XVEL,c1))

            du2dx3 = 0.5_RFREAL*(pGradCell(ZCOORD,GRC_MIXT_YVEL,c1) + &
                                 pGradCellOld(ZCOORD,GRC_MIXT_YVEL,c1))

            du3dx3 = 0.5_RFREAL*(pGradCell(ZCOORD,GRC_MIXT_ZVEL,c1) + &
                                 pGradCellOld(ZCOORD,GRC_MIXT_ZVEL,c1))

            t1N = 0.5_RFREAL*mu*(-pCv(CV_MIXT_XVEL,c1)-pCv(CV_MIXT_XVEL,c1) &
                     -pCvOld(CV_MIXT_XVEL,c1)-pCvOld(CV_MIXT_XVEL,c1))/deln
            t1D = ONE_THIRD*mu*(du1dx1 + du2dx2 + du3dx3)*nx
            t1T = mu*((du2dx1*ny-du2dx2*nx) + (du3dx1*nz-du3dx3*nx))
            t2N = 0.5_RFREAL*mu*(-pCv(CV_MIXT_YVEL,c1)-pCv(CV_MIXT_YVEL,c1) &
                     -pCvOld(CV_MIXT_YVEL,c1)-pCvOld(CV_MIXT_YVEL,c1))/deln
            t2D = ONE_THIRD*mu*(du1dx1 + du2dx2 + du3dx3)*ny
            t2T = mu*((du1dx2*nx-du1dx1*ny) + (du3dx2*nz-du3dx3*ny))
            t3N = 0.5_RFREAL*mu*(-pCv(CV_MIXT_ZVEL,c1)-pCv(CV_MIXT_ZVEL,c1) &
                     -pCvOld(CV_MIXT_ZVEL,c1)-pCvOld(CV_MIXT_ZVEL,c1))/deln
            t3D = ONE_THIRD*mu*(du1dx1 + du2dx2 + du3dx3)*nz
            t3T = mu*((du1dx3*nx-du1dx1*nz) + (du2dx3*ny-du2dx2*nz))

            t1 = t1N + t1D + t1T
            t2 = t2N + t2D + t2T
            t3 = t3N + t3D + t3T

! ------------------------------------------------------------------------------
!           Viscous flux
! ------------------------------------------------------------------------------

            value = (refGamma-1.0_RFREAL)*refM2* &
                      (t1*u1+t2*u2+t3*u3)*iRefReNum*nm 
            iRowGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE               
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                          iErrHypre)
#endif

! ------------------------------------------------------------------------------
!           Thermal flux
! ------------------------------------------------------------------------------

            SELECT CASE( pPatch%bcType )
            CASE (BC_NOSLIPWALL_HFLUX)
              qn = pPatch%mixt%vals(BCDAT_NOSLIP_Q,ifg*distrib)
              value = qn*nm*(refGamma-1.0_RFREAL) &
                      /(global%refDensity*global%refVelocity*aRef*aRef) 
            CASE (BC_NOSLIPWALL_TEMP)
              Twall    = pPatch%mixt%vals(BCDAT_NOSLIP_T,ifg*distrib) &
                         /global%refTemperature
              Temp1    = MixtPerf_HM_T_DGMP(pCv(CV_MIXT_DENS,c1),refGamma,refM, &
                         pCv(CV_MIXT_PRES,c1))
              Temp1Old = MixtPerf_HM_T_DGMP(pCvOld(CV_MIXT_DENS,c1),refGamma, &
                         refM,pCvOld(CV_MIXT_PRES,c1))
              value = (Twall+Twall-Temp1-Temp1Old)*mu*nm*iRefReNum/(Pr*deln) 
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END SELECT ! pPatch%bcType

#ifdef HYPRE   
            iRowGlobal = c1 + nCellsOffsetRegion
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                           iErrHypre)
#endif                                           

! --------- Friction coefficients ----------------------------------------------
            pPatch%cf(XCOORD,ifg) = -2.0_RFREAL*t1*iRefReNum 
            pPatch%cf(YCOORD,ifg) = -2.0_RFREAL*t2*iRefReNum
            pPatch%cf(ZCOORD,ifg) = -2.0_RFREAL*t3*iRefReNum
! --------- Heat transfer coefficients -----------------------------------------
            pPatch%ch(ifg) = 2.0_RFREAL*value/((refGamma-1.0_RFREAL)*refM2) 
          END DO ! ifg

        CASE (BC_OUTFLOW)
        CASE (BC_FARFIELD)
        CASE (BC_INFLOW_VELTEMP)
        CASE (BC_INJECTION)
        CASE (BC_PERIODIC,BC_SYMMETRY,BC_VIRTUAL)
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pPatch%bcType
        END DO ! iPatch

! ------------------------------------------------------------------------------
!     Compute coefficients due to viscous terms in axi-symmetry. Energy eqn.
! ------------------------------------------------------------------------------

        IF (pRegion%mixtInput%axiFlag .EQV. .TRUE.) THEN
          DO ifg = 1,pGrid%nFaces-pGrid%nFacesAV
            c1 = pGrid%f2c(1,ifg)
            c2 = pGrid%f2c(2,ifg)

            nx = pGrid%fn(XCOORD,ifg)
            ny = pGrid%fn(YCOORD,ifg)
            nz = pGrid%fn(ZCOORD,ifg)
            nm = pGrid%fn(XYZMAG,ifg)

            mu = 0.5_RFREAL*(pTv(TV_MIXT_MUEL,c1)+pTv(TV_MIXT_MUEL,c2)) 

            uysk = 0.25_RFREAL*(pCv(CV_MIXT_YVEL,c1)+ &
                                pCv(CV_MIXT_YVEL,c2)+ &
                                pCvOld(CV_MIXT_YVEL,c1)+ &
                                pCvOld(CV_MIXT_YVEL,c2))/ &
                               (pRegion%grid%fc(YCOORD,ifg)+EPS)

            v_nHalf = (pVfMixt(ifg)+pVfMixtOld(ifg))/2.0_RFREAL

#ifdef HYPRE             
            value = -(refGamma-1.0_RFREAL)*refM2*TWO_THIRD*mu*iRefREnum* &
                     uysk*v_nHalf*nm
            iRowGlobal = c1 + nCellsOffsetRegion
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal, &
                                           value,iErrHypre)

            value = (refGamma-1.0_RFREAL)*refM2*TWO_THIRD*mu*iRefREnum* &
                     uysk*v_nHalf*nm
            iRowGlobal = c2 + nCellsOffsetRegion
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal, &
                                           value,iErrHypre)
#endif                                         
          END DO ! ifg

          DO ifg = pGrid%nFaces-pGrid%nFacesAV+1,pGrid%nFaces
            c1 = pGrid%f2c(1,ifg)
            c2 = pGrid%f2c(2,ifg)

            direction = 1.0_RFREAL

            IF ( c1 .GT. pGrid%nCells ) THEN ! c1 is virtual
              dummy = c1
              c1    = c2
              c2    = dummy

              direction = -1.0_RFREAL
            END IF

            nx = direction*pGrid%fn(XCOORD,ifg)
            ny = direction*pGrid%fn(YCOORD,ifg)
            nz = direction*pGrid%fn(ZCOORD,ifg)
            nm = pGrid%fn(XYZMAG,ifg)

            mu = 0.5_RFREAL*(pTv(TV_MIXT_MUEL,c1)+pTv(TV_MIXT_MUEL,c2)) 

            uysk = 0.25_RFREAL*(pCv(CV_MIXT_YVEL,c1)+ &
                                pCv(CV_MIXT_YVEL,c2)+ &
                                pCvOld(CV_MIXT_YVEL,c1)+ &
                                pCvOld(CV_MIXT_YVEL,c2))/ &
                               (pRegion%grid%fc(YCOORD,ifg)+EPS)

            v_nHalf = (pVfMixt(ifg)+pVfMixtOld(ifg))/2.0_RFREAL

#ifdef HYPRE             
            value = -(refGamma-1.0_RFREAL)*refM2*TWO_THIRD*mu*iRefREnum* &
                     uysk*v_nHalf*nm
            iRowGlobal = c1 + nCellsOffsetRegion
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal, &
                                           value,iErrHypre)
#endif                                         
          END DO ! ifg

          DO iPatch = 1,pRegion%grid%nPatches
            pPatch => pRegion%patches(iPatch)
   
            distrib = pPatch%mixt%distrib
 
            SELECT CASE( pPatch%bcType )
            CASE (BC_SLIPWALL)
            CASE (BC_NOSLIPWALL_HFLUX,BC_NOSLIPWALL_TEMP)
            CASE (BC_OUTFLOW)
            CASE (BC_FARFIELD)
              DO ifg = 1,pPatch%nBFaces
                c1 = pPatch%bf2c(ifg)

                xc = pPatch%fc(XCOORD,ifg)
                yc = pPatch%fc(YCOORD,ifg)
                zc = pPatch%fc(ZCOORD,ifg)

                nx = pPatch%fn(XCOORD,ifg)
                ny = pPatch%fn(YCOORD,ifg)
                nz = pPatch%fn(ZCOORD,ifg)
                nm = pPatch%fn(XYZMAG,ifg)

                mu = pTv(TV_MIXT_MUEL,c1)

                uysk = 0.5_RFREAL*(pCv(CV_MIXT_YVEL,c1)+ &
                                   pCvOld(CV_MIXT_YVEL,c1))/ &
                                  (pPatch%fc(YCOORD,ifg)+EPS)

                v_nHalf = (pVfMixt(ifg)+pVfMixtOld(ifg))/2.0_RFREAL

#ifdef HYPRE             
                value = -(refGamma-1.0_RFREAL)*refM2*TWO_THIRD*mu*iRefREnum* &
                         uysk*v_nHalf*nm
                iRowGlobal = c1 + nCellsOffsetRegion
                CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal, &
                                               value,iErrHypre)
#endif                                           
              END DO ! ifg
            CASE (BC_INFLOW_VELTEMP)
            CASE (BC_INJECTION)
            CASE (BC_PERIODIC,BC_SYMMETRY,BC_VIRTUAL)
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END SELECT ! pPatch%bcType
          END DO ! iPatch
        END IF ! pRegion%mixtInput%axiFlag

      END IF ! pRegion%mixtInput%flowModel
    END DO ! iReg

! ==============================================================================
!   Assemble matrix and vector
! ==============================================================================

    CALL RFLU_HYPRE_AssembleMatrixVector(regions)

! ==============================================================================
!   Solve system of equation for density using HYPRE 
! ==============================================================================

    CALL RFLU_HYPRE_SolveMatrix(regions)

! ==============================================================================
!   extract solution from HYPRE solution vector parSolHypre
! ==============================================================================

! TEMPORARY: manoj: extract convergence data
! --- Extract convergence data -------------------------------------------------
!    CALL HYPRE_BoomerAMGGetNumIterations(regions(1)%hypreSolver,numIter, &
!                                         iErrHypre)
!    CALL HYPRE_BoomerAMGGetFinalRelativeResidualNorm(regions(1)%hypreSolver, &
!                                                     relRes,iErrHypre)
!    WRITE(*,*) 'Energy Eq: numIter=',numIter,' relative residual norm=',relRes
! END TEMPORARY

    global => regions(1)%global

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid

      iReg1 = pRegion%iRegionGlobal
      nCellsOffsetRegion = pRegion%nCellsOffset(iReg1) 

! --- Extract actual cell data -------------------------------------------------
      DO icg = 1,pGrid%nCells
        iRowGlobal = icg + nCellsOffsetRegion 
#ifdef HYPRE                  
        CALL HYPRE_IJVectorGetValues(regions(1)%SolHypre,1,iRowGlobal,value, &
                                     iErrHypre)
#endif                                     
        pRegion%mixt%delP(icg) = value
      END DO ! icg
    END DO ! iReg

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_SolveEnergyEq









! ******************************************************************************
!
! Purpose: Solve momentum equation in non-dissipative scheme by Hou and 
!          Mahesh.
!
! Description: None.
!
! Input:
!   regions            Pointer to region data
!   ICORD             I=X,Y,Z direction
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_HM_SolveMomentumEq(regions,ICORD)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

! *****************************************************************************
!   Declarations and definitions
! *****************************************************************************

! =============================================================================
!   Arguments
! =============================================================================

    INTEGER, INTENT(IN) :: ICORD
    TYPE(t_region), DIMENSION(:), POINTER :: regions
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: CV_MIXT_IVEL,CV_MIXT_JVEL,CV_MIXT_KVEL, &
               GRC_MIXT_IVEL,GRC_MIXT_JVEL,GRC_MIXT_KVEL,JCORD, &
               KCORD
    INTEGER :: c1,c2,distrib,dummy,errorFlag,icg,iColGlobal,iErrHypre,ifg, &
               iPatch,iRow,iRowGlobal,iReg,iReg1,nCellsOffsetRegion
    REAL(RFREAL), PARAMETER :: ONE_THIRD = 1.0_RFREAL/3.0_RFREAL
    REAL(RFREAL), PARAMETER :: TWO_THIRD = 2.0_RFREAL/3.0_RFREAL
    REAL(RFREAL), PARAMETER :: EPS = 1.0E-15_RFREAL
    REAL(RFREAL) :: deln,delT,direction,du1dx1,du1dx2,du1dx3,du2dx1,du2dx2, &
                    du2dx3,du3dx1,du3dx2,du3dx3,duidxi,dujdxi,dujdxj,dukdxi, &
                    dukdxk,iDt,iRefMu,iRefReNum,mu,nm,ni,nj,nk,nx,ny,nz, &
                    refGamma,refM2,rho_n,rhoAlpha_n,rhoAlpha_np1,rhoBeta_n, &
                    rhoBeta_np1,t1D,t1N,t1T,t2D,t2N,t2T,t3D,t3N,t3T,tiD,tiN, &
                    tiT,uysk,v_nHalf,value,vfnmi4,voliDt,volume,xc,y,yc,zc
    REAL(RFREAL) :: aoa,aos,aRef,mf,minj,SpeedOfsound,tf,tinj,nxi,uo,vo,wo, &
                    ro,po,sigma,VelX,VelY,VelZ
    REAL(RFREAL), DIMENSION(:), POINTER :: pVfMixt,pVfMixtOld
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pCvOld,pCvOld2,pGradP, &
                                             pGradPOld,pGradPOld2,pTv
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGradCell,pGradCellOld, &
                                               pGradCellOld2
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion
 
! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => regions(1)%global
 
    CALL RegisterFunction(global,'RFLU_HM_SolveMomentumEq',__FILE__)

    delT      = global%dtMin*global%refVelocity/global%refLength
    iDt       = 1.0_RFREAL/delT
    iRefMu    = 1.0_RFREAL/global%refVisc
    iRefReNum = 1.0_RFREAL/global%refReNum
    refGamma  = global%refGamma 
    aRef      = (refGamma*global%refPressure/global%refDensity)**0.5_RFREAL
    refM2     = (global%refVelocity/aRef)**2.0_RFREAL

    SELECT CASE( ICORD )
      CASE ( XCOORD )          
        JCORD = YCOORD
        KCORD = ZCOORD
        
        CV_MIXT_IVEL = CV_MIXT_XVEL
        CV_MIXT_JVEL = CV_MIXT_YVEL
        CV_MIXT_KVEL = CV_MIXT_ZVEL

        GRC_MIXT_IVEL = GRC_MIXT_XVEL
        GRC_MIXT_JVEL = GRC_MIXT_YVEL
        GRC_MIXT_KVEL = GRC_MIXT_ZVEL
      CASE ( YCOORD )          
        JCORD = ZCOORD
        KCORD = XCOORD
        
        CV_MIXT_IVEL = CV_MIXT_YVEL
        CV_MIXT_JVEL = CV_MIXT_ZVEL
        CV_MIXT_KVEL = CV_MIXT_XVEL

        GRC_MIXT_IVEL = GRC_MIXT_YVEL
        GRC_MIXT_JVEL = GRC_MIXT_ZVEL
        GRC_MIXT_KVEL = GRC_MIXT_XVEL
      CASE ( ZCOORD )          
        JCORD = XCOORD
        KCORD = YCOORD
        
        CV_MIXT_IVEL = CV_MIXT_ZVEL
        CV_MIXT_JVEL = CV_MIXT_XVEL
        CV_MIXT_KVEL = CV_MIXT_YVEL

        GRC_MIXT_IVEL = GRC_MIXT_ZVEL
        GRC_MIXT_JVEL = GRC_MIXT_XVEL
        GRC_MIXT_KVEL = GRC_MIXT_YVEL
      CASE DEFAULT       
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! ICORD 

! ******************************************************************************
!   Build system of equation for i-momentum
! ******************************************************************************

! ==============================================================================
!   Initialize Hypre matrix and vectors
! ==============================================================================

    CALL RFLU_HYPRE_InitializeMatrixVector(regions)

! ==============================================================================
!   loop over regions and assign pointers
! ==============================================================================

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid

      pVfMixt    => pRegion%mixt%vfMixt
      pVfMixtOld => pRegion%mixt%vfMixtOld
      pCv        => pRegion%mixt%cv
      pCvOld     => pRegion%mixt%cvOld
      pCvOld2    => pRegion%mixt%cvOld2
      pGradP     => pRegion%mixt%gradCell(:,GRC_MIXT_PRES,:)
      pGradPOld  => pRegion%mixt%gradCellOld(:,GRC_MIXT_PRES,:)
      pGradPOld2 => pRegion%mixt%gradCellOld2(:,GRC_MIXT_PRES,:)

      pGradCell     => pRegion%mixt%gradCell
      pGradCellOld  => pRegion%mixt%gradCellOld
      pGradCellOld2 => pRegion%mixt%gradCellOld2

      IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN
        pTv => pRegion%mixt%tv
      END IF ! pRegion%mixtInput%flowModel

      iReg1 = pRegion%iRegionGlobal
      nCellsOffsetRegion = pRegion%nCellsOffset(iReg1)

! ------------------------------------------------------------------------------
!     Loop over cells and initialize the coefficients
! ------------------------------------------------------------------------------

      DO icg = 1,pGrid%nCells
        volume = pGrid%vol(icg)

        rho_n = (pCvOld(CV_MIXT_DENS,icg)+pCvOld2(CV_OLD2_MIXT_DENS,icg)) &
                /2.0_RFREAL

        voliDt = volume*iDt

        value = voliDt
        iRowGlobal = icg + nCellsOffsetRegion
        iColGlobal = iRowGlobal
#ifdef HYPRE           
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value = rho_n*pCvOld(CV_MIXT_IVEL,icg)*voliDt &
                - (volume/4.0_RFREAL)*( pGradP(ICORD,icg) &
                                      + 2.0_RFREAL*pGradPOld(ICORD,icg) &
                                      + pGradPOld2(ICORD,icg) )
        iRowGlobal = icg + nCellsOffsetRegion
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif                                       
      END DO ! icg

! ------------------------------------------------------------------------------
!     Loop over interior faces and accumulate coefficients
! ------------------------------------------------------------------------------

      DO ifg = 1,pGrid%nFaces-pGrid%nFacesAV
        c1 = pGrid%f2c(1,ifg)
        c2 = pGrid%f2c(2,ifg)

        nm = pGrid%fn(XYZMAG,ifg)

        rhoAlpha_n = (pCvOld(CV_MIXT_DENS,c1)+pCvOld2(CV_OLD2_MIXT_DENS,c1)) &
                     /2.0_RFREAL
        rhoBeta_n  = (pCvOld(CV_MIXT_DENS,c2)+pCvOld2(CV_OLD2_MIXT_DENS,c2)) &
                     /2.0_RFREAL
        v_nHalf    = (pVfMixt(ifg)+pVfMixtOld(ifg))/2.0_RFREAL

        vfnmi4 = v_nHalf*nm/4.0_RFREAL

        value = vfnmi4
        iRowGlobal = c1 + pRegion%nCellsOffset(iReg1)
        iColGlobal = iRowGlobal
#ifdef HYPRE   
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value = - vfnmi4
        iRowGlobal = c2 + pRegion%nCellsOffset(iReg1)
        iColGlobal = iRowGlobal
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
 
        value = vfnmi4
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = c2 + nCellsOffsetRegion
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value = - vfnmi4
        iRowGlobal = c2 + nCellsOffsetRegion
        iColGlobal = c1 + nCellsOffsetRegion
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value = -( rhoAlpha_n*pCvOld(CV_MIXT_IVEL,c1) &
                 + rhoBeta_n*pCvOld(CV_MIXT_IVEL,c2) )*vfnmi4
        iRowGlobal = c1 + nCellsOffsetRegion
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)

        value = +( rhoAlpha_n*pCvOld(CV_MIXT_IVEL,c1) &
                 + rhoBeta_n*pCvOld(CV_MIXT_IVEL,c2) )*vfnmi4
        iRowGlobal = c2 + nCellsOffsetRegion
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif                                       
      END DO ! ifg

      DO ifg = pGrid%nFaces-pGrid%nFacesAV+1,pGrid%nFaces
        c1 = pGrid%f2c(1,ifg)
        c2 = pGrid%f2c(2,ifg)

        direction = 1.0_RFREAL

        IF ( c1 .GT. pGrid%nCells ) THEN ! c1 is virtual
          dummy = c1
          c1    = c2
          c2    = dummy

          direction = -1.0_RFREAL
        END IF

        nm = pGrid%fn(XYZMAG,ifg)

        rhoAlpha_n = (pCvOld(CV_MIXT_DENS,c1)+pCvOld2(CV_OLD2_MIXT_DENS,c1)) &
                     /2.0_RFREAL
        rhoBeta_n  = (pCvOld(CV_MIXT_DENS,c2)+pCvOld2(CV_OLD2_MIXT_DENS,c2)) &
                     /2.0_RFREAL
        v_nHalf    = direction*(pVfMixt(ifg)+pVfMixtOld(ifg))/2.0_RFREAL

        vfnmi4 = v_nHalf*nm/4.0_RFREAL

        value = vfnmi4
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = iRowGlobal
#ifdef HYPRE           
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value =  vfnmi4
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = pGrid%virt2GlobalIds(c2 - pGrid%nCells)
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)

        value = -( rhoAlpha_n*pCvOld(CV_MIXT_IVEL,c1) &
                 + rhoBeta_n*pCvOld(CV_MIXT_IVEL,c2) )*vfnmi4
        iRowGlobal = c1 + nCellsOffsetRegion
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif                                       
      END DO ! ifg

! ------------------------------------------------------------------------------
!     Loop over boundary faces. Inviscid. i-momentum eqn.
! ------------------------------------------------------------------------------

      DO iPatch = 1,pRegion%grid%nPatches
        pPatch => pRegion%patches(iPatch)

        distrib = pPatch%mixt%distrib
 
        SELECT CASE( pPatch%bcType )
        CASE (BC_SLIPWALL,BC_NOSLIPWALL_HFLUX,BC_NOSLIPWALL_TEMP)
        CASE (BC_OUTFLOW)
          DO ifg = 1,pPatch%nBFaces
            c1 = pPatch%bf2c(ifg)

            nx = pPatch%fn(XCOORD,ifg)
            ny = pPatch%fn(YCOORD,ifg)
            nz = pPatch%fn(ZCOORD,ifg)
            nm = pPatch%fn(XYZMAG,ifg)

            rhoAlpha_n = (pCvOld(CV_MIXT_DENS,c1)+pCvOld2(CV_OLD2_MIXT_DENS,c1)) &
                         /2.0_RFREAL

            v_nHalf = ((pCv(CV_MIXT_XVEL,c1)+pCvOld(CV_MIXT_XVEL,c1))*nx &
                      +(pCv(CV_MIXT_YVEL,c1)+pCvOld(CV_MIXT_YVEL,c1))*ny &
                      +(pCv(CV_MIXT_ZVEL,c1)+pCvOld(CV_MIXT_ZVEL,c1))*nz)/2.0_RFREAL

            vfnmi4 = v_nHalf*nm/4.0_RFREAL

            value = vfnmi4
            iRowGlobal = c1 + pRegion%nCellsOffset(iReg1)
            iColGlobal = iRowGlobal
#ifdef HYPRE   
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)
#endif                                       

            value = vfnmi4
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE   
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)
#endif                                       

            value = -( rhoAlpha_n*pCvOld(CV_MIXT_IVEL,c1) &
                     + rhoAlpha_n*pCvOld(CV_MIXT_IVEL,c1) )*vfnmi4
            iRowGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE   
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                           iErrHypre)
#endif                                       
          END DO ! ifg

        CASE (BC_FARFIELD)
          DO ifg = 1,pPatch%nBFaces
            c1 = pPatch%bf2c(ifg)

            xc = pPatch%fc(XCOORD,ifg)
            yc = pPatch%fc(YCOORD,ifg)
            zc = pPatch%fc(ZCOORD,ifg)
            
            nx = pPatch%fn(XCOORD,ifg)
            ny = pPatch%fn(YCOORD,ifg)
            nz = pPatch%fn(ZCOORD,ifg)
            nm = pPatch%fn(XYZMAG,ifg)
            
            deln = 2.0_RFREAL &
                   *( (xc-pGrid%cofg(XCOORD,c1))*(xc-pGrid%cofg(XCOORD,c1)) &
                    + (yc-pGrid%cofg(YCOORD,c1))*(yc-pGrid%cofg(YCOORD,c1)) &
                    + (zc-pGrid%cofg(ZCOORD,c1))*(zc-pGrid%cofg(ZCOORD,c1)) &
                    )**(0.5_RFREAL)

            
            rhoAlpha_n = (pCvOld(CV_MIXT_DENS,c1) &
                         +pCvOld2(CV_OLD2_MIXT_DENS,c1))/2.0_RFREAL
            rhoBeta_n  = rhoAlpha_n
            
!            v_nHalf = (pPatch%mixt%vfMixt(ifg) &
!                      +pPatch%mixt%vfMixtOld(ifg))/2.0_RFREAL

            mf  = pPatch%mixt%vals(BCDAT_FARF_MACH,distrib*ifg)
            aoa = pPatch%mixt%vals(BCDAT_FARF_ATTACK,distrib*ifg)
            aos = pPatch%mixt%vals(BCDAT_FARF_SLIP,distrib*ifg)
            tf  = pPatch%mixt%vals(BCDAT_FARF_TEMP,distrib*ifg)

            SpeedOfSound = SQRT(tf/global%refTemperature)

            VelX = mf*SpeedOfSound*COS(aoa)*COS(aos)
            VelY = mf*SpeedOfSound*SIN(aoa)*COS(aos)
            VelZ = mf*SpeedOfSound*SIN(aos)

            vfnmi4 = (nx*VelX+ny*VelY+nz*VelZ)*nm/4.0_RFREAL
 
#ifdef HYPRE   
            value = vfnmi4
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = iRowGlobal
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)
        
! Boundary rho*u is known, so it goes to RHS ------------------------------    
            value =  - vfnmi4*pRegion%mixt%cv(CV_MIXT_DENS,ifg) &
                             *pRegion%mixt%cv(CV_MIXT_IVEL,ifg)
            iRowGlobal = c1 + nCellsOffsetRegion
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal, &
                                           value,iErrHypre)

            value = -( rhoAlpha_n*pCvOld(CV_MIXT_IVEL,c1) &
                     + rhoBeta_n*pRegion%mixt%cvOld(CV_MIXT_IVEL,ifg) )*vfnmi4
            iRowGlobal = c1 + nCellsOffsetRegion
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal, &
                                           value,iErrHypre)
#endif                                       
          END DO ! ifg

        CASE (BC_INFLOW_VELTEMP)
          DO ifg = 1,pPatch%nBFaces
            c1 = pPatch%bf2c(ifg)

            nxi = pPatch%fn(ICORD,ifg)
            nm  = pPatch%fn(XYZMAG,ifg)

            minj = -pPatch%mixt%vals(BCDAT_INJECT_MFRATE,distrib*ifg)
            minj = minj/(global%refDensity*global%refVelocity)
            tinj = pPatch%mixt%vals(BCDAT_INJECT_TEMP,distrib*ifg)
            tinj = tinj/global%refTemperature

            v_nHalf = minj*tinj/(refGamma*refM2*pCvOld(CV_MIXT_PRES,c1)+1.0_RFREAL)

            value = -minj*(v_nHalf*nxi)*nm
            iRowGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE   
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif                                       
          END DO ! ifg

        CASE (BC_INJECTION)
          DO ifg = 1,pPatch%nBFaces
            c1 = pPatch%bf2c(ifg)

            nxi = pPatch%fn(ICORD,ifg)
            nm  = pPatch%fn(XYZMAG,ifg)

            minj = -pPatch%mixt%vals(BCDAT_INJECT_MFRATE,distrib*ifg)
            minj = minj/(global%refDensity*global%refVelocity)
            tinj = pPatch%mixt%vals(BCDAT_INJECT_TEMP,distrib*ifg)
            tinj = tinj/global%refTemperature

            v_nHalf = minj*tinj/(refGamma*refM2*pCvOld(CV_MIXT_PRES,c1)+1.0_RFREAL)

            value = -minj*(v_nHalf*nxi)*nm
            iRowGlobal = c1 + nCellsOffsetRegion
#ifdef HYPRE   
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif                                       
          END DO ! ifg

        CASE (BC_PERIODIC,BC_SYMMETRY,BC_VIRTUAL)
        CASE DEFAULT                
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pPatch%bcType
      END DO ! iPatch

! ------------------------------------------------------------------------------
!     Compute coefficients due to moving reference frame source terms. Mom eqn.
! ------------------------------------------------------------------------------
 
      IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
        DO icg = 1,pGrid%nCells
          volume = pGrid%vol(icg)

#ifdef HYPRE           
         value = -volume*pCvOld(CV_MIXT_DENS,icg)*pRegion%mvfAcc(ICORD)
          iRowGlobal = icg + nCellsOffsetRegion
          CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                         iErrHypre)
#endif                                       
        END DO ! icg
      END IF ! global%mvFrameFlag

! ------------------------------------------------------------------------------
!     Compute coefficients due to absorbing boundary condition. Momentum eqn. 
! ------------------------------------------------------------------------------

      IF ( (global%abcFlag .EQV. .TRUE.) .AND. (global%abcKind == 0) ) THEN
        DO icg = 1,pGrid%nCells
          distrib = global%abcDistrib
          ro = pRegion%mixt%cvRef(CV_MIXT_DENS,distrib*icg)
          uo = pRegion%mixt%cvRef(CV_MIXT_XVEL,distrib*icg)
          vo = pRegion%mixt%cvRef(CV_MIXT_YVEL,distrib*icg)
          wo = pRegion%mixt%cvRef(CV_MIXT_ZVEL,distrib*icg)
          po = pRegion%mixt%cvRef(CV_MIXT_PRES,distrib*icg)

          SELECT CASE ( ICORD )
            CASE ( XCOORD )
              uo = uo
            CASE ( YCOORD )
              uo = vo
            CASE ( ZCOORD )
              uo = wo
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! ICORD 

          volume = pGrid%vol(icg)

          sigma = pRegion%mixt%sigma(icg)*global%refLength/global%refVelocity 

          rhoAlpha_n = (pCvOld(CV_MIXT_DENS,icg) &
                       +pCvOld2(CV_OLD2_MIXT_DENS,icg))/2.0_RFREAL

          value = sigma*volume*0.5_RFREAL
          iRowGlobal = icg + nCellsOffsetRegion 
          iColGlobal = iRowGlobal
#ifdef HYPRE        
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)
#endif

          value = -(rhoAlpha_n*pCvOld(CV_MIXT_IVEL,icg)*0.5_RFREAL-ro*uo) &
                  *sigma*volume
          iRowGlobal = icg + nCellsOffsetRegion
#ifdef HYPRE           
          CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                         iErrHypre)
#endif                                       
        END DO ! icg
      END IF ! global%abcFlag

! ------------------------------------------------------------------------------
!     Compute viscous flux for i-momentum.
!     Loop over interior faces and accumulate coefficients
! ------------------------------------------------------------------------------

      IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN
        DO ifg = 1,pGrid%nFaces-pGrid%nFacesAV
          c1 = pGrid%f2c(1,ifg)
          c2 = pGrid%f2c(2,ifg)

          ni = pGrid%fn(ICORD,ifg)
          nj = pGrid%fn(JCORD,ifg)
          nk = pGrid%fn(KCORD,ifg)
          nm = pGrid%fn(XYZMAG,ifg)
 
          deln = ( (pGrid%cofg(ICORD,c2)-pGrid%cofg(ICORD,c1)) &
                  *(pGrid%cofg(ICORD,c2)-pGrid%cofg(ICORD,c1)) &
                 + (pGrid%cofg(JCORD,c2)-pGrid%cofg(JCORD,c1)) &
                  *(pGrid%cofg(JCORD,c2)-pGrid%cofg(JCORD,c1)) &
                 + (pGrid%cofg(KCORD,c2)-pGrid%cofg(KCORD,c1)) &
                  *(pGrid%cofg(KCORD,c2)-pGrid%cofg(KCORD,c1)) &
                 )**(0.5_RFREAL)

          mu = 0.5_RFREAL*(pTv(TV_MIXT_MUEL,c1)+pTv(TV_MIXT_MUEL,c2)) 

          rhoAlpha_np1 = (pCv(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c1)) &
                       /2.0_RFREAL
          rhoBeta_np1  = (pCv(CV_MIXT_DENS,c2)+pCvOld(CV_MIXT_DENS,c2)) &
                       /2.0_RFREAL

          duidxi = 0.25_RFREAL*(pGradCell(ICORD,GRC_MIXT_IVEL,c1) + &
                                pGradCell(ICORD,GRC_MIXT_IVEL,c2) + &
                                pGradCellOld(ICORD,GRC_MIXT_IVEL,c1) + &
                                pGradCellOld(ICORD,GRC_MIXT_IVEL,c2))

          dujdxi = 0.25_RFREAL*(pGradCell(ICORD,GRC_MIXT_JVEL,c1) + &
                                pGradCell(ICORD,GRC_MIXT_JVEL,c2) + &
                                pGradCellOld(ICORD,GRC_MIXT_JVEL,c1) + &
                                pGradCellOld(ICORD,GRC_MIXT_JVEL,c2))

          dukdxi = 0.25_RFREAL*(pGradCell(ICORD,GRC_MIXT_KVEL,c1) + &
                                pGradCell(ICORD,GRC_MIXT_KVEL,c2) + &
                                pGradCellOld(ICORD,GRC_MIXT_KVEL,c1) + &
                                pGradCellOld(ICORD,GRC_MIXT_KVEL,c2))

          dujdxj = 0.25_RFREAL*(pGradCell(JCORD,GRC_MIXT_JVEL,c1) + &
                                pGradCell(JCORD,GRC_MIXT_JVEL,c2) + &
                                pGradCellOld(JCORD,GRC_MIXT_JVEL,c1) + &
                                pGradCellOld(JCORD,GRC_MIXT_JVEL,c2))

          dukdxk = 0.25_RFREAL*(pGradCell(KCORD,GRC_MIXT_KVEL,c1) + &
                                pGradCell(KCORD,GRC_MIXT_KVEL,c2) + &
                                pGradCellOld(KCORD,GRC_MIXT_KVEL,c1) + &
                                pGradCellOld(KCORD,GRC_MIXT_KVEL,c2))

          tiN = mu*(pCv(CV_MIXT_IVEL,c2)-pCv(CV_MIXT_IVEL,c1) &
                   +pCvOld(CV_MIXT_IVEL,c2)-pCvOld(CV_MIXT_IVEL,c1)) &
                  /(2.0_RFREAL*deln)
          tiD = ONE_THIRD*mu*(duidxi + dujdxj + dukdxk)*ni
          tiT = mu*((dujdxi*nj-dujdxj*ni) + (dukdxi*nk-dukdxk*ni))

#ifdef HYPRE             
          value = nm*iRefREnum*mu/(2.0_RFREAL*deln*rhoAlpha_np1)
          iRowGlobal = c1 + pRegion%nCellsOffset(iReg1)
          iColGlobal = iRowGlobal
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)

          value = nm*iRefREnum*mu/(2.0_RFREAL*deln*rhoBeta_np1)
          iRowGlobal = c2 + pRegion%nCellsOffset(iReg1)
          iColGlobal = iRowGlobal
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)

          value = -nm*iRefREnum*mu/(2.0_RFREAL*deln*rhoBeta_np1)
          iRowGlobal = c1 + nCellsOffsetRegion
          iColGlobal = c2 + nCellsOffsetRegion
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)

          value = -nm*iRefREnum*mu/(2.0_RFREAL*deln*rhoAlpha_np1)
          iRowGlobal = c2 + nCellsOffsetRegion
          iColGlobal = c1 + nCellsOffsetRegion
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)

          value = iRefREnum*(tiD+tiT+mu*(pCvOld(CV_MIXT_IVEL,c2) &
                                         -pCvOld(CV_MIXT_IVEL,c1)) &
                                        /(2.0_RFREAL*deln))*nm
          iRowGlobal = c1 + nCellsOffsetRegion
          CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                         iErrHypre)

          value = -iRefREnum*(tiD+tiT+mu*(pCvOld(CV_MIXT_IVEL,c2) &
                                        -pCvOld(CV_MIXT_IVEL,c1)) &
                                       /(2.0_RFREAL*deln))*nm
          iRowGlobal = c2 + nCellsOffsetRegion
          CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                         iErrHypre)
#endif                                         
        END DO ! ifg

        DO ifg = pGrid%nFaces-pGrid%nFacesAV+1,pGrid%nFaces
          c1 = pGrid%f2c(1,ifg)
          c2 = pGrid%f2c(2,ifg)

          direction = 1.0_RFREAL

          IF ( c1 .GT. pGrid%nCells ) THEN ! c1 is virtual
            dummy = c1
            c1    = c2
            c2    = dummy

            direction = -1.0_RFREAL
          END IF

          ni = direction*pGrid%fn(ICORD,ifg)
          nj = direction*pGrid%fn(JCORD,ifg)
          nk = direction*pGrid%fn(KCORD,ifg)
          nm = pGrid%fn(XYZMAG,ifg)
 
          deln = ( (pGrid%cofg(ICORD,c2)-pGrid%cofg(ICORD,c1)) &
                  *(pGrid%cofg(ICORD,c2)-pGrid%cofg(ICORD,c1)) &
                 + (pGrid%cofg(JCORD,c2)-pGrid%cofg(JCORD,c1)) &
                  *(pGrid%cofg(JCORD,c2)-pGrid%cofg(JCORD,c1)) &
                 + (pGrid%cofg(KCORD,c2)-pGrid%cofg(KCORD,c1)) &
                  *(pGrid%cofg(KCORD,c2)-pGrid%cofg(KCORD,c1)) &
                 )**(0.5_RFREAL)

          mu = 0.5_RFREAL*(pTv(TV_MIXT_MUEL,c1)+pTv(TV_MIXT_MUEL,c2)) 

          rhoAlpha_np1 = (pCv(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c1)) &
                       /2.0_RFREAL
          rhoBeta_np1  = (pCv(CV_MIXT_DENS,c2)+pCvOld(CV_MIXT_DENS,c2)) &
                       /2.0_RFREAL

          duidxi = 0.25_RFREAL*(pGradCell(ICORD,GRC_MIXT_IVEL,c1) + &
                                pGradCell(ICORD,GRC_MIXT_IVEL,c2) + &
                                pGradCellOld(ICORD,GRC_MIXT_IVEL,c1) + &
                                pGradCellOld(ICORD,GRC_MIXT_IVEL,c2))

          dujdxi = 0.25_RFREAL*(pGradCell(ICORD,GRC_MIXT_JVEL,c1) + &
                                pGradCell(ICORD,GRC_MIXT_JVEL,c2) + &
                                pGradCellOld(ICORD,GRC_MIXT_JVEL,c1) + &
                                pGradCellOld(ICORD,GRC_MIXT_JVEL,c2))

          dukdxi = 0.25_RFREAL*(pGradCell(ICORD,GRC_MIXT_KVEL,c1) + &
                                pGradCell(ICORD,GRC_MIXT_KVEL,c2) + &
                                pGradCellOld(ICORD,GRC_MIXT_KVEL,c1) + &
                                pGradCellOld(ICORD,GRC_MIXT_KVEL,c2))

          dujdxj = 0.25_RFREAL*(pGradCell(JCORD,GRC_MIXT_JVEL,c1) + &
                                pGradCell(JCORD,GRC_MIXT_JVEL,c2) + &
                                pGradCellOld(JCORD,GRC_MIXT_JVEL,c1) + &
                                pGradCellOld(JCORD,GRC_MIXT_JVEL,c2))

          dukdxk = 0.25_RFREAL*(pGradCell(KCORD,GRC_MIXT_KVEL,c1) + &
                                pGradCell(KCORD,GRC_MIXT_KVEL,c2) + &
                                pGradCellOld(KCORD,GRC_MIXT_KVEL,c1) + &
                                pGradCellOld(KCORD,GRC_MIXT_KVEL,c2))

          tiN = mu*(pCv(CV_MIXT_IVEL,c2)-pCv(CV_MIXT_IVEL,c1) &
                   +pCvOld(CV_MIXT_IVEL,c2)-pCvOld(CV_MIXT_IVEL,c1)) &
                  /(2.0_RFREAL*deln)
          tiD = ONE_THIRD*mu*(duidxi + dujdxj + dukdxk)*ni
          tiT = mu*((dujdxi*nj-dujdxj*ni) + (dukdxi*nk-dukdxk*ni))

          value = nm*iRefREnum*mu/(2.0_RFREAL*deln*rhoAlpha_np1)
          iRowGlobal = c1 + pRegion%nCellsOffset(iReg1)
          iColGlobal = iRowGlobal
#ifdef HYPRE             
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)

          value = -nm*iRefREnum*mu/(2.0_RFREAL*deln*rhoBeta_np1)
          iRowGlobal = c1 + nCellsOffsetRegion
          iColGlobal = pGrid%virt2GlobalIds(c2 - pGrid%nCells)
          CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                         iColGlobal,value,iErrHypre)

          value = iRefREnum*(tiD+tiT+mu*(pCvOld(CV_MIXT_IVEL,c2) &
                                         -pCvOld(CV_MIXT_IVEL,c1)) &
                                        /(2.0_RFREAL*deln))*nm
          iRowGlobal = c1 + nCellsOffsetRegion
          CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                         iErrHypre)
#endif                                         
        END DO ! ifg

! ------------------------------------------------------------------------------
!       Compute viscous flux at patch boundaries for i-momentum equation.
! ------------------------------------------------------------------------------

        DO iPatch = 1,pRegion%grid%nPatches
          pPatch => pRegion%patches(iPatch)
   
          distrib = pPatch%mixt%distrib
 
        SELECT CASE( pPatch%bcType )
        CASE (BC_SLIPWALL)
        CASE (BC_NOSLIPWALL_HFLUX,BC_NOSLIPWALL_TEMP)
          DO ifg = 1,pPatch%nBFaces
            c1 = pPatch%bf2c(ifg)

            xc = pPatch%fc(XCOORD,ifg)
            yc = pPatch%fc(YCOORD,ifg)
            zc = pPatch%fc(ZCOORD,ifg)

            ni = pPatch%fn(ICORD,ifg)
            nj = pPatch%fn(JCORD,ifg)
            nk = pPatch%fn(KCORD,ifg)
            nm = pPatch%fn(XYZMAG,ifg)

            deln = 2.0_RFREAL &
                   *( (xc-pGrid%cofg(XCOORD,c1))*(xc-pGrid%cofg(XCOORD,c1)) &
                    + (yc-pGrid%cofg(YCOORD,c1))*(yc-pGrid%cofg(YCOORD,c1)) &
                    + (zc-pGrid%cofg(ZCOORD,c1))*(zc-pGrid%cofg(ZCOORD,c1)) &
                    )**(0.5_RFREAL)

            mu = pTv(TV_MIXT_MUEL,c1)

            rhoAlpha_np1 = (pCv(CV_MIXT_DENS,c1)+pCvOld(CV_MIXT_DENS,c1)) &
                           /2.0_RFREAL
            rhoBeta_np1  = rhoAlpha_np1 

            duidxi = 0.5_RFREAL*(pGradCell(ICORD,GRC_MIXT_IVEL,c1) + &
                                 pGradCellOld(ICORD,GRC_MIXT_IVEL,c1))

            dujdxi = 0.5_RFREAL*(pGradCell(ICORD,GRC_MIXT_JVEL,c1) + &
                                 pGradCellOld(ICORD,GRC_MIXT_JVEL,c1))

            dukdxi = 0.5_RFREAL*(pGradCell(ICORD,GRC_MIXT_KVEL,c1) + &
                                 pGradCellOld(ICORD,GRC_MIXT_KVEL,c1))

            dujdxj = 0.5_RFREAL*(pGradCell(JCORD,GRC_MIXT_JVEL,c1) + &
                                 pGradCellOld(JCORD,GRC_MIXT_JVEL,c1))

            dukdxk = 0.5_RFREAL*(pGradCell(KCORD,GRC_MIXT_KVEL,c1) + &
                                 pGradCellOld(KCORD,GRC_MIXT_KVEL,c1))

            tiN = mu*(-pCv(CV_MIXT_IVEL,c1)-pCv(CV_MIXT_IVEL,c1) &
                     -pCvOld(CV_MIXT_IVEL,c1)-pCvOld(CV_MIXT_IVEL,c1)) &
                    /(2.0_RFREAL*deln)
            tiD = ONE_THIRD*mu*(duidxi + dujdxj + dukdxk)*ni
            tiT = mu*((dujdxi*nj-dujdxj*ni) + (dukdxi*nk-dukdxk*ni))

#ifdef HYPRE               
            value = nm*iRefREnum*mu/(2.0_RFREAL*deln*rhoAlpha_np1)
            iRowGlobal = c1 + pRegion%nCellsOffset(iReg1)
            iColGlobal = iRowGlobal
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)

            value = nm*iRefREnum*mu/(2.0_RFREAL*deln*rhoBeta_np1)
            iRowGlobal = c1 + nCellsOffsetRegion
            iColGlobal = c1 + nCellsOffsetRegion
            CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                           iColGlobal,value,iErrHypre)

            value = iRefREnum*(tiD+tiT+mu*(-pCvOld(CV_MIXT_IVEL,c1) &
                                           -pCvOld(CV_MIXT_IVEL,c1)) &
                                          /(2.0_RFREAL*deln))*nm
            iRowGlobal = c1 + nCellsOffsetRegion
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                           iErrHypre)
#endif                                           
          END DO ! ifg
        CASE (BC_OUTFLOW)
        CASE (BC_FARFIELD)
        CASE (BC_INFLOW_VELTEMP)
        CASE (BC_INJECTION)
        CASE (BC_PERIODIC,BC_SYMMETRY,BC_VIRTUAL)
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pPatch%bcType
        END DO ! iPatch

! ------------------------------------------------------------------------------
!     Compute coefficients due to viscous terms in axi-symmetry. Momentum eqn.
! ------------------------------------------------------------------------------

        IF (pRegion%mixtInput%axiFlag .EQV. .TRUE.) THEN
          DO ifg = 1,pGrid%nFaces-pGrid%nFacesAV
            c1 = pGrid%f2c(1,ifg)
            c2 = pGrid%f2c(2,ifg)

            ni = pGrid%fn(ICORD,ifg)
            nj = pGrid%fn(JCORD,ifg)
            nk = pGrid%fn(KCORD,ifg)
            nm = pGrid%fn(XYZMAG,ifg)
 
            mu = 0.5_RFREAL*(pTv(TV_MIXT_MUEL,c1)+pTv(TV_MIXT_MUEL,c2)) 

            uysk = 0.25_RFREAL*(pCv(CV_MIXT_YVEL,c1)+ &
                                pCv(CV_MIXT_YVEL,c2)+ &
                                pCvOld(CV_MIXT_YVEL,c1)+ &
                                pCvOld(CV_MIXT_YVEL,c2))/ &
                               (pRegion%grid%fc(YCOORD,ifg)+EPS)

#ifdef HYPRE             
            value = -TWO_THIRD*mu*iRefREnum*uysk*ni*nm
            iRowGlobal = c1 + nCellsOffsetRegion
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal, &
                                           value,iErrHypre)

            value = TWO_THIRD*mu*iRefREnum*uysk*ni*nm
            iRowGlobal = c2 + nCellsOffsetRegion
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal, &
                                           value,iErrHypre)
#endif                                         
          END DO ! ifg

          DO ifg = pGrid%nFaces-pGrid%nFacesAV+1,pGrid%nFaces
            c1 = pGrid%f2c(1,ifg)
            c2 = pGrid%f2c(2,ifg)

            direction = 1.0_RFREAL

            IF ( c1 .GT. pGrid%nCells ) THEN ! c1 is virtual
              dummy = c1
              c1    = c2
              c2    = dummy

              direction = -1.0_RFREAL
            END IF

            ni = direction*pGrid%fn(ICORD,ifg)
            nj = direction*pGrid%fn(JCORD,ifg)
            nk = direction*pGrid%fn(KCORD,ifg)
            nm = pGrid%fn(XYZMAG,ifg)
 
            mu = 0.5_RFREAL*(pTv(TV_MIXT_MUEL,c1)+pTv(TV_MIXT_MUEL,c2)) 

            uysk = 0.25_RFREAL*(pCv(CV_MIXT_YVEL,c1)+ &
                                pCv(CV_MIXT_YVEL,c2)+ &
                                pCvOld(CV_MIXT_YVEL,c1)+ &
                                pCvOld(CV_MIXT_YVEL,c2))/ &
                               (pRegion%grid%fc(YCOORD,ifg)+EPS)

#ifdef HYPRE             
            value = -TWO_THIRD*mu*iRefREnum*uysk*ni*nm
            iRowGlobal = c1 + nCellsOffsetRegion
            CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal, &
                                           value,iErrHypre)
#endif                                         
          END DO ! ifg

          DO iPatch = 1,pRegion%grid%nPatches
            pPatch => pRegion%patches(iPatch)
   
            distrib = pPatch%mixt%distrib
 
            SELECT CASE( pPatch%bcType )
            CASE (BC_SLIPWALL)
            CASE (BC_NOSLIPWALL_HFLUX,BC_NOSLIPWALL_TEMP)
            CASE (BC_OUTFLOW)
            CASE (BC_FARFIELD)
              DO ifg = 1,pPatch%nBFaces
                c1 = pPatch%bf2c(ifg)

                xc = pPatch%fc(XCOORD,ifg)
                yc = pPatch%fc(YCOORD,ifg)
                zc = pPatch%fc(ZCOORD,ifg)

                ni = pPatch%fn(ICORD,ifg)
                nj = pPatch%fn(JCORD,ifg)
                nk = pPatch%fn(KCORD,ifg)
                nm = pPatch%fn(XYZMAG,ifg)

                mu = pTv(TV_MIXT_MUEL,c1)

                uysk = 0.5_RFREAL*(pCv(CV_MIXT_YVEL,c1)+ &
                                   pCvOld(CV_MIXT_YVEL,c1))/ &
                                  (pPatch%fc(YCOORD,ifg)+EPS)

#ifdef HYPRE             
                value = -TWO_THIRD*mu*iRefREnum*uysk*ni*nm
                iRowGlobal = c1 + nCellsOffsetRegion
                CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal, &
                                               value,iErrHypre)
#endif                                           
              END DO ! ifg
            CASE (BC_INFLOW_VELTEMP)
            CASE (BC_INJECTION)
            CASE (BC_PERIODIC,BC_SYMMETRY,BC_VIRTUAL)
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END SELECT ! pPatch%bcType
          END DO ! iPatch

          IF (ICORD == YCOORD) THEN
            DO icg = 1,pGrid%nCells
              volume = pGrid%vol(icg)

              uysk = 0.5_RFREAL*(pCv(CV_MIXT_YVEL,icg)+ &
                                 pCvOld(CV_MIXT_YVEL,icg))/ &
                              (pRegion%grid%cofg(YCOORD,icg)+EPS)**2.0_RFREAL

#ifdef HYPRE           
              value = -2.0_RFREAL*TWO_THIRD*mu*iRefREnum*uysk*volume
              iRowGlobal = icg + nCellsOffsetRegion
              CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal, &
                                             value,iErrHypre)
#endif                                       
            END DO ! icg

            DO ifg = 1,pGrid%nFaces-pGrid%nFacesAV
              c1 = pGrid%f2c(1,ifg)
              c2 = pGrid%f2c(2,ifg)

              ni = pGrid%fn(ICORD,ifg)
              nj = pGrid%fn(JCORD,ifg)
              nk = pGrid%fn(KCORD,ifg)
              nm = pGrid%fn(XYZMAG,ifg)
 
              mu = 0.5_RFREAL*(pTv(TV_MIXT_MUEL,c1)+pTv(TV_MIXT_MUEL,c2)) 

              uysk = 0.5_RFREAL*(pVfMixt(ifg)+pVfMixtOld(ifg))/ &
                                (pRegion%grid%fc(YCOORD,ifg)+EPS)

#ifdef HYPRE             
              value = TWO_THIRD*mu*iRefREnum*uysk*nm
              iRowGlobal = c1 + nCellsOffsetRegion
              CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal, &
                                             value,iErrHypre)

              value = -TWO_THIRD*mu*iRefREnum*uysk*nm
              iRowGlobal = c2 + nCellsOffsetRegion
              CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal, &
                                             value,iErrHypre)
#endif                                         
            END DO ! ifg

            DO ifg = pGrid%nFaces-pGrid%nFacesAV+1,pGrid%nFaces
              c1 = pGrid%f2c(1,ifg)
              c2 = pGrid%f2c(2,ifg)

              direction = 1.0_RFREAL

              IF ( c1 .GT. pGrid%nCells ) THEN ! c1 is virtual
                dummy = c1
                c1    = c2
                c2    = dummy

                direction = -1.0_RFREAL
              END IF

              ni = direction*pGrid%fn(ICORD,ifg)
              nj = direction*pGrid%fn(JCORD,ifg)
              nk = direction*pGrid%fn(KCORD,ifg)
              nm = pGrid%fn(XYZMAG,ifg)
 
              mu = 0.5_RFREAL*(pTv(TV_MIXT_MUEL,c1)+pTv(TV_MIXT_MUEL,c2)) 

              uysk = 0.5_RFREAL*direction*(pVfMixt(ifg)+ &
                                           pVfMixtOld(ifg))/ &
                                          (pRegion%grid%fc(YCOORD,ifg)+EPS)

#ifdef HYPRE             
              value = TWO_THIRD*mu*iRefREnum*uysk*nm
              iRowGlobal = c1 + nCellsOffsetRegion
              CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal, &
                                             value,iErrHypre)
#endif                                         
            END DO ! ifg

            DO iPatch = 1,pRegion%grid%nPatches
              pPatch => pRegion%patches(iPatch)
   
              distrib = pPatch%mixt%distrib
 
              SELECT CASE( pPatch%bcType )
              CASE (BC_SLIPWALL)
              CASE (BC_NOSLIPWALL_HFLUX,BC_NOSLIPWALL_TEMP)
              CASE (BC_OUTFLOW)
              CASE (BC_FARFIELD)
                DO ifg = 1,pPatch%nBFaces
                  c1 = pPatch%bf2c(ifg)

                  xc = pPatch%fc(XCOORD,ifg)
                  yc = pPatch%fc(YCOORD,ifg)
                  zc = pPatch%fc(ZCOORD,ifg)

                  ni = pPatch%fn(ICORD,ifg)
                  nj = pPatch%fn(JCORD,ifg)
                  nk = pPatch%fn(KCORD,ifg)
                  nm = pPatch%fn(XYZMAG,ifg)

                  mu = pTv(TV_MIXT_MUEL,c1)

                  uysk = 0.5_RFREAL*(pVfMixt(ifg)+pVfMixtOld(ifg))/ &
                                    (pPatch%fc(YCOORD,ifg)+EPS)

#ifdef HYPRE             
                  value = TWO_THIRD*mu*iRefREnum*uysk*nm
                  iRowGlobal = c1 + nCellsOffsetRegion
                  CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal, &
                                                 value,iErrHypre)
#endif                                           
                END DO ! ifg
              CASE (BC_INFLOW_VELTEMP)
              CASE (BC_INJECTION)
              CASE (BC_PERIODIC,BC_SYMMETRY,BC_VIRTUAL)
              CASE DEFAULT
                CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
              END SELECT ! pPatch%bcType
            END DO ! iPatch
          END IF ! ICORD
        END IF ! pRegion%mixtInput%axiFlag

      END IF ! flowModel
    END DO ! iReg

! ==============================================================================
!   Assemble matrix and vector
! ==============================================================================

    CALL RFLU_HYPRE_AssembleMatrixVector(regions)

! ==============================================================================
!   Solve system of equation for i-momentum using HYPRE 
! ==============================================================================

    CALL RFLU_HYPRE_SolveMatrix(regions)

! ==============================================================================
!   Extract i-vel solution from HYPRE solution vector parSolHypre
! ==============================================================================

    global => regions(1)%global

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid

      iReg1 = pRegion%iRegionGlobal
      nCellsOffsetRegion = pRegion%nCellsOffset(iReg1) 

      pCv     => pRegion%mixt%cv
      pCvOld  => pRegion%mixt%cvOld

! --- Extract actual cell data -------------------------------------------------
      DO icg = 1,pGrid%nCells
        rho_n = (pCv(CV_MIXT_DENS,icg)+pCvOld(CV_MIXT_DENS,icg))/2.0_RFREAL

        iRowGlobal = icg + nCellsOffsetRegion 
#ifdef HYPRE           
        CALL HYPRE_IJVectorGetValues(regions(1)%SolHypre,1,iRowGlobal,value, &
                                     iErrHypre)
#endif                                     
        pCv(CV_MIXT_IVEL,icg) = value/rho_n
      END DO ! icg
    END DO ! iReg

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_SolveMomentumEq












! ******************************************************************************
!
! Purpose: Wrapper routine to call momentum eq solver.
!
! Description: None.
!
! Input:
!   regions            Pointer to region data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_HM_SolveMomentumEqWrapper(regions)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

! *****************************************************************************
!   Declarations and definitions
! *****************************************************************************

! =============================================================================
!   Arguments
! =============================================================================

    TYPE(t_region), DIMENSION(:), POINTER :: regions
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: i 
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => regions(1)%global
 
    CALL RegisterFunction(global,'RFLU_HM_SolveMomentumEqWrapper',__FILE__)

    CALL RFLU_HM_SolveMomentumEq(regions,XCOORD)
    CALL RFLU_HM_SolveMomentumEq(regions,YCOORD)

    IF ( regions(1)%mixtInput%axiFlag .EQV. .FALSE. ) THEN
      CALL RFLU_HM_SolveMomentumEq(regions,ZCOORD)
    END IF ! pRegion%mixtInput%axiFlag

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_SolveMomentumEqWrapper











! ******************************************************************************
!
! Purpose: Test the Hypre solver.
!
! Description: None.
!
! Input:
!   regions            Pointer to region data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_HM_TestSolver(regions)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

! *****************************************************************************
!   Declarations and definitions
! *****************************************************************************

! =============================================================================
!   Arguments
! =============================================================================

    TYPE(t_region), DIMENSION(:), POINTER :: regions
    
! ==============================================================================
!   Locals
! ==============================================================================
   
    CHARACTER(CHRLEN) :: iFileName 
    INTEGER :: c1,c2,dummy,errorFlag,icg,iColGlobal,iErrHypre,ifg,iRow, &
               iRowGlobal,iReg,iReg1,nSizeHypre,nCellsOffsetRegion,procOffset
    REAL(RFREAL) :: value
    REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: soln
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
!   Start, assign pointers
! ******************************************************************************

    global => regions(1)%global
 
    CALL RegisterFunction(global,'RFLU_HM_TestSolver',__FILE__)

    nSizeHypre = 0
    procOffset = 0

    IF ( global%nRegions > 1 ) THEN
      procOffset = regions(1)%nCellsOffset(regions(1)%iRegionGlobal)
    END IF ! global%nRegions

! *****************************************************************************
!   Initialize Hypre matrix and vectors
! *****************************************************************************

    CALL RFLU_HYPRE_InitializeMatrixVector(regions)

! *****************************************************************************
!   Loop over regions and assign pointers
! *****************************************************************************

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid

      iReg1 = pRegion%iRegionGlobal
      nCellsOffsetRegion = pRegion%nCellsOffset(iReg1)
      nSizeHypre = nSizeHypre + pGrid%nCells

! ==============================================================================
!     Loop over cells and initialize the coefficients
! ==============================================================================

      DO icg = 1,pGrid%nCells

        value = 1.0_RFREAL
        iRowGlobal = icg + nCellsOffsetRegion
        iColGlobal = iRowGlobal
#ifdef HYPRE           
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
        value = iColGlobal
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif                                       
      END DO ! icg

! ==============================================================================
!     Loop over interior faces and accumulate coefficients
! ==============================================================================

      DO ifg = 1,pGrid%nFaces-pGrid%nFacesAV
        c1 = pGrid%f2c(1,ifg)
        c2 = pGrid%f2c(2,ifg)

        value = 1.0_RFREAL
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = iRowGlobal
#ifdef HYPRE           
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
        value = iColGlobal
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)

        value = 1.0_RFREAL
        iRowGlobal = c2 + nCellsOffsetRegion
        iColGlobal = iRowGlobal
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
        value = iColGlobal
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)

        value = 1.0_RFREAL
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = c2 + nCellsOffsetRegion
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
        value = iColGlobal
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)

        value = 1.0_RFREAL
        iRowGlobal = c2 + nCellsOffsetRegion
        iColGlobal = c1 + nCellsOffsetRegion
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
        value = iColGlobal
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif                                       
      END DO ! ifg

      DO ifg = pGrid%nFaces-pGrid%nFacesAV+1,pGrid%nFaces
        c1 = pGrid%f2c(1,ifg)
        c2 = pGrid%f2c(2,ifg)

        IF ( c1 .GT. pGrid%nCells ) THEN ! c1 is virtual
          dummy = c1
          c1    = c2
          c2    = dummy
        END IF ! c1

        value = 1.0_RFREAL
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = iRowGlobal
#ifdef HYPRE           
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
        value = iColGlobal
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)

        value = 1.0_RFREAL
        iRowGlobal = c1 + nCellsOffsetRegion
        iColGlobal = pGrid%virt2GlobalIds(c2 - pGrid%nCells)
        CALL HYPRE_IJMatrixAddToValues(pRegion%AHypre,1,1,iRowGlobal, &
                                       iColGlobal,value,iErrHypre)
        value = iColGlobal
        CALL HYPRE_IJVectorAddToValues(pRegion%RhsHypre,1,iRowGlobal,value, &
                                       iErrHypre)
#endif                                       
      END DO ! ifg
    END DO ! iReg

! ==============================================================================
!   Assemble matrix and vector
! ==============================================================================

    CALL RFLU_HYPRE_AssembleMatrixVector(regions)

! ==============================================================================
!   Solve system of equation for density using HYPRE 
! ==============================================================================

    CALL RFLU_HYPRE_SolveMatrix(regions)

! ==============================================================================
!   Extract solution from HYPRE solution vector parSolHypre
! ==============================================================================

    ALLOCATE(soln(nSizeHypre),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'soln')
    END IF ! global%error

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid

      iReg1 = pRegion%iRegionGlobal
      nCellsOffsetRegion = pRegion%nCellsOffset(iReg1) 

      DO icg = 1,pGrid%nCells
        iRowGlobal = icg + nCellsOffsetRegion 
#ifdef HYPRE           
        CALL HYPRE_IJVectorGetValues(regions(1)%SolHypre,1,iRowGlobal,value, &
                                     iErrHypre)
#endif                                     
        soln(iRowGlobal-procOffset) = value
      END DO ! icg
    END DO ! iReg

    IF ( global%myProcid == 0 ) THEN
      print *,"-------------------------------------------------------"
      print *,"nProcs=",global%nProcAlloc," nRegions=",global%nRegions
      print *,"-------------------------------------------------------"
    END IF ! global%myProcid

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid
      
      iReg1 = pRegion%iRegionGlobal
      nCellsOffsetRegion = pRegion%nCellsOffset(iReg1) 

      DO icg = 1,pGrid%nCells
        WRITE(*,'(A,I2,3X,A,I2,3X,A,I5,3X,A,I5,3X,A,E12.5)') &
                   "proc=",global%myProcid,"region=",iReg1, &
                   "local_cellno=",icg, &
                   "global_row=",icg+pRegion%nCellsOffset(iReg1), &
                   "soln=",soln(icg+pRegion%nCellsOffset(iReg1)-procOffset)
      END DO ! icg
    END DO ! iReg

    DEALLOCATE(soln,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'soln')
    END IF ! global%error

! ==============================================================================
!   Stop the run here.
! ==============================================================================

    CALL ErrorStop(global,ERR_HYPRE_STOP,__LINE__,'RFLU_HM_TestSolver')
 
! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_TestSolver










! ******************************************************************************
!
! Purpose: Update pressure by pressure-correction. 
!
! Description: None.
!
! Input:
!   pRegion            Pointer to region data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_HM_UpdatePressure(pRegion)

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
    
    INTEGER :: errorFlag,icg 
    TYPE(t_global), POINTER :: global  
 
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'RFLU_HM_UpdatePressure',__FILE__)

! ==============================================================================
!   Loop over cells and update pressure 
! ==============================================================================

    DO icg = 1,pRegion%grid%nCellsTot
      pRegion%mixt%cv(CV_MIXT_PRES,icg) = pRegion%mixt%cv(CV_MIXT_PRES,icg) &
                                        + pRegion%mixt%delP(icg)
    END DO ! icg

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_HM_UpdatePressure






END MODULE RFLU_ModHouMahesh

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModHouMahesh.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.9  2009/09/28 14:21:45  mparmar
! Few bug fixes, reordering predcorr cycle, added press, friction coeff etc.
!
! Revision 1.8  2009/08/13 01:25:38  mparmar
! Added RFLU_HM_SetPressureCoeff, bug fix in axisymm, removed patch%mixt
!
! Revision 1.7  2009/07/08 20:53:53  mparmar
! Removed RFLU_ModHouMaheshBoundCond and fixed minor bugs
!
! Revision 1.6  2009/07/08 19:11:56  mparmar
! Implemented boundary conditions, absorbing layer, moving reference frame, axisymmetric computations
!
! Revision 1.5  2008/12/06 08:43:41  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:16:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/12/04 20:02:00  haselbac
! Enclosed HYPRE calls within ifdef HYPRE sections
!
! Revision 1.2  2007/11/29 02:03:09  mparmar
! Removed ^M at the end of lines
!
! Revision 1.1  2007/11/28 23:04:47  mparmar
! Initial revision
!
!
! ******************************************************************************

