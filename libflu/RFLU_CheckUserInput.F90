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
! Purpose: Check parameters specified by the user.
!
! Description: None.
!
! Input:
!   regions        Region data
!
! Output: None.
!
! Notes:
!   1. Only check input for one region because input is the same for all
!      regions.
!
! ******************************************************************************
!
! $Id: RFLU_CheckUserInput.F90,v 1.3 2016/02/05 19:02:43 fred Exp $
!
! Copyright: (c) 2001-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_CheckUserInput(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModMixture, ONLY: t_mixt_input
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters
  USE ModMPI

  IMPLICIT NONE

! ******************************************************************************
! Arguments
! ******************************************************************************

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ******************************************************************************
! Locals
! ******************************************************************************

  INTEGER :: iReg
  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_CheckUserInput.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_CheckUserInput',__FILE__)

! ******************************************************************************
! Check global data
! ******************************************************************************

! ==============================================================================
! Solver type
! ==============================================================================

  IF ( global%solverType /= SOLV_EXPLICIT .AND. & 
       global%solverType /= SOLV_IMPLICIT_NK .AND. & 
       global%solverType /= SOLV_IMPLICIT_HM ) THEN 
    CALL ErrorStop(global,ERR_SOLVER_TYPE_INVALID,__LINE__)
  END IF ! global%solverType

! ==============================================================================
!   Check Hypre libraries for non-dissipative code 
! ==============================================================================

  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
#ifndef HYPRE
    CALL ErrorStop(global,ERR_HYPRE_LIBRARIES,__LINE__)
#endif 
  END IF ! global%solverType

! ==============================================================================
!   Check Hypre solvers for non-dissipative code 
! ==============================================================================

  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
    IF ( global%hypreSolver /= HYPRE_SOLV_BAMG .AND. & 
         global%hypreSolver /= HYPRE_SOLV_GMRES ) THEN 
      CALL ErrorStop(global,ERR_HYPRE_SOLVER_INVALID,__LINE__)
    END IF ! global%solverType
  END IF ! global%solverType

! ==============================================================================
! Check global data for Absorbing boundary Condition ABC
! ==============================================================================

  IF ( global%abcFlag .EQV. .TRUE. ) THEN
    IF (   global%abcDens  == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
      .OR. global%abcVelX  == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
      .OR. global%abcVelY  == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
      .OR. global%abcVelZ  == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
      .OR. global%abcPress == REAL(CRAZY_VALUE_INT,KIND=RFREAL) ) THEN
      CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__)
    END IF ! global%abcLeft

    IF (global%abcKind == 0 ) THEN

      IF ( global%abcSponge == REAL(CRAZY_VALUE_INT,KIND=RFREAL) ) THEN
        CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__)
      END IF ! global%abcSponge

      SELECT CASE( global%abcOrder )
        CASE (0)
          IF ( global%abcCoeff0 == REAL(CRAZY_VALUE_INT,KIND=RFREAL) ) THEN
            CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__)
          END IF ! global%abcCoeff0
        CASE (1)
          IF (   global%abcCoeff0 == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
            .OR. global%abcCoeff1 == REAL(CRAZY_VALUE_INT,KIND=RFREAL) ) THEN
            CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__)
          END IF ! global%abcCoeff0
        CASE (2)
          IF (   global%abcCoeff0 == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
            .OR. global%abcCoeff1 == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
            .OR. global%abcCoeff2 == REAL(CRAZY_VALUE_INT,KIND=RFREAL) ) THEN
            CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__)
          END IF ! global%abcCoeff0
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%abcOrder

      SELECT CASE( global%abcDomain )
        CASE (0)
          IF (   global%abcLeft  == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
            .OR. global%abcRight == REAL(CRAZY_VALUE_INT,KIND=RFREAL) ) THEN
            CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__)
          END IF ! global%abcLeft
        CASE (1)
          IF (   global%abcLeft   == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
            .OR. global%abcRight  == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
            .OR. global%abcBottom == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
            .OR. global%abcTop    == REAL(CRAZY_VALUE_INT,KIND=RFREAL) ) THEN
            CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__)
          END IF ! global%abcLeft
        CASE (2)
          IF (   global%abcRadius  == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
            .OR. global%abcCenterX == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
            .OR. global%abcCenterY == REAL(CRAZY_VALUE_INT,KIND=RFREAL) ) THEN
            CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__)
          END IF ! global%abcRadius
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%abcDomain

    END IF ! global%abcKind
  END IF ! global%abcFlag

! ==============================================================================
! Check global data for Ghost Fluid Method GFM
! ==============================================================================

  IF ( global%gfmFlag .EQV. .TRUE. ) THEN
    IF (   (global%gfmNRigids < 0) &
      .OR. (global%gfmNFluids < 0) &
      .OR. (global%gfmGridX1 == REAL(CRAZY_VALUE_INT,KIND=RFREAL)) &
      .OR. (global%gfmGridX2 == REAL(CRAZY_VALUE_INT,KIND=RFREAL)) &
      .OR. (global%gfmGridY1 == REAL(CRAZY_VALUE_INT,KIND=RFREAL)) &
      .OR. (global%gfmGridY2 == REAL(CRAZY_VALUE_INT,KIND=RFREAL)) &
      .OR. (global%gfmGridZ1 == REAL(CRAZY_VALUE_INT,KIND=RFREAL)) &
      .OR. (global%gfmGridZ2 == REAL(CRAZY_VALUE_INT,KIND=RFREAL)) &
      .OR. (global%gfmGridNx == CRAZY_VALUE_INT) &
      .OR. (global%gfmGridNy == CRAZY_VALUE_INT) &
      .OR. (global%gfmGridNz == CRAZY_VALUE_INT) ) THEN
      CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__)
    END IF ! global%gfmGrid

    IF (    global%gfmGridX1 >= global%gfmGridX2 &
      .OR.  global%gfmGridY1 >= global%gfmGridY2 &
      .OR.  global%gfmGridZ1 >= global%gfmGridZ2 ) THEN
      CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__,"Domain bounds invalid")
    END IF ! global%gfmGridX1
  END IF ! global%gfmFlag

! ==============================================================================
! Check global data for moving frame section MVFRAME
! ==============================================================================

  IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
    IF ( global%mvfMass <= 0 ) THEN
      CALL ErrorStop(global,ERR_MVF_MASS_INVALID,__LINE__)
    END IF ! global%mvfMass
 
    IF ( global%iPatchGlobalMvFrame == CRAZY_VALUE_INT &
         .OR. global%iPatchGlobalMvFrame < 1 ) THEN
      CALL ErrorStop(global,ERR_MVF_PATCH_INVALID,__LINE__)
    END IF ! global%iPatchGlobalMvFrame

    IF ( global%mvfAccType /= ACC_CONSTANT .AND. & 
         global%mvfAccType /= ACC_SINUSOIDAL ) THEN
      CALL ErrorStop(global,ERR_MVF_ACCTYPE_INVALID,__LINE__)
    END IF ! global%mvfAccType

    IF ( global%mvfAccFlag .EQV. .TRUE. ) THEN
      IF (    global%mvfAccTs == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
         .OR. global%mvfAccTe == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
         .OR. global%mvfAccX  == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
         .OR. global%mvfAccY  == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
         .OR. global%mvfAccZ  == REAL(CRAZY_VALUE_INT,KIND=RFREAL) ) THEN
        CALL ErrorStop(global,ERR_MVF_ACC_INVALID,__LINE__,"acc invalid")
      END IF ! global%mvfAccTs

      IF ( global%mvfAccType == ACC_SINUSOIDAL ) THEN
        IF (    global%mvfOmega == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
           .OR. global%mvfPhase == REAL(CRAZY_VALUE_INT,KIND=RFREAL) ) THEN
          CALL ErrorStop(global,ERR_MVF_ACC_INVALID,__LINE__,"acc invalid")
        END IF ! global%mvfOmega
      END IF ! global%mvfAccType
    END IF ! global%mvfAccFlag

    IF (    global%mvfVelInitX  == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
       .OR. global%mvfVelInitY  == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
       .OR. global%mvfVelInitZ  == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
       .OR. global%mvfLocInitX  == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
       .OR. global%mvfLocInitY  == REAL(CRAZY_VALUE_INT,KIND=RFREAL) &
       .OR. global%mvfLocInitZ  == REAL(CRAZY_VALUE_INT,KIND=RFREAL) ) THEN
      CALL ErrorStop(global,ERR_MVF_LOCVEL_INVALID,__LINE__,"acc invalid")
    END IF ! global%mvfVelX

    IF ( regions(1)%mixtInput%spaceOrder < 2 ) THEN
      CALL ErrorStop(global,ERR_MVF_ORDER_INVALID,__LINE__)
    END IF !  pMixtInput%spaceOrder
  END IF ! global%mvFrameFlag

! ==============================================================================
! Output format
! ==============================================================================
  
  IF ( global%postOutputFormat /= POST_OUTPUT_FORMAT_TECPLOT .AND. & 
       global%postOutputFormat /= POST_OUTPUT_FORMAT_ENSIGHT ) THEN 
    CALL ErrorStop(global,ERR_POST_OUTPUT_FORMAT_INVALID,__LINE__)
  END IF ! global%postOutputFormat

! ==============================================================================
! For ENSIGHT, check number of servers
! ==============================================================================

  IF ( global%postOutputFormat == POST_OUTPUT_FORMAT_ENSIGHT ) THEN 
    IF ( (global%postMergeFlag .EQV. .FALSE.) .AND. & 
         (global%postNServers > global%nRegions) ) THEN 
      CALL ErrorStop(global,ERR_POST_NSERVERS_INVALID,__LINE__)
    END IF ! global%postNServers

    IF ( (global%postMergeFlag .EQV. .TRUE.) .AND. & 
         (global%postNServers /= 1) ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
            '*** WARNING *** Invalid input for number of servers.'
      WRITE(STDOUT,'(A,20X,A)') SOLVER_NAME,'Setting number of servers to one.'             

      global%postNServers = 1
    END IF ! global%postNServers
  END IF ! global%postOutputFormat
  
#ifdef PLAG
! ==============================================================================
! Particle plotting fraction
! ==============================================================================

  IF ( global%plagUsed .EQV. .TRUE. ) THEN 
    IF ( (global%postPlagFrac <= 0.0_RFREAL) .OR. & 
         (global%postPlagFrac >  1.0_RFREAL) ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
            '*** WARNING *** Invalid input for particle plotting fraction.'
      WRITE(STDOUT,'(A,20X,A)') SOLVER_NAME, & 
            'Setting particle plotting fraction to one.'             

      global%postPlagFrac = 1.0_RFREAL
    END IF ! global%postPlagFrac
  END IF ! global%plagUsed
#endif 

! ******************************************************************************
! Check region related data
! ******************************************************************************

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
    pMixtInput => regions(iReg)%mixtInput
   
! ==============================================================================
!   Check nDv for non-dissipative code
! ==============================================================================

    IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
      IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN
        IF ( pMixtInput%nDv /= 1 ) THEN 
          CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
        END IF ! pMixtInput%nDv
      END IF ! pRegion%mixtInput%flowModel 
    END IF ! global%solverType
    
! ==============================================================================
!   Check dimensionality
! ==============================================================================

    IF ( pMixtInput%dimens < 1 .OR. pMixtInput%dimens > 3 ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! pMixtInput%dimens
    
! ==============================================================================
!   Check stencil dimensionality
! ==============================================================================

    IF ( global%solverType /= SOLV_IMPLICIT_HM ) THEN
      IF ( pMixtInput%stencilDimensCells < 1 .OR. &
           pMixtInput%stencilDimensCells > 3 ) THEN 
        CALL ErrorStop(global,ERR_STENCILDIMENS_INVALID,__LINE__)
      END IF ! pMixtInput%stencilDimensCells 
    
      IF ( pMixtInput%stencilDimensFaces < 1 .OR. &
           pMixtInput%stencilDimensFaces > 3 ) THEN 
        CALL ErrorStop(global,ERR_STENCILDIMENS_INVALID,__LINE__)
      END IF ! pMixtInput%stencilDimensFaces  
    
      IF ( pMixtInput%stencilDimensBFaces < 1 .OR. &
           pMixtInput%stencilDimensBFaces > 3 ) THEN 
        CALL ErrorStop(global,ERR_STENCILDIMENS_INVALID,__LINE__)
      END IF ! pMixtInput%stencilDimensBFaces          
    
      IF ( pMixtInput%stencilDimensCells > pMixtInput%dimens ) THEN 
        CALL ErrorStop(global,ERR_STENCILDIMENS_INVALID,__LINE__)
      END IF ! pMixtInput%stencilDimensCells   
    
      IF ( pMixtInput%stencilDimensFaces > pMixtInput%dimens ) THEN 
        CALL ErrorStop(global,ERR_STENCILDIMENS_INVALID,__LINE__)
      END IF ! pMixtInput%stencilDimensFaces    
    
      IF ( pMixtInput%stencilDimensBFaces > pMixtInput%dimens ) THEN 
        CALL ErrorStop(global,ERR_STENCILDIMENS_INVALID,__LINE__)
      END IF ! pMixtInput%stencilDimensBFaces                  
    END IF ! global%solverType

! ==============================================================================
!   Check axisymmetry
! ==============================================================================

    IF ( (pMixtInput%axiFlag .EQV. .TRUE.) .AND. (pMixtInput%dimens /= 2) ) THEN
      CALL ErrorStop(global,ERR_AXIFLAG_INCONSISTENT,__LINE__)
    END IF ! pMixtInput%axiFlag

! ==============================================================================
!   Check for valid value of viscosity for Navier-Stokes flow model. 
! ==============================================================================

    IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN
      IF ( pMixtInput%refVisc == REAL(CRAZY_VALUE_INT,RFREAL) ) THEN
        CALL ErrorStop(global,ERR_REFVISC_INVALID,__LINE__)
      END IF ! pMixtInput%refVisc
    END IF ! pMixtInput%flowModel

! ==============================================================================
!   Check order
! ==============================================================================

    IF ( global%solverType /= SOLV_IMPLICIT_HM ) THEN
      IF ( pMixtInput%stencilDimensCells > 1 ) THEN 
        IF ( pMixtInput%spaceOrder < 1 .OR. & 
             pMixtInput%spaceOrder > 2 ) THEN 
          CALL ErrorStop(global,ERR_ORDER_INVALID,__LINE__)
        END IF ! pMixtInput%spaceOrder
      END IF ! pMixtInput%stencilDimensCells
   
      IF ( pMixtInput%spaceOrderBFaces < 1 .OR. & 
           pMixtInput%spaceOrderBFaces > 2 ) THEN 
        CALL ErrorStop(global,ERR_ORDER_INVALID,__LINE__)
      END IF ! pMixtInput%spaceOrderBFaces
    END IF ! global%solverType
   
! ==============================================================================
!   Check reconstruction
! ==============================================================================

    IF ( global%solverType /= SOLV_IMPLICIT_HM ) THEN
      IF ( pMixtInput%reconst /= RECONST_NONE          .AND. & 
           pMixtInput%reconst /= RECONST_WENO_SIMPLE   .AND. & 
           pMixtInput%reconst /= RECONST_WENO_XYZ      .AND. & 
           pMixtInput%reconst /= RECONST_LIM_BARTHJESP .AND. & 
           pMixtInput%reconst /= RECONST_LIM_VENKAT    ) THEN 
        CALL ErrorStop(global,ERR_RECONST_INVALID,__LINE__)
      END IF ! pMixtInput%reconst
    END IF ! global%solverType
    
! ==============================================================================
!   Check constraints
! ==============================================================================

    IF ( global%solverType /= SOLV_IMPLICIT_HM ) THEN
      IF ( pMixtInput%cReconstCells /= CONSTR_NONE .AND. & 
           pMixtInput%cReconstCells /= CONSTR_WEIGHTED ) THEN 
        CALL ErrorStop(global,ERR_CONSTR_INVALID,__LINE__)
      END IF ! pMixtInput%cReconstCells

      IF ( pMixtInput%cReconstFaces /= CONSTR_NONE .AND. & 
           pMixtInput%cReconstFaces /= CONSTR_WEIGHTED ) THEN 
        CALL ErrorStop(global,ERR_CONSTR_INVALID,__LINE__)
      END IF ! pMixtInput%cReconstFaces
    END IF ! global%solverType

! ==============================================================================
!   Check discretization
! ==============================================================================

    IF ( global%solverType /= SOLV_IMPLICIT_HM ) THEN
      IF ( pMixtInput%spaceDiscr /= DISCR_UPW_ROE  .AND. & 
           pMixtInput%spaceDiscr /= DISCR_UPW_HLLC .AND. & 
           pMixtInput%spaceDiscr /= DISCR_UPW_AUSMPLUS .AND. & 
           pMixtInput%spaceDiscr /= DISCR_UPW_AUSMPLUSUP ) THEN 
        CALL ErrorStop(global,ERR_DISCR_INVALID,__LINE__)
      END IF ! pMixtInput%spaceDiscr
    END IF ! global%solverType

! ==============================================================================
!   Check gas model and compatibility non-dissipative solver and flux function
! ==============================================================================

    IF ( global%solverType /= SOLV_IMPLICIT_HM ) THEN
      IF ( pMixtInput%gasModel /= GAS_MODEL_TCPERF      .AND. & 
           pMixtInput%gasModel /= GAS_MODEL_MIXT_TCPERF .AND. & 
           pMixtInput%gasModel /= GAS_MODEL_MIXT_PSEUDO .AND. &
           pMixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ .AND. &
           pMixtInput%gasModel /= GAS_MODEL_MIXT_JWL ) THEN
        CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__)
      END IF ! pMixtInput%gasModel
    
      IF ( pMixtInput%gasModel == GAS_MODEL_MIXT_TCPERF .OR. & 
           pMixtInput%gasModel == GAS_MODEL_MIXT_PSEUDO ) THEN 
        IF ( pMixtInput%spaceDiscr /= DISCR_UPW_HLLC .AND. &
             pMixtInput%spaceDiscr /= DISCR_UPW_AUSMPLUS .AND. & 
             pMixtInput%spaceDiscr /= DISCR_UPW_AUSMPLUSUP ) THEN 
          CALL ErrorStop(global,ERR_GASMODEL_DISCR_MISMATCH,__LINE__)
        END IF ! pMixtInput%spaceDiscr
      END IF ! pMixtInput%gasModel

      IF ( pMixtInput%gasModel == GAS_MODEL_MIXT_JWL ) THEN 
        IF ( pMixtInput%spaceDiscr /= DISCR_UPW_AUSMPLUS .AND. &
             pMixtInput%spaceDiscr /= DISCR_UPW_AUSMPLUSUP  ) THEN 
          CALL ErrorStop(global,ERR_GASMODEL_DISCR_MISMATCH,__LINE__)
        END IF ! pMixtInput%spaceDiscr
      END IF ! pMixtInput%gasModel

      IF ( pMixtInput%gasModel == GAS_MODEL_MIXT_GASLIQ ) THEN 
        IF ( pMixtInput%spaceDiscr /= DISCR_UPW_ROE .AND. &
             pMixtInput%spaceDiscr /= DISCR_UPW_HLLC ) THEN 
          CALL ErrorStop(global,ERR_GASMODEL_DISCR_MISMATCH,__LINE__)
        END IF ! pMixtInput%spaceDiscr
      END IF ! pMixtInput%gasModel
    ELSE
      IF ( pMixtInput%gasModel /= GAS_MODEL_TCPERF ) THEN 
        CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__)
      END IF ! pMixtInput%gasModel
    END IF ! global%solverType

! ==============================================================================
!   Check in-cell test tolerance
! ==============================================================================

    IF ( pMixtInput%tolerICT < 0.0_RFREAL ) THEN 
      CALL ErrorStop(global,ERR_TOLERICT_INVALID,__LINE__)
    END IF ! pMixtInput%tolerICT

! ==============================================================================
!   Check for valid input for grid motion type
! ==============================================================================

    IF ( pMixtInput%moveGrid .EQV. .TRUE. ) THEN
      IF ( pMixtInput%moveGridType /= MOVEGRID_TYPE_DISP .AND. &
           pMixtInput%moveGridType /= MOVEGRID_TYPE_XYZ .AND. & 
           pMixtInput%moveGridType /= MOVEGRID_TYPE_GENX ) THEN
        global%warnCounter = global%warnCounter + 1

        IF ( iReg == LBOUND(regions,1) .AND.  &
             global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_NONE ) THEN
          WRITE(STDOUT,'(A,3X,A,A)') SOLVER_NAME, &
                                     '*** WARNING *** Invalid input for ', &
                                     'grid motion type. Overriding user input.'
        END IF ! iReg

#ifndef GENX
        pMixtInput%moveGridType = MOVEGRID_TYPE_DISP
#else
        pMixtInput%moveGridType = MOVEGRID_TYPE_GENX
#endif
      END IF ! pMixtInput%moveGridType

#ifdef GENX
      IF ( pMixtInput%moveGridType /= MOVEGRID_TYPE_GENX ) THEN
        global%warnCounter = global%warnCounter + 1

        IF ( iReg == LBOUND(regions,1) .AND.  &
             global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_NONE ) THEN
          WRITE(STDOUT,'(A,3X,A,A)') SOLVER_NAME, &
                                     '*** WARNING *** Invalid input for ', &
                                     'grid motion type. Overriding user input.'
        END IF ! iReg

        pMixtInput%moveGridType = MOVEGRID_TYPE_GENX
      END IF ! pMixtInput%moveGridType
#endif
    END IF ! pMixtInput%moveGrid

#ifdef GENX
! ==============================================================================
!   Check for valid input for grid motion within coupled calculations
! ==============================================================================

    IF ( pMixtInput%moveGrid .EQV. .FALSE. ) THEN
      global%warnCounter = global%warnCounter + 1

      IF ( iReg == LBOUND(regions,1) .AND. &
           global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,A)') SOLVER_NAME, &
                                   '*** WARNING *** Invalid input for ', &
                                   'grid motion. Overriding user input.'
      END IF ! global%myProcid

      pMixtInput%moveGrid = .TRUE.
    END IF ! pMixtInput
#endif

! ==============================================================================
!   Check for valid input for OLES computations
! ==============================================================================

    IF ( pMixtInput%spaceDiscr == DISCR_OPT_LES .AND. &
         pMixtInput%flowModel  /= FLOW_NAVST ) THEN
      CALL ErrorStop(global,ERR_OLES_FLOWMODEL,__LINE__)
    END IF ! pMixtInput
  END DO ! iReg

! ==============================================================================
!   Check for valid input for random number generator
! ==============================================================================

  IF ( global%randSeedType < RAND_SEED_TYPE_FIXED .OR. & 
       global%randSeedType > RAND_SEED_TYPE_CLOCK      ) THEN 
    CALL ErrorStop(global,ERR_RAND_SEED_TYPE_INVALID,__LINE__)
  END IF ! global%randSeedType

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CheckUserInput

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CheckUserInput.F90,v $
! Revision 1.3  2016/02/05 19:02:43  fred
! Adding iterative JWL EOS solver for cylindrical detonation problem
!
! Revision 1.2  2015/12/18 22:53:37  rahul
! Added AUSM+up in the list of flux schemes.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.10  2010/03/14 23:40:33  mparmar
! Added checks for Hypre solver and encapsulated dissipative solver specific checks
!
! Revision 1.9  2009/07/08 19:11:28  mparmar
! Added check for absorbing layer related variables and refVisc
!
! Revision 1.8  2008/12/06 08:43:34  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:16:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2008/05/29 01:35:16  mparmar
! Added checks for non-dissipative solver & moving reference frame
!
! Revision 1.5  2008/03/27 12:08:47  haselbac
! Added check for axiFlag
!
! Revision 1.4  2007/11/28 23:04:56  mparmar
! Added consistency check for SOLV_IMPLICIT_HM
!
! Revision 1.3  2007/06/18 17:44:17  mparmar
! Added checking of moving reference frame data
!
! Revision 1.2  2007/04/12 17:53:25  haselbac
! Added checking of postPlagFrac
!
! Revision 1.1  2007/04/09 18:48:52  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:48  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.30  2007/02/27 13:01:33  haselbac
! Enabled 1d computations
!
! Revision 1.29  2006/10/20 21:23:07  mparmar
! Added upper limit check of spaceOrderBFaces
!
! Revision 1.28  2006/08/19 15:38:43  mparmar
! Added checking of mixtInput%spaceOrderBFaces
!
! Revision 1.27  2006/05/01 20:59:30  haselbac
! Changed check for spaceDiscr with MTCP gas model
!
! Revision 1.26  2006/04/15 16:54:19  haselbac
! Expanded check for reconst
!
! Revision 1.25  2006/04/13 18:06:15  haselbac
! Fixed check for gas model, added check for GASLIQ gas model and discr
!
! Revision 1.24  2006/04/07 14:41:21  haselbac
! Added checks for stencilDimens params, shuffled checks around
!
! Revision 1.23  2006/03/26 20:21:28  haselbac
! Extended check of gas model
!
! Revision 1.22  2006/01/06 22:05:46  haselbac
! Added check for stencil dimensionality and order
!
! Revision 1.21  2005/12/25 15:22:16  haselbac
! Added checks for constrained reconstruction
!
! Revision 1.20  2005/12/24 21:25:11  haselbac
! Added check for ICT tolerance
!
! Revision 1.19  2005/12/01 17:09:53  fnajjar
! Added appropriate initialization, checking and printing of random seed type
!
! Revision 1.18  2005/11/14 16:54:04  haselbac
! Added error checking for pseudo-gas model
!
! Revision 1.17  2005/11/10 01:59:24  haselbac
! Added checks for gas-model/discr consistency
!
! Revision 1.16  2005/10/31 19:25:35  haselbac
! Added checking of gasModel
!
! Revision 1.15  2005/10/27 18:55:11  haselbac
! Added check for constraints, clean-up
!
! Revision 1.14  2005/10/05 20:01:46  haselbac
! Added checks for pot output format and nservers
!
! Revision 1.13  2005/08/03 18:53:48  hdewey2
! Added check for solverType
!
! Revision 1.12  2005/07/14 21:58:30  haselbac
! Added checking for flux function and order of discretization
!
! Revision 1.11  2005/07/11 19:23:02  mparmar
! Added check of reconst option
!
! Revision 1.10  2005/06/29 18:07:28  haselbac
! Added check for gm type within GENX
!
! Revision 1.9  2005/03/09 14:52:46  haselbac
! Added check for dimensionality
!
! Revision 1.8  2004/10/19 19:36:58  haselbac
! Updated for GEN3
!
! Revision 1.7  2004/07/30 22:47:35  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.6  2004/03/05 22:09:01  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.5  2004/01/09 21:02:41  haselbac
! Cosmetics only
!
! Revision 1.4  2003/07/22 01:55:21  haselbac
! Added global%warnCounter
!
! Revision 1.3  2003/04/10 23:29:20  fnajjar
! Checking consistency of viscosity models
!
! Revision 1.2  2003/03/31 16:10:55  haselbac
! Cosmetics, added grid-motion type check
!
! Revision 1.1  2003/01/28 15:53:31  haselbac
! Initial revision, moved from rocflu
!
! Revision 1.7  2003/01/08 21:08:31  haselbac
! Fixed problems with double slashes in ifdefs (Absoft 8.0)
!
! Revision 1.6  2002/12/20 23:19:22  haselbac
! Fixed output bug: increased indentation
!
! Revision 1.5  2002/10/27 19:11:27  haselbac
! Added check for GENX calculations
!
! Revision 1.4  2002/09/09 15:37:29  haselbac
! Cleaned up routine and moved global under regions
!
! Revision 1.3  2002/09/02 23:44:42  wasistho
! Removed TURB compilation check
!
! Revision 1.2  2002/08/24 03:18:35  wasistho
! modify TURB error msg
!
! Revision 1.1  2002/08/18 02:34:09  wasistho
! Added to check TURB module activation and other input data
!
! ******************************************************************************

