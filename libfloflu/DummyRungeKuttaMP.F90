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
! Purpose: calculate solution at a new time level.
!
! Description: the governing equations are integrated in time using
!              the classical 4-stage Runge-Kutta method (4th-order in
!              time) in low-storage formulation.
!
! Input: regions = data of all regions.
!
! Output: regions%levels%mixt = new solution after one time step.
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: DummyRungeKuttaMP.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE DummyRungeKuttaMP( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters

  USE RFLU_ModBoundXvUtils
  USE RFLU_ModGridSpeedUtils, ONLY: RFLU_DescaleGridSpeeds, &
                                    RFLU_ScaleGridSpeeds, &
                                    RFLU_SetGridSpeedScaleFactor
  USE RFLU_ModMovingFrame, ONLY: RFLU_MVF_ComputeAcceleration, &
                                 RFLU_MVF_SetVelocity, &
                                 RFLU_MVF_UpdateBC
  USE RFLU_ModTimeZoom, ONLY: RFLU_TimeZoomDriver
  USE RFLU_ModMPI
  USE RFLU_ModNSCBC
  USE RFLU_ModGFM
  USE RFLU_ModRelatedPatches, ONLY: RFLU_RELP_TransformWrapper
  USE RFLU_ModTime, ONLY: RFLU_SetTimeRK

  USE ModInterfaces, ONLY : AfterUpdateMP, &
                            CellGradientsMP, ConvectiveFluxesMP, &
                            GlobalCommunicationMP, InitCommunicationMP, &
                            RKInitMP, RKUpdateMP, &
                            SourceTermsMP, UpdateBoundaryConditionsMP, &
                            UpdateDependentVarsMP, ViscousFluxesMP, &
                            ZeroDummyCellsMP, ZeroResidualsMP
  USE ModInterfaces, ONLY : RFLU_EquilibriumEulerian, &
                            RFLU_SetVarsContWrapper, &
                            RFLU_SetVarsDiscWrapper, & 
                            RFLU_UpdateBoundaryValues
  USE ModInterfaces, ONLY : RFLU_MinimumTimeStepImpulse, &
                            RFLU_TimeStepImpulse

#ifdef GENX
  USE RFLU_ModGENXTools, ONLY: RFLU_GENX_InitBFLAG
#endif
#ifdef PLAG
  USE PLAG_RFLU_ModComm
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iRegLocal, istage

! ... local variables
  INTEGER :: flowModel

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'DummyRungeKuttaMP',__FILE__ )

! loop over stages and regions ================================================

!  DO istage=1,global%nrkSteps
  DO istage=1,1
    DO iRegLocal=1,global%nRegionsLocal
      iReg = iRegLocal

! ----- set pointer and get models --------------------------------------------

        pRegion => regions(iReg)

        flowModel  = regions(iReg)%mixtInput%flowModel
        regions(iReg)%irkStep = istage
        regions(iReg)%dummyStep = .TRUE.

! ----- Set RK time -----------------------------------------------------------

        CALL RFLU_SetTimeRK(pRegion,iStage)

! ----- store previous solution; set dissipation to zero ----------------------

        CALL RKInitMP( regions(iReg),istage )

! ----- compute cell-gradients for higher-order scheme ------------------------

        CALL CellGradientsMP( regions(iReg) )

! ----- compute viscous fluxes ------------------------------------------------

        IF ( flowModel == FLOW_NAVST ) THEN
          CALL ViscousFluxesMP( regions(iReg) )
        END IF ! flowModel

! ----- compute convective fluxes; form residual ------------------------------

        CALL ConvectiveFluxesMP( regions(iReg) )

! ----- add source terms ------------------------------------------------------

        CALL SourceTermsMP( regions(iReg) )
        
! ----- Compute dt for each cell due to particle momentum exchange ------------

        CALL RFLU_TimeStepImpulse( pRegion )

    ENDDO    ! iReg
  END DO ! istage

! compute minimum dt due to particle momentum exchange ========================

  CALL RFLU_MinimumTimeStepImpulse(regions)

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE DummyRungeKuttaMP

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: DummyRungeKuttaMP.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
!
! ******************************************************************************
