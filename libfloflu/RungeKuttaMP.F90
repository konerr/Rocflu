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
! $Id: RungeKuttaMP.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RungeKuttaMP( regions )

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

  USE ModInterfaces, ONLY: AfterUpdateMP, &
                           CellGradientsMP, ConvectiveFluxesMP, &
                           GlobalCommunicationMP, InitCommunicationMP, &
                           RKInitMP, RKUpdateMP, &
                           SourceTermsMP, UpdateBoundaryConditionsMP, &
                           UpdateDependentVarsMP, ViscousFluxesMP, &
                           ZeroDummyCellsMP, ZeroResidualsMP
  USE ModInterfaces, ONLY: RFLU_EquilibriumEulerian, &
                           RFLU_SetVarsContWrapper, &
                           RFLU_SetVarsDiscWrapper, & 
                           RFLU_UpdateBoundaryValues

#ifdef GENX
  USE RFLU_ModGENXTools, ONLY: RFLU_GENX_InitBFLAG
#endif
#ifdef SPEC
  USE SPEC_RFLU_ModPBA, ONLY: SPEC_RFLU_PBA_ProgramBurn, &
                              SPEC_RFLU_PBA_ReactionZone
#endif
#ifdef PLAG
  USE PLAG_RFLU_ModComm
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iRegLocal, istage, icg

! ... local variables
  INTEGER :: flowModel

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RungeKuttaMP',__FILE__ )

! loop over stages and regions ================================================

  DO istage=1,global%nrkSteps

! ----- compute particle accelerations and velocities  ------------------------
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

    DO iRegLocal=1,global%nRegionsLocal
      iReg = iRegLocal

! ----- set pointer and get models --------------------------------------------

        pRegion => regions(iReg)

        flowModel  = regions(iReg)%mixtInput%flowModel
        regions(iReg)%irkStep = istage
        regions(iReg)%dummyStep = .FALSE.

! ----- Set RK time -----------------------------------------------------------

        CALL RFLU_SetTimeRK(pRegion,iStage)


! ----- RFLU fill GENX incoming buffers ---------------------------------------

#ifdef GENX
	CALL RFLU_GENX_InitBFLAG(pRegion)
#endif
        CALL RFLU_UpdateBoundaryValues(regions(iReg),istage)

! ----- Scale grid speeds -----------------------------------------------------

        CALL RFLU_SetGridSpeedScaleFactor(pRegion)
        CALL RFLU_ScaleGridSpeeds(pRegion)

! ----- program burn, before RKInitMP -----------------------------------------
#ifdef SPEC        
! DEBUG: Manoj-PBA1D, Notes: To match Jianghui's implementation
!                         1. Move call to ProgramBurn to before update dependent vars
!                    changing it back to old implementation of Manoj
!                         2. Remove call to ReactionZone
! END DEBUG
        IF ( (global%specUsed .EQV. .TRUE.) .AND. &
             (global%pbaFlag .EQV. .TRUE.) ) THEN
! DEBUG: Manoj-PBA1D
          CALL SPEC_RFLU_PBA_ProgramBurn( pRegion )
! END DEBUG
        END IF ! pbaFlag
#endif

! ----- set ghost fluid solution ----------------------------------------------

! TEMPORARY: Manoj: GFM, need to put more thinking
        CALL RFLU_GFM_SetGhostFluid( pRegion )

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

! ----- zero residuals --------------------------------------------------------

        CALL ZeroResidualsMP(regions(iReg))

! ----- add source terms ------------------------------------------------------

        CALL SourceTermsMP( regions(iReg) )

! ----- zero residuals --------------------------------------------------------

        CALL ZeroResidualsMP(regions(iReg))

! ----- add Equilibrium Eulerian corrections ----------------------------------

        CALL RFLU_EquilibriumEulerian( pRegion )

! ----- zero out residuals in dummy cells -------------------------------------

        CALL ZeroDummyCellsMP( regions(iReg) )

! ----- Descale grid speeds -----------------------------------------------------
        CALL RFLU_DescaleGridSpeeds(pRegion)

    ENDDO    ! iReg

    IF(global%zoomFactor > 1) THEN
       CALL RFLU_TimeZoomDriver(regions)
    ENDIF ! global%zoomFactor    

    DO iRegLocal=1,global%nRegionsLocal
      iReg = iRegLocal

! ----- set Region pointer

        pRegion => regions(iReg)

! ----- update solution; sum up residuals -------------------------------------

        CALL RKUpdateMP( regions(iReg),iReg,istage )

        IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_BXV_ComputeVarsCv(pRegion)
        END IF ! RFLU_NSCBC_DecideHaveNSCBC

! ----- program burn, after RKUpdateMP ----------------------------------------
#ifdef SPEC        
        IF ( (global%specUsed .EQV. .TRUE.) .AND. &
             (global%pbaFlag .EQV. .TRUE.) ) THEN
! DEBUG: Manoj-PBA1D
!          CALL SPEC_RFLU_PBA_ReactionZone( pRegion )
! END DEBUG
        END IF ! pbaFlag
#endif
        
! ----- set ghost fluid solution ----------------------------------------------

! TEMPORARY: Manoj: GFM, need to put more thinking
        CALL RFLU_GFM_SetGhostFluid( pRegion )

! ----- perform checks and enforce after-update conditions --------------------

!        CALL AfterUpdateMP( pRegion,istage )

! ----- Descale grid speeds -----------------------------------------------------

        CALL RFLU_DescaleGridSpeeds(pRegion)

! ----- program burn, before SetDependentVars ---------------------------------
#ifdef SPEC        
        IF ( (global%specUsed .EQV. .TRUE.) .AND. &
             (global%pbaFlag .EQV. .TRUE.) ) THEN
! DEBUG: Manoj-PBA1D
!          CALL SPEC_RFLU_PBA_ProgramBurn( pRegion )
! END DEBUG
        END IF ! pbaFlag
#endif

! ----- update dependent variables --------------------------------------------

        CALL RFLU_MPI_ISendWrapper(pRegion)
        CALL RFLU_SetVarsContWrapper(pRegion,1,pRegion%grid%nCells)

! Subbu - Perform check after computing pressure & temperature
! ----- perform checks and enforce after-update conditions --------------------

        CALL AfterUpdateMP( pRegion,istage )
! Subbu - End Perform check

! ----- update dependent variables on boundary faces --------------------------

        IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_BXV_SetDependentVars(pRegion)
        END IF !
    END DO ! iReg

    CALL RFLU_MPI_CopyWrapper(regions)

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)

      CALL RFLU_MPI_RecvWrapper(pRegion)
      CALL RFLU_SetVarsContWrapper(pRegion,pRegion%grid%nCells+1, & 
                                   pRegion%grid%nCellsTot)
      CALL RFLU_RELP_TransformWrapper(pRegion)                                   
    END DO ! iReg

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)

      CALL RFLU_MPI_ClearRequestWrapper(pRegion)
    END DO ! iReg

#ifdef PLAG
    IF ( global%plagUsed .EQV. .TRUE. ) THEN
      CALL PLAG_RFLU_CommDriver(regions)

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_SetVarsDiscWrapper(pRegion)
      END DO ! iReg

! --- update particle volume fraction in virtual cells ------------------------      
! TEMPORARY: Manoj: 2012-05-29: checking effect of not communicating vFracE
!IF (1==2) THEN
      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_MPI_PLAG_ISendWrapper(pRegion)
      END DO ! iReg 

      CALL RFLU_MPI_PLAG_CopyWrapper(regions)
    
      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_MPI_PLAG_RecvWrapper(pRegion)
      END DO ! iReg 

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_MPI_ClearRequestWrapper(pRegion)
      END DO ! iReg 
!END IF ! 1==2
! END TEMPORARY
    END IF ! global%plagUsed
 
#endif
   
  END DO ! istage

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE RungeKuttaMP

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RungeKuttaMP.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.8  2009/08/28 18:29:39  mtcampbe
! RocfluMP integration with Rocstar and some makefile tweaks.  To build
! Rocstar with new Rocflu:
! make ROCFLU=RocfluMP
! To build Rocstar with the new RocfluND:
! make ROCFLU=RocfluMP HYPRE=/the/hypre/install/path
!
! Revision 1.7  2008/12/06 08:43:32  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:16:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2008/05/29 01:35:11  mparmar
! Added Setting of reference frame velocity and update of BC here
!
! Revision 1.4  2007/12/04 13:36:59  haselbac
! Bug fix: Removed USE RFLU_ModPatchVelocity
!
! Revision 1.3  2007/12/03 16:34:03  mparmar
! Removed RFLU_SetPatchVelocity
!
! Revision 1.2  2007/06/18 17:42:13  mparmar
! Added calls for moving reference frame implementation
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.16  2007/03/27 00:39:50  haselbac
! Removed call to PLAG_CalcnPclsTotGlobal, now in RFLU_TimeStepping
!
! Revision 1.15  2007/03/20 22:02:29  fnajjar
! Included call to PLAG_CalcnPclsTotGlobal
!
! Revision 1.14  2006/08/21 16:10:01  haselbac
! Adapted to name change
!
! Revision 1.13  2006/08/19 15:48:25  mparmar
! Added computations of boundary Cv and Dv for NSCBC implementation
!
! Revision 1.12  2006/08/18 21:09:27  fnajjar
! Removed IF around PLAG_RFLU_CommDriver for serial periodic cases
!
! Revision 1.11  2006/03/25 21:40:03  haselbac
! Added call to transforming data on related patches, cosmetics
!
! Revision 1.10  2005/12/03 19:44:54  haselbac
! Apparent bug fix: Separated call to RFLU_MPI_ClearRequestWrapper into separate loop
!
! Revision 1.9  2005/12/01 21:52:14  fnajjar
! Added IF statement around PLAG_RFLU_CommDriver, only active for more than one nRegions
!
! Revision 1.8  2005/11/10 22:21:07  fnajjar
! ACH: Proper fix for updating PLAG dv
!
! Revision 1.7  2005/11/10 16:51:28  fnajjar
! Added plagUsed IF statement around PLAG routines
!
! Revision 1.6  2005/11/02 14:53:24  haselbac
! Fady: Temporary fix so comm particles get non-cv vars updated properly
!
! Revision 1.5  2005/05/18 22:04:41  fnajjar
! Added PLAG communication routines; only initial implementation
!
! Revision 1.4  2005/04/29 00:06:09  haselbac
! Added routines to clear send requests
!
! Revision 1.3  2005/04/15 15:06:06  haselbac
! Converted to MPI
!
! Revision 1.2  2005/03/31 16:31:02  haselbac
! Added call to RFLU_SetTimeRK
!
! Revision 1.1  2004/12/01 16:51:13  haselbac
! Initial revision after changing case
!
! Revision 1.18  2004/11/14 19:36:23  haselbac
! Replaced call to UpdateDependentVarsMP by RFLU_SetVarsWrapper
!
! Revision 1.17  2004/07/30 22:47:33  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.16  2004/04/14 02:07:02  haselbac
! Added grid-speed scaling calls for RFLU
!
! Revision 1.15  2004/03/25 21:14:20  jferry
! changed AfterUpdate to call most subroutines only after final RK stage
!
! Revision 1.14  2004/03/02 21:47:28  jferry
! Added After Update interactions
!
! Revision 1.13  2004/02/26 21:11:58  wasistho
! added globalCommunication
!
! Revision 1.12  2004/02/26 21:01:46  haselbac
! Enclosed updateBoundaryConditionsMP within ifdef RFLO
!
! Revision 1.11  2004/01/29 22:52:47  haselbac
! Added calls to RFLU_EnforceBoundsWrapper and updateDependentVarsMP
!
! Revision 1.10  2003/12/04 03:23:06  haselbac
! Added call to CellGradientsMP and validity check
!
! Revision 1.9  2003/11/25 21:01:45  haselbac
! Added calls to RFLU_UpdateDummyCells and ZeroResidualsMP
!
! Revision 1.8  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/10/03 20:42:07  haselbac
! Added Rocflu calls
!
! Revision 1.4  2003/10/01 23:52:09  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.3  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.2  2003/04/10 01:22:41  jblazek
! Got rid of pRegion in ViscousFluxesMP.
!
! Revision 1.1  2003/03/28 19:42:55  fnajjar
! Initial import for RocfluidMP
!
! ******************************************************************************

