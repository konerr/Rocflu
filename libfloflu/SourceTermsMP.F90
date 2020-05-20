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
! Purpose: add source terms to the residual.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: pRegion%levels%mixt%rhs = complete right-hand side (residual).
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: SourceTermsMP.F90,v 1.2 2016/05/06 00:49:37 rahul Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SourceTermsMP( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : SourceTerms, &
                            SourceTermsMovingFrame
  USE ModError
  USE ModParameters

  USE RFLU_ModABC, ONLY: RFLU_ABC_SourceTerms
  USE RFLU_ModAxisymmetry, ONLY: RFLU_AXI_SourceTerms

  USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, &
                               RFLU_ConvertCvPrim2Cons

#ifdef INRT
  USE ModInterfacesInteract, ONLY : INRT_SourceTerms
#endif
#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY: PLAG_RFLU_CorrectMixtProperties, &
                                     PLAG_RFLU_ComputeMaxImpulse, &
                                     PLAG_RFLU_GetMixtPG, &
                                     PLAG_RFLU_InterParticleForce
#endif
#ifdef SPEC
  USE SPEC_RFLU_ModChemistry, ONLY: SPEC_RFLU_IntegrateChemSrcTerm
  USE ModInterfacesSpecies, ONLY : SPEC_RFLU_SourceTerms_GL 
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), TARGET :: region

! ... local variables
  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! ******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'SourceTermsMP',__FILE__ )

! add source terms
  IF ( global%accelOn .EQV. .TRUE. ) THEN  
    CALL SourceTerms( region )
  END IF ! global%accelOn

  IF ( global%abcFlag .EQV. .TRUE. ) THEN
    IF ( global%abcKind == 0 ) THEN
      CALL RFLU_ABC_SourceTerms( region )
    END IF ! global%abcKind
  END IF ! global%abcFlag

  IF ( region%mixtInput%axiFlag .EQV. .TRUE. ) THEN
    CALL RFLU_AXI_SourceTerms( region )
  END IF ! region%mixtInput%axiFlag

  IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
    CALL SourceTermsMovingFrame( region )
  END IF ! global%mvFrameFlag

  pRegion => region
  CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)

#ifdef PLAG
    IF ( global%plagUsed .EQV. .TRUE. ) THEN
    CALL PLAG_RFLU_CorrectMixtProperties(pRegion)

    IF ( pRegion%plagInput%flagCollision .EQV. .TRUE. ) THEN
      CALL PLAG_RFLU_InterParticleForce(pRegion)
    END IF ! pRegion%plagInput%flagCollision
    CALL PLAG_RFLU_GetMixtPG(pRegion)

!   Compute max Impulse transfer possible
!   Which can be used by force computation routines to limit F
!   and can be used by timeStepImpulse computation routine to limit dt
    CALL PLAG_RFLU_ComputeMaxImpulse( pRegion )

! **********************************************************************
! Rahul-
! 1. Add force due to hyperbolicity preserving interface pressure to
!    partilce momentum equation.
! 2. Compute non-cnservative terms for gas phase

    IF ( (region%mixtInput%spaceDiscr .EQ. DISCR_UPW_AUSMPLUSUP) &
        .AND. (region%plag%nPcls .GE. 1) ) THEN
!      CALL SourceTermsPint( region )
      CALL PLAG_RFLU_ComputeSourcePint( region )
      CALL PLAG_RFLU_CalcForceInterface( region )
    END IF ! region%mixtInput%spcaeDiscr

! Rahul-End
! **********************************************************************

  END IF ! global%plagUsed
#endif  

#ifdef INRT
  IF (global%inrtUsed) THEN
! TEMPORARY: Manoj: 2012-05-29: Only allowing inter-particle collision      
!WRITE(*,*) "===== Commented out  INRT_SourceTerms"
    CALL INRT_SourceTerms( region )
! END TEMPORARY
  ENDIF
#endif

  pRegion => region
  CALL RFLU_ConvertCvPrim2Cons(pRegion,CV_MIXT_STATE_CONS)

#ifdef SPEC
! ==============================================================================
! Species
! ==============================================================================

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    IF ( pRegion%specInput%sourceFlag .EQV. .TRUE. ) THEN   
   
! ------------------------------------------------------------------------------
!     Cavitation source term for gas-liquid-vapor mixture model
! ------------------------------------------------------------------------------   
   
      IF ( pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_GASLIQ ) THEN
        CALL SPEC_RFLU_SourceTerms_GL(pRegion)
      END IF ! pRegion%mixtInput%gasModel

! ------------------------------------------------------------------------------
!     Combustion source term for burning-crack simulations
! ------------------------------------------------------------------------------   

! TEMPORARY - At present, do not use operator-split integration of chemistry
!             source terms based on modifications by Luca. This means that no
!             longer call source term integration function directly in 
!             RFLU_TimeStepping.F90.
!      CALL SPEC_RFLU_IntegrateChemSrcTerm(pRegion,0)
! END TEMPORARY
    END IF ! pRegion%specInput%sourceFlag
  END IF ! global%specUsed
#endif

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE SourceTermsMP

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SourceTermsMP.F90,v $
! Revision 1.2  2016/05/06 00:49:37  rahul
! 1. Added global%accelOn check for SourceTerms.
! 2. Added calls to PLAG_RFLU_CalcForceInterface and SourceTermsPint. These
!    are related to E-L AUSM+up.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.6  2009/07/08 19:11:27  mparmar
! Added call for RFLU_ABC_SourceTerms
!
! Revision 1.5  2008/12/06 08:43:32  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:16:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2008/03/27 12:07:52  haselbac
! Added call to axisymm source term routine
!
! Revision 1.2  2007/06/18 17:42:41  mparmar
! Added calls to compute source terms due to moving reference frame
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.6  2006/03/30 20:47:22  haselbac
! Clean-up of source terms for species
!
! Revision 1.5  2006/03/26 20:21:22  haselbac
! Added call to GL src term routine
!
! Revision 1.4  2005/11/30 22:05:25  fnajjar
! Added call to PLAG_RFLU_CorrectMixtProperties
!
! Revision 1.3  2005/10/05 13:48:57  haselbac
! Bug fix: Enclosed call to chem src term within IF
!
! Revision 1.2  2005/06/06 14:23:03  haselbac
! Adapted to Lucas changes
!
! Revision 1.1  2004/12/01 16:51:28  haselbac
! Initial revision after changing case
!
! Revision 1.15  2004/05/03 15:09:41  jferry
! added equilibrium Eulerian capability for smoke
!
! Revision 1.14  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.13  2004/03/02 21:49:21  jferry
! Added inrtUsed flag to mixture data structure
!
! Revision 1.12  2004/01/31 03:56:19  haselbac
! Added RFLU state vector conversion routines
!
! Revision 1.11  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.8  2003/10/03 20:13:02  wasistho
! initial installation of turbModel SA and DES
!
! Revision 1.7  2003/10/01 23:52:09  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.6  2003/09/26 21:45:34  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.5  2003/08/28 20:33:00  wasistho
! excluced ModInterfacesTurbulence,Radiation,Periodic from ModInterfaces
!
! Revision 1.4  2003/07/17 00:55:20  wasistho
! initial activation rocrad
!
! Revision 1.3  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.2  2003/04/05 01:59:33  wasistho
! install ROCPERI
!
! Revision 1.1  2003/03/28 19:47:43  fnajjar
! Initial import for RocfluidMP
!
! ******************************************************************************

