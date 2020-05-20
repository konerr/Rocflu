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
! *****************************************************************************
!
! Purpose: Check parameters specified by the user:
!          These checks require knowledge of input to all MP modules
!
! Description: None.
!
! Input:
!   regions        Region data
!
! Output: None.
!
! Notes:
!
! *****************************************************************************
!
! $Id: RFLU_CheckDerivedUserInput.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! *****************************************************************************

SUBROUTINE RFLU_CheckDerivedUserInput(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModMixture, ONLY: t_mixt_input
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters
  USE ModMPI

  IMPLICIT NONE

! *****************************************************************************
! Arguments
! *****************************************************************************

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! *****************************************************************************
! Locals
! *****************************************************************************

  INTEGER :: iReg
  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = & 
    '$RCSfile: RFLU_CheckDerivedUserInput.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_CheckDerivedUserInput',__FILE__)

! *****************************************************************************
! Check region independent data
! *****************************************************************************

! =============================================================================
! Program burn algorithm
! =============================================================================

  IF ( global%pbaFlag .EQV. .TRUE. ) THEN 
#ifdef SPEC
    IF ( global%specUsed .EQV. .FALSE. ) THEN 
      CALL ErrorStop(global,ERR_PBA_STOP,__LINE__, &
                     'Species module must be active.')
    END IF ! global%specUsed
#else    
    CALL ErrorStop(global,ERR_PBA_STOP,__LINE__, & 
                   'Can only be used with species module.')
#endif      
  END IF ! global%pbaFlag

! *****************************************************************************
! Check region related data
! *****************************************************************************

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
    pMixtInput => regions(iReg)%mixtInput

#ifdef PLAG
! =============================================================================
!   Check for consistency between RK scheme and PLAG module
! =============================================================================

    IF ( (global%plagUsed .EQV. .TRUE.) .AND. &
         ( (global%rkScheme /= RK_SCHEME_3_WRAY) .AND. &
           (global%rkScheme /= RK_SCHEME_1_EULER) ) ) THEN
      CALL ErrorStop(global,ERR_RK_SCHEME_INVALID,__LINE__, &
                     'PLAG requires RK3.')
    END IF ! global

! =============================================================================
!   Check for consistency between PLAG module and flag for converting 
!   Lagrangian to Eulerian field
! =============================================================================
 
    IF ( global%plagUsed .EQV. .FALSE. ) THEN
      IF ( global%postLag2EulFlag .EQV. .TRUE. ) THEN 
        IF ( iReg == LBOUND(regions,1) ) THEN
          global%warnCounter = global%warnCounter + 1

          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'*** WARNING *** Invalid '// &
                                   'input for Eulerian postprocessing flag.' 
          WRITE(STDOUT,'(A,20X,A)') SOLVER_NAME,'Setting flag to false.'
        END IF ! iReg

        global%postLag2EulFlag = .FALSE.
      END IF ! global%postLag2EulFlag
    END IF ! global%plagUsed
#endif

! =============================================================================
!   Check for valid input for viscosity model
! =============================================================================

    IF (  pMixtInput%computeTv .AND. &
        ( pMixtInput%viscModel <  VISC_SUTHR .OR.  &
          ( ( pMixtInput%viscModel >  VISC_ANTIB ) .AND. &
            ( pMixtInput%viscModel /=  VISC_WILKE_SUTH ) ) ) ) THEN
      CALL ErrorStop(global,ERR_UNKNOWN_VISCMODEL,__LINE__)
    END IF ! pMixtInput
  END DO ! iReg

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CheckDerivedUserInput

! *****************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CheckDerivedUserInput.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2009/07/09 20:43:58  mparmar
! Added check for VISC_WILKE_SUTH
!
! Revision 1.3  2008/12/06 08:43:34  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:52  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:48  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2005/12/10 16:54:03  haselbac
! Added check for postLag2EulFlag
!
! Revision 1.3  2004/12/04 02:33:12  haselbac
! Bug fix: Missing brackets in logical expression
!
! Revision 1.2  2004/11/30 20:07:48  fnajjar
! Added error trap for running PLAG with none RK3 scheme
!
! Revision 1.1  2004/07/30 22:47:35  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! *****************************************************************************

