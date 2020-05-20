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
! Purpose: Read and process user input for the core solver and all physical
!   modules.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes:
!   1. timeStampPrep is written into timeStamp and currentTime only if this
!      routine is called from the preprocessor. This is indicated by the
!      optional argument isPrep. This construct is required because the
!      routine reading the timeStamp from the input file only sets timeStamp
!      and currentTime if not compiled within GENX. This is done because
!      RFLU_InitFlowSolver gets the actual time passed in as an argument when
!      running within GENX, and hence it should not be overwritten. For this
!      to work with the preprocessor, the little detour has to be taken...
!
!******************************************************************************
!
! $Id: RFLU_GetUserInput.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002,2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_GetUserInput(regions,inPrep)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModMPI
  USE ModParameters

  USE ModInterfaces, ONLY: RFLU_CheckDerivedUserInput, &
                           RFLU_PrintUserInput, &
                           RFLU_SetDerivedUserInput, &
                           RFLU_UserInput

#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY: PLAG_UserInput
#endif

#ifdef PERI
  USE ModInterfacesPeriodic, ONLY: PERI_UserInput
#endif

#ifdef SPEC
  USE ModInterfacesSpecies, ONLY: SPEC_UserInput
#endif

#ifdef TURB
  USE ModInterfacesTurbulence, ONLY: TURB_UserInput
#endif

#ifdef INRT
  USE ModInterfacesInteract, ONLY: INRT_UserInput
#endif

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  LOGICAL, OPTIONAL :: inPrep
  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_GetUserInput.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_GetUserInput',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading user input file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Input for mixture & base solver
! ******************************************************************************

  CALL RFLU_UserInput(regions)

#ifdef GENX
! write timeStampPrep into time variables: See note above

  IF ( PRESENT(inPrep) .EQV. .TRUE. ) THEN
    IF ( inPrep .EQV. .TRUE. ) THEN
      global%timeStamp   = global%timeStampPrep
      global%currentTime = global%timeStampPrep
    END IF ! inPrep
  END IF ! PRESENT
#endif

! ******************************************************************************
! Input for physical modules
! ******************************************************************************

#ifdef PLAG
  CALL PLAG_UserInput(regions)
#endif

#ifdef PERI
  CALL PERI_UserInput(regions)
#endif

#ifdef RADI
  CALL RADI_UserInput(regions)
#endif

#ifdef SPEC
  CALL SPEC_UserInput(regions)
#endif

#ifdef TURB
  CALL TURB_UserInput(regions)
#endif

#ifdef INRT
  CALL INRT_UserInput(regions) ! must follow all other *_UserInput routines
#endif

! ******************************************************************************
! Derived values which require knowledge of input to all MP modules
! ******************************************************************************

  CALL RFLU_SetDerivedUserInput(regions)

! ******************************************************************************
! Checks which require knowledge of input to all MP modules
! ******************************************************************************

  CALL RFLU_CheckDerivedUserInput(regions)

! ******************************************************************************
! Write out user input for base solver and physical modules
! ******************************************************************************

  IF ( (global%myProcid == MASTERPROC) .AND. &
       (global%verbLevel > VERBOSE_LOW) ) THEN
    CALL RFLU_PrintUserInput(regions)
  END IF ! global

! ******************************************************************************
! End
! ******************************************************************************

  IF ( (global%myProcid == MASTERPROC) .AND. &
       (global%verbLevel > VERBOSE_NONE) ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading user input files done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_GetUserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_GetUserInput.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:35  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:16:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2008/05/23 21:45:53  mtcampbe
! Set printing of user input to fire if verblevel > verb_low instead of
! high.
!
! Revision 1.1  2007/04/09 18:48:52  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.10  2004/07/30 22:47:35  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.9  2004/07/23 22:43:15  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.8  2004/06/17 23:04:30  wasistho
! added PERI_UserInput
!
! Revision 1.7  2004/02/02 22:49:22  haselbac
! Added interface for PLAG_UserInput
!
! Revision 1.6  2004/01/29 22:56:31  haselbac
! Added call to RFLU_SetDerivedUserInput
!
! Revision 1.5  2003/11/25 21:02:46  haselbac
! Added SPEC support and made cosmetic changes
!
! Revision 1.4  2003/03/04 22:12:34  jferry
! Initial import of Rocinteract
!
! Revision 1.3  2003/02/02 21:15:36  haselbac
! Bug fix: Need to check for presence and value separately...
!
! Revision 1.2  2003/01/30 19:07:34  haselbac
! Added optional argument, rfluprep-GENX problem solved
!
! Revision 1.1  2003/01/28 15:53:32  haselbac
! Initial revision, moved from rocflu
!
! Revision 1.7  2002/09/17 22:51:23  jferry
! Removed Fast Eulerian particle type
!
! Revision 1.6  2002/09/09 15:50:47  haselbac
! global and mixtInput now under region
!
! Revision 1.5  2002/08/24 03:20:56  wasistho
! put safety within #ifdef TURB
!
! Revision 1.4  2002/07/25 14:27:55  haselbac
! Added MASTERPROC distinction for output
!
! Revision 1.3  2002/06/17 13:34:12  haselbac
! Prefixed SOLVER_NAME to all screen output
!
! Revision 1.2  2002/05/04 17:09:28  haselbac
! Cosmetic changes
!
! Revision 1.1  2002/03/26 19:24:07  haselbac
! Initial revision
!
!******************************************************************************

