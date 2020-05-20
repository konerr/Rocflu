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
! Purpose: read in user input related to time stepping.
!
! Description: none.
!
! Input: user input file.
!
! Output: global = flow type and global time-stepping parameters.
!
! Notes: 
!   1. RFLU: Must not overwrite global%currentTime if running within GENX 
!      because solver gets actual time from Roccom. In preprocessor, however,
!      do not get time from Roccom, so need to read timeStamp into separate
!      variable. This will be used in RFLU_GetUserInput to set the variable
!      global%currentTime iff RFLU_GetUserInput is called from within the
!      preprocessing module. 
!
!******************************************************************************
!
! $Id: ReadTimestepSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadTimestepSection( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
!  INTEGER, PARAMETER :: NVALS_MAX = 26
! BBR - new value for NVALS_MAX
  INTEGER, PARAMETER :: NVALS_MAX = 28

  CHARACTER(10) :: keys(NVALS_MAX)

  LOGICAL :: defined(NVALS_MAX)

  REAL(RFREAL) :: vals(NVALS_MAX)

!******************************************************************************

  CALL RegisterFunction( global,'ReadTimestepSection',__FILE__ )

! specify keywords and search for them

  keys( 1) = 'FLOWTYPE'
  keys( 2) = 'TIMESTEP'
  keys( 3) = 'MAXTIME'
  keys( 4) = 'WRITIME'
  keys( 5) = 'PRNTIME'
  keys( 6) = 'MAXITER'
  keys( 7) = 'RESTOL'
  keys( 8) = 'WRIITER'
  keys( 9) = 'PRNITER'
  keys(10) = 'STARTTIME'
  keys(11) = 'STARTITER'
  keys(12) = 'SOLVERTYPE'
  keys(13) = 'ORDER'
  keys(14) = 'MAXSUBITER'
  keys(15) = 'TOLSUBITER'
  keys(16) = 'PREDICTSOL'
  keys(17) = 'DTMINLIMIT'
  keys(18) = 'RKSCHEME'
  keys(19) = 'RUNTIME'
  keys(20) = 'FOMWRITIME'
  keys(21) = 'TOLPC'
  keys(22) = 'TOLHYPRE'
  keys(23) = 'MAXITPC'
  keys(24) = 'MAXITHYPRE'
  keys(25) = 'FOMWRIITER'
  keys(26) = 'HYPRESOLVER'
  keys(27) = 'WALLTIME'
  keys(28) = 'SAFEWRTIME'
  
  CALL ReadSection( global,IF_INPUT,NVALS_MAX,keys,vals,defined )

  IF (defined(1)) THEN
                              global%flowType = FLOW_STEADY
    IF (vals(1) > 0.9_RFREAL) global%flowType = FLOW_UNSTEADY
  ENDIF
  IF (defined( 2)) global%dtImposed     = ABS(vals(2))
  IF (defined( 3)) global%maxTime       = ABS(vals(3))
  IF (defined( 4)) global%writeTime     = ABS(vals(4))
  IF (defined( 5)) global%printTime     = ABS(vals(5))
  IF (defined( 6)) global%maxIter       = INT(ABS(vals(6))+0.5_RFREAL)
  IF (defined( 7)) global%resTol        = ABS(vals(7))
  IF (defined( 8)) global%writeIter     = MAX(1,INT(ABS(vals(8))+0.5_RFREAL))
  IF (defined( 9)) global%printIter     = MAX(1,INT(ABS(vals(9))+0.5_RFREAL))

#ifndef GENX
  IF (defined(10)) global%timeStamp     = ABS(vals(10))
#else
  IF (defined(10)) global%timeStampPrep = ABS(vals(10))
#endif

#ifndef GENX
  IF ( defined(10) .EQV. .TRUE. ) THEN
    global%currentTime = global%timeStamp
  ELSE 
    global%currentTime = 0.0_RFREAL
  END IF ! defined
#else
  IF ( defined(10) .EQV. .FALSE. ) THEN 
    global%timeStampPrep = 0.0_RFREAL
  END IF ! defined
#endif

  IF ( defined(11) .EQV. .TRUE. ) THEN 
    global%currentIter = INT(ABS(vals(11))+0.5_RFREAL)
  END IF ! defined(11)

  IF ( defined(12) .EQV. .TRUE. ) THEN 
    IF ( NINT(vals(12)) == SOLV_EXPLICIT ) THEN
      global%solverType = SOLV_EXPLICIT
    ELSE IF ( NINT(vals(12)) == SOLV_IMPLICIT_NK ) THEN
      global%solverType = SOLV_IMPLICIT_NK
    ELSE
      global%solverType = SOLV_IMPLICIT_HM
    END IF ! NINT
  ENDIF

  IF ( defined(13) .EQV. .TRUE. ) THEN 
    global%tstepOrder = MAX(2,INT(ABS(vals(13))+0.5_RFREAL))
  END IF ! defined 

  IF ( defined(14) .EQV. .TRUE. ) THEN 
    global%maxSubIter = MAX(1,INT(ABS(vals(14))+0.5_RFREAL))
  END IF ! defined 

  IF ( defined(15) .EQV. .TRUE. ) THEN 
    global%tolSubIter = ABS(vals(15))
  END IF ! defined 

  IF ( defined(16) .EQV. .TRUE. ) THEN
    IF ( vals(16) < 0.9_RFREAL ) THEN
      global%predictSol = .false.
    ELSE
      global%predictSol = .true.
    END IF ! vals
  END IF ! defined

  IF ( defined(17) .EQV. .TRUE. ) THEN 
    global%dtMinLimit = ABS(vals(17))
  END IF ! defined 

  IF ( defined(18) ) THEN 
    global%rkScheme = INT(ABS(vals(18)) + 0.5_RFREAL)
  END IF ! defined

  IF ( defined(19) .EQV. .TRUE. ) THEN 
    global%runTime = ABS(vals(19))
  END IF ! defined 

  IF (defined(20)) THEN
    global%fomWriteTime = ABS(vals(20)) 
  ELSE
    global%fomWriteTime = global%writeTime
  END IF ! defined(20)

  IF ( defined(21) .EQV. .TRUE. ) THEN 
    global%tolPC = ABS(vals(21))
  END IF ! defined 

  IF ( defined(22) .EQV. .TRUE. ) THEN 
    global%tolHypre = ABS(vals(22))
  END IF ! defined 

  IF ( defined(23) .EQV. .TRUE. ) THEN 
    global%maxIterPC = INT(vals(23))
  END IF ! defined 

  IF ( defined(24) .EQV. .TRUE. ) THEN 
    global%maxIterHypre = INT(vals(24))
  END IF ! defined 

  IF ( defined(25) .EQV. .TRUE. ) THEN 
    global%fomWriteIter = MAX(1,INT(ABS(vals(25))+0.5_RFREAL))
  ELSE
    global%fomWriteIter = global%writeIter 
  END IF ! defined 

  IF ( defined(26) .EQV. .TRUE. ) THEN 
    IF ( NINT(vals(26)) == HYPRE_SOLV_BAMG ) THEN
      global%hypreSolver = HYPRE_SOLV_BAMG
    ELSE IF ( NINT(vals(26)) == HYPRE_SOLV_GMRES ) THEN
      global%hypreSolver = HYPRE_SOLV_GMRES
    ELSE
      global%hypreSolver = HYPRE_SOLV_BAMG
    END IF ! NINT
  ENDIF

  IF ( defined(27) .EQV. .TRUE. ) THEN
    global%wallTime = ABS(vals(27))
  END IF

  IF ( defined(28) .EQV. .TRUE. ) THEN
    global%safeWriteTime = ABS(vals(28))
  END IF

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadTimestepSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadTimestepSection.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.8  2010/03/14 23:20:59  mparmar
! Added HYPRESOLVER
!
! Revision 1.7  2009/07/08 19:11:25  mparmar
! Added FOMWRIITER
!
! Revision 1.6  2008/12/06 08:43:32  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:16:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2007/11/28 23:17:45  mparmar
! Added reading of parameters for SOLV_IMPLICIT_HM
!
! Revision 1.3  2007/06/18 17:40:16  mparmar
! Added reading of FOMWRITIME
!
! Revision 1.2  2007/06/14 01:47:38  haselbac
! Added new keyword, clean-up
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2006/05/09 23:36:17  wasistho
! added DTFIXED for implicit
!
! Revision 1.2  2005/08/03 18:14:07  hdewey2
! Added reading of solverType
!
! Revision 1.1  2004/12/01 16:50:55  haselbac
! Initial revision after changing case
!
! Revision 1.19  2004/11/17 16:23:21  haselbac
! Added RKSCHEME as input parameter
!
! Revision 1.18  2004/08/10 00:22:24  wasistho
! added RFREAL to real number 0.9 in IF statements
!
! Revision 1.17  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.14  2003/10/15 02:38:58  haselbac
! Added new key and parameter NVALS_MAX
!
! Revision 1.13  2003/07/03 21:48:44  jblazek
! Implemented dual-time stepping.
!
! Revision 1.12  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.11  2003/03/25 19:21:00  haselbac
! Removed old DEBUG variable
!
! Revision 1.10  2003/01/30 19:06:40  haselbac
! Added timeStampPrep variable, see note
!
! Revision 1.9  2002/10/16 21:09:59  haselbac
! Fixed bug in RFLU code segment
!
! Revision 1.8  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.7  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.6  2002/05/28 13:46:55  haselbac
! Set currentTime to timeStamp for RFLU
!
! Revision 1.5  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.4  2002/02/01 00:00:24  jblazek
! Edge and corner cells defined for each level.
!
! Revision 1.3  2002/01/11 17:18:31  jblazek
! Updated description of I/O variables.
!
! Revision 1.2  2001/12/22 00:09:38  jblazek
! Added routines to store grid and solution.
!
! Revision 1.1  2001/12/07 16:54:32  jblazek
! Added files to read user input.
!
!******************************************************************************

