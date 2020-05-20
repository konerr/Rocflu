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
! Purpose: Read in user input related to time-dependent boundary conditions
!   (done on all processors).
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Reads TBCs for all modules
!
! ******************************************************************************
!
! $Id: RFLU_ReadTbcInputFile.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ReadTbcInputFile(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  USE ModMPI

  USE ModInterfaces, ONLY : RFLU_ReadTbcSection

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN+9) :: fname
  CHARACTER(256)      :: line

  INTEGER :: errorFlag,loopCounter

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction( global,'RFLU_ReadTbcInputFile',__FILE__ )

! ******************************************************************************
! Open file
! ******************************************************************************

  WRITE(fname,'(A)') TRIM(global%inDir)//TRIM(global%casename)//'.bc'
  OPEN(IF_INPUT,file=fname,form='formatted',status='old',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! read file looking for keywords

  loopCounter = 0

  DO
    READ(IF_INPUT,'(A256)',ERR=10) line
    
    SELECT CASE( TRIM(line) )

    CASE ('# TBC_PIECEWISE')
      CALL RFLU_ReadTbcSection( pRegion, TBC_PIECEWISE )

    CASE ('# TBC_SINUSOIDAL')
      CALL RFLU_ReadTbcSection( pRegion, TBC_SINUSOIDAL )

    CASE ('# TBC_STOCHASTIC')
      CALL RFLU_ReadTbcSection( pRegion, TBC_STOCHASTIC )

    CASE ('# TBC_WHITENOISE')
      CALL RFLU_ReadTbcSection( pRegion, TBC_WHITENOISE )

    CASE ('# END')
      EXIT
    END SELECT ! TRIM(line)

    loopCounter = loopCounter + 1 ! Prevent infinite loop
    IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
      CALL ErrorStop(global,ERR_INFINITE_LOOP ,__LINE__)
    END IF ! loopCounter
  ENDDO ! <empty>

! close file

  CLOSE(IF_INPUT,IOSTAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname))
  END IF ! global%error

! finalization & error handling -----------------------------------------------

  GOTO 999

! error handling

10  CONTINUE
  CALL ErrorStop(global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname))

999 CONTINUE

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ROCFLU time-dependent '// & 
                             'boundary condition file done.'
  END IF ! global%verbLevel

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ReadTbcInputFile

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ReadTbcInputFile.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:36  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:53  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2005/04/27 02:08:02  haselbac
! Cosmetics only
!
! Revision 1.2  2003/06/10 22:54:42  jferry
! Added Piecewise TBC
!
! Revision 1.1  2003/06/04 20:05:54  jferry
! re-worked implementation of TBCs in unstructured code
!
! ******************************************************************************

