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
! Purpose: read in user input (done on all processors).
!
! Description: none.
!
! Input: user input file.
!
! Output: regions = reference values, numerics, probe`s position.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_ReadInputFile.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_ReadInputFile( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY      : t_region
  USE ModGlobal, ONLY          : t_global
  USE PLAG_ModInterfaces, ONLY : PLAG_ReadDisPartSection, &
                                 PLAG_ReadDisPartnContSection, &
                                 PLAG_ReadDisPartInitSection

  USE ModError
  USE ModParameters
  USE ModMPI
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  CHARACTER(CHRLEN)   :: RCSIdentString
  CHARACTER(CHRLEN+4) :: fname
  CHARACTER(256)      :: line

  LOGICAL :: usedSomewhere, unusedSomewhere

  INTEGER :: errorFlag, readStatus

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ReadInputFile.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction( global, 'PLAG_ReadInputFile',__FILE__ )

! Open file -------------------------------------------------------------------

  fname = TRIM(global%inDir)//TRIM(global%casename)//'.inp'
  OPEN(IF_INPUT,file=fname,form='formatted',status='old',iostat=errorFlag )
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop( global, ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )
  END IF ! global%error

! Read file looking for keywords ----------------------------------------------

  DO
    READ(IF_INPUT,'(A256)',err=10,end=86) line

    SELECT CASE(TRIM(line))
      CASE ('# DISPART')
        IF ( global%myProcid == MASTERPROC .AND. &
          global%verbLevel > VERBOSE_NONE ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
           'Reading PLAG_ReadDisPartSection...'
        END IF ! global%verbLevel
        CALL PLAG_ReadDisPartSection( regions )

      CASE ('# DISPART_NCONT')
        IF ( global%myProcid == MASTERPROC .AND. &
          global%verbLevel > VERBOSE_NONE ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
           'Reading PLAG_ReadDisPartnContSection...'
        END IF ! global%verbLevel
        CALL PLAG_ReadDisPartnContSection( regions )

      CASE ('# DISPART_INIT')
        IF ( global%myProcid == MASTERPROC .AND. &
          global%verbLevel > VERBOSE_NONE ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
           'Reading PLAG_ReadDisPartInitSection...'
        END IF ! global%verbLevel
        CALL PLAG_ReadDisPartInitSection( regions )
    END SELECT
  ENDDO

86   CONTINUE

! close file ------------------------------------------------------------------

  CLOSE( IF_INPUT,iostat=errorFlag )
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop( global, ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )
  END IF ! global%error

! set global%plagUsed ----------------------------------------------------------

  usedSomewhere   = .FALSE.
  unusedSomewhere = .FALSE.

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
    readStatus = regions(iReg)%plagInput%readStatus
    usedSomewhere   = usedSomewhere  .OR.(readStatus == 1)
    unusedSomewhere = unusedSomewhere.OR.(readStatus /= 1) ! == 0 or -1
  END DO ! iReg

  IF (usedSomewhere.AND.unusedSomewhere) THEN
    CALL ErrorStop( global,ERR_MP_ALLORNONE,__LINE__ )
  END IF ! usedSomewhere.AND.unusedSomewhere

  global%plagUsed = usedSomewhere

! finalization & error handling -----------------------------------------------

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global, ERR_FILE_READ,__LINE__,'File: '//TRIM(fname) )

999  CONTINUE

END SUBROUTINE PLAG_ReadInputFile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ReadInputFile.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:52  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:50:27  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:34  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 20:58:06  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2004/08/20 23:27:13  fnajjar
! Added Infrastructure for Plag prep tool
!
! Revision 1.4  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.3  2003/03/17 18:42:45  jblazek
! Added inDir to path of the input file.
!
! Revision 1.2  2003/02/26 23:38:30  jferry
! eliminated end=999 convention to ensure that files get closed
!
! Revision 1.1  2002/10/25 14:19:16  f-najjar
! Initial Import of Rocpart
!
!******************************************************************************

