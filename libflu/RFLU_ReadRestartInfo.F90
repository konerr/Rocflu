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
! Purpose: Read restart info file to get last output iteration or time.
!
! Description: None.
!
! Input: 
!   global              Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ReadRestartInfo.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ReadRestartInfo(global)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError
  USE ModMPI
  
  USE ModInterfaces, ONLY: RFLU_CloseRestartInfo, & 
                           RFLU_OpenRestartInfo

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: fileExists
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: dummyInteger,errorFlag
  REAL(RFREAL) :: dummyRFReal

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ReadRestartInfo.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction(global,'RFLU_ReadRestartInfo',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Reading restart info file...'
  END IF ! global%verbLevel

! ==============================================================================
! Open file
! ==============================================================================

  CALL RFLU_OpenRestartInfo(global,FILE_POSITION_START,fileExists)

! ==============================================================================
! Set or read restart iteration or time
! ==============================================================================

  IF ( global%flowType == FLOW_STEADY ) THEN ! steady flow
    global%currentIter = 0

    IF ( fileExists .EQV. .TRUE. ) THEN
      DO 
        READ(IF_RESTINFO,*,IOSTAT=errorFlag) dummyInteger
        
        IF ( errorFlag /= ERR_NONE ) THEN 
          EXIT
        ELSE 
          global%currentIter = dummyInteger
        END IF ! errorFlag
      END DO ! <empty>
    END IF ! fileExists
  ELSE ! unsteady flow
    ! BBR - Begin - initializing previousTime
    global%previousTime = 0.0_RFREAL
    ! BBR - End
    global%currentTime = 0.0_RFREAL
  
    IF ( fileExists .EQV. .TRUE. ) THEN
      DO 
        READ(IF_RESTINFO,*,IOSTAT=errorFlag) dummyRFReal
        
        IF ( errorFlag /= ERR_NONE ) THEN 
          EXIT
        ELSE
          ! BBR - Begin - Record dump time before last
          ! in case safe data dump out of sync with frequency dump
          global%previousTime = global%currentTime
          ! BBR - end 
          global%currentTime = dummyRFReal
        END IF ! errorFlag
      END DO ! <empty>
    END IF ! fileExists
  END IF ! global%flowType

! ==============================================================================
! Write info
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    IF ( global%flowType == FLOW_STEADY ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME, & 
                                       'Restart iteration:',global%currentIter
    ELSE 
      WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME, & 
                                          'Restart time:',global%currentTime      
    END IF ! global%flowType
  END IF ! global%myProcid

! ==============================================================================
! Close file
! ==============================================================================

  CALL RFLU_CloseRestartInfo(global)        

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Reading restart info file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ReadRestartInfo


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ReadRestartInfo.F90,v $
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
! Revision 1.4  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.3  2004/06/25 20:04:51  haselbac
! Removed some code and put into RFLU_SetRestartTimeFlag
!
! Revision 1.2  2004/06/16 20:00:33  haselbac
! Added variables to make file handling easier, cosmetics
!
! Revision 1.1  2003/06/20 22:32:30  haselbac
! Initial revision
!
! ******************************************************************************

