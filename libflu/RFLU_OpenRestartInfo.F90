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
! Purpose: Open restart file.
!
! Description: None.
!
! Input: 
!   global              Pointer to global data
!   filePosition        Position at which file to be opened (start or end)
!
! Output: 
!   fileExists          Logical indicating whether file exists
!
! Notes: 
!   1. The filePosition parameter is needed because the restart info file
!      is opened with two goals. The first is to open it with the goal of
!      getting the last output iteration or time. The second is to be able
!      to write to it by appending additional lines.
!   2. The fileExists parameter is needed because on if the file exists
!      on restarting the code, it will need to be read to determine the 
!      last output iteration or time.
!
! ******************************************************************************
!
! $Id: RFLU_OpenRestartInfo.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_OpenRestartInfo(global,filePosition,fileExists)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError
  USE ModMPI
  
  USE ModBuildFileNames, ONLY: BuildFileNamePlain

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: filePosition
  TYPE(t_global), POINTER :: global

  LOGICAL, INTENT(OUT) :: fileExists

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: iFileName,RCSIdentString
  INTEGER :: errorFlag

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_OpenRestartInfo.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction(global,'RFLU_OpenRestartInfo',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Opening restart info file...'
  END IF ! global%verbLevel

! =============================================================================
! Open file
! =============================================================================

  CALL BuildFileNamePlain(global,FILEDEST_OUTDIR,'.rin',iFileName)

  INQUIRE(FILE=iFileName,EXIST=fileExists)
    
  IF ( fileExists .EQV. .TRUE. ) THEN
    IF ( filePosition == FILE_POSITION_START ) THEN 
      OPEN(IF_RESTINFO,FILE=iFileName,FORM='FORMATTED',STATUS='OLD', &
           IOSTAT=errorFlag) 
    ELSE IF ( filePosition == FILE_POSITION_END ) THEN 
      OPEN(IF_RESTINFO,FILE=iFileName,FORM='FORMATTED',STATUS='OLD', &
           POSITION='APPEND',IOSTAT=errorFlag)     
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! filePosition
  ELSE
    OPEN(IF_RESTINFO,FILE=iFileName,FORM='FORMATTED',STATUS='NEW', &
         IOSTAT=errorFlag)    
  END IF ! file

  global%error = errorFlag
  IF ( global%error /= 0 ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! *****************************************************************************
! End
! *****************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Opening restart info file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_OpenRestartInfo

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_OpenRestartInfo.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:35  mtcampbe
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
! Revision 1.3  2004/06/16 20:00:21  haselbac
! Added use of ModBuildFileNames, cosmetics
!
! Revision 1.2  2003/08/20 02:09:58  haselbac
! Changed verbosity conditions to reduce solver output in GENx runs
!
! Revision 1.1  2003/06/20 22:32:30  haselbac
! Initial revision
!
! ******************************************************************************

