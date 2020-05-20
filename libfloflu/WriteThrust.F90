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
! Purpose: write thrust history into a file.
!
! Description: none.
!
! Input: global%thrustTotal = total thrust
!
! Output: to file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: WriteThrust.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE WriteThrust( global )

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY : t_global
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  CHARACTER(CHRLEN+4) :: fname

#ifdef MPI
  REAL(RFREAL) :: localThrust(2), globalThrust(2)
#endif

!******************************************************************************

  CALL RegisterFunction( global,'WriteThrust',__FILE__ )

! sum up data from other processors

#ifdef MPI
  localThrust(1) = global%thrustMom
  localThrust(2) = global%thrustPress

  CALL MPI_Reduce( localThrust,globalThrust,2,MPI_RFREAL,MPI_SUM,MASTERPROC, &
                   global%mpiComm,global%mpierr )
  IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

  global%thrustMom   = globalThrust(1)
  global%thrustPress = globalThrust(2)
#endif
  global%thrustTotal = global%thrustMom + global%thrustPress

! steady flow

  IF (global%flowType==FLOW_STEADY .AND. global%myProcid==MASTERPROC) THEN
    WRITE(IF_THRUST,1000,err=10) global%currentIter,global%thrustMom, &
                                 global%thrustPress,global%thrustTotal

! unsteady flow

  ELSE IF (global%flowType==FLOW_UNSTEADY .AND. global%myProcid==MASTERPROC) THEN
    WRITE(IF_THRUST,2000,err=10) global%currentTime,global%thrustMom, &
                                 global%thrustPress,global%thrustTotal
  ENDIF

! close and open file (instead of fflush)

  IF (global%thrustOpenClose .AND. global%myProcid==MASTERPROC) THEN
    WRITE(fname,'(A)') TRIM(global%outDir)//TRIM(global%casename)//'.thr'
    CLOSE(IF_THRUST)
    OPEN(IF_THRUST,FILE=fname,FORM='FORMATTED',STATUS='OLD',POSITION='APPEND')
  ENDIF

! finalize

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_WRITE,__LINE__,'Thrust history file.' )

1000 FORMAT(I6,1PE13.4,2E13.4)
2000 FORMAT(1PE12.5,3E13.4)

999  CONTINUE

END SUBROUTINE WriteThrust

!******************************************************************************
!
! RCS Revision history:
!
! $Log: WriteThrust.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 16:52:27  haselbac
! Initial revision after changing case
!
! Revision 1.1  2003/06/02 17:12:00  jblazek
! Added computation of thrust.
!
!******************************************************************************

