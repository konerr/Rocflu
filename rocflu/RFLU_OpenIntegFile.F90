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
! Purpose: Open file for integral quantities history.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
!
! Copyright: (c) 2002-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_OpenIntegFile(global)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY : t_global
  USE ModMPI
  USE ModParameters
  
  USE ModBuildFileNames, ONLY: BuildFileNamePlain  
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: fname,RCSIdentString
  INTEGER :: errorFlag

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_OpenIntegFile.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction(global,'RFLU_OpenIntegFile',__FILE__)

! ==============================================================================
! Open file
! ==============================================================================

  IF ( global%myProcid == MASTERPROC ) THEN
    CALL BuildFileNamePlain(global,FILEDEST_OUTDIR,'.int',fname)

! ==============================================================================
!   Append to existing file (restart) or create new file
! ==============================================================================

    IF ( ( global%flowType == FLOW_UNSTEADY .AND. & 
           global%currentTime > 0.0_RFREAL ) .OR. &
         ( global%flowType == FLOW_STEADY   .AND. & 
           global%currentIter>1 ) ) THEN
      OPEN(IF_INTEG,FILE=fname,FORM='FORMATTED',STATUS='OLD', &
                     POSITION='APPEND',IOSTAT=errorFlag)                     
    ELSE
      OPEN(IF_INTEG,FILE=fname,FORM='FORMATTED',STATUS='UNKNOWN', &
                     IOSTAT=errorFlag)
    END IF ! global
   
    global%error = errorFlag
    IF ( global%error /= 0 ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname))
    END IF ! global%error
  END IF ! global%myProcid

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_OpenIntegFile

! ******************************************************************************
!
! RCS Revision history:
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:48  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:57  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:01  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.6  2004/06/16 20:01:10  haselbac
! Added use of ModBuildFileNames, cosmetics
!
! Revision 1.5  2003/01/28 14:44:07  haselbac
! Use common building of file name
!
! Revision 1.4  2002/10/08 15:49:29  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.3  2002/10/05 19:24:07  haselbac
! GENX integration: added outDir to file name
!
! Revision 1.2  2002/09/09 15:51:56  haselbac
! global now under region
!
! Revision 1.1  2002/05/04 17:01:59  haselbac
! Initial revision
!
! ******************************************************************************

