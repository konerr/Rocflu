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
! Purpose: Build basic file name, that is, file name consisting only of 
!   directory, case name, extension, and domain id.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!   id          Region index
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
!******************************************************************************
!
! $Id: BuildFileNameBasic.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE BuildFileNameBasic(global,dest,ext,id,fileName)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: destString,RCSIdentString
  
! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(*), INTENT(IN) :: ext
  INTEGER, INTENT(IN) :: dest,id
  TYPE(t_global), POINTER :: global 
  
  CHARACTER(CHRLEN), INTENT(OUT) :: fileName  
  
! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: BuildFileNameBasic.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction(global,'BuildFileNameBasic',__FILE__)
   
  IF ( dest == FILEDEST_INDIR ) THEN 
    WRITE(destString,'(A)') global%inDir
  ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
    WRITE(destString,'(A)') global%outDir
  ELSE 
    CALL ErrorStop(global,ERR_FILEDEST_INVALID,__LINE__)
  END IF ! dest
           
  WRITE(fileName,'(A,A,I5.5)') TRIM(destString)//TRIM(global%caseName)// & 
                               TRIM(ext),'_',id
   
  CALL DeregisterFunction(global)     
    
! ******************************************************************************
! End
! ******************************************************************************
  
END SUBROUTINE BuildFileNameBasic


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: BuildFileNameBasic.F90,v $
!   Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
!   merged rocflu micro and macro
!
!   Revision 1.1.1.1  2014/07/15 14:31:37  brollin
!   New Stable version
!
!   Revision 1.3  2008/12/06 08:43:31  mtcampbe
!   Updated license.
!
!   Revision 1.2  2008/11/19 22:16:46  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.1  2007/04/09 18:48:32  haselbac
!   Initial revision after split from RocfloMP
!
!   Revision 1.1  2007/04/09 17:59:24  haselbac
!   Initial revision after split from RocfloMP
!
!   Revision 1.2  2006/04/07 15:19:15  haselbac
!   Removed tabs
!
!   Revision 1.1  2004/12/01 16:48:06  haselbac
!   Initial revision after changing case
!
!   Revision 1.1  2003/01/28 16:12:15  haselbac
!   Initial revision
!
! ******************************************************************************

