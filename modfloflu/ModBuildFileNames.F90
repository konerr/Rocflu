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
! Purpose: Collection of utility routines for building file names.
!
! Description: None
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ModBuildFileNames.F90,v 1.4 2016/01/31 05:59:03 rahul Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE ModBuildFileNames

  USE ModParameters
  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModError

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Public data
! ==============================================================================

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: ModBuildFileNames.F90,v $ $Revision: 1.4 $'        

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: BuildFileNameBasic, &
            BuildFileNamePlain, & 
            BuildFileNamePlainSteady, &
            BuildFileNamePlainUnsteady, &
            BuildFileNameSteady, &
            BuildFileNameUnsteady, & 
            BuildRegionIdString, &
            BuildFileNameSteadyVTK, &
            BuildFileNameUnsteadyVTK, &
            BuildFileNameSteadyPVTK, &
            BuildFileNameUnsteadyPVTK            

! ==============================================================================
! Private functions
! ==============================================================================
                
! ******************************************************************************
! Routines
! ******************************************************************************
                
  CONTAINS




! ******************************************************************************
!
! Purpose: Build basic file name, that is, file name consisting only of 
!   directory, case name, extension, and region id.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!   regId       Region index
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE BuildFileNameBasic(global,dest,ext,regId,fileName)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: destString,regIdString
  
! ==============================================================================
!   Arguments
! ==============================================================================

    CHARACTER(*), INTENT(IN) :: ext
    INTEGER, INTENT(IN) :: dest,regId
    TYPE(t_global), POINTER :: global 

    CHARACTER(CHRLEN), INTENT(OUT) :: fileName  
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildFileNameBasic',__FILE__)

! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( dest == FILEDEST_INDIR ) THEN 
      WRITE(destString,'(A)') global%inDir
    ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
      WRITE(destString,'(A)') global%outDir
    ELSE 
      CALL ErrorStop(global,ERR_FILEDEST_INVALID,__LINE__)
    END IF ! dest

    CALL BuildRegionIdString(global,regId,regIdString)
#ifdef FOLDER 
!write(*,*) 'In Folder call'
!stop
    WRITE(destString,'(A)') 'GRIDDATA/'
    WRITE(fileName,'(A)') TRIM(destString)//TRIM(global%caseName)// & 
                          TRIM(ext)//'_'//TRIM(regIdString)
#else
    WRITE(fileName,'(A)') TRIM(destString)//TRIM(global%caseName)// & 
                          TRIM(ext)//'_'//TRIM(regIdString)
#endif    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)     
  
  END SUBROUTINE BuildFileNameBasic








! ******************************************************************************
!
! Purpose: Build plain file name, that is, file name consisting only of 
!   directory, case name, and extension.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE BuildFileNamePlain(global,dest,ext,fileName)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: destString
  
! ==============================================================================
!   Arguments
! ==============================================================================

    CHARACTER(*), INTENT(IN) :: ext
    INTEGER, INTENT(IN) :: dest
    TYPE(t_global), POINTER :: global 

    CHARACTER(CHRLEN), INTENT(OUT) :: fileName  
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildFileNamePlain',__FILE__)

! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( dest == FILEDEST_INDIR ) THEN 
      WRITE(destString,'(A)') global%inDir
    ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
      WRITE(destString,'(A)') global%outDir
    ELSE 
      CALL ErrorStop(global,ERR_FILEDEST_INVALID,__LINE__)
    END IF ! dest

    WRITE(fileName,'(A)') TRIM(destString)//TRIM(global%caseName)//TRIM(ext)   
    
! ******************************************************************************
!   End
! ******************************************************************************
  
    CALL DeregisterFunction(global)     
  
  END SUBROUTINE BuildFileNamePlain








! ******************************************************************************
!
! Purpose: Build plain file name with stamp, that is, file name consisting 
!   only of directory, case name, iteration stamp, and extension.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!   it          Iteration counter
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE BuildFileNamePlainSteady(global,dest,ext,it,fileName)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: destString
  
! ==============================================================================
!   Arguments
! ==============================================================================

    CHARACTER(*), INTENT(IN) :: ext
    INTEGER, INTENT(IN) :: dest,it
    TYPE(t_global), POINTER :: global 

    CHARACTER(CHRLEN), INTENT(OUT) :: fileName  
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildFileNamePlainSteady',__FILE__)

! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( dest == FILEDEST_INDIR ) THEN 
      WRITE(destString,'(A)') global%inDir
    ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
      WRITE(destString,'(A)') global%outDir
    ELSE 
      CALL ErrorStop(global,ERR_FILEDEST_INVALID,__LINE__)
    END IF ! dest
#ifdef FOLDER
    ! BBR - begin
    WRITE(destString,'(A,I6.6,A)') 'SOL_',it,'/'
    ! BBR -end
#endif
    WRITE(fileName,'(A,A,I6.6,A)') & 
      TRIM(destString)//TRIM(global%caseName),'_',it,TRIM(ext) 
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)     
  
  END SUBROUTINE BuildFileNamePlainSteady








! ******************************************************************************
!
! Purpose: Build plain file name with stamp, that is, file name consisting 
!   only of directory, case name, time stamp, and extension.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!   tm          Time stamp
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE BuildFileNamePlainUnsteady(global,dest,ext,tm,fileName)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: destString
  
! ==============================================================================
!   Arguments
! ==============================================================================

    CHARACTER(*), INTENT(IN) :: ext
    INTEGER, INTENT(IN) :: dest
    REAL(RFREAL), INTENT(IN) :: tm
    TYPE(t_global), POINTER :: global 

    CHARACTER(CHRLEN), INTENT(OUT) :: fileName  

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildFileNamePlainUnsteady',__FILE__)

! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( dest == FILEDEST_INDIR ) THEN 
      WRITE(destString,'(A)') global%inDir
    ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
      WRITE(destString,'(A)') global%outDir
    ELSE 
      CALL ErrorStop(global,ERR_FILEDEST_INVALID,__LINE__)
    END IF ! dest
#ifdef FOLDER
    ! BBR - begin
    WRITE(destString,'(A,1PE11.5,A)') 'PARAVIEW_',tm,'/'
    ! BBR -end
#endif

    WRITE(fileName,'(A,A,1PE11.5,A)') & 
      TRIM(destString)//TRIM(global%caseName),'_',tm,TRIM(ext) 

! ******************************************************************************
! End
! ******************************************************************************

    CALL DeregisterFunction(global)     
  
  END SUBROUTINE BuildFileNamePlainUnsteady








! ******************************************************************************
!
! Purpose: Build file name for steady flow, that is, file name consisting of 
!   directory, case name, extension, region id, and iteration counter.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!   regId       Region index
!   it          Iteration counter
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE BuildFileNameSteady(global,dest,ext,regId,it,fileName)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: destString,regIdString
  
! ==============================================================================
!   Arguments
! ==============================================================================

    CHARACTER(*), INTENT(IN) :: ext
    INTEGER, INTENT(IN) :: dest,regId,it
    TYPE(t_global), POINTER :: global 

    CHARACTER(CHRLEN), INTENT(OUT) :: fileName  
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildFileNameSteady',__FILE__)

! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( dest == FILEDEST_INDIR ) THEN 
      WRITE(destString,'(A)') global%inDir
    ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
      WRITE(destString,'(A)') global%outDir
    ELSE 
      CALL ErrorStop(global,ERR_FILEDEST_INVALID,__LINE__)
    END IF ! dest

    CALL BuildRegionIdString(global,regId,regIdString)
#ifdef FOLDER
    ! BBR - begin
    WRITE(destString,'(A,I6.6,A)') 'SOL_',it,'/'
    ! BBR - end
#endif
    WRITE(fileName,'(A,I6.6)') TRIM(destString)//TRIM(global%caseName)// & 
                               TRIM(ext)//'_'//TRIM(regIdString)//'_',it
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)     
  
  END SUBROUTINE BuildFileNameSteady



!begin BBR
! ******************************************************************************
!
! Purpose: Build file name for "region" parallel Paraview visualization of 
!          steady flow, that is, file name consisting of directory, case name,
!           extension, region id, and iteration counter.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!   regId       Region index
!   it          Iteration counter
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE BuildFileNameSteadyVTK(global,dest,ext,regId,it,fileName)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: destString,regIdString

! ==============================================================================
!   Arguments
! ==============================================================================

    CHARACTER(*), INTENT(IN) :: ext
    INTEGER, INTENT(IN) :: dest,regId,it
    TYPE(t_global), POINTER :: global

    CHARACTER(CHRLEN), INTENT(OUT) :: fileName
    CHARACTER(LEN=6) :: iter

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildFileNameSteadyVTK',__FILE__)

! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( dest == FILEDEST_INDIR ) THEN
      WRITE(destString,'(A)') global%inDir
    ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
      WRITE(destString,'(A)') global%outDir
    ELSE
      CALL ErrorStop(global,ERR_FILEDEST_INVALID,__LINE__)
    END IF ! dest

    CALL BuildRegionIdString(global,regId,regIdString)

!    WRITE(fileName,'(A,I6.6)') TRIM(destString)//TRIM(global%caseName)// &
!                               TRIM(ext)//'_'//TRIM(regIdString)//'_',it
#ifdef FOLDER
    ! BBR - begin
    WRITE(destString,'(A,I6.6,A)') 'PARAVIEW_',it,'/' 
    ! BBR - end
#endif
    WRITE(fileName,'(A,I6.6,A)') TRIM(destString)//TRIM(global%caseName)//'_', &
                               it,'_'//TRIM(regIdString)//TRIM(ext)
!#else
!    WRITE(iter,'(I6.6)') it
!    WRITE(fileName,'(A)') TRIM(destString)//TRIM(global%caseName)//'_'// &
!                               iter//'_'//TRIM(regIdString)//TRIM(ext)
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE BuildFileNameSteadyVTK




! ******************************************************************************
!
! Purpose: Build file name for grouping regions for parallel Paraview 
!          visualization of steady flow, that is, 
!          file name consisting of directory, case name, extension, region id, 
!          and iteration counter.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!   it          Iteration counter
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE BuildFileNameSteadyPVTK(global,dest,ext,it,fileName)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: destString,regIdString

! ==============================================================================
!   Arguments
! ==============================================================================

    CHARACTER(*), INTENT(IN) :: ext
    INTEGER, INTENT(IN) :: dest,it
    TYPE(t_global), POINTER :: global

    CHARACTER(CHRLEN), INTENT(OUT) :: fileName
    CHARACTER(LEN=6) :: iter

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildFileNameSteadyPVTK',__FILE__)

! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( dest == FILEDEST_INDIR ) THEN
      WRITE(destString,'(A)') global%inDir
    ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
      WRITE(destString,'(A)') global%outDir
    ELSE
      CALL ErrorStop(global,ERR_FILEDEST_INVALID,__LINE__)
    END IF ! dest

!    WRITE(fileName,'(A,I6.6)') TRIM(destString)//TRIM(global%caseName)// &
!                               TRIM(ext)//'_',it
#ifdef FOLDER
    ! BBR - begin
    WRITE(destString,'(A,I6.6,A)') 'PARAVIEW_',it,'/'
    ! BBR -end
#endif
    WRITE(fileName,'(A,I6.6,A)') TRIM(destString)//TRIM(global%caseName)// &
                               '_',it,TRIM(ext)

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE BuildFileNameSteadyPVTK
!end BBR - line 521




! ******************************************************************************
!
! Purpose: Build file name for unsteady flow, that is, file name consisting of 
!   directory, case name, extension, region id, and time stamp.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!   regId       Region index
!   tm          Time stamp
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
! ******************************************************************************
 
  SUBROUTINE BuildFileNameUnsteady(global,dest,ext,regId,tm,fileName)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: destString,regIdString
  
! ==============================================================================
!   Arguments
! ==============================================================================

    CHARACTER(*), INTENT(IN) :: ext
    INTEGER, INTENT(IN) :: dest,regId
    REAL(RFREAL), INTENT(IN) :: tm
    TYPE(t_global), POINTER :: global 

    CHARACTER(CHRLEN), INTENT(OUT) :: fileName  
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildFileNameUnsteady',__FILE__)

! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( dest == FILEDEST_INDIR ) THEN 
      WRITE(destString,'(A)') global%inDir
    ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
      WRITE(destString,'(A)') global%outDir
    ELSE 
      CALL ErrorStop(global,ERR_FILEDEST_INVALID,__LINE__)
    END IF ! dest

    CALL BuildRegionIdString(global,regId,regIdString)
#ifdef FOLDER
    ! BBR - begin
    WRITE(destString,'(A,1PE11.5,A)') 'SOL_',tm,'/'
    ! BBR - end
#endif
    WRITE(fileName,'(A,1PE11.5)') TRIM(destString)//TRIM(global%caseName)// & 
                                  TRIM(ext)//'_'//TRIM(regIdString)//'_',tm
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)     
  
  END SUBROUTINE BuildFileNameUnsteady




! begin BBR
! ******************************************************************************
!
! Purpose: Build file name for "region" parallel PARAVIEW visualization of 
!   unsteady flow, that is, file name consisting of directory, case name,
!   extension, region id, and time stamp.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!   regId       Region index
!   tm          Time stamp
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE BuildFileNameUnsteadyVTK(global,dest,ext,regId,tm,fileName)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: destString,regIdString

! ==============================================================================
!   Arguments
! ==============================================================================

    CHARACTER(*), INTENT(IN) :: ext
    INTEGER, INTENT(IN) :: dest,regId
    REAL(RFREAL), INTENT(IN) :: tm
    TYPE(t_global), POINTER :: global

    CHARACTER(CHRLEN), INTENT(OUT) :: fileName
    CHARACTER(LEN=11) :: atm

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildFileNameUnsteadyVTK',__FILE__)

! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( dest == FILEDEST_INDIR ) THEN
      WRITE(destString,'(A)') global%inDir
    ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
      WRITE(destString,'(A)') global%outDir
    ELSE
      CALL ErrorStop(global,ERR_FILEDEST_INVALID,__LINE__)
    END IF ! dest

    CALL BuildRegionIdString(global,regId,regIdString)

!    WRITE(fileName,'(A,1PE11.5)') TRIM(destString)//TRIM(global%caseName)// &
!                                  TRIM(ext)//'_'//TRIM(regIdString)//'_',tm
#ifdef FOLDER
   ! BBR - begin
    WRITE(destString,'(A,1PE11.5,A)') 'PARAVIEW_',tm,'/'
   ! BBR -end
#endif
    WRITE(fileName,'(A,1PE11.5,A)') TRIM(destString)//TRIM(global%caseName)//'_', &
                          tm,'_'//TRIM(regIdString)//TRIM(ext)
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE BuildFileNameUnsteadyVTK





! ******************************************************************************
!
! Purpose: Build file name for grouping parallel PARAVIEW visualization of 
!   unsteady flow, that is, file name consisting of directory, case name,
!   extension, region id, and time stamp.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!   tm          Time stamp
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE BuildFileNameUnsteadyPVTK(global,dest,ext,tm,fileName)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: destString,regIdString

! ==============================================================================
!   Arguments
! ==============================================================================

    CHARACTER(*), INTENT(IN) :: ext
    INTEGER, INTENT(IN) :: dest
    REAL(RFREAL), INTENT(IN) :: tm
    TYPE(t_global), POINTER :: global

    CHARACTER(CHRLEN), INTENT(OUT) :: fileName
    CHARACTER(LEN=11) :: atm

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildFileNameUnsteadyPVTK',__FILE__)

! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( dest == FILEDEST_INDIR ) THEN
      WRITE(destString,'(A)') global%inDir
    ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
      WRITE(destString,'(A)') global%outDir
    ELSE
      CALL ErrorStop(global,ERR_FILEDEST_INVALID,__LINE__)
    END IF ! dest

!    WRITE(fileName,'(A,1PE11.5)') TRIM(destString)//TRIM(global%caseName)// &
!                                  TRIM(ext)//'_'//TRIM(regIdString)//'_',tm
#ifdef FOLDER
   ! BBR - begin
    WRITE(destString,'(A,1PE11.5,A)') 'PARAVIEW_',tm,'/'
   ! BBR -end
#endif
    WRITE(fileName,'(A,1PE11.5,A)') TRIM(destString)//TRIM(global%caseName)// &
                          '_',tm,TRIM(ext)

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE BuildFileNameUnsteadyPVTK
!end BBR -line 771


! ******************************************************************************
!
! Purpose: Build region id string.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   regId       Region index
!
! Output: 
!   regIdString Region string
!
! Notes: None.
!
! ******************************************************************************
 
  SUBROUTINE BuildRegionIdString(global,regId,regIdString)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: regId
    TYPE(t_global), POINTER :: global 

    CHARACTER(CHRLEN), INTENT(OUT) :: regIdString 
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildRegionIdString',__FILE__)

! ******************************************************************************
!   Write region id into region string
! ******************************************************************************

    WRITE(regIdString,'(I5.5)') regId
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)     
  
  END SUBROUTINE BuildRegionIdString







END MODULE ModBuildFileNames

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModBuildFileNames.F90,v $
! Revision 1.4  2016/01/31 05:59:03  rahul
! Fixed a bug pertaining to FOLDER definition.
!
! Revision 1.3  2015/12/19 00:38:29  rahul
! Added compiler flag FOLDER.
!
! Revision 1.2  2015/07/23 23:11:18  brollin
! 1) The pressure coefficient of the  collision model has been changed back to its original form
! 2) New options in the format of the grid and solutions have been added. Now the user can choose the endianness, and convert from one to the over in rfluconv.
! 3) The solutions are now stored in folders named by timestamp or iteration number.
! 4) The address enty in the hashtable has been changed to an integer(8) for cases when the grid becomes very large.
! 5) RFLU_WritePM can now compute PM2 on the fly for the Macroscale problem
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:37  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:10  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:16  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2006/04/07 15:19:18  haselbac
! Removed tabs
!
! Revision 1.2  2004/10/19 19:28:38  haselbac
! Added BuildRegionIdString bcos needed in GENX modules, cosmetics
!
! Revision 1.1  2004/06/16 20:00:43  haselbac
! Initial revision
!
! ******************************************************************************

