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
! Purpose: Collection of routines to read boundary condition data file.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModReadWriteBcDataFile.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModReadWriteBcDataFile

  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global 
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY: t_grid
  USE ModMPI
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_DecideReadWriteBcDataFile, & 
            RFLU_ReadBcDataFile, & 
            RFLU_WriteBcDataFile
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModReadWriteBcDataFile.F90,v $ $Revision: 1.1.1.1 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS







! ******************************************************************************
!
! Purpose: Decide whether need to read and/or write boundary-condition data 
!   file.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  LOGICAL FUNCTION RFLU_DecideReadWriteBcDataFile(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DecideReadWriteBcDataFile',__FILE__)
    
! ******************************************************************************
!   Set pointers and initialize
! ******************************************************************************    
    
    pGrid => pRegion%grid  
    
    RFLU_DecideReadWriteBcDataFile = .FALSE.  
    
! ******************************************************************************
!   Loop over patches
! ******************************************************************************

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
  
! ==============================================================================
!     If not coupled and have distribution, write data
! ==============================================================================    
  
      IF ( (pPatch%bcCoupled == BC_NOT_COUPLED) .AND. & 
           (pPatch%mixt%distrib == BCDAT_DISTRIB) ) THEN
! TEMPORARY
        RFLU_DecideReadWriteBcDataFile = .TRUE.
!        RFLU_DecideReadWriteBcDataFile = .FALSE.
! END TEMPORARY

        EXIT
      END IF ! pPatch%bcCoupled
    END DO ! iPatch 

! ******************************************************************************
!   End
! ******************************************************************************
  
    CALL DeregisterFunction(global)

  END FUNCTION RFLU_DecideReadWriteBcDataFile








! ******************************************************************************
!
! Purpose: Read boundary-condition data file.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: 
!   1. At present restricted to mixture. Best way to extend to MP not clear.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadBcDataFile(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: errorString,iFileName,sectionString
    INTEGER :: dummyInteger,errorFlag,iData,iFile,ifl,iPatch,iPatchGlobal, & 
               loopCounter,nBFaces,nData
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadBcDataFile',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading boundary-condition data...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid
    
! ******************************************************************************
!   Open file
! ******************************************************************************

    iFile = IF_DISTR

    CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.bcd', & 
                            pRegion%iRegionGlobal,iFileName) 

    OPEN(iFile,FILE=iFileName,FORM='FORMATTED',STATUS='OLD',IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error  

! ******************************************************************************
!   Read header
! ******************************************************************************

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU boundary-condition data file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM

! ******************************************************************************
!   Read data
! ******************************************************************************

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile,'(A)') sectionString

      SELECT CASE ( TRIM(sectionString) ) 

! ==============================================================================
!       Patch data
! ==============================================================================
        
        CASE ( '# Patch data' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Patch data...'
          END IF ! global%verbLevel       

          READ(iFile,*) iPatch,iPatchGlobal,nBFaces,nData  
            
! ------------------------------------------------------------------------------
!         Check that input correct
! ------------------------------------------------------------------------------            
                    
          IF ( iPatch > pGrid%nPatches ) THEN
            WRITE(errorString,'(A,1X,I3)') 'Patch index invalid:',iPatch 
            CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, & 
                           TRIM(errorString))
          END IF ! iPatch                                                  

          pPatch => pRegion%patches(iPatch)

          IF ( nBFaces /= pPatch%nBFaces ) THEN 
            WRITE(errorString,'(A,1X,I6)') 'Number of faces invalid:',nBFaces
            CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, & 
                           TRIM(errorString))
          END IF ! nBFaces          

          IF ( iPatchGlobal /= pPatch%iPatchGlobal ) THEN 
            WRITE(errorString,'(A,1X,I3)') 'Global patch index invalid:', & 
                                           iPatchGlobal
            CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, & 
                           TRIM(errorString))
          END IF ! iPatchGlobal           

          IF ( nData /= pPatch%mixt%nData ) THEN 
            WRITE(errorString,'(A,1X,I3)') & 
              'Number of pieces of data invalid:',nData
            CALL ErrorStop(global,ERR_BCDATA_VALUE_INVALID,__LINE__, & 
                           TRIM(errorString))
          END IF ! iPatchGlobal    
    
! ------------------------------------------------------------------------------
!         Read data
! ------------------------------------------------------------------------------    
    
          DO ifl = 1,pPatch%nBFaces
            DO iData = 1,pPatch%mixt%nData
              READ(iFile,'(1X,I6,1X,I2,1X,E23.16)') & 
                dummyInteger,dummyInteger,pPatch%mixt%vals(iData,ifl)
            END DO ! iData  
          END DO ! ifl

! ==============================================================================
!       End marker
! ==============================================================================
      
        CASE ( '# End' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel           

          EXIT

! ==============================================================================
!       Invalid section string
! ==============================================================================

        CASE DEFAULT
          IF ( global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
          END IF ! verbosityLevel           

          CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)      

      END SELECT ! TRIM

! ==============================================================================
!     Guard against infinite loop - might be unnecessary because of read errors?
! ==============================================================================
  
      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter

    END DO ! <empty> 

! ******************************************************************************
!   Close file
! ******************************************************************************
    
    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error  

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading boundary-condition data done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadBcDataFile








! ******************************************************************************
!
! Purpose: Write boundary-condition data file.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: 
!   1. At present restricted to mixture. Best way to extend to MP not clear.
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteBcDataFile(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: dummyInteger,errorFlag,iData,iFile,ifl,iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteBcDataFile',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing boundary-condition data...'
    END IF ! global%verbLevel
  
! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid
  
! ******************************************************************************
!   Open file
! ******************************************************************************

    iFile = IF_DISTR

    CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.bcd', & 
                            pRegion%iRegionGlobal,iFileName) 

    OPEN(iFile,FILE=iFileName,FORM='FORMATTED',STATUS='UNKNOWN', & 
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error  

! ******************************************************************************
!   Write header
! ******************************************************************************
 
    sectionString = '# ROCFLU boundary-condition data file'   
    WRITE(iFile,'(A)') sectionString

! ******************************************************************************
!   Loop over patches
! ******************************************************************************

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
  
! ==============================================================================
!     If not coupled and have distribution, write data
! ==============================================================================    
  
      IF ( (pPatch%bcCoupled == BC_NOT_COUPLED) .AND. & 
           (pPatch%mixt%distrib == BCDAT_DISTRIB) ) THEN
        WRITE(iFile,'(A)') '# Patch data' 
        WRITE(iFile,'(2(1X,I3),1X,I6,1X,I2)') iPatch,pPatch%iPatchGlobal, & 
                                              pPatch%nBFaces, & 
                                              pPatch%mixt%nData

        DO ifl = 1,pPatch%nBFaces
          DO iData = 1,pPatch%mixt%nData
            WRITE(iFile,'(1X,I6,1X,I2,1X,E23.16)') ifl,iData, &
              pPatch%mixt%vals(iData,ifl)
          END DO ! iData  
        END DO ! ifl
      END IF ! pPatch%bcCoupled
    END DO ! iPatch

! ******************************************************************************
!   Write footer
! ******************************************************************************
 
    sectionString = '# End'   
    WRITE(iFile,'(A)') sectionString

! ******************************************************************************
!   Close file
! ******************************************************************************

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error  

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
        'Writing boundary-condition data done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteBcDataFile









! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModReadWriteBcDataFile


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModReadWriteBcDataFile.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.5  2009/07/08 19:12:06  mparmar
! Uncommented RFLU_DecideReadWriteBcDataFile = .TRUE.
!
! Revision 1.4  2008/12/06 08:43:44  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:16:57  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/12/05 13:23:45  haselbac
! Rm unnecessary RCSIdentString decl, ifort compiler on vonkarman complained
!
! Revision 1.1  2007/04/09 18:49:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:41  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/08/19 15:39:14  mparmar
! Renamed patch variables
!
! Revision 1.3  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.2  2005/12/23 13:47:01  haselbac
! Temporarily disable bc data files to avoid hardcoding for runs with tbc
!
! Revision 1.1  2004/07/06 15:14:30  haselbac
! Initial revision
!
! ******************************************************************************

