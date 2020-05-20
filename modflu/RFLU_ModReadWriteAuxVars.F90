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
! Purpose: Suite of routines to read and write auxiliary vars files.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModReadWriteAuxVars.F90,v 1.2 2015/07/23 23:11:18 brollin Exp $
!
! Copyright: (c) 2006-2007 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModReadWriteAuxVars

  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  USE ModBuildFileNames, ONLY: BuildFileNameSteady, &
                               BuildFileNameUnsteady

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_ReadAuxVarsWrapper, &
            RFLU_WriteAuxVarsWrapper

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

  CHARACTER(CHRLEN), PRIVATE :: &
    RCSIdentString = '$RCSfile: RFLU_ModReadWriteAuxVars.F90,v $'

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS






! ******************************************************************************
!
! Purpose: Read auxiliary variable file in ASCII ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadAuxVarsASCII(pRegion)

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

    LOGICAL :: fileExists
    CHARACTER(CHRLEN) :: errorString,iFileName,sectionString,timeString1, &
                         timeString2
    INTEGER :: errorFlag,i,iDataSet,iFile,iLoc,iVars,j,loopCounter,nCellsTot, &
               nCellsExpected,nVars,nVarsExpected,precActual,precExpected, &
               rangeActual,rangeExpected
    REAL(RFREAL) :: currentTime
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvOld
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, open file
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadAuxVarsASCII',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading auxiliary ASCII vars '// &
                               'file...'
    END IF ! global%verbLevel

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.mixt.auxa', & 
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName) 
                                
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime
      END IF ! global%verbLevel
    ELSE
      CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.mixt.auxa', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)      

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                         'number:',global%currentIter
      END IF ! global%verbLevel
    ENDIF ! global%flowType

    iFile = IF_SOLUT
       
    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD", &
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU auxiliary vars file' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END IF ! TRIM

! -----------------------------------------------------------------------------
!   Precision and range
! -----------------------------------------------------------------------------

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# Precision and range' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END IF ! TRIM

    precExpected  = PRECISION(1.0_RFREAL)
    rangeExpected = RANGE(1.0_RFREAL)

    READ(iFile,'(2(I8))') precActual,rangeActual
    IF ( precActual < precExpected .OR. rangeActual < rangeExpected ) THEN
      CALL ErrorStop(global,ERR_PREC_RANGE,__LINE__)
    END IF ! precActual

! -----------------------------------------------------------------------------
!   Initial residual and physical time
! -----------------------------------------------------------------------------

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# Initial residual' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
    END IF ! TRIM

    READ(iFile,'(E23.16)') global%resInit

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# Physical time' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
    END IF ! TRIM

    READ(iFile,'(E23.16)') currentTime

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      IF ( global%currentTime < 0.0_RFREAL ) THEN
        global%currentTime = currentTime
      ELSE
        WRITE(timeString1,'(1PE11.5)') global%currentTime
        WRITE(timeString2,'(1PE11.5)') currentTime
        
        IF ( TRIM(timeString1) /= TRIM(timeString2) ) THEN
          WRITE(STDOUT,'(A,4(1X,A))') SOLVER_NAME, & 
                                      '*** WARNING *** Time mismatch:', & 
                                      TRIM(timeString1),'vs.',TRIM(timeString2)
                                      
          global%warnCounter = global%warnCounter + 1
        END IF ! global%currentTime
      END IF ! global%currentTime
    END IF ! global%flowType  

! ==============================================================================
!   Dimensions
! ==============================================================================

    pGrid => pRegion%grid

    nVarsExpected  = pRegion%mixtInput%nCvOld2
    nCellsExpected = pGrid%nCellsTot

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# Dimensions' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END IF ! TRIM
  
    READ(iFile,'(2(I16))') nCellsTot,nVars
    IF ( nCellsTot /= nCellsExpected ) THEN
      WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nCellsTot, &
                                                'but expected:',nCellsExpected
      CALL ErrorStop(global,ERR_INVALID_NCELLS,__LINE__,errorString)
    END IF ! nCellsExpected

    IF ( nVars /= nVarsExpected ) THEN
      WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nVars, &
                                                'but expected:',nVarsExpected
      CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
    END IF ! nVarsExpected

! ==============================================================================
!   Rest of file
! ==============================================================================

    iVars       = 0
    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile,'(A)') sectionString

! ------------------------------------------------------------------------------
!     Invalid section string
! ------------------------------------------------------------------------------

      IF ( sectionString(1:2) /= '# ' ) THEN
        IF ( global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
        END IF ! verbosityLevel

        CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
      END IF ! sectionString(1:2)

      sectionString = sectionString(3:)

      SELECT CASE ( TRIM(sectionString) )

! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------

        CASE ( 'End' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel

          EXIT

! ------------------------------------------------------------------------------
!       Extract dataset location in cv from section string
! ------------------------------------------------------------------------------

        CASE DEFAULT

          READ(sectionString,'(I2)') iLoc

          IF ( iLoc < 1 .AND. iLoc > pRegion%mixtInput%nCvOld2 ) THEN
            IF ( global%verbLevel > VERBOSE_LOW ) THEN
              sectionString = '# ' // sectionString
              WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
            END IF ! verbosityLevel

            CALL ErrorStop(global,ERR_INVALID_CVILOC,__LINE__)
          END IF ! iLoc

          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A,I2,A)') SOLVER_NAME,'CvOld location ', &
                                          iLoc,'...'
          END IF ! global%verbLevel

          pCvOld => pRegion%mixt%cvOld

          iVars = iVars + 1
          
          READ(iFile,'(5(E23.16))') (pCvOld(iLoc,j),j=1,pGrid%nCellsTot)

      END SELECT ! TRIM

! ------------------------------------------------------------------------------
!     Guard against infinite loop - might be unnecessary because of read errors?
! ------------------------------------------------------------------------------

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter
    END DO ! <empty>

! ==============================================================================
!   Check and information about number of variables read
! ==============================================================================

    IF ( iVars /= nVars ) THEN
      CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
    END IF ! iVar

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading auxiliary ASCII vars '// &
                               'file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadAuxVarsASCII







! ******************************************************************************
!
! Purpose: Read auxiliary vars file in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadAuxVarsBinary(pRegion)

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

    LOGICAL :: fileExists
    CHARACTER(CHRLEN) :: errorString,iFileName,sectionString,timeString1, &
                         timeString2
    INTEGER :: errorFlag,i,iDataSet,iFile,iLoc,iVars,j,loopCounter,nCellsTot, &
               nCellsExpected,nVars,nVarsExpected,precActual,precExpected, &
               rangeActual,rangeExpected
    REAL(RFREAL) :: currentTime
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvOld
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, open file
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadAuxVarsBinary',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary auxiliary vars ' // &
                               'file...'
    END IF ! global%verbLevel

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.mixt.aux', & 
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime
      END IF ! global%verbLevel
    ELSE
      CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.mixt.aux', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                         'number:',global%currentIter
      END IF ! global%verbLevel
    ENDIF ! global%flowType

    iFile = IF_SOLUT
       
!    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
!         IOSTAT=errorFlag)
! BBR - begin
    IF( global%solutFormat .EQ. FORMAT_BINARY )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
         IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_L )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
         ACCESS="SEQUENTIAL",CONVERT="LITTLE_ENDIAN",IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_B )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
         ACCESS="SEQUENTIAL",CONVERT="BIG_ENDIAN",IOSTAT=errorFlag)
    END IF 
! BBR - end 

    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
! Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU auxiliary vars file' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END IF ! TRIM

! -----------------------------------------------------------------------------
!   Precision and range
! -----------------------------------------------------------------------------

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Precision and range' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END IF ! TRIM

    precExpected  = PRECISION(1.0_RFREAL)
    rangeExpected = RANGE(1.0_RFREAL)

    READ(iFile) precActual,rangeActual
    IF ( precActual < precExpected .OR. rangeActual < rangeExpected ) THEN
      CALL ErrorStop(global,ERR_PREC_RANGE,__LINE__)
    END IF ! precActual

! -----------------------------------------------------------------------------
!   Initial residual and physical time
! -----------------------------------------------------------------------------

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Initial residual' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
    END IF ! TRIM

    READ(iFile) global%resInit

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Physical time' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
    END IF ! TRIM

    READ(iFile) currentTime
      
    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      IF ( global%currentTime < 0.0_RFREAL ) THEN
        global%currentTime = currentTime
      ELSE
        WRITE(timeString1,'(1PE11.5)') global%currentTime
        WRITE(timeString2,'(1PE11.5)') currentTime

        IF ( TRIM(timeString1) /= TRIM(timeString2) ) THEN
          WRITE(STDOUT,'(A,4(1X,A))') SOLVER_NAME, & 
                                      '*** WARNING *** Time mismatch:', & 
                                      TRIM(timeString1),'vs.',TRIM(timeString2)
                                      
          global%warnCounter = global%warnCounter + 1
        END IF ! global%currentTime
      END IF ! global%currentTime
    END IF ! global%flowType  

! ==============================================================================
!   Dimensions
! ==============================================================================

    pGrid => pRegion%grid

    nVarsExpected  = pRegion%mixtInput%nCvOld2
    nCellsExpected = pGrid%nCellsTot

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Dimensions' ) THEN
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END IF ! TRIM

    READ(iFile) nCellsTot,nVars

    IF ( nCellsTot /= nCellsExpected ) THEN
      WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nCellsTot, &
                                                'but expected:',nCellsExpected
      CALL ErrorStop(global,ERR_INVALID_NCELLS,__LINE__,errorString)
    END IF ! nCellsExpected

    IF ( nVars /= nVarsExpected ) THEN
      WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nVars, &
                                                'but expected:',nVarsExpected
      CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
    END IF ! nVarsExpected

! ==============================================================================
!   Rest of file
! ==============================================================================

    iVars       = 0
    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile) sectionString

! ------------------------------------------------------------------------------
!     Invalid section string
! ------------------------------------------------------------------------------

      IF ( sectionString(1:2) /= '# ' ) THEN
        IF ( global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
        END IF ! verbosityLevel

        CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
      END IF ! sectionString(1:2)

      sectionString = sectionString(3:)

      SELECT CASE ( TRIM(sectionString) )

! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------

        CASE ( 'End' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel

          EXIT

! ------------------------------------------------------------------------------
!       Extract dataset location in cv from section string
! ------------------------------------------------------------------------------

        CASE DEFAULT

          READ(sectionString,'(I2)') iLoc

          IF ( iLoc < 1 .AND. iLoc > pRegion%mixtInput%nCvOld2 ) THEN
            IF ( global%verbLevel > VERBOSE_LOW ) THEN
              sectionString = '# ' // sectionString
              WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
            END IF ! verbosityLevel

            CALL ErrorStop(global,ERR_INVALID_CVILOC,__LINE__)
          END IF ! iLoc

          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A,I2,A)') SOLVER_NAME,'CvOld location ', &
                                          iLoc,'...'
          END IF ! global%verbLevel

          pCvOld => pRegion%mixt%cvOld

          iVars = iVars + 1

          READ(iFile) (pCvOld(iLoc,j),j=1,pGrid%nCellsTot)

      END SELECT ! TRIM

! ------------------------------------------------------------------------------
!     Guard against infinite loop - might be unnecessary because of read errors?
! ------------------------------------------------------------------------------

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter

    END DO ! <empty>

! ==============================================================================
!   Check and information about number of variables read
! ==============================================================================

    IF ( iVars /= nVars ) THEN
      CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
    END IF ! iVar

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary auxiliary vars ' // &
                               'file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadAuxVarsBinary








! ******************************************************************************
!
! Purpose: Wrapper for reading of auxiliary vars files in ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadAuxVarsWrapper(pRegion)

#ifdef GENX
    USE RFLU_ModGENXIO, ONLY: RFLU_GENX_DecideReadFile, & 
                              RFLU_GENX_GetDataFlow
#endif

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

    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadAuxVarsWrapper',__FILE__)

! ******************************************************************************
!   Read mixture solution files
! ******************************************************************************

#ifdef GENX
    IF ( RFLU_GENX_DecideReadFile(global) .EQV. .FALSE. ) THEN 
#endif
      IF ( global%solutFormat == FORMAT_ASCII ) THEN
        CALL RFLU_ReadAuxVarsASCII(pRegion)
!      ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
! BBR - begin
       ELSE IF ( global%solutFormat == FORMAT_BINARY .OR. & 
                global%solutFormat == FORMAT_BINARY_L .OR. & 
                global%solutFormat == FORMAT_BINARY_B ) THEN 
! BBR - end
        CALL RFLU_ReadAuxVarsBinary(pRegion)
      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! global%solutFormat
#ifdef GENX
    ELSE 
      CALL RFLU_GENX_GetDataFlow(pRegion) 
    END IF ! RFLU_GENX_DecideReadFile         
#endif

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadAuxVarsWrapper








! ******************************************************************************
!
! Purpose: Write auxiliary vars file in ASCII ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteAuxVarsASCII(pRegion)

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
    INTEGER :: errorFlag,i,iFile,iLoc,j
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvOld
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, open file
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteAuxVarsASCII',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing auxiliary ASCII vars ' // &
                               'file...'
    END IF ! global%verbLevel

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.mixt.auxa', &
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime
      END IF ! global%verbLevel
    ELSE
      CALL BuildFileNameSteady(global,FILEDEST_OUTDIR,'.mixt.auxa', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                         'number:',global%currentIter
      END IF ! global%verbLevel
    ENDIF ! global%flowType

    iFile = IF_SOLUT      
    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", &
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    sectionString = '# ROCFLU auxiliary vars file'
    WRITE(iFile,'(A)') sectionString

    sectionString = '# Precision and range'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(2(I8))') PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)

    sectionString = '# Initial residual'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(E23.16)') global%resInit

    sectionString = '# Physical time'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(E23.16)') global%currentTime
    
! ==============================================================================
!   Dimensions
! ==============================================================================

    pGrid => pRegion%grid

    sectionString = '# Dimensions'
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(2(I16))') pGrid%nCellsTot,pRegion%mixtInput%nCvOld2

! ==============================================================================
!   Data
! ==============================================================================

    pCvOld => pRegion%mixt%cvOld

! ------------------------------------------------------------------------------
!   Writing density
! ------------------------------------------------------------------------------

    iLoc = CV_MIXT_DENS

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A,I2,A)') SOLVER_NAME,'CvOld location ',iLoc,'...'
    END IF ! global%verbLevel

    WRITE(sectionString,'(A,I2)') '# ',iLoc
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(5(E23.16))') (pCvOld(iLoc,j),j=1,pGrid%nCellsTot)
                
! ------------------------------------------------------------------------------
!   Writing pressure
! ------------------------------------------------------------------------------

    iLoc = CV_MIXT_PRES

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A,I2,A)') SOLVER_NAME,'CvOld location ',iLoc,'...'
    END IF ! global%verbLevel

    WRITE(sectionString,'(A,I2)') '# ',iLoc
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(5(E23.16))') (pCvOld(iLoc,j),j=1,pGrid%nCellsTot)
                
! ==============================================================================
!   End marker
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile,'(A)') sectionString

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII auxiliary vars ' // &
                               'file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteAuxVarsASCII







! ******************************************************************************
!
! Purpose: Write auxiliary vars file in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteAuxVarsBinary(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,i,iFile,iLoc,j
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvOld
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
!   Start, open file
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteAuxVarsBinary',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary auxiliary ' // &
                               'vars file...'
    END IF ! global%verbLevel

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.mixt.aux', &
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName)    

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime
      END IF ! global%verbLevel
    ELSE
      CALL BuildFileNameSteady(global,FILEDEST_OUTDIR,'.mixt.aux', &
                               pRegion%iRegionGlobal,global%currentIter, &
                               iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                         'number:',global%currentIter
      END IF ! global%verbLevel
    ENDIF ! global%flowType

    iFile = IF_SOLUT
!    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
!         IOSTAT=errorFlag)
! BBR - begin
    IF ( global%solutFormat .EQ. FORMAT_BINARY )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
         IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_L )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
         ACCESS="SEQUENTIAL",CONVERT="LITTLE_ENDIAN",IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_B )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
         ACCESS="SEQUENTIAL",CONVERT="BIG_ENDIAN",IOSTAT=errorFlag)
    END IF
! BBR - end 

    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    sectionString = '# ROCFLU auxiliary vars file'
    WRITE(iFile) sectionString

    sectionString = '# Precision and range'
    WRITE(iFile) sectionString
    WRITE(iFile) PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)

    sectionString = '# Initial residual'
    WRITE(iFile) sectionString
    WRITE(iFile) global%resInit

    sectionString = '# Physical time'
    WRITE(iFile) sectionString
    WRITE(iFile) global%currentTime
    
! ==============================================================================
!   Dimensions
! ==============================================================================

    pGrid => pRegion%grid

    sectionString = '# Dimensions'
    WRITE(iFile) sectionString
    WRITE(iFile) pGrid%nCellsTot,pRegion%mixtInput%nCvOld2

! ==============================================================================
!   Data
! ==============================================================================

    pCvOld => pRegion%mixt%cvOld

! ------------------------------------------------------------------------------
!   Writing density
! ------------------------------------------------------------------------------

    iLoc = CV_MIXT_DENS

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A,I2,A)') SOLVER_NAME,'CvOld location ',iLoc,'...'
    END IF ! global%verbLevel

    WRITE(sectionString,'(A,I2)') '# ',iLoc
    WRITE(iFile) sectionString
    WRITE(iFile) (pCvOld(iLoc,j),j=1,pGrid%nCellsTot)

! ------------------------------------------------------------------------------
!   Writing pressure
! ------------------------------------------------------------------------------

    iLoc = CV_MIXT_PRES

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A,I2,A)') SOLVER_NAME,'CvOld location ',iLoc,'...'
    END IF ! global%verbLevel

    WRITE(sectionString,'(A,I2)') '# ',iLoc
    WRITE(iFile) sectionString
    WRITE(iFile) (pCvOld(iLoc,j),j=1,pGrid%nCellsTot)

! ==============================================================================
!   End marker
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile) sectionString

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary auxiliary vars ' // &
                               'file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteAuxVarsBinary








! ******************************************************************************
!
! Purpose: Wrapper for writing of auxiliary vars files in ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_WriteAuxVarsWrapper(pRegion)

#ifdef GENX
    USE RFLU_ModGENXIO, ONLY: RFLU_GENX_DecideWriteFile, & 
                              RFLU_GENX_PutDataFlow
#endif

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

    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteAuxVarsWrapper',__FILE__)

! ******************************************************************************
!   Write mixture solution files
! ******************************************************************************

#ifdef GENX
    IF ( RFLU_GENX_DecideWriteFile(global) .EQV. .FALSE. ) THEN 
#endif
      IF ( global%solutFormat == FORMAT_ASCII ) THEN
        CALL RFLU_WriteAuxVarsASCII(pRegion)
!      ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
! BBR - begin
       ELSE IF ( global%solutFormat == FORMAT_BINARY .OR. & 
                global%solutFormat == FORMAT_BINARY_L .OR. &
                global%solutFormat == FORMAT_BINARY_B ) THEN
! BBR - end
        CALL RFLU_WriteAuxVarsBinary(pRegion)
      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! global%solutFormat
#ifdef GENX
    ELSE 
      CALL RFLU_GENX_PutDataFlow(pRegion) 
    END IF ! RFLU_GENX_DecideReadFile         
#endif

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteAuxVarsWrapper







! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLU_ModReadWriteAuxVars


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModReadWriteAuxVars.F90,v $
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
! Revision 1.1.1.1  2014/07/15 14:31:36  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:44  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:16:57  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2008/01/12 14:33:23  haselbac
! Bug fix: Superfluous closing parenthesis in WRITE statements
!
! Revision 1.1  2007/11/28 23:04:41  mparmar
! Initial revision
!
!
! ******************************************************************************

