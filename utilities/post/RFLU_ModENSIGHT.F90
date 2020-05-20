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
! Purpose: Collection of routines to write ENSIGHT files.
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
! $Id: RFLU_ModENSIGHT.F90,v 1.2 2015/07/23 23:11:19 brollin Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModENSIGHT

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region

  USE RFLU_ModENSIGHTUtils
  USE RFLU_ModPlottingVars, ONLY: RFLU_CountPlottingVars

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PARAMETER, PRIVATE :: &
    RCSIdentString = '$RCSfile: RFLU_ModENSIGHT.F90,v $ $Revision: 1.2 $'

  INTEGER, PRIVATE :: nScalars,nScalars1,nVectors,partNumber,partNumberSave

  TYPE t_var_info
    PRIVATE
    CHARACTER(CHRLEN) :: iFileName,varName,varNameShort
    INTEGER :: iFile
  END TYPE t_var_info
  
  TYPE(t_var_info) :: geoInfo  
  TYPE(t_var_info), DIMENSION(:), ALLOCATABLE :: scalarInfo
  TYPE(t_var_info), DIMENSION(:), ALLOCATABLE :: vectorInfo   

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: RFLU_ENS_BuildDataInfo, & 
            RFLU_ENS_CloseFileGeometry, & 
            RFLU_ENS_DestroyDataInfo, & 
            RFLU_ENS_InitPartNumber, &
            RFLU_ENS_MapRegion2Server, &   
            RFLU_ENS_OpenFileGeometry, & 
            RFLU_ENS_WriteFlowWrapper, &  
            RFLU_ENS_WriteGridWrapper

! ******************************************************************************
! Subroutines and functions
! ******************************************************************************

  CONTAINS






! ******************************************************************************
!
! Purpose: Build scalar information.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region type
!   iServer             Server index
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ENS_BuildDataInfo(pRegion,iServer)

  USE ModBuildFileNames, ONLY: BuildFileNameSteady, & 
                               BuildFileNameUnsteady

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iServer
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: dummyString
  INTEGER :: errorFlag,iFileOffs,iPv,iPv2,iScalar,iScalarOffs,iVector
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ENS_BuildDataInfo', &
                        __FILE__)

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building ENSIGHT data info...'
  END IF ! global%verbLevel

! ******************************************************************************
! Count number of scalars
! ******************************************************************************

! ==============================================================================
! Mixture
! ==============================================================================

  nScalars = 5
  nVectors = 1
  
! ==============================================================================
! Plotting variables
! ==============================================================================

  IF ( pRegion%plot%nPv > 0 ) THEN
    nScalars = nScalars + pRegion%plot%nPv
  END IF ! nPv 

! ==============================================================================
! Physical modules
! ==============================================================================
  
#ifdef SPEC
! ------------------------------------------------------------------------------
! Species
! ------------------------------------------------------------------------------

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    nScalars = nScalars + global%nSpecies
  END IF ! global%specUsed
#endif
  
! ******************************************************************************
! Allocate memory
! ******************************************************************************

  ALLOCATE(scalarInfo(nScalars),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'scalarInfo')
  END IF ! global%error

  ALLOCATE(vectorInfo(nVectors),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vectorInfo')
  END IF ! global%error

! ******************************************************************************
! Set names and file indices
! ******************************************************************************

! ==============================================================================
! Mixture
! ==============================================================================

  nScalars = 5
  nVectors = 1

! ==============================================================================
! Plotting variables
! ==============================================================================

  IF ( pRegion%plot%nPv > 0 ) THEN
    nScalars = nScalars + pRegion%plot%nPv
  END IF ! nPv 

! ------------------------------------------------------------------------------
! Scalars
! ------------------------------------------------------------------------------

  scalarInfo(1)%varName = 'Density'
  scalarInfo(2)%varName = 'Energy'
  scalarInfo(3)%varName = 'Pressure'
  scalarInfo(4)%varName = 'Temperature'
  scalarInfo(5)%varName = 'Soundspeed'
  
  scalarInfo(1)%varNameShort = 'r'
  scalarInfo(2)%varNameShort = 'rE'
  scalarInfo(3)%varNameShort = 'p'
  scalarInfo(4)%varNameShort = 'T'
  scalarInfo(5)%varNameShort = 'a'

  nScalars1 = 5
  iScalar   = 0

  IF ( pRegion%plot%nPv > 0 ) THEN
    DO iPv = 1,pRegion%plot%nPv
      iPv2 = pRegion%plot%pvi2pv(iPv)
      
      iScalar = iScalar+1
      
      scalarInfo(nScalars1+iScalar)%varName      = pRegion%plot%pvNameLong(iPv2)
      scalarInfo(nScalars1+iScalar)%varNameShort = pRegion%plot%pvNameShort(iPv2)      
    END DO ! iPv
  END IF ! nPv 

  scalarInfo(1)%iFile = IF_ENS_SCALAR
  scalarInfo(2)%iFile = IF_ENS_SCALAR + 1
  scalarInfo(3)%iFile = IF_ENS_SCALAR + 2
  scalarInfo(4)%iFile = IF_ENS_SCALAR + 3
  scalarInfo(5)%iFile = IF_ENS_SCALAR + 4 

  nScalars1 = 5

  IF ( pRegion%plot%nPv > 0 ) THEN
    DO iScalar=1,pRegion%plot%nPv
      scalarInfo(5+iScalar)%iFile = IF_ENS_SCALAR + nScalars1 + iScalar - 1 
    END DO ! iScalar  
  END IF ! nPv 

  iScalarOffs = 5

  IF ( pRegion%plot%nPv > 0 ) THEN
    iScalarOffs = iScalarOffs + pRegion%plot%nPv
  END IF ! nPv 

! ------------------------------------------------------------------------------
! Vectors
! ------------------------------------------------------------------------------

  vectorInfo(1)%varName = 'Momentum'
  
  vectorInfo(1)%varNameShort = 'rv'       

  vectorInfo(1)%iFile = IF_ENS_VECTOR
    
! ==============================================================================
! Physical modules
! ==============================================================================

  iFileOffs = scalarInfo(iScalarOffs)%iFile 

#ifdef SPEC
! ------------------------------------------------------------------------------
! Species
! ------------------------------------------------------------------------------

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    DO iScalar = 1,global%nSpecies
      WRITE(dummyString,'(A,I2.2)') 'Density',iScalar  
      scalarInfo(iScalar+iScalarOffs)%varName = TRIM(dummyString)

      WRITE(dummyString,'(A,I2.2)') 'rY',iScalar 
      scalarInfo(iScalar+iScalarOffs)%varNameShort = TRIM(dummyString)

      scalarInfo(iScalar+iScalarOffs)%iFile = iScalar + iFileOffs
    END DO ! iScalar
  
    iScalarOffs = iScalarOffs + global%nSpecies
    nScalars    = nScalars    + global%nSpecies
  END IF ! global%specUsed
#endif

! ******************************************************************************
! Write info
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN 
    WRITE(STDOUT,'(A,3X,A,1X,I2)') SOLVER_NAME,'Server:',iServer
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Scalar variables:'   
    WRITE(STDOUT,'(A,5X,A,2X,A,2X,A,26X,A)') SOLVER_NAME,'#','Index', &
                                             'Long name','Short name'

    DO iScalar = 1,nScalars
      WRITE(STDOUT,'(A,4X,I2,3X,I3,3X,A32,3X,A5)') & 
            SOLVER_NAME,iScalar,scalarInfo(iScalar)%iFile, & 
            ADJUSTL(scalarInfo(iScalar)%varName), & 
            ADJUSTL(scalarInfo(iScalar)%varNameShort)
    END DO ! iScalar

    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Vector variables:'   
    WRITE(STDOUT,'(A,5X,A,2X,A,2X,A,26X,A)') SOLVER_NAME,'#','Index', &
                                             'Long name','Short name'

    DO iVector = 1,nVectors
      WRITE(STDOUT,'(A,4X,I2,3X,I3,3X,A32,3X,A5)') &
            SOLVER_NAME,iVector,vectorInfo(iVector)%iFile, &      
            ADJUSTL(vectorInfo(iVector)%varName), & 
            ADJUSTL(vectorInfo(iVector)%varNameShort)
    END DO ! iVector
  END IF ! global%verbLevel

! ******************************************************************************
! Set file names
! ******************************************************************************

! ==============================================================================
! Geometry
! ==============================================================================

  CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.geo',iServer, &
                           global%currentIter,geoInfo%iFileName)

! ==============================================================================
! Scalars
! ==============================================================================

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    DO iScalar = 1,nScalars
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR, & 
                                 '.'//TRIM(scalarInfo(iScalar)%varNameShort), & 
                                 iServer,global%currentTime, &
                                 scalarInfo(iScalar)%iFileName)
    END DO ! iScalar
  ELSE IF ( global%flowType == FLOW_STEADY ) THEN
    DO iScalar = 1,nScalars
      CALL BuildFileNameSteady(global,FILEDEST_INDIR, & 
                               '.'//TRIM(scalarInfo(iScalar)%varNameShort), &
                               iServer,global%currentIter, &
                               scalarInfo(iScalar)%iFileName)
    END DO ! iScalar
  ELSE ! defensive coding
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END IF ! global%flowType 
  
! ==============================================================================
! Vectors
! ==============================================================================
  
  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    DO iVector = 1,nVectors
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR, & 
                                 '.'//TRIM(vectorInfo(iVector)%varNameShort), & 
                                 iServer,global%currentTime, &
                                 vectorInfo(iVector)%iFileName)
    END DO ! iScalar
  ELSE IF ( global%flowType == FLOW_STEADY ) THEN
    DO iVector = 1,nVectors
      CALL BuildFileNameSteady(global,FILEDEST_INDIR, & 
                               '.'//TRIM(vectorInfo(iVector)%varNameShort), &
                               iServer,global%currentIter, &
                               vectorInfo(iVector)%iFileName)
    END DO ! iVector
  ELSE ! defensive coding
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END IF ! global%flowType      
         
! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building ENSIGHT data info done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ENS_BuildDataInfo








! ******************************************************************************
!
! Purpose: Open ENSIGHT scalar file.
!
! Description: None.
!
! Input:
!   global                Pointer to global type
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ENS_CloseFileFlowWrapper(global)

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

  CHARACTER(80) :: dummyString
  CHARACTER(CHRLEN) :: iFileName 
  INTEGER :: errorFlag,iScalar,iVector

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_ENS_CloseFileFlowWrapper', &
                        __FILE__)

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing ENSIGHT flow files...'
  END IF ! global%verbLevel

! ******************************************************************************
! Close files
! ******************************************************************************

  DO iScalar = 1,nScalars
    CALL RFLU_ENS_CloseFileScalarVector(global,scalarInfo(iScalar)%iFile)
  END DO ! iScalar
  
  DO iVector = 1,nVectors
    CALL RFLU_ENS_CloseFileScalarVector(global,vectorInfo(iVector)%iFile)
  END DO ! iVector
         
! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing ENSIGHT flow files done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ENS_CloseFileFlowWrapper








! ******************************************************************************
!
! Purpose: Close ENSIGHT geometry file.
!
! Description: None.
!
! Input:
!   global                Pointer to global type
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ENS_CloseFileGeometry(global)

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

  INTEGER :: errorFlag

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_ENS_CloseFileGeometry', &
                        __FILE__)

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing ENSIGHT geometry file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_ENS_GEOMETRY,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__)
  END IF ! global%error    

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing ENSIGHT geometry file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ENS_CloseFileGeometry






! ******************************************************************************
!
! Purpose: Close ENSIGHT scalar/vector file.
!
! Description: None.
!
! Input:
!   global                Pointer to global type
!   iFile                 File index
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ENS_CloseFileScalarVector(global,iFile)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iFile
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_ENS_CloseFileScalarVector', &
                        __FILE__)

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(iFile,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__)
  END IF ! global%error    

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ENS_CloseFileScalarVector







! ******************************************************************************
!
! Purpose: Destroy data information.
!
! Description: None.
!
! Input:
!   global                Pointer to global type
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ENS_DestroyDataInfo(global)

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

  INTEGER :: errorFlag

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_ENS_DestroyDataInfo', &
                        __FILE__)

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying ENSIGHT data info...'
  END IF ! global%verbLevel

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DEALLOCATE(scalarInfo,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'scalarInfo')
  END IF ! global%error
  
  DEALLOCATE(vectorInfo,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vectorInfo')
  END IF ! global%error  

! ******************************************************************************
! Reset number of data items
! ******************************************************************************

  nScalars = 0
  nVectors = 0

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying ENSIGHT data info done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ENS_DestroyDataInfo







! ******************************************************************************
!
! Purpose: Initialize part number.
!
! Description: None.
!
! Input:
!   global                Pointer to global type
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ENS_InitPartNumber(global)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global%postPartNumber = 0
         
! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE RFLU_ENS_InitPartNumber







! ******************************************************************************
!
! Purpose: Map regions to servers.
!
! Description: None.
!
! Input:
!   iReg        Region index
!   nRegions    Number of regions
!   nServers    Number of servers
!
! Output: None.
!
! Notes: 
!   1. Note use of integer division.
!
! ******************************************************************************

INTEGER FUNCTION RFLU_ENS_MapRegion2Server(iReg,nRegions,nServers)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iReg,nRegions,nServers

! ******************************************************************************
! Start
! ******************************************************************************

  RFLU_ENS_MapRegion2Server = MIN((iReg-1)/nServers + 1,nServers)
         
! ******************************************************************************
! End
! ******************************************************************************

END FUNCTION RFLU_ENS_MapRegion2Server






! ******************************************************************************
!
! Purpose: Open ENSIGHT scalar file.
!
! Description: None.
!
! Input:
!   global              Pointer to global type
!   iServer             Index of server
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ENS_OpenFileFlowWrapper(global,iServer)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iServer
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(80) :: dummyString
  CHARACTER(CHRLEN) :: iFileName 
  INTEGER :: errorFlag,iScalar,iVector

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_ENS_OpenFileFlowWrapper', &
                        __FILE__)

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Opening ENSIGHT flow files...'
  END IF ! global%verbLevel

! ******************************************************************************
! Open files
! ******************************************************************************

  DO iScalar = 1,nScalars
    CALL RFLU_ENS_OpenFileScalarVector(global,scalarInfo(iScalar)%iFile, & 
                                       scalarInfo(iScalar)%iFileName)
  END DO ! iScalar                                                                              
                  
  DO iVector = 1,nVectors                                                                          
    CALL RFLU_ENS_OpenFileScalarVector(global,vectorInfo(iVector)%iFile, &
                                       vectorInfo(iVector)%iFileName)
  END DO ! iVector       
         
! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Opening ENSIGHT flow files done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ENS_OpenFileFlowWrapper






! ******************************************************************************
!
! Purpose: Open ENSIGHT geometry file.
!
! Description: None.
!
! Input:
!   global              Pointer to global type
!   iServer             Index of server
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ENS_OpenFileGeometry(global,iServer)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iServer
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(80) :: dummyString
  CHARACTER(CHRLEN) :: iFileName 
  INTEGER :: errorFlag

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_ENS_OpenFileGeometry', &
                        __FILE__)

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Opening ENSIGHT geometry file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Open file
! ******************************************************************************

!  OPEN(IF_ENS_GEOMETRY,FILE=TRIM(geoInfo%iFileName),FORM="UNFORMATTED", &
!       STATUS="UNKNOWN",IOSTAT=errorFlag)
! BBR - begin
    IF( global%solutFormat .EQ. FORMAT_BINARY )THEN
    OPEN(IF_ENS_GEOMETRY,FILE=TRIM(geoInfo%iFileName),FORM="UNFORMATTED", &
        STATUS="UNKNOWN",IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_L )THEN
    OPEN(IF_ENS_GEOMETRY,FILE=TRIM(geoInfo%iFileName),FORM="UNFORMATTED", &
       ACCESS="SEQUENTIAL",CONVERT="LITTLE_ENDIAN",STATUS="UNKNOWN", &
       IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_B )THEN
    OPEN(IF_ENS_GEOMETRY,FILE=TRIM(geoInfo%iFileName),FORM="UNFORMATTED", &
    ACCESS="SEQUENTIAL",CONVERT="BIG_ENDIAN",STATUS="UNKNOWN", &
    IOSTAT=errorFlag)
    END IF
! BBR - end 
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
  END IF ! global%error  

! ******************************************************************************
! Write header
! ******************************************************************************

  dummyString = 'Fortran Binary'
  WRITE(IF_ENS_GEOMETRY) dummyString
  
  dummyString = TRIM(global%casename)
  WRITE(IF_ENS_GEOMETRY) dummyString
  
  dummyString = TRIM(global%casename)
  WRITE(IF_ENS_GEOMETRY) dummyString
  
  dummyString = 'node id assign'
  WRITE(IF_ENS_GEOMETRY) dummyString
  
  dummyString = 'element id assign'
  WRITE(IF_ENS_GEOMETRY) dummyString        

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Opening ENSIGHT geometry file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ENS_OpenFileGeometry








! ******************************************************************************
!
! Purpose: Open ENSIGHT scalar/vector file.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   iFile               File index 
!   iFileName           File name
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ENS_OpenFileScalarVector(global,iFile,iFileName)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName
  INTEGER, INTENT(IN) :: iFile
  TYPE(t_global), POINTER :: global
  
! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(80) :: dummyString
  INTEGER :: errorFlag

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_ENS_OpenFileScalarVector', &
                        __FILE__)

! ******************************************************************************
! Open file
! ******************************************************************************

!  OPEN(iFile,FILE=TRIM(iFileName),FORM="UNFORMATTED",STATUS="UNKNOWN", & 
!       IOSTAT=errorFlag)
! BBR - begin
    IF( global%solutFormat .EQ. FORMAT_BINARY )THEN
    OPEN(iFile,FILE=TRIM(iFileName),FORM="UNFORMATTED",STATUS="UNKNOWN", & 
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
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,TRIM(iFileName))
  END IF ! global%error  

! ******************************************************************************
! Write header
! ******************************************************************************
  
  dummyString = TRIM(global%casename)
  WRITE(iFile) dummyString
         
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ENS_OpenFileScalarVector









! ******************************************************************************
!
! Purpose: Restore part number.
!
! Description: None.
!
! Input:
!   global                Pointer to global type
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ENS_RestorePartNumber(global)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global%postPartNumber = global%postPartNumberSave
         
! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE RFLU_ENS_RestorePartNumber







! ******************************************************************************
!
! Purpose: Store part number.
!
! Description: None.
!
! Input:
!   global                Pointer to global type
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ENS_StorePartNumber(global)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global%postPartNumberSave = global%postPartNumber 
         
! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE RFLU_ENS_StorePartNumber









! ******************************************************************************
!
! Purpose: Write case file.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   iServer             Index of server
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ENS_WriteFileCase(global,iServer)

  USE ModBuildFileNames, ONLY: BuildFileNameBasic

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iServer
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(80) :: dummyString
  CHARACTER(CHRLEN) :: iFileName
  INTEGER :: errorFlag,iScalar,iVector

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_ENS_WriteFileCase', &
                        __FILE__)

! ******************************************************************************
! Write 
! ******************************************************************************

  CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.case',iServer,iFileName)
  
  OPEN(IF_ENS_CASE,FILE=TRIM(iFileName),FORM="FORMATTED",STATUS="UNKNOWN", & 
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,TRIM(iFileName))
  END IF ! global%error 

  WRITE(IF_ENS_CASE,'(A)') 'FORMAT'
  WRITE(IF_ENS_CASE,'(A)') 'type: ensight gold' 
  WRITE(IF_ENS_CASE,'(A)') 'GEOMETRY'
  WRITE(IF_ENS_CASE,'(A,1X,A)') 'model:',TRIM(geoInfo%iFileName)
  WRITE(IF_ENS_CASE,'(A)') 'VARIABLE'
        
  DO iScalar = 1,nScalars
    WRITE(IF_ENS_CASE,'(A,2(1X,A))') 'scalar per element:', & 
                                     TRIM(scalarInfo(iScalar)%varName), & 
                                     TRIM(scalarInfo(iScalar)%iFileName)    
  END DO ! iScalar
  
  DO iVector = 1,nVectors
    WRITE(IF_ENS_CASE,'(A,2(1X,A))') 'vector per element:', & 
                                     TRIM(vectorInfo(iVector)%varName), & 
                                     TRIM(vectorInfo(iVector)%iFileName)   
  END DO ! iVector         

  CLOSE(IF_ENS_CASE,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ENS_WriteFileCase







! ******************************************************************************
!
! Purpose: Write solution to ENSIGHT files.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   emptyPartFlagIn     Flag indicating whether part should be empty
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ENS_WriteFlowWrapper(pRegion,emptyPartFlagIn)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  LOGICAL, OPTIONAL :: emptyPartFlagIn
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: emptyPartFlag
  CHARACTER(80) :: dummyString
  INTEGER :: errorFlag,icg,icl,ifl,iPart,iPatch,iPv,iScalar,iScalarOffs,offs
  REAL(RFREAL), DIMENSION(:), POINTER :: pScalar
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pVector  
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ENS_WriteFlowWrapper', &
                        __FILE__)

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing solution to ENSIGHT files...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                     pRegion%iRegionGlobal    
  END IF ! global%verbLevel

  pGrid => pRegion%grid

! ******************************************************************************
! Set emptyPartFlag if missing
! ******************************************************************************

  IF ( PRESENT(emptyPartFlagIn) .EQV. .FALSE. ) THEN 
    emptyPartFlag = .FALSE.
  ELSE 
    emptyPartFlag = emptyPartFlagIn 
  END IF ! PRESENT(emptyPartFlagIn)
  
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A,1X,L1)') SOLVER_NAME,'Writing empty parts:', &
                                   emptyPartFlag 
  END IF ! global%verbLevel  
  
! ******************************************************************************
! Write scalars
! ******************************************************************************

! ==============================================================================
! Mixture
! ==============================================================================

! ------------------------------------------------------------------------------
! Scalars
! ------------------------------------------------------------------------------

  IF ( emptyPartFlag .EQV. .FALSE. ) THEN
    pScalar => pRegion%mixt%cv(CV_MIXT_DENS,:)
  ELSE 
    NULLIFY(pScalar)
  END IF ! pScalar
  
  CALL RFLU_ENS_RestorePartNumber(global)
  CALL RFLU_ENS_WriteScalar(pRegion,pScalar,scalarInfo(1)%iFile,emptyPartFlag)
  
  
  IF ( emptyPartFlag .EQV. .FALSE. ) THEN
    pScalar => pRegion%mixt%cv(CV_MIXT_ENER,:)
  ELSE 
    NULLIFY(pScalar)
  END IF ! pScalar

  CALL RFLU_ENS_RestorePartNumber(global)  
  CALL RFLU_ENS_WriteScalar(pRegion,pScalar,scalarInfo(2)%iFile,emptyPartFlag)  
  
  
  IF ( emptyPartFlag .EQV. .FALSE. ) THEN
    pScalar => pRegion%mixt%dv(DV_MIXT_PRES,:)
  ELSE 
    NULLIFY(pScalar)
  END IF ! pScalar  
  
  CALL RFLU_ENS_RestorePartNumber(global)  
  CALL RFLU_ENS_WriteScalar(pRegion,pScalar,scalarInfo(3)%iFile,emptyPartFlag)  

  
  IF ( emptyPartFlag .EQV. .FALSE. ) THEN
    pScalar => pRegion%mixt%dv(DV_MIXT_TEMP,:)
  ELSE 
    NULLIFY(pScalar)
  END IF ! pScalar  
  
  CALL RFLU_ENS_RestorePartNumber(global)  
  CALL RFLU_ENS_WriteScalar(pRegion,pScalar,scalarInfo(4)%iFile,emptyPartFlag)  
  
  
  IF ( emptyPartFlag .EQV. .FALSE. ) THEN
    pScalar => pRegion%mixt%dv(DV_MIXT_SOUN,:)
  ELSE 
    NULLIFY(pScalar)
  END IF ! pScalar 
  
  CALL RFLU_ENS_RestorePartNumber(global) 
  CALL RFLU_ENS_WriteScalar(pRegion,pScalar,scalarInfo(5)%iFile,emptyPartFlag)  
  
  iScalarOffs = 5
 
! ------------------------------------------------------------------------------
! Plotting variables 
! ------------------------------------------------------------------------------

  IF ( pRegion%plot%nPv > 0 ) THEN
    iScalar = 0
    
    DO iPv = 1,pRegion%plot%nPv
      iScalar = iScalar+1

      IF ( emptyPartFlag .EQV. .FALSE. ) THEN
        pScalar => pRegion%plot%pv(iPv,:)
      ELSE 
        NULLIFY(pScalar)
      END IF ! pScalar 
    
      CALL RFLU_ENS_RestorePartNumber(global) 
      CALL RFLU_ENS_WriteScalar(pRegion,pScalar, &
                                scalarInfo(iScalarOffs+iScalar)%iFile, & 
                                emptyPartFlag)  
    END DO ! iPv

    iScalarOffs = iScalarOffs + pRegion%plot%nPv
  END IF ! nPv 
  
! ------------------------------------------------------------------------------
! Vectors
! ------------------------------------------------------------------------------
  
  IF ( emptyPartFlag .EQV. .FALSE. ) THEN
    pVector => pRegion%mixt%cv(CV_MIXT_XMOM:CV_MIXT_ZMOM,:)
  ELSE 
    NULLIFY(pVector)
  END IF ! pScalar 

  CALL RFLU_ENS_RestorePartNumber(global)
  CALL RFLU_ENS_WriteVector(pRegion,pVector,vectorInfo(1)%iFile,emptyPartFlag)          

! ==============================================================================
! Physical modules
! ==============================================================================

#ifdef SPEC
! ------------------------------------------------------------------------------
! Species
! ------------------------------------------------------------------------------

  IF ( global%specUsed .EQV. .TRUE. ) THEN 
    DO iScalar = 1,pRegion%specInput%nSpecies
      IF ( emptyPartFlag .EQV. .FALSE. ) THEN
        pScalar => pRegion%spec%cv(iScalar,:)
      ELSE 
        NULLIFY(pScalar)
      END IF ! pScalar
      
      CALL RFLU_ENS_RestorePartNumber(global)  
      CALL RFLU_ENS_WriteScalar(pRegion,pScalar, &
                                scalarInfo(iScalar+iScalarOffs)%iFile, &
                                emptyPartFlag)        
    END DO ! iScalar
  END IF ! global%specUsed
#endif

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Writing solution to ENSIGHT files done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ENS_WriteFlowWrapper







! ******************************************************************************
!
! Purpose: Write grid to ENSIGHT geometry file.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region
!   emptyPartFlagIn     Flag indicating whether part should be empty 
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ENS_WriteGridWrapper(pRegion,emptyPartFlagIn)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  LOGICAL, OPTIONAL :: emptyPartFlagIn
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: emptyPartFlag
  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ENS_WriteGrid', &
                        __FILE__)

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Writing grid to ENSIGHT geometry file...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                     pRegion%iRegionGlobal  
  END IF ! global%verbLevel

! ******************************************************************************
! Set emptyPartFlag if missing
! ******************************************************************************

  IF ( PRESENT(emptyPartFlagIn) .EQV. .FALSE. ) THEN 
    emptyPartFlag = .FALSE.
  ELSE 
    emptyPartFlag = emptyPartFlagIn 
  END IF ! PRESENT(emptyPartFlagIn)
  
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A,1X,L1)') SOLVER_NAME,'Writing empty parts:', &
                                   emptyPartFlag 
  END IF ! global%verbLevel   
  
! ******************************************************************************
! Write geometry
! ******************************************************************************

  CALL RFLU_ENS_WriteGrid(pRegion,emptyPartFlag)

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Writing grid to ENSIGHT geometry file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ENS_WriteGridWrapper







END MODULE RFLU_ModENSIGHT

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModENSIGHT.F90,v $
! Revision 1.2  2015/07/23 23:11:19  brollin
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
! Revision 1.3  2008/12/06 08:43:58  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:12  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:58:09  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2007/03/19 21:42:36  haselbac
! Adapted to changes related to plotting variables
!
! Revision 1.3  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.2  2006/01/24 21:20:38  mparmar
! Added plotting variables
!
! Revision 1.1  2005/10/05 20:23:34  haselbac
! Initial revision
!
! ******************************************************************************

