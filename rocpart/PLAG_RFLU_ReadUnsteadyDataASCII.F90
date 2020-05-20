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
! Purpose: Read flow file for particles in ASCII ROCFLU format.
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
!
! $Id: PLAG_RFLU_ReadUnsteadyDataASCII.F90,v 1.2 2015/08/12 20:18:19 rahul Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_ReadUnsteadyDataASCII(pRegion)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region 
  USE ModPartLag, ONLY: t_plag,t_tile_plag     
  USE ModMPI

  USE PLAG_ModParameters

  USE ModBuildFileNames, ONLY: BuildFileNameUnsteady

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
! Local variables
! ==============================================================================

  LOGICAL :: timeDataRead
  CHARACTER(CHRLEN) :: errorString,iFileName,sectionString,RCSIdentString, & 
                       timeString1,timeString2
  INTEGER :: errorFlag,i,iCont,iFile,ifl,iMass,iPatch,iVars,j,loopCounter, &
             nCont,nData,nDataAllocate,nDataDummy,nDataExpected,nDataMin, &
             nPcls,nPclsAllocate,nPclsExpected,nVars,nVarsExpected,precActual, &
             precExpected,rangeActual,rangeExpected
  INTEGER, DIMENSION(:,:), POINTER :: pAiv
  REAL(RFREAL) :: currentTime
  REAL(RFREAL), DIMENSION(:,:), POINTER :: dummyVal
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pDudtMixt,pDudtPlag
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_plag), POINTER :: pPlag
  TYPE(t_tile_plag), POINTER :: pTilePlag  
  
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion  
  
! ******************************************************************************
! Start, open file
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_ReadUnsteadyDataASCII.F90,v $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_ReadUnsteadyDataASCII',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII particle unsteady data file...'
  END IF ! global%verbLevel

  CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.plag_unsda', & 
                             pRegion%iRegionGlobal,global%currentTime, & 
                             iFileName)

  iFile = IF_SOLUT
  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
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

  READ(iFile,'(A)') sectionString
  IF ( TRIM(sectionString) /= '# ROCFLU particle unsteady data file' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM  
  
! -----------------------------------------------------------------------------
! Precision and range
! -----------------------------------------------------------------------------
    
  READ(iFile,'(A)') sectionString
  IF ( TRIM(sectionString) /= '# Precision and range' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM  
  
  precExpected  = PRECISION(1.0_RFREAL)
  rangeExpected = RANGE(1.0_RFREAL)
  
  READ(iFile,'(2(I16))') precActual,rangeActual
  IF ( precActual < precExpected .OR. rangeActual < rangeExpected ) THEN 
    CALL ErrorStop(global,ERR_PREC_RANGE,__LINE__)
  END IF ! precActual
  
! -----------------------------------------------------------------------------
! Initial residual and physical time
! -----------------------------------------------------------------------------
    
  READ(iFile,'(A)') sectionString
  IF ( TRIM(sectionString) /= '# Physical time' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
  END IF ! TRIM    
   
  READ(iFile,'(E23.16)') currentTime 

! ==============================================================================
! Dimensions
! ==============================================================================
  
  pGrid => pRegion%grid  
  pPlag => pRegion%plag
   
  pDudtMixt => pPlag%dudtMixt
  pDudtPlag => pPlag%dudtPlag
      
  nCont = pRegion%plagInput%nCont 
   
  nVarsExpected = 7
  nDataExpected = pRegion%plagInput%nUnsteadyData
  nPclsExpected = pPlag%nPcls    
  
  nDataMin   = nDataExpected
  nDataDummy = 0

  READ(iFile,'(A)') sectionString
  IF ( TRIM(sectionString) /= '# Dimensions' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM

  READ(iFile,'(3(I16))') nPcls,nVars,nData
  
  IF ( nPcls /= nPclsExpected ) THEN 
    WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nPcls, & 
                                              'but expected:',nPclsExpected
    CALL ErrorStop(global,ERR_PLAG_INVALID_NPCLS,__LINE__,errorString)
  END IF ! nCellsExpected     
  
  IF ( nVars /= nVarsExpected ) THEN 
    WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nVars, & 
                                              'but expected:',nVarsExpected  
    CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
  END IF ! nVarsExpected    

  IF ( nData /= nDataExpected ) THEN 
    WRITE(STDOUT,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nData, & 
                                         'but expected:',nDataExpected  
    nDataMin = MIN(nData,nDataExpected)
    nDataDummy = MAX(0,(nData-nDataExpected))
    WRITE(STDOUT,'(A,1X,I6,1X,A)') 'Using only ',nDataMin,  ' data.'
  END IF ! nVarsExpected    

  nPclsAllocate = MAX(1,nPcls)
  nDataAllocate = MAX(1,nData)

  ALLOCATE(dummyVal(nDataAllocate,nPclsAllocate),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'dummyVal')
  END IF ! global%error

  dummyVal = 0.0_RFREAL

  timeDataRead = .FALSE.

! ==============================================================================
! Rest of file
! ==============================================================================

  iCont       = 0
  iVars       = 0
  loopCounter = 0

  DO ! set up infinite loop
    loopCounter = loopCounter + 1
  
    READ(iFile,'(A)') sectionString

    SELECT CASE ( TRIM(sectionString) ) 

      CASE ( '# Relevant time history' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Relevant time history...'
        END IF ! global%verbLevel    
     
        timeDataRead = .TRUE. 
        iVars = iVars + 1
        READ(iFile,'(5(E23.16))') (dummyVal(i,1),i=1,nData)

        pRegion%plagInput%nTimeBH = nDataMin
        DO i=1,nDataMin
          pPlag%timeBH(i) = dummyVal(i,1)
        END DO ! i
        pPlag%timeBH(1) = 0.0_RFREAL ! Index 1 is for current time

      CASE ( '# Particle x-velocity substantial derivative' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Particle x-vel sd...'
        END IF ! global%verbLevel    
      
        iVars = iVars + 1
        READ(iFile,'(5(E23.16))') ((dummyVal(i,j),i=1,nData),j=1,nPcls)

        DO j=1,nPcls
          DO i=1,nDataMin
            pDudtPlag(XCOORD,i,j) = dummyVal(i,j)
          END DO ! i
        END DO ! j

      CASE ( '# Particle y-velocity substantial derivative' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Particle y-vel sd...'
        END IF ! global%verbLevel    
      
        iVars = iVars + 1
        READ(iFile,'(5(E23.16))') ((dummyVal(i,j),i=1,nData),j=1,nPcls)

        DO j=1,nPcls
          DO i=1,nDataMin
            pDudtPlag(YCOORD,i,j) = dummyVal(i,j)
          END DO ! i
        END DO ! j

      CASE ( '# Particle z-velocity substantial derivative' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Particle z-vel sd...'
        END IF ! global%verbLevel    
      
        iVars = iVars + 1
        READ(iFile,'(5(E23.16))') ((dummyVal(i,j),i=1,nData),j=1,nPcls)

        DO j=1,nPcls
          DO i=1,nDataMin
            pDudtPlag(ZCOORD,i,j) = dummyVal(i,j)
          END DO ! i
        END DO ! j

      CASE ( '# Fluid x-velocity substantial derivative' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Fluid x-vel sd...'
        END IF ! global%verbLevel    
      
        iVars = iVars + 1
        READ(iFile,'(5(E23.16))') ((dummyVal(i,j),i=1,nData),j=1,nPcls)

        DO j=1,nPcls
          DO i=1,nDataMin
            pDudtMixt(XCOORD,i,j) = dummyVal(i,j)
          END DO ! i
        END DO ! j

      CASE ( '# Fluid y-velocity substantial derivative' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Fluid y-vel sd...'
        END IF ! global%verbLevel    
      
        iVars = iVars + 1
        READ(iFile,'(5(E23.16))') ((dummyVal(i,j),i=1,nData),j=1,nPcls)

        DO j=1,nPcls
          DO i=1,nDataMin
            pDudtMixt(YCOORD,i,j) = dummyVal(i,j)
          END DO ! i
        END DO ! j

      CASE ( '# Fluid z-velocity substantial derivative' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Fluid z-vel sd...'
        END IF ! global%verbLevel    
      
        iVars = iVars + 1
        READ(iFile,'(5(E23.16))') ((dummyVal(i,j),i=1,nData),j=1,nPcls)

        DO j=1,nPcls
          DO i=1,nDataMin
            pDudtMixt(ZCOORD,i,j) = dummyVal(i,j)
          END DO ! i
        END DO ! j

! ------------------------------------------------------------------------------
!     End marker
! ------------------------------------------------------------------------------ 
      
      CASE ( '# End' ) 
        IF ( global%myProcid == MASTERPROC .AND. & 
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
        END IF ! global%verbLevel           
      
        EXIT
      
! ------------------------------------------------------------------------------
!     Invalid section string
! ------------------------------------------------------------------------------ 
      
      CASE DEFAULT
        IF ( global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
        END IF ! verbosityLevel           
      
        CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)        
             
    END SELECT ! TRIM
  
! ------------------------------------------------------------------------------
!   Guard against infinite loop - might be unnecessary because of read errors?
! ------------------------------------------------------------------------------  
  
    IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
      CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
    END IF ! loopCounter
  
  END DO ! <empty>

! ==============================================================================
! Check and information about number of variables read
! ==============================================================================

  IF ( iVars /= nVars ) THEN 
    CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
  END IF ! iVar

! - Time history data is not available then construct it          
!  WRITE(*,*) "timeDataRead=",timeDataRead
!  WRITE(*,*) "iVars=",iVars
!  IF ( (timeDataRead .EQV. .FALSE.) .AND. (iVars==nVars-1) ) THEN
!    pRegion%plagInput%nTimeBH = nDataMin
!    DO i=1,nDataMin
!      pPlag%timeBH(i) = (i-1)*global%dtMin
!    END DO ! i
!    pPlag%timeBH(1) = 0.0_RFREAL ! Index 1 is for current time
!  ELSE
!    CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
!  END IF ! timeDataRead

  DEALLOCATE(dummyVal,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'dummyVal')
  END IF ! global%error

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(iFile,IOSTAT=errorFlag)
  global%error = errorFlag      
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
  END IF ! global%error
   
  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII particle unsteady data file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)
    
! ******************************************************************************
! End
! ******************************************************************************
 
END SUBROUTINE PLAG_RFLU_ReadUnsteadyDataASCII


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_ReadUnsteadyDataASCII.F90,v $
! Revision 1.2  2015/08/12 20:18:19  rahul
! Changed integer format from I8 to I16 in dimensions and precision & range
! in PLAG_RFLU_ReadUnsteadyDataASCII.F90
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
! ******************************************************************************

