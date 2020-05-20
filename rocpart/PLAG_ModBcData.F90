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
! Purpose: Collection of routines to initialize, check and read 
!          boundary condition input file.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: PLAG_ModBcData.F90,v 1.3 2016/02/24 06:10:38 rahul Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
! ******************************************************************************

MODULE PLAG_ModBcData

  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch, t_bcvalues_plag
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  USE PLAG_ModParameters

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: PLAG_CheckBcData, &
            PLAG_CreateBcData, &
	    PLAG_DestroyBcData, &
            PLAG_InitBcData, &
            PLAG_ReadBcInputFile, &
            PLAG_PrintBcData

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

  CHARACTER(CHRLEN) :: RCSIdentString = &
    '$RCSfile: PLAG_ModBcData.F90,v $ $Revision: 1.3 $'

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS






!******************************************************************************
!
! Purpose: Check user boundary conditions parameters for 
!          Lagrangians particles to default values.
!
! Description: none.
!
! Input:
!   pRegion     Pointer to region
!
! Notes: none.
!
!******************************************************************************

  SUBROUTINE PLAG_CheckBcData( pRegion )

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: bcType,iPatch, iReg
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_bcvalues_plag), POINTER :: pPlagBc
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_CheckBcData',__FILE__)

! ******************************************************************************
!   Set variables
! ******************************************************************************

    iReg = pRegion%iRegionGlobal

! *****************************************************************************
!   Initialize variables on patch inflow-based data
! *****************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch  => pRegion%patches(iPatch)
      pPlagBc => pPatch%plag
      bcType = pPatch%bcType

      SELECT CASE(bcType)
        CASE( BC_INFLOW,         &
              BC_INFLOW_VELTEMP, &
              BC_INJECTION       )
	  IF ( pPlagBc%inflowDiamDist < PLAG_BC_INFLOW_LOGNORM .OR. &
               pPlagBc%inflowDiamDist > PLAG_BC_INFLOW_PDF          ) THEN
            WRITE(STDOUT,1050) iReg,iPatch,pPlagBc%inflowDiamDist
            CALL ErrorStop( global,ERR_PLAG_INFLOWDIAMDIST,__LINE__ )
          END IF ! inflowDiamDist

	  IF ( pPlagBc%inflowDiamDist == PLAG_BC_INFLOW_LOGSKWD .AND. &
               pPlagBc%inflowDiamMean/&
               pPlagBc%inflowDiamMax > 0.8_RFREAL ) THEN
	    WRITE(STDOUT,1060) iReg,iPatch, &
	                       pPlagBc%inflowDiamMean/&
                               pPlagBc%inflowDiamMax
	    CALL ErrorStop( global,ERR_PLAG_INFLOWDIAM,__LINE__ )
	  END IF ! inflowDiamDist

	  IF (  pPlagBc%inflowDiamDist == PLAG_BC_INFLOW_LOGSKWD .AND. &
                pPlagBc%inflowDiamMin/&
                pPlagBc%inflowDiamMean > 0.8_RFREAL ) THEN
	    WRITE(STDOUT,1070) iReg,iPatch, &
	                       pPlagBc%inflowDiamMin/&
			       pPlagBc%inflowDiamMean
	    CALL ErrorStop( global,ERR_PLAG_INFLOWDIAM,__LINE__ )
	  END IF ! inflowDiamDist
! Rahul-temp
! changed .OR. to .AND. to accommodate new inflowModel i.e PLAG_BC_NOINFLOW=0
! This fixed an error in PLAG_ModInflow. I don't know if this is the 
! right way to do it.
!          IF ( pPlagBc%inflowModel < PLAG_BC_INFLOW_MODEL1 .AND. &
!               pPlagBc%inflowModel > PLAG_BC_INFLOW_CRE  ) THEN
          IF ( pPlagBc%inflowModel < PLAG_BC_NOINFLOW .OR. &
               pPlagBc%inflowModel > PLAG_BC_INFLOW_CRE  ) THEN
            WRITE(STDOUT,1080) iReg,iPatch,pPlagBc%inflowModel
            CALL ErrorStop( global,ERR_PLAG_EJECMODEL,__LINE__ )
          END IF ! inflowModel
! End-temp
  
      END SELECT  ! bcType  
    END DO !  iPatch

! *****************************************************************************
!   Formats
! *****************************************************************************

1050 FORMAT('Region ',I5,3X,'Patch ',I5,', diameter Distribution ',I1)
1060 FORMAT('Region ',I5,3X,'Patch ',I5,', inflowDiamMean/inflowDiamMax  = ',F12.5)
1070 FORMAT('Region ',I5,3X,'Patch ',I5,', inflowDiamMin /inflowDiamMean = ',F12.5)
1080 FORMAT('Region ',I5,3X,'Patch ',I5,', inflowModel ',I1)

! *****************************************************************************
!   End
! *****************************************************************************

  CALL DeregisterFunction( global )

  END SUBROUTINE PLAG_CheckBcData

  
  
  
  
  
  
  
  
!******************************************************************************
!
! Purpose: Create arrays specific for boundary conditions.
!
! Description: none.
!
! Input:
!   pRegion     Pointer to region
!
! Notes: none.
!
!******************************************************************************

  SUBROUTINE PLAG_CreateBcData( pRegion )

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: bcType,errorFlag,iPatch,nCont
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_bcvalues_plag), POINTER :: pPlagBc
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_CreateBcData',__FILE__)

! *****************************************************************************
!   Set values
! *****************************************************************************
    
    nCont   = pRegion%plagInput%nCont 

! *****************************************************************************
!   Allocate memory
! *****************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch  => pRegion%patches(iPatch)
      pPlagBc => pPatch%plag
      bcType = pPatch%bcType

      SELECT CASE(bcType)
        CASE( BC_INFLOW,         &
              BC_INFLOW_VELTEMP, &
              BC_INJECTION       )
	  ALLOCATE(pPlagBc%inflowMassFluxRatio(nCont),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'inflowMassFluxRatio')
          END IF ! global%error

	  ALLOCATE(pPlagBc%inflowTemp(nCont),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'inflowTemp')
          END IF ! global%error	  

        CASE DEFAULT 
	  NULLIFY(pPlagBc%inflowMassFluxRatio)
	  NULLIFY(pPlagBc%inflowTemp)

      END SELECT  ! bcType  
    END DO !  iPatch

! *****************************************************************************
!   End
! *****************************************************************************
  
  CALL DeregisterFunction( global )

  END SUBROUTINE PLAG_CreateBcData
  
  
  
  
  
  
  
  
!******************************************************************************
!
! Purpose: Destroy arrays specific for boundary conditions.
!
! Description: none.
!
! Input:
!   pRegion     Pointer to region
!
! Notes: none.
!
!******************************************************************************

  SUBROUTINE PLAG_DestroyBcData( pRegion )

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: bcType,errorFlag,iPatch,nCont
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_bcvalues_plag), POINTER :: pPlagBc
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_DestroyBcData',__FILE__)

! *****************************************************************************
!   Deallocate memory
! *****************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch  => pRegion%patches(iPatch)
      pPlagBc => pPatch%plag
      bcType = pPatch%bcType

      SELECT CASE(bcType)
        CASE( BC_INFLOW,         &
              BC_INFLOW_VELTEMP, &
              BC_INJECTION       )
	  DEALLOCATE(pPlagBc%inflowMassFluxRatio,STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'inflowMassFluxRatio')
          END IF ! global%error

	  DEALLOCATE(pPlagBc%inflowTemp,STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'inflowTemp')
          END IF ! global%error	  
      END SELECT  ! bcType  
    END DO !  iPatch

! *****************************************************************************
!   End
! *****************************************************************************
  
  CALL DeregisterFunction( global )

  END SUBROUTINE PLAG_DestroyBcData








!******************************************************************************
!
! Purpose: Initialize user boundary conditions parameters for 
!          Lagrangians particles to default values.
!
! Description: none.
!
! Input:
!   pRegion     Pointer to region
!
! Notes: none.
!
!******************************************************************************

  SUBROUTINE PLAG_InitBcData( pRegion )

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: bcType,iPatch
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_bcvalues_plag), POINTER :: pPlagBc
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_InitBcData',__FILE__)

! *****************************************************************************
!   Initialize variables on patch inflow-based data
! *****************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch  => pRegion%patches(iPatch)
      pPlagBc => pPatch%plag
      bcType = pPatch%bcType

      SELECT CASE(bcType)
        CASE( BC_INFLOW,         &
              BC_INFLOW_VELTEMP, &
              BC_INJECTION       )
	  pPlagBc%inflowModel     = PLAG_BC_INFLOW_MODEL1
          pPlagBc%inflowDiamDist  = PLAG_BC_INFLOW_LOGNORM
          pPlagBc%inflowVelRatio  = 0.0_RFREAL
          pPlagBc%inflowSpLoad    = 1.0_RFREAL
          pPlagBc%inflowDiamMean  = 10.0E-06_RFREAL
          pPlagBc%inflowDiamMin   = 1.0E-06_RFREAL
          pPlagBc%inflowDiamMax   = 100.0E-06_RFREAL
          pPlagBc%inflowStdDev    = 0.0_RFREAL
          pPlagBc%inflowBeta      = 1.0_RFREAL
      END SELECT  ! bcType  
    END DO !  iPatch

! *****************************************************************************
!   End
! *****************************************************************************
  
  CALL DeregisterFunction( global )

  END SUBROUTINE PLAG_InitBcData









! ******************************************************************************
!
! Purpose: Read in user input related to boundary conditions for Rocpart 
!   (done on all processors).
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

  SUBROUTINE PLAG_ReadBcInputFile(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNamePlain
  
    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName
    CHARACTER(256) :: line
    INTEGER :: errorFlag,iPatch,loopCounter
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_ReadBcInputFile',__FILE__)

! ******************************************************************************
!   Open file
! ******************************************************************************

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.bc',iFileName)
  
    OPEN(IF_INPUT,FILE=iFileName,FORM='FORMATTED',STATUS='OLD',IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error

! ******************************************************************************
!   Read file looking for keywords
! ******************************************************************************

    loopCounter = 0

    KeyWordLoop: DO
      READ(IF_INPUT,'(A256)',IOSTAT=errorFlag) line

      IF ( errorFlag > 0 ) THEN ! Error occurred
	CALL ErrorStop(global,ERR_FILE_READ,__LINE__,'File: '//TRIM(iFileName))
      ELSE IF ( errorFlag < 0 ) THEN ! Encountered end of file
	EXIT KeyWordLoop
      END IF ! errorFlag    

      SELECT CASE( TRIM(line) )    
  ! TEMPORARY - Keep this for backward compatibility
	CASE ('# BC_INFLOW')
          CALL PLAG_ReadBcInflowSection(pRegion)
  ! END TEMPORARY    
	CASE ('# BC_INFLOW_VELTEMP')
          CALL PLAG_ReadBcInflowSection(pRegion)             
	CASE ('# BC_INJECT')
          CALL PLAG_ReadBcInflowSection(pRegion)
	CASE ('# END')
          EXIT
      END SELECT ! TRIM(line)

      loopCounter = loopCounter + 1 ! Prevent infinite loop
      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
	CALL ErrorStop(global,ERR_INFINITE_LOOP ,__LINE__)
      END IF ! loopCounter
    END DO KeyWordLoop

! ******************************************************************************
!   Close file
! ******************************************************************************

    CLOSE(IF_INPUT,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_ReadBcInputFile







!******************************************************************************
!
! Purpose: Read in user input related to ejection boundary condition for 
!            Rocpart.
!
! Description: None.
!
! Input: 
!   pRegion    Region pointer
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************

  SUBROUTINE PLAG_ReadBcInflowSection(pRegion)
 
    USE ModInterfaces, ONLY: MakeNumberedKeys,ReadPatchSection 

    IMPLICIT NONE

! *****************************************************************************
!   Definitions and declarations
! *****************************************************************************

! =============================================================================
!   Arguments
! =============================================================================

    TYPE(t_region), POINTER :: pRegion

! =============================================================================
!   Locals
! =============================================================================

    INTEGER, PARAMETER :: NVALS_MAX = 9

    CHARACTER(CHRLEN) :: bcName
    CHARACTER(15), DIMENSION(:), ALLOCATABLE :: keys
    CHARACTER(256) :: fileName
    LOGICAL, DIMENSION(:), ALLOCATABLE :: defined
    INTEGER :: checkSum,checkSumCont,distrib,errorFlag,iCont,iKey,iPatch,iPatchBeg, &
               iPatchEnd,ival,nBegin,nBegin1,nBegin2,nCont,nEnd,nEnd1,nEnd2,&
               nKeys,nKeysShift1,nKeysShift2
    REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: vals
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_bcvalues_plag), POINTER :: pPlagBc
    TYPE(t_global), POINTER :: global

! *****************************************************************************
!   Start
! *****************************************************************************
 
    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_ReadBcInflowSection',__FILE__)

! *****************************************************************************
!   Set values
! *****************************************************************************
    
    nCont   = pRegion%plagInput%nCont 
    
    nKeys       = NVALS_MAX   +2 *nCont
    nKeysShift1 = NVALS_MAX   +1
    nKeysShift2 = nKeysShift1 +nCont
    
    nBegin1 = nKeysShift1 
    nEnd1   = nBegin1 +nCont -1
    nBegin2 = nKeysShift2 
    nEnd2   = nBegin2 +nCont -1

    nBegin = 1 
    nEnd   = nCont

! *****************************************************************************
!   Allocate memory 
! *****************************************************************************

    ALLOCATE(keys(nKeys),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'keys')
    END IF ! global%error

    ALLOCATE(vals(nKeys),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vals')
    END IF ! global%error  

    ALLOCATE(defined(nKeys),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'defined')
    END IF ! global%error 

! *****************************************************************************
!   Generate keys 
! *****************************************************************************

    keys(1) = 'PLAG_MODEL'
    keys(2) = 'PLAG_VELRATIO'
    keys(3) = 'PLAG_SPLOAD'
    keys(4) = 'PLAG_BETA'
    keys(5) = 'PLAG_DIAMDIST'
    keys(6) = 'PLAG_DIAMMEAN'
    keys(7) = 'PLAG_DIAMMIN'
    keys(8) = 'PLAG_DIAMMAX'
    keys(9) = 'PLAG_STDDEV'

    CALL MakeNumberedKeys(keys,nKeysShift1,'PLAG_MASSRATIO',nBegin,nEnd,1)
    CALL MakeNumberedKeys(keys,nKeysShift2,'PLAG_TEMP',nBegin,nEnd,1)

    defined(:) = .FALSE.

! *****************************************************************************
!   Read section
! *****************************************************************************

    CALL ReadPatchSection(global,IF_INPUT,nKeys,keys,vals,iPatchBeg,iPatchEnd, &
                          distrib,fileName,bcName,defined)

! Saptarshi - begin
   defined(:) = .TRUE. !BBR - temporary - WHY??
! Saptarshi - end
! Rahul - Resolved this issue. Use PLAG_MASSRATIO1 and PLAG_TEMP1 instead of 
!         PLAG_MASSRATIO and PLAG_TEMP respectively in the .bc file. 
!         No need to comment this.
! *****************************************************************************
!   Check if specified number of patches exceeds available ones
! *****************************************************************************

    IF ( iPatchEnd > global%nPatches ) THEN 
      CALL ErrorStop(global,ERR_PATCH_RANGE,__LINE__)
    END IF ! iPatchEnd

! *****************************************************************************
!   Copy values
! *****************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch  => pRegion%patches(iPatch)
      pPlagBc => pPatch%plag

! =============================================================================
!     Check whether this global patch exists in this region
! =============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. & 
           pPatch%iPatchGlobal <= iPatchEnd       ) THEN
        IF (defined(1)) THEN
          pPlagBc%inflowModel     = NINT(ABS(vals(1)))
          IF ( vals(1) > 0.9 .AND. vals(3) < 1.1 ) &
            pPlagBc%inflowModel  = PLAG_BC_INFLOW_MODEL1
          IF ( vals(1) > 1.9 .AND. vals(3) < 2.1 ) &
            pPlagBc%inflowModel  = PLAG_BC_INFLOW_CRE
! Rahul - Added PLAG_BC_NOINFLOW to fix the injection BC issue in shktb case
          IF ( vals(1) > -0.9 .AND. vals(3) < 0.1 ) &
            pPlagBc%inflowModel  = PLAG_BC_NOINFLOW
! Rahul - end
        ENDIF ! defined
       
        IF (defined(2)) &
          pPlagBc%inflowVelRatio = ABS(vals(2))
	
	IF (defined(3)) &
	  pPlagBc%inflowSpload   = ABS(vals(3))
        
	IF (defined(4)) &
	  pPlagBc%inflowBeta     = ABS(vals(4))
        
	IF (defined(5)) THEN
	  pPlagBc%inflowDiamDist = ABS(vals(5))
          IF ( vals(5) > 0.9 .AND. vals(5) < 1.1 ) &
            pPlagBc%inflowDiamDist  = PLAG_BC_INFLOW_LOGNORM   
          IF ( vals(5) > 1.9 .AND. vals(5) < 2.1 ) &
            pPlagBc%inflowDiamDist  = PLAG_BC_INFLOW_LOGSKWD   
          IF ( vals(5) > 2.9 .AND. vals(5) < 3.1 ) THEN
            pPlagBc%inflowDiamDist  = PLAG_BC_INFLOW_PDF   
            CALL PLAG_RFLU_ReadPdfFromFile( pRegion)  
          END IF ! vals(5)
        END IF ! defined(5)
	
	IF (defined(6)) & 
	  pPlagBc%inflowDiamMean = ABS(vals(6))
	
	IF (defined(7)) & 
	  pPlagBc%inflowDiamMin  = ABS(vals(7))
	
	IF (defined(8)) & 
	  pPlagBc%inflowDiamMax  = ABS(vals(8))
	
	IF (defined(9)) & 
	  pPlagBc%inflowStdDev   = ABS(vals(9))

! -----------------------------------------------------------------------------
!       Check whether all values defined
!       Note: Not really required but strongly suggested
!             since initial values are set.
! -----------------------------------------------------------------------------
        
	checkSum = 0

        DO iKey = 1, NVALS_MAX
          IF ( defined(iKey) .EQV. .TRUE. ) THEN 
            checkSum = checkSum + 1
          END IF ! defined
        END DO ! iKey
        
        IF ( checkSum /= NVALS_MAX ) THEN 
          CALL ErrorStop(global,ERR_BCVAL_MISSING,__LINE__)
        END IF ! checkSum 

! -----------------------------------------------------------------------------
!       Copy vaues pertinent to number of constituents
! -----------------------------------------------------------------------------
                
	ival= NVALS_MAX
        DO iCont = 1, nCont
          ival= ival+1
	  IF (defined(ival) ) &
            pPlagBc%inflowMassFluxRatio(iCont) = ABS(vals(ival))
        END DO ! iCont	

        DO iCont = 1, nCont
          ival= ival+1
	  IF (defined(ival) ) &
            pPlagBc%inflowTemp(iCont) = ABS(vals(ival))
        END DO ! iCont 	

! -----------------------------------------------------------------------------
!       Check whether all values defined for constituents
! -----------------------------------------------------------------------------
        
        checkSumCont = 0

        DO iKey = nBegin1, nEnd2
          IF ( defined(iKey) .EQV. .TRUE. ) THEN 
            checkSumCont = checkSumCont + 1
          END IF ! defined
        END DO ! iKey
        
        IF ( checkSumCont /= 2*nCont ) THEN 
          CALL ErrorStop(global,ERR_BCVAL_MISSING,__LINE__)
        END IF ! checkSumCont

      END IF ! pPatch%iPatchGlobal
    END DO ! iPatch

! *****************************************************************************
!   Deallocate memory 
! *****************************************************************************

    DEALLOCATE(keys,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'keys')
    END IF ! global%error

    DEALLOCATE(vals,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vals')
    END IF ! global%error  

    DEALLOCATE(defined,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'defined')
    END IF ! global%error 
  
! *****************************************************************************
!   End
! *****************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_ReadBcInflowSection








!******************************************************************************
!
! Purpose: read in user defined porbability density function (PDF) of
!          the mass injection process.
!
! Description: none
!
! Input: user input file.
!
! Output: structure with information relative to the imposed pdf
!
! Notes:  none
!
!******************************************************************************

  SUBROUTINE PLAG_RFLU_ReadPdfFromFile(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
  

    CHARACTER(CHRLEN)       :: iFileName
    INTEGER                 :: bcType,errorFlag,iPatch,k,nbins,nrow,tmpint
    REAL(RFREAL)            :: tmpscal
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: tmpvec
    
    TYPE(t_patch), POINTER   :: pPatch
    TYPE(t_bcvalues_plag), POINTER :: pPlagBc
    TYPE(t_global), POINTER  :: global
    
! *****************************************************************************
!   Start
! *****************************************************************************
 
    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_RFLU_ReaPdfFromFle',__FILE__)

! *****************************************************************************
!   Open File
! *****************************************************************************

    WRITE(iFileName,'(A,A,A)') &
    TRIM(global%inDir),TRIM(global%casename),'.plag_injcpdf'
    
    OPEN(IF_PLAG_INFLOWPDF,FILE=iFileName,FORM='formatted',STATUS='old',&
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) & 
      CALL ErrorStop( global, ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName) )
    
    READ(IF_PLAG_INFLOWPDF,*,err=10,end=10) nbins

! *****************************************************************************
!   Allocate memory
! *****************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch  => pRegion%patches(iPatch)
      pPlagBc => pPatch%plag
      bcType = pPatch%bcType

      SELECT CASE(bcType)
        CASE( BC_INFLOW,         &
              BC_INFLOW_VELTEMP, &
              BC_INJECTION       )
          ALLOCATE(pPlagBc%PDFBc%pdfvalues(nbins+1,3), STAT=errorFlag )
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) &
	    CALL ErrorStop( global, ERR_ALLOCATE,__LINE__, 'pdfvalues' )  

        CASE DEFAULT 
	  NULLIFY(pPlagBc%PDFBc%pdfvalues)
      END SELECT ! bcType
    END DO ! iPatch         
 
    ALLOCATE(tmpvec(nbins+1,5), STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) & 
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__, 'tmpvec' )

! *****************************************************************************
!   Read in the discrete PDF
! *****************************************************************************

    DO k = 1,nbins
      READ(IF_PLAG_INFLOWPDF,*,err=10,end=10) nrow,tmpvec(k,1:2)
    END DO ! k

    tmpvec(1,3) = (tmpvec(2,1)-tmpvec(1,1))        
    DO k = 2,nbins
      tmpvec(k,3) = 2_RFREAL*(tmpvec(k,1)-tmpvec(k-1,1)) - tmpvec(k-1,3)
    END DO ! k

! *****************************************************************************
!  Evaluate the two following cumulative sums 
! *****************************************************************************
  
   tmpvec(1,4:5) = 0_RFREAL
   
   DO k = 1,nbins
     tmpvec(k+1,4) = tmpvec(k,4) + tmpvec(k,2)*tmpvec(k,3)
     tmpvec(k+1,5) = tmpvec(k,5) + tmpvec(k,3)
   END DO ! k
   
   tmpvec(:,4) = tmpvec(:,4) / tmpvec(nbins+1,4)
   tmpvec(:,5) = tmpvec(:,5) + tmpvec(1,1) - tmpvec(1,3)/2_RFREAL

! *****************************************************************************
!  Assign the tmp vector entries to the region pointer
! *****************************************************************************
    
    DO iPatch = 1,pRegion%grid%nPatches
      pPatch  => pRegion%patches(iPatch)
      pPlagBc => pPatch%plag
      bcType = pPatch%bcType

      SELECT CASE(bcType)
        CASE( BC_INFLOW,         &
              BC_INFLOW_VELTEMP, &
              BC_INJECTION       )
          pPlagBc%PDFBc%nbins = nbins        
          pPlagBc%PDFBc%pdfvalues(:,1:3) = tmpvec(:,3:5)
      END SELECT ! bcType
    END DO ! iPatch  

! *****************************************************************************
!   Find the maximum in the PDF curve. This maximum represents 
!    the most probable diameter value, from which the search will start.  
! *****************************************************************************                                                      

    tmpscal = 0_RFREAL
    
    DO k = 1,nbins
      IF ( tmpvec(k+1,4) - tmpvec(k,4) > tmpscal ) THEN
        tmpint = k+1
        tmpscal = tmpvec(k+1,4) - tmpvec(k,4)
      END IF ! tmpvec
    END DO ! k 
    
    tmpscal = tmpvec(tmpint,4)
    
    DO iPatch = 1,pRegion%grid%nPatches
      pPatch  => pRegion%patches(iPatch)
      pPlagBc => pPatch%plag
      bcType = pPatch%bcType

      SELECT CASE(bcType)
        CASE( BC_INFLOW,         &
	      BC_INFLOW_VELTEMP, &
              BC_INJECTION       )     
           pPlagBc%PDFBc%locmax = tmpint
           pPlagBc%PDFBc%valmax = tmpscal
      END SELECT ! bcType
    END DO ! iPatch  

! *****************************************************************************
!   Deallocate memory and close file
! *****************************************************************************
  
    DEALLOCATE(tmpvec, STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) & 
      CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__, 'tmpvec' )
  
    CLOSE(IF_PLAG_INFLOWPDF,IOSTAT=errorFlag)
    global%error = errorFlag
      IF (global%error /= ERR_NONE) &
       CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__, &
                       'File: '//TRIM(iFileName)       )

    GOTO 999

10  CONTINUE
    CALL ErrorStop( global, ERR_FILE_READ,__LINE__,'File: '//TRIM(iFileName) )

! *****************************************************************************
!   End
! *****************************************************************************

999 CONTINUE

    CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLU_ReadPdfFromFile




!******************************************************************************
!
! Purpose: Print user boundary conditions parameters for 
!          Lagrangians particles to default values.
!
! Description: none.
!
! Input:
!   pRegion     Pointer to region
!
! Notes: none.
!
!******************************************************************************

  SUBROUTINE PLAG_PrintBcData( pRegion )

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: bcType,iCont,iPatch,iReg
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_bcvalues_plag), POINTER :: pPlagBc
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_PrintBcData',__FILE__)

! ******************************************************************************
!   Set variables
! ******************************************************************************

    iReg = pRegion%iRegionGlobal

! *****************************************************************************
!   Print variables on patch inflow-based data
! *****************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch  => pRegion%patches(iPatch)
      pPlagBc => pPatch%plag
      
      bcType  = pPatch%bcType

      SELECT CASE(bcType)
        CASE( BC_INFLOW,         &
              BC_INFLOW_VELTEMP, &
              BC_INJECTION       )
	  WRITE(STDOUT,1010) SOLVER_NAME//'        inflowModel      ', &
                	     pPlagBc%inflowModel
	  WRITE(STDOUT,1020) SOLVER_NAME//'        inflowVelRatio   ', &
                	     pPlagBc%inflowVelRatio
	  WRITE(STDOUT,1020) SOLVER_NAME//'        spLoad           ', &
                	     pPlagBc%inflowspLoad
	  WRITE(STDOUT,1020) SOLVER_NAME//'        inflowBeta       ', &
                	     pPlagBc%inflowBeta
	  WRITE(STDOUT,1010) SOLVER_NAME//'        inflowDiamDist   ', &
                	     pPlagBc%inflowDiamDist
	  WRITE(STDOUT,1020) SOLVER_NAME//'        inflowDiamMean   ', &
                	     pPlagBc%inflowDiamMean
	  WRITE(STDOUT,1020) SOLVER_NAME//'        inflowDiamMin    ', &
                	     pPlagBc%inflowDiamMin
	  WRITE(STDOUT,1020) SOLVER_NAME//'        inflowDiamMax    ', &
                	     pPlagBc%inflowDiamMax
	  WRITE(STDOUT,1020) SOLVER_NAME//'        inflowStdDev     ', &
                	     pPlagBc%inflowStdDev
	  DO iCont = 1, pRegion%plagInput%nCont
	    WRITE(STDOUT,2015) SOLVER_NAME//'          inflowMassFluxRatio(',iCont,')=', &
                	       pPlagBc%inflowMassFluxRatio(iCont)
	    WRITE(STDOUT,2015) SOLVER_NAME//'          inflowTemp(',iCont,')         =', &
                	       pPlagBc%inflowTemp(iCont)
	  ENDDO ! iCont 	  
	  
      END SELECT  ! bcType  
    END DO !  iPatch

! *****************************************************************************
!   Formats
! *****************************************************************************

1000 FORMAT(/,A,1X,80('-'))
1005 FORMAT(/,A)
1010 FORMAT(A,' = ',I2)
1015 FORMAT(A,' = ',I8)
1020 FORMAT(A,' = ',E12.5)

2010 FORMAT(A,I2,A,ES15.5)
2015 FORMAT(A,I2,A,EN15.5)
2020 FORMAT(A,I2,A)
2025 FORMAT(A,I2,A,I4)

! *****************************************************************************
!   End
! *****************************************************************************

  CALL DeregisterFunction( global )

  END SUBROUTINE PLAG_PrintBcData

  
  
  
  
  
  
  
  

!******************************************************************************


END MODULE PLAG_ModBcData

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModBcData.F90,v $
! Revision 1.3  2016/02/24 06:10:38  rahul
! Made changes to accomodate the PLAG_BC_NOINFLOW case in the .bc file.
! This is added to bypass particle injection/ejection issue in shktb case.
! Set PLAG_MODEL to 0 in the .bc file to enable this case.
!
! Revision 1.2  2015/07/27 04:45:42  brollin
! 1) Corrected bug in RFLUCONV where global%gridFormat was used instead of global%gridSrcFormat
! 2) Implemented new subroutine for shock tube problems (Shktb)
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:51  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/08/07 21:54:01  fnajjar
! Removed statements with BC_RANGE since obsolete for RocfluMP
!
! Revision 1.1  2007/05/16 22:44:50  fnajjar
! Initial import
!
! ******************************************************************************

