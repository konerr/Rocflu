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
! Purpose: Collection of routines for patch coefficients.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModPatchCoeffs.F90,v 1.2 2015/07/23 23:11:18 brollin Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModPatchCoeffs

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModMPI

  USE ModBuildFileNames, ONLY: BuildFileNameSteady, &
                               BuildFileNameUnsteady

  IMPLICIT NONE
  
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModPatchCoeffs.F90,v $ $Revision: 1.2 $'        

! ==============================================================================
! Public functions
! ==============================================================================
  
  PUBLIC :: RFLU_CreatePatchCoeffs, & 
            RFLU_DestroyPatchCoeffs, & 
            RFLU_ReadPatchCoeffsWrapper, &
            RFLU_WritePatchCoeffsWrapper

! ==============================================================================
! Private functions
! ==============================================================================

  PRIVATE :: RFLU_ClosePatchCoeffs, & 
             RFLU_InitPatchCoeffs, &
             RFLU_NullifyPatchCoeffs, &
             RFLU_OpenPatchCoeffsASCII, &
             RFLU_OpenPatchCoeffsBinary, &
             RFLU_ReadPatchCoeffsASCII, & 
             RFLU_ReadPatchCoeffsBinary, & 
             RFLU_WritePatchCoeffsASCII, & 
             RFLU_WritePatchCoeffsBinary

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  
  
  
! *******************************************************************************
!
! Purpose: Close patch-coefficients file.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ClosePatchCoeffs(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag
    TYPE(t_global), POINTER :: global    
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_ClosePatchCoeffs',__FILE__)
             
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Closing patch-coefficients file...'
    END IF ! global%verbLevel             
             
! ******************************************************************************
!   Close file
! ******************************************************************************

    CLOSE(IF_PATCH_COEF,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__)
    END IF ! global%error

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Closing patch-coefficients file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ClosePatchCoeffs  
  
  
  
  
  

! *******************************************************************************
!
! Purpose: Create patch coefficients.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_CreatePatchCoeffs(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_CreatePatchCoeffs',__FILE__)
    
    pGrid => pRegion%grid
        
! ******************************************************************************
!   Allocate memory
! ******************************************************************************
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      ALLOCATE(pPatch%cp(pPatch%nBFaces),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%cp')
      END IF ! global%error

      ALLOCATE(pPatch%cf(XCOORD:ZCOORD,pPatch%nBFaces),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%cf')
      END IF ! global%error
      
      ALLOCATE(pPatch%ch(pPatch%nBFaces),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%ch')
      END IF ! global%error                              
      
      ALLOCATE(pPatch%cmass(pPatch%nBFaces),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%cmass')
      END IF ! global%error
      
      ALLOCATE(pPatch%cmom(XCOORD:ZCOORD,pPatch%nBFaces),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%cmom')
      END IF ! global%error
    END DO ! iPatch

! ******************************************************************************
!   Initialize memory
! ******************************************************************************

    CALL RFLU_InitPatchCoeffs(pRegion)

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_CreatePatchCoeffs
 








! *******************************************************************************
!
! Purpose: Destroy patch coefficients.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_DestroyPatchCoeffs(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,iPatch
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_DestroyPatchCoeffs',__FILE__)
    
    pGrid => pRegion%grid
        
! ******************************************************************************
!   Deallocate memory
! ******************************************************************************
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      DEALLOCATE(pPatch%cp,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%cp')
      END IF ! global%error

      DEALLOCATE(pPatch%cf,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%cf')
      END IF ! global%error
      
      DEALLOCATE(pPatch%ch,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%ch')
      END IF ! global%error              
      
      DEALLOCATE(pPatch%cmass,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%cmass')
      END IF ! global%error
      
      DEALLOCATE(pPatch%cmom,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%cmom')
      END IF ! global%error
    END DO ! iPatch

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyPatchCoeffs(pRegion)

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_DestroyPatchCoeffs





! *******************************************************************************
!
! Purpose: Initialize patch coefficients.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_InitPatchCoeffs(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,ifl,iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_InitPatchCoeffs',__FILE__)
    
    pGrid => pRegion%grid
        
! ******************************************************************************
!   Allocate memory
! ******************************************************************************
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      DO ifl = 1,pPatch%nBFaces
        pPatch%cp(ifl)        = 0.0_RFREAL
        pPatch%cf(XCOORD,ifl) = 0.0_RFREAL
        pPatch%cf(YCOORD,ifl) = 0.0_RFREAL
        pPatch%cf(ZCOORD,ifl) = 0.0_RFREAL
        pPatch%ch(ifl)        = 0.0_RFREAL                
        pPatch%cmass(ifl)     = 0.0_RFREAL                
        pPatch%cmom(XCOORD,ifl) = 0.0_RFREAL                
        pPatch%cmom(YCOORD,ifl) = 0.0_RFREAL                
        pPatch%cmom(ZCOORD,ifl) = 0.0_RFREAL                
      END DO ! ifl                      
    END DO ! iPatch

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_InitPatchCoeffs







! *******************************************************************************
!
! Purpose: Nullify patch coefficients.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_NullifyPatchCoeffs(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: iPatch
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_NullifyPatchCoeffs',__FILE__)
    
    pGrid => pRegion%grid
        
! ******************************************************************************
!   Deallocate memory
! ******************************************************************************
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      NULLIFY(pPatch%cp)
      NULLIFY(pPatch%cf)
      NULLIFY(pPatch%ch)           
      NULLIFY(pPatch%cmass)           
      NULLIFY(pPatch%cmom)           
    END DO ! iPatch

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_NullifyPatchCoeffs







! *******************************************************************************
!
! Purpose: Open patch-coefficient file in ASCII format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   fileStatus          File status
!
! Output: 
!   fileExists          Flag indicating existence of file
!
! Notes:
!   1. File may not exist when it is opened in postprocessor before solver 
!      was run, so need to deal with this gracefully.
!
! ******************************************************************************

  SUBROUTINE RFLU_OpenPatchCoeffsASCII(pRegion,fileStatus,fileExists)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    INTEGER, INTENT(IN) :: fileStatus
    LOGICAL, INTENT(OUT), OPTIONAL :: fileExists
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName
    INTEGER :: errorFlag,iFile
    TYPE(t_global), POINTER :: global 
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_OpenPatchCoeffsASCII',__FILE__)
         
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Opening ASCII patch-coefficients file...'
    END IF ! global%verbLevel     
     
    iFile = IF_PATCH_COEF
     
! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.pcoa', & 
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
      CALL BuildFileNameSteady(global,FILEDEST_OUTDIR,'.pcoa', & 
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

! ******************************************************************************
!   Open file
! ******************************************************************************

    IF ( fileStatus == FILE_STATUS_OLD ) THEN 
      INQUIRE(FILE=iFileName,EXIST=fileExists)
    
      IF ( fileExists .EQV. .TRUE. ) THEN     
        OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD", &
             IOSTAT=errorFlag)
        global%error = errorFlag          
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
        END IF ! global%error
      END IF ! fileExists
    ELSE IF ( fileStatus == FILE_STATUS_UNKNOWN ) THEN 
      OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", &
           IOSTAT=errorFlag)
      global%error = errorFlag          
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
      END IF ! global%error      
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)      
    END IF ! fileStatus       

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Opening ASCII patch-coefficients file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_OpenPatchCoeffsASCII







! *******************************************************************************
!
! Purpose: Open patch-coefficient file in binary format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   fileStatus          File status
!
! Output: 
!   fileExists          Flag indicating existence of file
!
! Notes: 
!   1. File may not exist when it is opened in postprocessor before solver 
!      was run, so need to deal with this gracefully.
!
! ******************************************************************************

  SUBROUTINE RFLU_OpenPatchCoeffsBinary(pRegion,fileStatus,fileExists)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    INTEGER, INTENT(IN) :: fileStatus
    LOGICAL, INTENT(OUT), OPTIONAL :: fileExists
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName
    INTEGER :: errorFlag,iFile
    TYPE(t_global), POINTER :: global 
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_OpenPatchCoeffsBinary',__FILE__)
         
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Opening binary patch-coefficients file...'
    END IF ! global%verbLevel     
     
    iFile = IF_PATCH_COEF     
     
! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.pco', & 
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
      CALL BuildFileNameSteady(global,FILEDEST_OUTDIR,'.pco', & 
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

! ******************************************************************************
!   Open file
! ******************************************************************************

    IF ( fileStatus == FILE_STATUS_OLD ) THEN 
      INQUIRE(FILE=iFileName,EXIST=fileExists)
    
      IF ( fileExists .EQV. .TRUE. ) THEN     
!        OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
!             IOSTAT=errorFlag)
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
      END IF ! fileExists
    ELSE IF ( fileStatus == FILE_STATUS_UNKNOWN ) THEN 
!      OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
!           IOSTAT=errorFlag)
! BBR - begin
    IF( global%solutFormat .EQ. FORMAT_BINARY )THEN
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
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)      
    END IF ! fileStatus         

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Opening binary patch-coefficients file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_OpenPatchCoeffsBinary







! *******************************************************************************
!
! Purpose: Read patch coefficients in ASCII format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadPatchCoeffsASCII(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,iFile,ifl,iPatch,loopCounter
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_ReadPatchCoeffsASCII',__FILE__)
        
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Reading ASCII patch-coefficients file...'
    END IF ! global%verbLevel        
   
    pGrid => pRegion%grid
    
    iFile = IF_PATCH_COEF   
        
! ******************************************************************************
!   Header and general information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU patch-coefficients file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM

! ******************************************************************************
!   Read rest of file
! ******************************************************************************

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile,'(A)') sectionString

      SELECT CASE ( TRIM(sectionString) ) 

! ==============================================================================
!       Pressure coefficient
! ==============================================================================

        CASE ( '# Pressure coefficient' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Pressure coefficient...'
          END IF ! global%verbLevel    

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)

            READ(iFile,'(5(E23.16))') (pPatch%cp(ifl),ifl=1,pPatch%nBFaces)       
          END DO ! iPatch
        
! ==============================================================================
!       Skin-friction coefficient
! ==============================================================================

        CASE ( '# Skin-friction coefficient' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Skin-friction coefficient...'
          END IF ! global%verbLevel    

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)

            READ(iFile,'(5(E23.16))') (pPatch%cf(XCOORD,ifl), & 
                                       ifl=1,pPatch%nBFaces)       
            READ(iFile,'(5(E23.16))') (pPatch%cf(YCOORD,ifl), &
                                       ifl=1,pPatch%nBFaces)       
            READ(iFile,'(5(E23.16))') (pPatch%cf(ZCOORD,ifl), &
                                       ifl=1,pPatch%nBFaces) 
          END DO ! iPatch

! ==============================================================================
!       Heat-transfer coefficient
! ==============================================================================

        CASE ( '# Heat-transfer coefficient' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Heat-transfer coefficient...'
          END IF ! global%verbLevel    

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)

            READ(iFile,'(5(E23.16))') (pPatch%ch(ifl),ifl=1,pPatch%nBFaces)       
          END DO ! iPatch

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
!     Guard against infinite loop - unnecessary because of read errors?
! ==============================================================================
  
      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter
    END DO ! <empty>       

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Reading ASCII patch-coefficients file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadPatchCoeffsASCII






! *******************************************************************************
!
! Purpose: Read patch coefficients in binary format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadPatchCoeffsBinary(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,iFile,ifl,iPatch,loopCounter
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_ReadPatchCoeffsBinary',__FILE__)
   
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Reading binary patch-coefficients file...'
    END IF ! global%verbLevel
   
    pGrid => pRegion%grid   
     
    iFile = IF_PATCH_COEF
         
! ******************************************************************************
!   Header and general information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU patch-coefficients file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM

! ******************************************************************************
!   Read rest of file
! ******************************************************************************

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile) sectionString

      SELECT CASE ( TRIM(sectionString) ) 

! ==============================================================================
!       Pressure coefficient
! ==============================================================================

        CASE ( '# Pressure coefficient' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Pressure coefficient...'
          END IF ! global%verbLevel    

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)

            READ(iFile) (pPatch%cp(ifl),ifl=1,pPatch%nBFaces)       
          END DO ! iPatch
        
! ==============================================================================
!       Skin-friction coefficient
! ==============================================================================

        CASE ( '# Skin-friction coefficient' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Skin-friction coefficient...'
          END IF ! global%verbLevel    

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)

            READ(iFile) (pPatch%cf(XCOORD,ifl),ifl=1,pPatch%nBFaces)       
            READ(iFile) (pPatch%cf(YCOORD,ifl),ifl=1,pPatch%nBFaces)       
            READ(iFile) (pPatch%cf(ZCOORD,ifl),ifl=1,pPatch%nBFaces) 
          END DO ! iPatch

! ==============================================================================
!       Heat-transfer coefficient
! ==============================================================================

        CASE ( '# Heat-transfer coefficient' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Heat-transfer coefficient...'
          END IF ! global%verbLevel    

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)

            READ(iFile) (pPatch%ch(ifl),ifl=1,pPatch%nBFaces)       
          END DO ! iPatch

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
!     Guard against infinite loop - unnecessary because of read errors?
! ==============================================================================
  
      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter
    END DO ! <empty>       

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Reading binary patch-coefficients file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadPatchCoeffsBinary





! *******************************************************************************
!
! Purpose: Wrapper for reading patch coefficients.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   fileStatus          File status
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadPatchCoeffsWrapper(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    LOGICAL :: fileExists
    TYPE(t_global), POINTER :: global  
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_ReadPatchCoeffsWrapper',__FILE__)

! ******************************************************************************
!   Read file (if it exists)
! ******************************************************************************
    
    IF ( global%solutFormat == FORMAT_ASCII ) THEN
      CALL RFLU_OpenPatchCoeffsASCII(pRegion,FILE_STATUS_OLD,fileExists)
      
      IF ( fileExists .EQV. .TRUE. ) THEN
        CALL RFLU_ReadPatchCoeffsASCII(pRegion)
        CALL RFLU_ClosePatchCoeffs(pRegion)
      END IF ! fileExists
!    ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
! BBR - begin
       ELSE IF ( global%solutFormat == FORMAT_BINARY .OR. & 
                global%solutFormat == FORMAT_BINARY_L .OR. &
                global%solutFormat == FORMAT_BINARY_B ) THEN
! BBR - end
      CALL RFLU_OpenPatchCoeffsBinary(pRegion,FILE_STATUS_OLD,fileExists)
      
      IF ( fileExists .EQV. .TRUE. ) THEN 
        CALL RFLU_ReadPatchCoeffsBinary(pRegion)
        CALL RFLU_ClosePatchCoeffs(pRegion)        
      END IF ! fileExists
    ELSE
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! global%solutFormat   

! ******************************************************************************
!   Write warning if file is missing
! ******************************************************************************
    
    IF ( fileExists .EQV. .FALSE. ) THEN 
      global%warnCounter = global%warnCounter + 1
    
      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,2(1X,A))') SOLVER_NAME,'*** WARNING ***', &
                                    'Patch coefficient file missing, not read.'        
      END IF ! global%myProcid
    END IF ! fileExists 
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadPatchCoeffsWrapper






! *******************************************************************************
!
! Purpose: Write patch coefficients in ASCII format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_WritePatchCoeffsASCII(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,iFile,ifl,iPatch
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_WritePatchCoeffsASCII',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing ASCII patch-coefficients file...'
    END IF ! global%verbLevel
    
    pGrid => pRegion%grid

    iFile = IF_PATCH_COEF
            
! ******************************************************************************
!   Header and general information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    sectionString = '# ROCFLU patch-coefficients file'
    WRITE(iFile,'(A)') sectionString
              
! ******************************************************************************
!   Write patch coefficients to file
! ******************************************************************************
    
! ==============================================================================
!   Pressure coefficient
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Pressure coefficient...'
    END IF ! global%verbLevel

    sectionString = '# Pressure coefficient'
    WRITE(iFile,'(A)') sectionString
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      WRITE(iFile,'(5(E23.16))') (pPatch%cp(ifl),ifl=1,pPatch%nBFaces)       
    END DO ! iPatch
    
! ==============================================================================
!   Skin-friction coefficient
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Skin-friction coefficient...'
    END IF ! global%verbLevel

    sectionString = '# Skin-friction coefficient'
    WRITE(iFile,'(A)') sectionString
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      WRITE(iFile,'(5(E23.16))') (pPatch%cf(XCOORD,ifl),ifl=1,pPatch%nBFaces)
      WRITE(iFile,'(5(E23.16))') (pPatch%cf(YCOORD,ifl),ifl=1,pPatch%nBFaces) 
      WRITE(iFile,'(5(E23.16))') (pPatch%cf(ZCOORD,ifl),ifl=1,pPatch%nBFaces)                    
    END DO ! iPatch    
    
! ==============================================================================
!   Heat-transfer coefficient
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Heat-transfer coefficient...'
    END IF ! global%verbLevel

    sectionString = '# Heat-transfer coefficient'
    WRITE(iFile,'(A)') sectionString
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      WRITE(iFile,'(5(E23.16))') (pPatch%ch(ifl),ifl=1,pPatch%nBFaces)       
    END DO ! iPatch    

! ******************************************************************************
!   End marker
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile,'(A)') sectionString  

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Writing ASCII patch-coefficients file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WritePatchCoeffsASCII






! *******************************************************************************
!
! Purpose: Write patch coefficients in binary format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_WritePatchCoeffsBinary(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,iFile,ifl,iPatch
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_WritePatchCoeffsBinary',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing binary patch-coefficients file...'
    END IF ! global%verbLevel
    
    pGrid => pRegion%grid

    iFile = IF_PATCH_COEF
             
! ******************************************************************************
!   Header and general information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    sectionString = '# ROCFLU patch-coefficients file'
    WRITE(iFile) sectionString
              
! ******************************************************************************
!   Write patch coefficients to file
! ******************************************************************************
    
! ==============================================================================
!   Pressure coefficient
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Pressure coefficient...'
    END IF ! global%verbLevel

    sectionString = '# Pressure coefficient'
    WRITE(iFile) sectionString
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      WRITE(iFile) (pPatch%cp(ifl),ifl=1,pPatch%nBFaces)       
    END DO ! iPatch
    
! ==============================================================================
!   Skin-friction coefficient
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Skin-friction coefficient...'
    END IF ! global%verbLevel

    sectionString = '# Skin-friction coefficient'
    WRITE(iFile) sectionString
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      WRITE(iFile) (pPatch%cf(XCOORD,ifl),ifl=1,pPatch%nBFaces)
      WRITE(iFile) (pPatch%cf(YCOORD,ifl),ifl=1,pPatch%nBFaces) 
      WRITE(iFile) (pPatch%cf(ZCOORD,ifl),ifl=1,pPatch%nBFaces)                    
    END DO ! iPatch    
    
! ==============================================================================
!   Heat-transfer coefficient
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Heat-transfer coefficient...'
    END IF ! global%verbLevel

    sectionString = '# Heat-transfer coefficient'
    WRITE(iFile) sectionString
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      WRITE(iFile) (pPatch%ch(ifl),ifl=1,pPatch%nBFaces)       
    END DO ! iPatch    

! ******************************************************************************
!   End marker
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile) sectionString  

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Writing binary patch-coefficients file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WritePatchCoeffsBinary






! *******************************************************************************
!
! Purpose: Wrapper for writing patch coefficients.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_WritePatchCoeffsWrapper(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    TYPE(t_global), POINTER :: global  
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_WritePatchCoeffsWrapper',__FILE__)
    
    IF ( global%solutFormat == FORMAT_ASCII ) THEN
      CALL RFLU_OpenPatchCoeffsASCII(pRegion,FILE_STATUS_UNKNOWN)
      CALL RFLU_WritePatchCoeffsASCII(pRegion)
      CALL RFLU_ClosePatchCoeffs(pRegion)
!    ELSE IF ( global%solutFormat == FORMAT_BINARY ) THEN
! BBR - begin
       ELSE IF ( global%solutFormat == FORMAT_BINARY .OR. & 
                global%solutFormat == FORMAT_BINARY_L .OR. &
                global%solutFormat == FORMAT_BINARY_B ) THEN
! BBR - end 
      CALL RFLU_OpenPatchCoeffsBinary(pRegion,FILE_STATUS_UNKNOWN)
      CALL RFLU_WritePatchCoeffsBinary(pRegion)
      CALL RFLU_ClosePatchCoeffs(pRegion)      
    ELSE
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! global%solutFormat   
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WritePatchCoeffsWrapper






END MODULE RFLU_ModPatchCoeffs

! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLU_ModPatchCoeffs.F90,v $
!   Revision 1.2  2015/07/23 23:11:18  brollin
!   1) The pressure coefficient of the  collision model has been changed back to its original form
!   2) New options in the format of the grid and solutions have been added. Now the user can choose the endianness, and convert from one to the over in rfluconv.
!   3) The solutions are now stored in folders named by timestamp or iteration number.
!   4) The address enty in the hashtable has been changed to an integer(8) for cases when the grid becomes very large.
!   5) RFLU_WritePM can now compute PM2 on the fly for the Macroscale problem
!
!   Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
!   merged rocflu micro and macro
!
!   Revision 1.1.1.1  2014/07/15 14:31:37  brollin
!   New Stable version
!
!   Revision 1.3  2008/12/06 08:43:44  mtcampbe
!   Updated license.
!
!   Revision 1.2  2008/11/19 22:16:56  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.1  2007/04/09 18:49:25  haselbac
!   Initial revision after split from RocfloMP
!
!   Revision 1.1  2007/04/09 18:00:41  haselbac
!   Initial revision after split from RocfloMP
!
!   Revision 1.4  2006/10/20 21:18:56  mparmar
!   Added allocation/deallocation of cmass and cmom
!
!   Revision 1.3  2006/04/07 15:19:20  haselbac
!   Removed tabs
!
!   Revision 1.2  2004/12/07 19:39:14  haselbac
!   Bug fix: Incorrect WRITE statements when writing binary data
!
!   Revision 1.1  2004/06/16 20:00:59  haselbac
!   Initial revision
!
! ******************************************************************************

