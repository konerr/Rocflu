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
! Purpose: Collection of routines to read and convert STARCD grids.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModSTARCD.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2008 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModSTARCD

  USE ModParameters
  USE ModDataTypes  
  USE ModError  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  USE RFLU_ModGrid

  USE ModBuildFileNames, ONLY: BuildFileNamePlain

  USE RFLU_ModDimensions, ONLY: RFLU_SetMaxDimension

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

  PRIVATE

! ==============================================================================
! Data
! ==============================================================================

  TYPE t_gridSTARCD
    INTEGER :: nMappings,nPatches
    INTEGER, DIMENSION(:), POINTER :: bf2p
    INTEGER, DIMENSION(:,:), POINTER :: bf2v,patch2bc
  END TYPE t_gridSTARCD

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModSTARCD.F90,v $ $Revision: 1.1.1.1 $'
            
  TYPE(t_gridSTARCD) :: gridSTARCD

! ==============================================================================
! Public procedures
! ==============================================================================

  PUBLIC :: RFLU_ConvSTARCD2ROCFLU, & 
            RFLU_ReadGridSTARCD

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS








! ******************************************************************************
!
! Purpose: Check connectivity.
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

  SUBROUTINE RFLU_CheckGridSTARCD(pRegion)
     
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

    INTEGER :: errorFlag,icg,ifg,ivgMax,ivgMin
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CheckGridSTARCD', &
                          __FILE__)

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Checking connectivity arrays...'    
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Volume grid
! ******************************************************************************

    DO icg = 1,pGrid%nHexsTot
      ivgMin = MINVAL(pGrid%hex2v(1:8,icg))
      ivgMax = MINVAL(pGrid%hex2v(1:8,icg))
            
      IF ( ivgMin < 1 .OR. ivgMax > pGrid%nVertTot ) THEN 
        global%error = ERR_VERTEX_NUMBER      
      END IF ! ivgMin

      IF ( global%error /= ERR_NONE ) THEN
        IF ( global%verbLevel > VERBOSE_NONE ) THEN          
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
        END IF ! global%verbLevel           
        CALL ErrorStop(global,global%error,__LINE__)
      END IF ! global%error
    END DO ! icg

! ******************************************************************************
!   Surface grid
! ******************************************************************************

    DO ifg = 1,pGrid%nBFaces
      ivgMin = MINVAL(gridSTARCD%bf2v(1:4,ifg))
      ivgMax = MAXVAL(gridSTARCD%bf2v(1:4,ifg))
            
      IF ( ivgMin < 1 .OR. ivgMax > pGrid%nVertTot ) THEN 
        global%error = ERR_VERTEX_NUMBER      
      END IF ! ivgMin

      IF ( global%error /= ERR_NONE ) THEN
        IF ( global%verbLevel > VERBOSE_NONE ) THEN          
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
        END IF ! global%verbLevel           
        CALL ErrorStop(global,global%error,__LINE__)
      END IF ! global%error
    END DO ! ifg

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Checking connectivity arrays done.'
    END IF ! global%verbLevel       

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_CheckGridSTARCD







! ******************************************************************************
!
! Purpose: Convert grid format from STARCD to ROCFLU.
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

  SUBROUTINE RFLU_ConvSTARCD2ROCFLU(pRegion)

    USE RFLU_ModCellMapping
    USE RFLU_ModFaceList
    USE RFLU_ModGeometry
     
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

    CHARACTER(CHRLEN) :: iFileName
    INTEGER :: c1,errorFlag,iBegMax,iBegMin,iBeg1,iBeg2,iEndMax,iEndMin,iEnd1, &
               iEnd2,icg,icl,ifg,iFile,ifl,iMap,iMap2,iPatch,iPatch2,ivg,ivl, &
               j,v1,v2,v3,v4
    REAL(RFREAL) :: dotProd,fCenX,fCenY,fCenZ,fVecX,fVecY,fVecZ
    REAL(RFREAL) :: xyzNodes(XCOORD:ZCOORD,4)
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ConvSTARCD2ROCFLU', &
                          __FILE__)

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Converting from STARCD to ROCFLU format...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer and initialize variables
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Convert patch data structure
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Converting patch data structure...'  
    END IF ! global%verbLevel

! ==============================================================================
!   Read patch mapping file
! ==============================================================================

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Reading patch mapping file...'  
    END IF ! global%verbLevel

! ------------------------------------------------------------------------------
!   Open file
! ------------------------------------------------------------------------------

    iFile = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.sgi',iFileName)

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error  

! ------------------------------------------------------------------------------
!   Read file
! ------------------------------------------------------------------------------
  
    READ(iFile,*) pGrid%nPatches
    READ(iFile,*) gridSTARCD%nMappings

    ALLOCATE(gridSTARCD%patch2bc(3,gridSTARCD%nMappings),STAT=errorFlag)
    global%error = errorFlag 
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridSTARCD%patch2bc')
    END IF ! global%error   

    DO iMap = 1,gridSTARCD%nMappings
      READ(iFile,*) (gridSTARCD%patch2bc(j,iMap),j=1,3)
    END DO ! iMap

! ------------------------------------------------------------------------------
!   Close file
! ------------------------------------------------------------------------------

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Reading patch mapping file done.'  
    END IF ! global%verbLevel    

! ==============================================================================
!   Check for consistent input - somewhat complicated...
! ==============================================================================

    IF ( global%checkLevel > CHECK_NONE ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Checking patch mapping entries...'
      END IF ! global%verbLevel

      DO iMap = 1,gridSTARCD%nMappings
        IF ( gridSTARCD%patch2bc(2,iMap) < gridSTARCD%patch2bc(1,iMap) ) THEN
          IF ( global%verbLevel > VERBOSE_NONE ) THEN  
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.' 
          END IF ! global%verbLevel   
          CALL ErrorStop(global,ERR_PATCH_NUMBERING,__LINE__)
        END IF ! gridSTARCD
      END DO ! iMap   

      IF ( MINVAL(gridSTARCD%patch2bc(3,:)) /= 1 .OR. & 
           MAXVAL(gridSTARCD%patch2bc(3,:)) /= pGrid%nPatches ) THEN
        IF ( global%verbLevel > VERBOSE_NONE ) THEN              
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
        END IF ! global%verbLevel               
        CALL ErrorStop(global,ERR_PATCH_NUMBERING,__LINE__)
      END IF ! gridSTARCD

      DO iMap = 1,gridSTARCD%nMappings
        DO iMap2 = 1,gridSTARCD%nMappings

          IF ( iMap /= iMap2 ) THEN 
            iBeg1 = gridSTARCD%patch2bc(1,iMap)
            iEnd1 = gridSTARCD%patch2bc(2,iMap)

            iBeg2 = gridSTARCD%patch2bc(1,iMap2)
            iEnd2 = gridSTARCD%patch2bc(2,iMap2)        

            IF ( iBeg1 < iBeg2 ) THEN 
              iBegMin = iBeg1
              iEndMin = iEnd1
              iBegMax = iBeg2
              iEndMax = iEnd2
            ELSE IF ( iBeg1 > iBeg2 ) THEN 
              iBegMin = iBeg2
              iEndMin = iEnd2
              iBegMax = iBeg1
              iEndMax = iEnd1         
            ELSE ! iBeg1 and iBeg2 have the same value
              IF ( global%verbLevel > VERBOSE_NONE ) THEN
                WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
              END IF ! global%verbLevel    
              CALL ErrorStop(global,ERR_PATCH_NUMBERING,__LINE__)
            END IF ! iBeg1

            IF ( iEndMin >= iBegMax ) THEN
              IF ( global%verbLevel > VERBOSE_NONE ) THEN          
                WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
              END IF ! global%verbLevel    
              CALL ErrorStop(global,ERR_PATCH_NUMBERING,__LINE__)
            END IF ! iEndMin
          END IF ! iMap

        END DO ! iMap2
      END DO ! iMap

      IF ( global%verbLevel > VERBOSE_NONE ) THEN  
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, & 
                                 'Checking patch mapping entries done.'
      END IF ! global%verbLevel  
    END IF ! global%checkLevel

! ==============================================================================
!   Allocate patch memory and initialize patch structure
! ==============================================================================

    ALLOCATE(pRegion%patches(pGrid%nPatches),STAT=errorFlag)
    global%error = errorFlag 
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%patches')
    END IF ! global%error       

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      pPatch%nBTris  = 0
      pPatch%nBQuads = 0 
      pPatch%nBVert  = 0  
        
      pPatch%iPatchGlobal = iPatch
      pPatch%iBorder      = PATCH_IBORDER_DEFAULT      
      pPatch%renumFlag    = .FALSE.
    END DO ! iPatch

    global%nPatches = pGrid%nPatches

! ==============================================================================
!   Determine number of faces on each patch. NOTE assume only quads for now.
! ==============================================================================

    DO ifg = 1,pGrid%nBFaces
      iPatch = gridSTARCD%bf2p(ifg)

      DO iMap = 1,gridSTARCD%nMappings
        IF ( iPatch >= gridSTARCD%patch2bc(1,iMap) .AND. & 
             iPatch <= gridSTARCD%patch2bc(2,iMap) ) THEN 
          iPatch2 = gridSTARCD%patch2bc(3,iMap)
        END IF ! iPatch
      END DO ! iMap
    
      pPatch => pRegion%patches(iPatch2)
    
      pPatch%nBQuads = pPatch%nBQuads + 1
    END DO ! iPatch
          
! ==============================================================================
!   Set total boundary patch quantities and number of boundary faces
! ==============================================================================

    pGrid%nBFaces = 0 
    
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      pPatch%nBFaces = pPatch%nBTris + pPatch%nBQuads
      pGrid%nBFaces  = pGrid%nBFaces + pPatch%nBFaces

      pPatch%nBFacesTot = pPatch%nBFaces
      pPatch%nBQuadsTot = pPatch%nBQuads  
      pPatch%nBTrisTot  = pPatch%nBTris    
      pPatch%nBVertTot  = pPatch%nBVert  
      
      pPatch%nBTrisMax  = RFLU_SetMaxDimension(global,pPatch%nBTrisTot)
      pPatch%nBQuadsMax = RFLU_SetMaxDimension(global,pPatch%nBQuadsTot) 
      pPatch%nBFacesMax = RFLU_SetMaxDimension(global,pPatch%nBFacesTot)             
      pPatch%nBVertMax  = RFLU_SetMaxDimension(global,pPatch%nBVertTot)
      
      pPatch%nBCellsVirt = 0                    
    END DO ! iPatch  
    
    pGrid%nBFacesTot = pGrid%nBFaces    
    
! ==============================================================================
!   Allocate memory for boundary-face connectivity
! ==============================================================================

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      IF ( pPatch%nBTrisMax > 0 ) THEN 
        ALLOCATE(pPatch%bTri2v(3,pPatch%nBTrisMax),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bTri2v')
        END IF ! global%error 
      ELSE 
        NULLIFY(pPatch%bTri2v)
      END IF ! pPatch%nBTrisMax
      
      IF ( pPatch%nBQuadsMax > 0 ) THEN      
        ALLOCATE(pPatch%bQuad2v(4,pPatch%nBQuadsMax),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bQuad2v')
        END IF ! global%error                       
      ELSE 
        NULLIFY(pPatch%bQuad2v)
      END IF ! pPatch%nBQuadsMax
    END DO ! iPatch    
    
! ==============================================================================
!   Build boundary face connectivity
! ==============================================================================
    
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      pPatch%nBTris  = 0
      pPatch%nBQuads = 0              
    END DO ! iPatch    
    
    DO ifg = 1,pGrid%nBFaces
      iPatch = gridSTARCD%bf2p(ifg)

      DO iMap = 1,gridSTARCD%nMappings
        IF ( iPatch >= gridSTARCD%patch2bc(1,iMap) .AND. & 
             iPatch <= gridSTARCD%patch2bc(2,iMap) ) THEN 
          iPatch2 = gridSTARCD%patch2bc(3,iMap)
        END IF ! iPatch
      END DO ! iMap
            
      pPatch => pRegion%patches(iPatch2)

      pPatch%nBQuads = pPatch%nBQuads + 1
            
      DO ivl = 1,4
        pPatch%bQuad2v(ivl,pPatch%nBQuads) = gridSTARCD%bf2v(ivl,ifg)
      END DO ! ivl            
    END DO ! ifg

! ******************************************************************************
!   Check that boundary connectivity is correct, i.e., that it leads to outward-
!   pointing face normals. If that is not the case, reverse ordering of faces.
! ******************************************************************************

! ==============================================================================
!   Build cell mapping and face list
! ==============================================================================
    
    CALL RFLU_CreateCellMapping(pRegion)
    CALL RFLU_BuildLoc2GlobCellMapping(pRegion)
    CALL RFLU_BuildGlob2LocCellMapping(pRegion)

    CALL RFLU_CreateFaceList(pRegion)
    CALL RFLU_BuildFaceList(pRegion)

! ==============================================================================
!   Compute approximate centroids and flip numbering if necessary. NOTE cell 
!   mapping (built above) required for computation of approximate centroids.
! ==============================================================================

    CALL RFLU_CreateApproxCentroids(pRegion)
    CALL RFLU_ComputeApproxCentroids(pRegion)
    
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      DO ifl = 1,pPatch%nBQuads
        ifg = ifl + pPatch%nBTris 

        c1 = pPatch%bf2c(ifg)

        v1 = pPatch%bQuad2v(1,ifl) 
        v2 = pPatch%bQuad2v(2,ifl) 
        v3 = pPatch%bQuad2v(3,ifl) 
        v4 = pPatch%bQuad2v(4,ifl) 

        xyzNodes(XCOORD:ZCOORD,1) = pGrid%xyz(XCOORD:ZCOORD,v1)
        xyzNodes(XCOORD:ZCOORD,2) = pGrid%xyz(XCOORD:ZCOORD,v2)
        xyzNodes(XCOORD:ZCOORD,3) = pGrid%xyz(XCOORD:ZCOORD,v3)
        xyzNodes(XCOORD:ZCOORD,4) = pGrid%xyz(XCOORD:ZCOORD,v4)

        CALL FaceVectorQuad(xyzNodes(XCOORD:ZCOORD,1:4), &
                            fVecX,fVecY,fVecZ)
        CALL FaceCentroidQuad(xyzNodes(XCOORD:ZCOORD,1:4), &
                              fCenX,fCenY,fCenZ)

        dotProd = fVecX*(fCenX - pGrid%cofgApp(XCOORD,c1)) &
                + fVecY*(fCenY - pGrid%cofgApp(YCOORD,c1)) &
                + fVecZ*(fCenZ - pGrid%cofgApp(ZCOORD,c1))

        IF ( dotProd < 0.0_RFREAL ) THEN 
          pPatch%bQuad2v(1,ifl) = v1
          pPatch%bQuad2v(2,ifl) = v4
          pPatch%bQuad2v(3,ifl) = v3
          pPatch%bQuad2v(4,ifl) = v2
        END IF ! dotProd
      END DO ! ifl
    END DO ! iPatch
    
! ******************************************************************************
!   Dellocate temporary memory 
! ******************************************************************************

    CALL RFLU_DestroyApproxCentroids(pRegion)
    CALL RFLU_DestroyFaceList(pRegion)
    CALL RFLU_DestroyCellMapping(pRegion)   
 
! ******************************************************************************
!   Allocate memory for other Rocflu data structures 
! ******************************************************************************
    
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      ALLOCATE(pPatch%bf2c(pPatch%nBFacesMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2c')
      END IF ! global%error

      ALLOCATE(pPatch%bf2v(4,pPatch%nBFacesMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2v')
      END IF ! global%error  

      DO ifl = 1,pPatch%nBFacesMax
        pPatch%bf2v(1,ifl) = VERT_NONE
        pPatch%bf2v(2,ifl) = VERT_NONE 
        pPatch%bf2v(3,ifl) = VERT_NONE 
        pPatch%bf2v(4,ifl) = VERT_NONE                                    
      END DO ! ifl
    END DO ! iPatch    
        
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Converting from STARCD to ROCFLU format done.'
    END IF ! global%verbLevel       

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_ConvSTARCD2ROCFLU









! *******************************************************************************
!
! Purpose: Print STARCD grid information.
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
! *******************************************************************************

  SUBROUTINE RFLU_PrintGridSTARCDInfo(pRegion)

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
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start, set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid
    
! ******************************************************************************
!   Write information
! ******************************************************************************

    WRITE(STDOUT,'(A,3X,A)')       SOLVER_NAME,'Grid Statistics:'
    WRITE(STDOUT,'(A,5X,A,2X,I9)') SOLVER_NAME,'Vertices:       ', & 
                                   pGrid%nVertTot
    WRITE(STDOUT,'(A,5X,A,2X,I9)') SOLVER_NAME,'Cells:          ', & 
                                   pGrid%nCellsTot                   
    WRITE(STDOUT,'(A,5X,A,2X,I9)') SOLVER_NAME,'Patches:        ', & 
                                   gridSTARCD%nPatches   

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_PrintGridSTARCDInfo










! *******************************************************************************
!
! Purpose: Read grid file from STARCD format.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
! 
! Output: None.
!
! Notes:
!   1. CENTAUR cell and node pointers are not read in - read into the 
!      dummy integer idum.
!
! *******************************************************************************

  SUBROUTINE RFLU_ReadGridSTARCD(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNamePlain
  
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

    CHARACTER(CHRLEN) :: iFileName,line
    INTEGER :: dummyInteger,errorFlag,icg,iFile,ifg,ivg
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadGridSTARCD', &
                          __FILE__)

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading STARCD grid file...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer and initialize variables
! ******************************************************************************

    pGrid => pRegion%grid
    
    pGrid%nVert = 0
    pGrid%nTets = 0        
    pGrid%nHexs = 0        
    pGrid%nPris = 0        
    pGrid%nPyrs = 0     
                  
    pGrid%nPatches = 0

! ******************************************************************************
!   Read and set dimensions
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading dimensions...'
    END IF ! global%verbLevel

    iFile = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.scd',iFileName)

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF

    emptyLoop: DO
      READ(iFile,'(A)',IOSTAT=errorFlag) line

      IF ( errorFlag == 0 ) THEN
        IF ( line(1:4) == 'RNAM' ) THEN
          gridSTARCD%nPatches = gridSTARCD%nPatches + 1
        ELSE
          SELECT CASE ( line(1:5) )
            CASE ( 'VREAD' )
              Call RFLU_GetNumber(line,pGrid%nVert)
            CASE ( 'CREAD' )
              Call RFLU_GetNumber(line,pGrid%nCells)
            CASE ( 'BREAD' )
              Call RFLU_GetNumber(line,pGrid%nBFaces)
            CASE DEFAULT
          END SELECT
        END IF ! line
      ELSE
        EXIT emptyLoop
      END IF ! errorFlag
    END DO emptyLoop

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF

    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading dimensions done.'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set total and max dimensions
! ******************************************************************************

    pGrid%nHexs = pGrid%nCells ! TEMPORARY restriction

    pGrid%nVertTot  = pGrid%nVert
    pGrid%nCellsTot = pGrid%nCells
    pGrid%nTetsTot  = pGrid%nTets
    pGrid%nHexsTot  = pGrid%nHexs
    pGrid%nPrisTot  = pGrid%nPris
    pGrid%nPyrsTot  = pGrid%nPyrs

    pGrid%nVertMax  = RFLU_SetMaxDimension(global,pGrid%nVertTot)
    pGrid%nCellsMax = RFLU_SetMaxDimension(global,pGrid%nCellsTot)
    pGrid%nTetsMax  = RFLU_SetMaxDimension(global,pGrid%nTetsTot)
    pGrid%nHexsMax  = RFLU_SetMaxDimension(global,pGrid%nHexsTot)
    pGrid%nPrisMax  = RFLU_SetMaxDimension(global,pGrid%nPrisTot)
    pGrid%nPyrsMax  = RFLU_SetMaxDimension(global,pGrid%nPyrsTot)

! ******************************************************************************
!   Read coordinates
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading coordinates...'
    END IF ! global%verbLevel

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.vrt',iFileName)

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF

    ALLOCATE(pGrid%xyz(XCOORD:ZCOORD,pGrid%nVertMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'grid%xyz')
    END IF ! global%error

    DO ivg = 1,pGrid%nVert
      READ(iFile,*) dummyInteger,pGrid%xyz(XCOORD:ZCOORD,ivg)
    END DO ! ivg

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF

    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading coordinates done.'
    END IF ! global%verbLevel

! ******************************************************************************
!   Read connectivity
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading connectivity...'
    END IF ! global%verbLevel

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.cel',iFileName)

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF

    ALLOCATE(pGrid%hex2v(8,pGrid%nHexsMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%hex2v')
    END IF ! global%error

    DO icg = 1,pGrid%nHexs
      READ(iFile,*) dummyInteger,pGrid%hex2v(1:8,icg)
    END DO ! icg

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF

    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading connectivity done.'
    END IF ! global%verbLevel

! ******************************************************************************
!   Read boundary faces
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading boundary faces...'
    END IF ! global%verbLevel

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.bnd',iFileName)

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF

    ALLOCATE(gridSTARCD%bf2v(4,pGrid%nBFaces),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridSTARCD%bf2v')
    END IF ! global%error

    ALLOCATE(gridSTARCD%bf2p(pGrid%nBFaces),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridSTARCD%bf2p')
    END IF ! global%error

    DO ifg = 1,pGrid%nBFaces
      READ(iFile,*) dummyInteger,gridSTARCD%bf2v(1:4,ifg),gridSTARCD%bf2p(ifg)
    END DO ! ifg

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF

    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading boundary faces done.'
    END IF ! global%verbLevel

! ******************************************************************************
!   Check validity of connectivity arrays
! ******************************************************************************

    IF ( global%checkLevel > CHECK_NONE ) THEN
      CALL RFLU_CheckGridSTARCD(pRegion)
    END IF ! global%checkLevel

! ******************************************************************************
!   Print grid statistics
! ******************************************************************************
 
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      CALL RFLU_PrintGridSTARCDInfo(pRegion) 
    END IF ! global%verbLevel   

! ******************************************************************************
!   Close file
! ******************************************************************************

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error  

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN    
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Reading STARCD grid file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)
   
  END SUBROUTINE RFLU_ReadGridSTARCD





! *******************************************************************************
!
! Purpose: Get number from string.
!
! Description: None.
!
! Input: 
!   line	String containing number with commas
! 
! Output: 
!   n		Number
!
! Notes:
!   1. Assume number is located between 4th and 5th comma in string.
!
! *******************************************************************************

  SUBROUTINE RFLU_GetNumber(line,n)

    CHARACTER(80), INTENT(IN) :: line
    INTEGER, INTENT(OUT) :: n

    INTEGER :: i
    INTEGER :: cp(0:5)

    cp(0) = 0

    DO i = 1,5
      cp(i) = INDEX(line(cp(i-1)+1:LEN_TRIM(line)),",")
      cp(i) = cp(i-1) + cp(i)
    END DO ! i

    READ(line(cp(4)+1:cp(5)-1),*) n

  END SUBROUTINE RFLU_GetNumber



! ******************************************************************************
! End
! ******************************************************************************


END MODULE RFLU_ModSTARCD

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModSTARCD.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:56  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:10  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2008/02/09 23:30:20  haselbac
! Cosmetics only
!
! ******************************************************************************

