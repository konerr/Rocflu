! ******************************************************************************
!
! Purpose: Collection of routines to read and convert OMG grids.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModOMG.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2010 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModOMG

  USE ModParameters
  USE ModDataTypes  
  USE ModError  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  USE RFLU_ModDimensions, ONLY: RFLU_SetMaxDimension

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

  PRIVATE

! ==============================================================================
! Data
! ==============================================================================

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModOMG.F90,v $ $Revision: 1.1.1.1 $'        

  INTEGER, DIMENSION(4,4), PARAMETER :: f2vTetOMG = &
    RESHAPE((/1,2,3,VERT_NONE,1,4,2,VERT_NONE,2,4,3,VERT_NONE,1,3,4, &
              VERT_NONE/), (/4,4/))

! ==============================================================================
! Public procedures
! ==============================================================================

  PUBLIC :: RFLU_ConvOMG2ROCFLU, & 
            RFLU_ReadGridOMGASCII

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

  SUBROUTINE RFLU_CheckGridOMG(pRegion)
     
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

    INTEGER :: errorFlag,cvmax,cvmin
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CheckGridOMG', &
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

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, &
                               'Volume grid...'    
    END IF ! global%verbLevel

    IF ( pGrid%nTetsTot > 0 ) THEN 
      cvmin = MINVAL(pGrid%tet2v(1:4,1:pGrid%nTetsTot))
      cvmax = MAXVAL(pGrid%tet2v(1:4,1:pGrid%nTetsTot))

      IF ( pGrid%nTetsTot == pGrid%nCellsTot ) THEN       
        IF ( cvmin /= 1 .OR. cvmax /= pGrid%nVertTot ) THEN 
          global%error = ERR_VERTEX_NUMBER
        END IF ! cvmin
      ELSE  
        IF ( cvmin < 1 .OR. cvmax > pGrid%nVertTot ) THEN 
          global%error = ERR_VERTEX_NUMBER
        END IF ! vmin
      END IF ! cvmin

      IF ( global%error /= ERR_NONE ) THEN
        IF ( global%verbLevel > VERBOSE_NONE ) THEN          
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
        END IF ! global%verbLevel           
        CALL ErrorStop(global,global%error,__LINE__)
      END IF ! global%error
    END IF ! pGrid

! ******************************************************************************
!   Surface grid
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN     
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Surface grid...'    
    END IF ! global%verbLevel

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                               'Checking connectivity arrays done.'
    END IF ! global%verbLevel       

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_CheckGridOMG





! ******************************************************************************
!
! Purpose: Convert grid format from OMG to ROCFLU.
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

  SUBROUTINE RFLU_ConvOMG2ROCFLU(pRegion)

    USE RFLU_ModCellMapping
    USE RFLU_ModGeometry, ONLY: RFLU_ComputeApproxCentroids, & 
                                RFLU_CreateApproxCentroids, &
                                RFLU_DestroyApproxCentroids
     
    USE ModInterfaces, ONLY: FaceCentroidTria, &
                             FaceVectorTria

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

    INTEGER :: cntr,errorFlag,icg,ifg,ifl,iPatch,nBQuads,nBTris,v1,v2,v3
    REAL(RFREAL) :: cofgAppX,cofgAppY,cofgAppZ,dotProd,fCenX,fCenY,fCenZ, &
                    fVecX,fVecY,fVecZ
    REAL(RFREAL) :: xyzNodes(XCOORD:ZCOORD,4)
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ConvOMG2ROCFLU', &
                          __FILE__)

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Converting from OMG to ROCFLU format...'
    END IF ! global%verbLevel

! ==============================================================================
!   Set grid pointer and initialize variables
! ==============================================================================

    pGrid => pRegion%grid

    pGrid%nEdges    = 0
    pGrid%nEdgesTot = 0

    pGrid%nFaces    = 0
    pGrid%nFacesTot = 0

! ==============================================================================
!   Check that number of triangular and quadrilateral faces correct and set 
!   number of boundary faces
! ==============================================================================

    nBTris  = 0
    nBQuads = 0

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      pPatch%iPatchGlobal = iPatch
      pPatch%iBorder      = PATCH_IBORDER_DEFAULT
      pPatch%renumFlag    = .FALSE.

      nBTris  = nBTris  + pPatch%nBTris
      nBQuads = nBQuads + pPatch%nBQuads
    END DO ! iPatch 

    pGrid%nBFaces    = nBTris + nBQuads
    pGrid%nBFacesTot = pGrid%nBFaces    

    global%nPatches = pGrid%nPatches

! *****************************************************************************
!   Compute approximate centroids. NOTE cell mapping (built above) required for
!   computation of approximate centroids.
! *****************************************************************************

    CALL RFLU_CreateCellMapping(pRegion)
    CALL RFLU_BuildLoc2GlobCellMapping(pRegion)
    CALL RFLU_BuildGlob2LocCellMapping(pRegion)

    CALL RFLU_CreateApproxCentroids(pRegion)
    CALL RFLU_ComputeApproxCentroids(pRegion)

! *****************************************************************************
!   Convert bTri2v list
! *****************************************************************************

    cntr = 0

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      DO ifg = 1,pPatch%nBTris
        icg = pPatch%bTri2v(1,ifg)
        ifl = pPatch%bTri2v(2,ifg)

        cofgAppX = pGrid%cofgApp(XCOORD,icg)
        cofgAppY = pGrid%cofgApp(YCOORD,icg)
        cofgAppZ = pGrid%cofgApp(ZCOORD,icg)

        v1 = pGrid%tet2v(f2vTetOMG(1,ifl),icg)
        v2 = pGrid%tet2v(f2vTetOMG(2,ifl),icg)
        v3 = pGrid%tet2v(f2vTetOMG(3,ifl),icg)

        xyzNodes(XCOORD:ZCOORD,1) = pGrid%xyz(XCOORD:ZCOORD,v1)
        xyzNodes(XCOORD:ZCOORD,2) = pGrid%xyz(XCOORD:ZCOORD,v2)
        xyzNodes(XCOORD:ZCOORD,3) = pGrid%xyz(XCOORD:ZCOORD,v3)

        CALL FaceVectorTria(xyzNodes(XCOORD:ZCOORD,1:3),fVecX,fVecY,fVecZ)
        CALL FaceCentroidTria(xyzNodes(XCOORD:ZCOORD,1:3),fCenX,fCenY,fCenZ)

        dotProd = fVecX*(fCenX - cofgAppX) &
                + fVecY*(fCenY - cofgAppY) &
                + fVecZ*(fCenZ - cofgAppZ)

        IF ( dotProd < 0.0_RFREAL ) THEN
          cntr = cntr + 1
          
          pPatch%bTri2v(1,ifg) = v1
          pPatch%bTri2v(2,ifg) = v3
          pPatch%bTri2v(3,ifg) = v2
        ELSE 
          pPatch%bTri2v(1,ifg) = v1
          pPatch%bTri2v(2,ifg) = v2
          pPatch%bTri2v(3,ifg) = v3
        END IF ! dotProd
      END DO ! ifg
    END DO ! iPatch

    CALL RFLU_DestroyApproxCentroids(pRegion)
    CALL RFLU_DestroyCellMapping(pRegion)

! *****************************************************************************
!   Allocate memory for boundary face lists bf2c and bf2v
! *****************************************************************************

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
                               'Converting from OMG to ROCFLU format done.'
    END IF ! global%verbLevel       

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_ConvOMG2ROCFLU







! *******************************************************************************
!
! Purpose: Print CENTAUR grid information.
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

  SUBROUTINE RFLU_PrintGridOMGInfo(pRegion)

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
    WRITE(STDOUT,'(A,7X,A,I9)')    SOLVER_NAME,'Tetrahedra:     ', & 
                                   pGrid%nTetsTot
    WRITE(STDOUT,'(A,7X,A,I9)')    SOLVER_NAME,'Hexahedra:      ', & 
                                   pGrid%nHexsTot
    WRITE(STDOUT,'(A,7X,A,I9)')    SOLVER_NAME,'Prisms:         ', & 
                                   pGrid%nPrisTot
    WRITE(STDOUT,'(A,7X,A,I9)')    SOLVER_NAME,'Pyramids:       ', & 
                                   pGrid%nPyrsTot
    WRITE(STDOUT,'(A,5X,A,2X,I9)') SOLVER_NAME,'Patches:        ', & 
                                   pGrid%nPatches   

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_PrintGridOMGInfo






! *******************************************************************************
!
! Purpose: Read grid file from CENTAUR in ASCII format.
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

  SUBROUTINE RFLU_ReadGridOMGASCII(pRegion)

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

    CHARACTER(CHRLEN) :: dummyString,iFileName
    CHARACTER(12), DIMENSION(:), ALLOCATABLE :: dummyStringArray
    INTEGER :: dummyInteger,errorFlag,i,iFile,iPatch,j,nCellsInDomain,nDomain, &
               nHex,nWedge,nOct
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!   Start, open file and read title
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadGridOMGASCII', &
                          __FILE__)

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII OMG grid file...'
    END IF ! global%verbLevel

    iFile  = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.msh',iFileName)

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF

! ******************************************************************************
!   Read header. NOTE assume have only tets.
! ******************************************************************************

    pGrid => pRegion%grid

    READ(iFile,*) dummyString
    READ(iFile,*) dummyString

    READ(iFile,*) pGrid%nVertTot,pGrid%nTetsTot,nHex,nWedge,nOct,nDomain, &
                  pGrid%nPatches

    IF ( (nHex > 0) .OR. (nWedge > 0) .OR. (nOct > 0) .OR. (nDomain > 1) ) THEN 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! nHex
  
! ******************************************************************************
!   Coordinates
! ******************************************************************************

    pGrid%nVertMax = RFLU_SetMaxDimension(global,pGrid%nVertTot)

    ALLOCATE(pGrid%xyz(3,pGrid%nVertMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'grid%xyz')
    END IF ! global%error

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Coordinates...'
    END IF ! global%verbLevel

    DO j = 1,pGrid%nVertTot 
      READ(iFile,*) (pGrid%xyz(i,j),i=XCOORD,ZCOORD)  
    END DO ! j

! ******************************************************************************
!   Cell connectivity
! ******************************************************************************

    pGrid%nTetsMax = RFLU_SetMaxDimension(global,pGrid%nTetsTot)

    IF ( pGrid%nTetsMax > 0 ) THEN 
      ALLOCATE(pGrid%tet2v(4,pGrid%nTetsMax),STAT=errorFlag)
      global%error = errorFlag 
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%tet2v')
      END IF ! global%error
    ELSE 
      NULLIFY(pGrid%tet2v)
    END IF ! pGrid%nTetsMax

    IF ( pGrid%nTetsTot > 0 ) THEN 
      IF ( global%verbLevel > VERBOSE_NONE ) THEN  
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Tetrahedra...'
      END IF ! global%verbLevel 

      DO j = 1,pGrid%nTetsTot         
        READ(iFile,*) (pGrid%tet2v(i,j),i=1,4)
      END DO ! j
    END IF ! pGrid%nTetsTot

    pGrid%nHexsTot = 0 
    pGrid%nPrisTot = 0
    pGrid%nPyrsTot = 0

! ******************************************************************************
!   Set grid size variables
! ******************************************************************************

    pGrid%nCellsTot = pGrid%nTetsTot + pGrid%nHexsTot + pGrid%nPrisTot & 
                    + pGrid%nPyrsTot
    pGrid%nCellsMax = RFLU_SetMaxDimension(global,pGrid%nCellsTot)

    pGrid%nVert  = pGrid%nVertTot
    pGrid%nCells = pGrid%nCellsTot
    pGrid%nTets  = pGrid%nTetsTot
    pGrid%nHexs  = pGrid%nHexsTot
    pGrid%nPris  = pGrid%nPrisTot
    pGrid%nPyrs  = pGrid%nPyrsTot

! ******************************************************************************
!   Domain information. NOTE for the moment, assume have single domain and that 
!   therefore the cell-in-domain array simply has entries starting with 1 and 
!   ending with nCells, so not necessary to read array
! ******************************************************************************

    READ(iFile,'(I8,A)') nCellsInDomain,dummyString

    IF ( nCellsInDomain /= pGrid%nCells ) THEN 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! nCellsInDomain

    READ(iFile,'(10I8)') (dummyInteger,i=1,nCellsInDomain)

! ******************************************************************************
!   Boundary information. NOTE somewhat cumbersome way of reading bTri2v info
!   caused by IBM compiler, which did not like READ statement with implied DO 
!   loop. 
! ******************************************************************************
  
    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary information...'   
    END IF ! global%verbLevel

    ALLOCATE(pRegion%patches(pGrid%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%patches')
    END IF ! global%errorFlag

    DO iPatch = 1,pGrid%nPatches 
      pPatch => pRegion%patches(iPatch)

      READ(iFile,'(I8,A)') pPatch%nBTrisTot,dummyString

      pPatch%nBQuadsTot = 0

      pPatch%nBTris  = pPatch%nBTrisTot 
      pPatch%nBQuads = pPatch%nBQuadsTot
                
      pPatch%nBFacesTot = pPatch%nBTrisTot + pPatch%nBQuadsTot
      pPatch%nBFaces    = pPatch%nBFacesTot

      pPatch%nBTrisMax  = RFLU_SetMaxDimension(global,pPatch%nBTrisTot)
      pPatch%nBQuadsMax = RFLU_SetMaxDimension(global,pPatch%nBQuadsTot)
      pPatch%nBFacesMax = RFLU_SetMaxDimension(global,pPatch%nBFacesTot)

      ALLOCATE(pPatch%bTri2v(3,pPatch%nBTrisMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bTri2v')
      END IF ! global%errorFlag

      ALLOCATE(dummyStringArray(pPatch%nBTris),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'dummyStringArray')
      END IF ! global%errorFlag

      NULLIFY(pPatch%bQuad2v)

      READ(iFile,'(5(A12))') (dummyStringArray(i),i=1,pPatch%nBTris)

      DO i = 1,pPatch%nBTris
        READ(dummyStringArray(i),'(I10,I2)') pPatch%bTri2v(1,i), &
                                             pPatch%bTri2v(2,i)
      END DO ! i

      DEALLOCATE(dummyStringArray,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'dummyStringArray')
      END IF ! global%errorFlag

    END DO ! iPatch 

! ******************************************************************************
!   Check validity of connectivity arrays
! ******************************************************************************

    IF ( global%checkLevel > CHECK_NONE ) THEN
      CALL RFLU_CheckGridOMG(pRegion)
    END IF ! global%checkLevel

! ******************************************************************************
!   Print grid statistics
! ******************************************************************************
 
    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      CALL RFLU_PrintGridOMGInfo(pRegion) 
    END IF ! global%verbLevel  

! ******************************************************************************
!   Close file
! ******************************************************************************

    CLOSE(iFile, IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error  

! ******************************************************************************
! End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN    
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Reading ASCII OMG grid file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)
   
  END SUBROUTINE RFLU_ReadGridOMGASCII






! ******************************************************************************
! End
! ******************************************************************************


END MODULE RFLU_ModOMG

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModOMG.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.1  2010/05/24 16:07:38  haselbac
! Initial revision
!
! ******************************************************************************

