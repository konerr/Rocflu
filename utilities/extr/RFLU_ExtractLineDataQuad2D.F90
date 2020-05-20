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
! Purpose: Extract data along line in quadrilateral grid. 
!
! Description: None.
!
! Input: 
!   regions		Region data
!   iRegStart		Starting region index
!   icgStart  		Starting cell index
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_ExtractLineDataQuad2D.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_ExtractLineDataQuad2D(regions,iRegStart,icgStart)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModMPI
  USE ModParameters
  
  USE ModSortSearch
  
  USE RFLU_ModBoundLists
  USE RFLU_ModDimensions
  USE RFLU_ModCellFaceEdgeInfo
  USE RFLU_ModCellMapping
  USE RFLU_ModFaceList
  USE RFLU_ModGeometry
  USE RFLU_ModGeometryTools
  USE RFLU_ModReadBcInputFile
  USE RFLU_ModReadWriteGrid
  USE RFLU_ModRenumberings

  USE ModInterfaces, ONLY: RFLU_BuildStencilInterp, & 
                           RFLU_ComputeBilinearInterpWghts, & 
                           RFLU_CreateGrid, & 
                           RFLU_DestroyGrid
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: icgStart,iRegStart
  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Local variables
! ==============================================================================

  LOGICAL :: exitRegionLoopFlag 
  CHARACTER :: infoType
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: c1,c1k,c2,c2k,icg,icgs,ifg,iloc,iReg,iRegNew,trajLoopCounter
  REAL(RFREAL) :: dist,distTot,eps,iMagTraj,xLoc,xTraj,yLoc,yTraj,zLoc,zTraj
  REAL(RFREAL) :: wghts(4)
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid,pGridNew,pGridSerial
  TYPE(t_region), POINTER :: pRegion,pRegionNew,pRegionSerial

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ExtractLineDataQuad2D.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_ExtractLineDataQuad2D',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extracting line...'
  END IF ! global%verbLevel

! ******************************************************************************
! Some preparatory stuff
! ******************************************************************************

  pRegionSerial => regions(0)
  pGridSerial   => pRegionSerial%grid

  trajLoopCounter = 0

  eps = EPSILON(1.0_RFREAL)

  iReg = iRegStart
  icg  = icgStart

  xLoc = global%extrXCoordBeg
  yLoc = global%extrYCoordBeg 
  zLoc = global%extrZCoordBeg

  xTraj = global%extrXCoordEnd - xLoc
  yTraj = global%extrYCoordEnd - yLoc
  zTraj = global%extrZCoordEnd - zLoc

  distTot  = SQRT(xTraj*xTraj + yTraj*yTraj + zTraj*zTraj)
  iMagTraj = 1.0_RFREAL/distTot 

  xTraj = iMagTraj*xTraj
  yTraj = iMagTraj*yTraj
  zTraj = iMagTraj*zTraj    

! ******************************************************************************
! Loop over regions
! ******************************************************************************

  regionLoop: DO
    exitRegionLoopFlag = .FALSE.

    pRegion => regions(iReg)
    pGrid   => pRegion%grid

    CALL RFLU_ReadDimensions(pRegion)
    CALL RFLU_CreateGrid(pRegion)

    IF ( pRegion%grid%nPatches > 0 ) THEN
      CALL RFLU_ReadBCInputFileWrapper(pRegion)
    END IF ! pRegion%grid%nPatches

    CALL RFLU_ReadGridWrapper(pRegion)

    CALL RFLU_CreateCellMapping(pRegion)
    CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
    CALL RFLU_BuildGlob2LocCellMapping(pRegion)

    CALL RFLU_CreateBVertexLists(pRegion)
    CALL RFLU_BuildBVertexLists(pRegion)

    CALL RFLU_CreateFaceList(pRegion)
    CALL RFLU_BuildFaceList(pRegion)
    CALL RFLU_RenumberBFaceLists(pRegion)

    CALL RFLU_CreateCell2FaceList(pRegion)
    CALL RFLU_BuildCell2FaceList(pRegion)

    CALL RFLU_CreateGeometry(pRegion)
    CALL RFLU_BuildGeometry(pRegion)

! ******************************************************************************
!   Loop over trajectory 
! ******************************************************************************

    trajLoop: DO 
      trajLoopCounter = trajLoopCounter + 1
      
! ==============================================================================
!     Find appropriate intersection and associated face  
! ==============================================================================

      CALL RFLU_ComputeLineCellXSectSafe(pRegion,xLoc,yLoc,zLoc,xTraj,yTraj, &
                                         zTraj,icg,dist,iloc,ifg)

! ==============================================================================
!     Update total distance travelled                             
! ==============================================================================
                                   
      distTot = distTot - dist

! ==============================================================================
!     Check whether have remaining total distance 
! ==============================================================================

! ------------------------------------------------------------------------------
!     No distance remaining
! ------------------------------------------------------------------------------

      IF ( distTot <= 0.0_RFREAL ) THEN 
        exitRegionLoopFlag = .TRUE.

        EXIT trajLoop 
        
! ------------------------------------------------------------------------------
!     Distance remaining: Keep searching 
! ------------------------------------------------------------------------------
        
      ELSE 

! ----- Trajectory intersects interior face ------------------------------------

        IF ( iloc == 0 ) THEN           
          c1 = pGrid%f2c(1,ifg)
          c2 = pGrid%f2c(2,ifg)
          
          c1k = RFLU_GetGlobalCellKind(global,pGrid,c1)
          c2k = RFLU_GetGlobalCellKind(global,pGrid,c2)

! ------- If intersection distance is greater than epsilon, add to list

          IF ( dist > eps ) THEN 
            IF ( pGridSerial%nXSect < pGridSerial%nXSectMax ) THEN 
              pGridSerial%nXSect = pGridSerial%nXSect + 1

              pGridSerial%xSectList(1,pGridSerial%nXSect) = iReg
              pGridSerial%xSectGeom(1,pGridSerial%nXSect) = xLoc
              pGridSerial%xSectGeom(2,pGridSerial%nXSect) = yLoc
              pGridSerial%xSectGeom(3,pGridSerial%nXSect) = zLoc

              CALL RFLU_BuildStencilInterp(pRegion,ifg,xLoc,yLoc,zLoc, & 
                   pGridSerial%xSectList(2:5,pGridSerial%nXSect))
              CALL RFLU_ComputeBilinearInterpWghts(pRegion,xLoc,yLoc,zLoc, &
                   pGridSerial%xSectList(2:5,pGridSerial%nXSect), & 
                   pGridSerial%xSectGeom(4:7,pGridSerial%nXSect))
            ELSE 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END IF ! pGridSerial%nXSectMax 
          END IF ! dist

! ------- Distinguish between AA and AV intersections 

          SELECT CASE ( RFLU_GetFaceKind(global,c1k,c2k) ) 
            CASE ( FACE_KIND_AA ) ! Actual-actual face
              IF ( c1 == icg ) THEN 
                icg = c2
              ELSE 
                icg = c1
              END IF ! c1                                    
            CASE ( FACE_KIND_AV ) ! Actual-virtual face
              IF ( c1 == icg ) THEN 
                icg = c2
              ELSE 
                icg = c1
              END IF ! c1

              CALL RFLU_ReadDimensions(pRegion)

              CALL RFLU_RNMB_CreatePC2SCMap(pRegion)
              CALL RFLU_RNMB_CreatePV2SVMap(pRegion)    
              CALL RFLU_RNMB_CreatePBF2SBFMap(pRegion) 

              CALL RFLU_RNMB_ReadPxx2SxxMaps(pRegion)
              
              icgs = pGrid%pc2sc(icg)

              iRegNew = pGridSerial%sc2r(icgs)
              
              CALL RFLU_RNMB_DestroyPC2SCMap(pRegion)
              CALL RFLU_RNMB_DestroyPV2SVMap(pRegion)    
              CALL RFLU_RNMB_DestroyPBF2SBFMap(pRegion)   

              pRegionNew => regions(iRegNew)
              pGridNew   => pRegionNew%grid
              
              CALL RFLU_ReadDimensions(pRegionNew)
              CALL RFLU_CreateGrid(pRegionNew)

              IF ( pRegion%grid%nPatches > 0 ) THEN
                CALL RFLU_ReadBCInputFileWrapper(pRegionNew)
              END IF ! pRegion%grid%nPatches

              CALL RFLU_RNMB_CreatePC2SCMap(pRegionNew)
              CALL RFLU_RNMB_CreatePV2SVMap(pRegionNew)    
              CALL RFLU_RNMB_CreatePBF2SBFMap(pRegionNew)
                                          
              CALL RFLU_RNMB_ReadPxx2SxxMaps(pRegionNew)
              
              CALL RFLU_RNMB_BuildSC2PCMap(pRegionNew)

              CALL BinarySearchInteger(pGridNew%sc2pc(1:1,1:pGridNew%nCellsTot), & 
                                       pGridNew%nCellsTot,icgs,iloc)

              IF ( iloc /= ELEMENT_NOT_FOUND ) THEN 
                icg = pGridNew%sc2pc(2,iloc)
              ELSE 
                CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
              END IF ! iloc                          
              
              CALL RFLU_RNMB_DestroySC2PCMap(pRegionNew)
              CALL RFLU_RNMB_DestroyPC2SCMap(pRegionNew)
              CALL RFLU_RNMB_DestroyPV2SVMap(pRegionNew)    
              CALL RFLU_RNMB_DestroyPBF2SBFMap(pRegionNew)              
	      
              EXIT trajLoop                  
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! RFLU_GetFaceKind                    
                  
! ------------------------------------------------------------------------------  
!       Trajectory intersects boundary face 
! ------------------------------------------------------------------------------  
          
        ELSE ! Boundary face
          exitRegionLoopFlag = .TRUE.

          EXIT trajLoop 
        END IF ! iLoc 
      END IF ! distTot      
      
! ==============================================================================
!     Guard against infinite loop 
! ==============================================================================
  
      IF ( trajLoopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! trajLoopCounter      
    END DO trajLoop

! ******************************************************************************
!   Destroy temporary memory
! ******************************************************************************

    CALL RFLU_DestroyCellMapping(pRegion)
    CALL RFLU_DestroyBVertexLists(pRegion)
    CALL RFLU_DestroyFaceList(pRegion)
    CALL RFLU_DestroyCell2FaceList(pRegion)
    CALL RFLU_DestroyGeometry(pRegion)
    CALL RFLU_DestroyGrid(pRegion)
  
! ******************************************************************************
!   Exit out of region loop if necessary otherwise set new region index and 
!   continue
! ******************************************************************************

    IF ( exitRegionLoopFlag .EQV. .TRUE. ) THEN 
      EXIT regionLoop
    END IF ! exitRegionLoopFlag

    iReg = iRegNew
  END DO regionLoop
  
! ******************************************************************************
! End
! ******************************************************************************

  IF ( (global%myProcid == MASTERPROC) .AND. & 
       (global%verbLevel > VERBOSE_NONE) ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Extracting line done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractLineDataQuad2D

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ExtractLineDataQuad2D.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:55  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:07  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/11/27 13:17:26  haselbac
! Initial revision
!
!******************************************************************************

