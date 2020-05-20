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
! Purpose: Suite of utility routines to partition region.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModPartitionRegionUtils.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005-2007 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModPartitionRegionUtils

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  USE ModSortSearch

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: RFLU_PART_AddVirtualCellsInv1, & 
            RFLU_PART_AddVirtualCellsInv2, & 
            RFLU_PART_MarkVirtualCells, & 
            RFLU_PART_SetVirtualCellRegionIndex

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModPartitionRegionUtils.F90,v $ $Revision: 1.1.1.1 $'


! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS






! ******************************************************************************
!
! Purpose: Add virtual cells for inviscid fluxes based on first-order scheme.
!
! Description: None.
!
! Input:
!   pRegion	        Pointer to region
!   pRegionSerial	Pointer to serial region	
!   vc			List of virtual cells (empty)
!   nCellsVirtMax	Maximum allowable number of virtual cells  
!
! Output: 
!   vc			List of virtual cells 
!   nCellsVirt   	Number of virtual cells
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_AddVirtualCellsInv1(pRegion,pRegionSerial,vc, &
                                           nCellsVirtMax,nCellsVirt)                                   

    IMPLICIT NONE
      
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: nCellsVirtMax
    INTEGER, INTENT(OUT) :: nCellsVirt
    INTEGER, INTENT(INOUT) :: vc(nCellsVirtMax)    
    TYPE(t_region), POINTER :: pRegion,pRegionSerial

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: c1,c2,errorFlag,icg,icl,ifg,ifl,iflBeg,iflEnd,iLoc
    TYPE(t_grid), POINTER :: pGrid,pGridSerial
    TYPE(t_global), POINTER :: global     
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegionSerial%global

    CALL RegisterFunction(global,'RFLU_PART_AddVirtualCellsInv1', &
         __FILE__)

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
            'Adding virtual cells for inviscid first-order stencil...' 
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers and initialize variables
! ******************************************************************************

    pGrid       => pRegion%grid
    pGridSerial => pRegionSerial%grid    

    nCellsVirt = 0

    iflBeg = pGridSerial%avfCSRInfo(pRegion%iRegionGlobal)

    IF ( pRegion%iRegionGlobal /= global%nRegionsLocal ) THEN
      iflEnd = pGridSerial%avfCSRInfo(pRegion%iRegionGlobal+1)-1
    ELSE 
      iflEnd = 2*pGridSerial%nFacesCut
    END IF ! pRegion%iRegionGlobal
   
! ******************************************************************************
!  Loop over actual-virtual faces
! ******************************************************************************
    
    DO ifl = iflBeg,iflEnd        
      ifg = pGridSerial%avfCSR(ifl)

      DO icl = 1,2
        icg = pGridSerial%f2c(icl,ifg)

        IF ( pGridSerial%sc2r(icg) /= pRegion%iRegionGlobal ) THEN
          IF ( nCellsVirt > 0 ) THEN 
            CALL BinarySearchInteger(vc(1:nCellsVirt),nCellsVirt,icg,iLoc)
          ELSE 
            iLoc = ELEMENT_NOT_FOUND
          END IF ! nCellsVirt

          IF ( iLoc == ELEMENT_NOT_FOUND ) THEN 
            nCellsVirt = nCellsVirt + 1

            IF ( nCellsVirt <= nCellsVirtMax ) THEN 
              vc(nCellsVirt) = icg
            ELSE 
              CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'vc')
            END IF ! nCellsVirt
              
            IF ( nCellsVirt > 1 ) THEN    
              CALL QuickSortInteger(vc(1:nCellsVirt),nCellsVirt)              
            END IF ! nCellsVirt 
          END IF ! iLoc
        END IF ! pGridSerial%sc2r 
      END DO ! icl     
    END DO ! ifl 
                          
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
            'Adding virtual cells for inviscid first-order stencil done.' 
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_AddVirtualCellsInv1








! ******************************************************************************
!
! Purpose: Add virtual cells for inviscid fluxes based on higher-order scheme.
!
! Description: Building list of virtual cells proceeds in several steps:
!   1. Build list of vertices from list of actual-virtual faces
!   2. Build list of source cells adjacent to actual-virtual faces which are 
!      in same region as that for which virtual cells are constructed. One 
!      layer of source cells is sufficient. If more layers of source cells 
!      are specified here, this should not change the number of virtual cells.
!   3. Build stencils for source cells and add stencil members to list of 
!      virtual cells if not in same region
!   4. Loop over existing virtual cells and add layers.
!
! Input:
!   pRegion	        Pointer to region
!   pRegionSerial	Pointer to serial region	
!   vc			List of virtual cells (empty)
!   nCellsVirtMax	Maximum allowable number of virtual cells  
!
! Output: 
!   vc			List of virtual cells 
!   nCellsVirt   	Number of virtual cells
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_AddVirtualCellsInv2(pRegion,pRegionSerial,vc, &
                                           nCellsVirtMax,nCellsVirt)                                   

    USE RFLU_ModStencilsCells
    USE RFLU_ModTopologyUtils, ONLY: RFLU_BuildFaceVertList, & 
                                     RFLU_BuildVertCellNghbList

    IMPLICIT NONE
      
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: nCellsVirtMax
    INTEGER, INTENT(OUT) :: nCellsVirt
    INTEGER, INTENT(INOUT) :: vc(nCellsVirtMax)    
    TYPE(t_region), POINTER :: pRegion,pRegionSerial

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER, PARAMETER :: LOOP_COUNTER_LIMIT = 5
    INTEGER :: errorFlag,i,icg,icg2,iDimens,iDir,iflBeg,iflEnd,iLayer,iLoc,j, &
              loopCounter,nLayers,nVert,nVertEst,scDim,vcNewDim,vcNewDimMax, &
              vcOldDim
    INTEGER, DIMENSION(:), ALLOCATABLE :: avv,sc,vcNew,vcOld
    TYPE(t_grid), POINTER :: pGrid,pGridSerial
    TYPE(t_global), POINTER :: global     
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegionSerial%global

    CALL RegisterFunction(global,'RFLU_PART_AddVirtualCellsInv2', &
         __FILE__)

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
            'Adding virtual cells for inviscid higher-order stencil...' 
    END IF ! global%myProcid

! ******************************************************************************
!   Set pointers and initialize variables
! ******************************************************************************

    pGrid       => pRegion%grid
    pGridSerial => pRegionSerial%grid    

    nCellsVirt = 0

    loopCounter = 0

! ******************************************************************************
!   Build list of interface vertices
! ******************************************************************************

! ==============================================================================
!   Estimate number of interface vertices 
! ==============================================================================

    iflBeg = pGridSerial%avfCSRInfo(pRegion%iRegionGlobal)

    IF ( pRegion%iRegionGlobal /= global%nRegionsLocal ) THEN
      iflEnd = pGridSerial%avfCSRInfo(pRegion%iRegionGlobal+1)-1
    ELSE 
      iflEnd = 2*pGridSerial%nFacesCut
    END IF ! pRegion%iRegionGlobal
      
    nVertEst = 2*(iflEnd - iflBeg + 1)  

    IF ( nVertEst < 100 ) THEN ! Kludge
      nVertEst = 100 + 4*nVertEst
    ELSE IF ( nVertEst < 1000 ) THEN 
      nVertEst = 2*nVertEst
    END IF ! nVertEst      

! ==============================================================================
!   Build list of interface vertices. NOTE this sometimes failed for some cases
!   because the estimated size was too small, so added capability of adaptively 
!   increasing size within some limits to prevent code from aborting.
! ==============================================================================

    emptyLoop: DO 
      loopCounter = loopCounter + 1

      ALLOCATE(avv(nVertEst),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'avv')
      END IF ! global%error      
      
      CALL RFLU_BuildFaceVertList(global,pGridSerial, & 
                                  pGridSerial%avfCSR(iflBeg:iflEnd), & 
                                  iflEnd-iflBeg+1,avv,nVertEst,nVert,errorFlag)
     
      IF ( errorFlag /= ERR_NONE ) THEN 
        IF ( loopCounter <= LOOP_COUNTER_LIMIT ) THEN 
          DEALLOCATE(avv,STAT=errorFlag)
          global%error = errorFlag   
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'avv')
          END IF ! global%error      
  
          nVertEst = 2*nVertEst

          global%warnCounter = global%warnCounter + 1

          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A,1X,I2,1X,A)') SOLVER_NAME, & 
                         '*** WARNING *** Attempt ',loopCounter, & 
                         'to build vertex list failed because array too small.' 
            WRITE(STDOUT,'(A,19X,A,1X,I7)') SOLVER_NAME, & 
                         'Attempting again with array size:',nVertEst
          END IF ! global%myProcid
        ELSE
          CALL ErrorStop(global,ERR_ALLOCATE_ADAPTIVE,__LINE__) 
        END IF ! loopCounter
      ELSE 
        EXIT emptyLoop
      END IF ! errorFlag
    END DO emptyLoop

! ******************************************************************************
!   Based on list of interface vertices, build list of source cells adjacent to 
!   actual-virtual faces for construction of virtual cells. One layer of source 
!   cells is sufficient. If more layers of source cells are specified here, this
!   should not change the number of virtual cells.
! ******************************************************************************
                               
    ALLOCATE(sc(nCellsVirtMax),STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'sc')
    END IF ! global%error             

    nLayers = 1    

    CALL RFLU_BuildVertCellNghbList(global,pGridSerial,avv,nVert,nLayers, &
                                    pRegion%iRegionGlobal,sc,nCellsVirtMax, &
                                    scDim)

    DEALLOCATE(avv,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'avv')
    END IF ! global%error 
                                
! ******************************************************************************
!   Loop over list of source cells, build stencil for each cell, and add cells
!   in stencil as virtual cells if not in same region and not already added
! ******************************************************************************

    CALL RFLU_SetInfoC2CStencilWrapper(pRegionSerial, & 
                                       pRegionSerial%mixtInput%spaceOrder-1)        
    CALL RFLU_CreateC2CStencilWrapper(pRegionSerial)   

    DO i = 1,scDim
      icg = sc(i)
      
      CALL RFLU_BuildC2CStencilWrapper(pRegionSerial,icg,CONSTR_NONE)
      
      SELECT CASE ( pRegionSerial%mixtInput%stencilDimensCells ) 
        CASE ( 1 ) 
          DO iDimens = 1,pRegionSerial%mixtInput%dimens
            iDir = XCOORD - 1 + iDimens 

            DO j = 1,pGridSerial%c2cs1D(iDir,icg)%nCellMembs
              icg2 = pGridSerial%c2cs1D(iDir,icg)%cellMembs(j)

              IF ( pGridSerial%sc2r(icg2) /= pRegion%iRegionGlobal ) THEN
                IF ( nCellsVirt > 0 ) THEN
                  CALL BinarySearchInteger(vc(1:nCellsVirt),nCellsVirt,icg2,iLoc)
                ELSE
                  iLoc = ELEMENT_NOT_FOUND
                END IF ! nCellsVirt

                IF ( iLoc == ELEMENT_NOT_FOUND ) THEN
                  nCellsVirt = nCellsVirt + 1
  
                  vc(nCellsVirt) = icg2
  
                  IF ( nCellsVirt > 1 ) THEN
                    CALL QuickSortInteger(vc(1:nCellsVirt),nCellsVirt)
                  END IF ! nCellsVirt 
                END IF ! iLoc
              END IF ! pGridSerial%sc2r                      
            END DO ! j
          END DO ! iDimens
        CASE ( 2,3 ) 
          DO j = 1,pGridSerial%c2cs(icg)%nCellMembs
            icg2 = pGridSerial%c2cs(icg)%cellMembs(j)
                                
            IF ( pGridSerial%sc2r(icg2) /= pRegion%iRegionGlobal ) THEN
              IF ( nCellsVirt > 0 ) THEN 
                CALL BinarySearchInteger(vc(1:nCellsVirt),nCellsVirt,icg2,iLoc)
              ELSE 
                iLoc = ELEMENT_NOT_FOUND
              END IF ! nCellsVirt

              IF ( iLoc == ELEMENT_NOT_FOUND ) THEN 
                nCellsVirt = nCellsVirt + 1
              
                vc(nCellsVirt) = icg2
                  
                IF ( nCellsVirt > 1 ) THEN    
                  CALL QuickSortInteger(vc(1:nCellsVirt),nCellsVirt)              
                END IF ! nCellsVirt 
              END IF ! iLoc
            END IF ! pGridSerial%sc2r                      
          END DO ! j      
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pRegionSerial%mixtInput%stencilDimensCells
    END DO ! i      
          
    DEALLOCATE(sc,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'sc')
    END IF ! global%error             

! ******************************************************************************
!   Loop over current list of virtual cells and add layers of cells. For the 
!   current scheme, two additional layers are required.
! ******************************************************************************

! ==============================================================================
!   Allocate temporary memory 
! ==============================================================================

    vcOldDim    = nCellsVirt
    vcNewDimMax = nCellsVirtMax

    ALLOCATE(vcOld(vcOldDim),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vcOld')
    END IF ! global%error

    DO i = 1,vcOldDim
      vcOld(i) = vc(i)
    END DO ! i

    ALLOCATE(vcNew(vcNewDimMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vcNew')
    END IF ! global%error

! ==============================================================================
!   Loop over layers
! ==============================================================================

    nLayers = 2

    DO iLayer = 1,nLayers
      vcNewDim = 0            
         
      DO i = 1,vcOldDim
        icg = vcOld(i)

        CALL RFLU_BuildC2CStencilWrapper(pRegionSerial,icg,CONSTR_NONE)
        
        SELECT CASE ( pRegionSerial%mixtInput%stencilDimensCells ) 
          CASE ( 1 ) 
            DO iDimens = 1,pRegionSerial%mixtInput%dimens
              iDir = XCOORD - 1 + iDimens 
  
              DO j = 1,pGridSerial%c2cs1D(iDir,icg)%nCellMembs
                icg2 = pGridSerial%c2cs1D(iDir,icg)%cellMembs(j)
                
                IF ( pGridSerial%sc2r(icg2) /= pRegion%iRegionGlobal ) THEN
                  CALL BinarySearchInteger(vc(1:nCellsVirt),nCellsVirt,icg2,iLoc)

                  IF ( iLoc == ELEMENT_NOT_FOUND ) THEN
                    IF ( nCellsVirt < nCellsVirtMax ) THEN  
                      nCellsVirt = nCellsVirt + 1              
                      vc(nCellsVirt) = icg2
                    ELSE 
                      CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'vc')          
                    END IF ! nCellsVirt

                    IF ( vcNewDim < vcNewDimMax ) THEN 
                      vcNewDim = vcNewDim + 1              
                      vcNew(vcNewDim) = icg2
                    ELSE
                      CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'vcNew')         
                    END IF ! vcNewDim 

                    CALL QuickSortInteger(vc(1:nCellsVirt),nCellsVirt) 
                  END IF ! iLoc
                END IF ! pGridSerial%sc2r                
              END DO ! j
            END DO ! iDimens
          CASE ( 2,3 )       
            DO j = 1,pGridSerial%c2cs(icg)%nCellMembs
              icg2 = pGridSerial%c2cs(icg)%cellMembs(j)

              IF ( pGridSerial%sc2r(icg2) /= pRegion%iRegionGlobal ) THEN
                CALL BinarySearchInteger(vc(1:nCellsVirt),nCellsVirt,icg2,iLoc)

                IF ( iLoc == ELEMENT_NOT_FOUND ) THEN
                  IF ( nCellsVirt < nCellsVirtMax ) THEN  
                    nCellsVirt = nCellsVirt + 1              
                    vc(nCellsVirt) = icg2
                  ELSE 
                    CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'vc')          
                  END IF ! nCellsVirt

                  IF ( vcNewDim < vcNewDimMax ) THEN 
                    vcNewDim = vcNewDim + 1              
                    vcNew(vcNewDim) = icg2
                  ELSE
                    CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'vcNew')         
                  END IF ! vcNewDim 

                  CALL QuickSortInteger(vc(1:nCellsVirt),nCellsVirt) 
                END IF ! iLoc
              END IF ! pGridSerial%sc2r
            END DO ! j
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pRegionSerial%mixtInput%stencilDimensCells                      
      END DO ! i
             
      DEALLOCATE(vcOld,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vcOld')
      END IF ! global%error
      
      IF ( iLayer /= nLayers ) THEN 
        vcOldDim = vcNewDim
      
        ALLOCATE(vcOld(vcOldDim),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vcOld')
        END IF ! global%error    
      
        DO i = 1,vcOldDim
          vcOld(i) = vcNew(i)
        END DO ! i
      END IF ! iLayer      
    END DO ! iLayer   

! ==============================================================================
!   Deallocate temporary memory 
! ==============================================================================

    DEALLOCATE(vcNew,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vcNew')
    END IF ! global%error
           
    CALL RFLU_DestroyC2CStencilWrapper(pRegionSerial)          
                          
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
            'Adding virtual cells for inviscid higher-order stencil done.' 
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_AddVirtualCellsInv2





! ******************************************************************************
!
! Purpose: Mark virtual cells in serial region (arising from SYPE boundaries)
!   that are to be assigned region index.
!
! Description:
!   pRegion		Pointer to region
!   c2cCSR		CSR data structure containing AV and VV face neighbours of 
!               each cell
!   c2cCSRInfo	Access data structure for c2cCSR
!   cellStatus	Cell status array
!
! Input:   
!
! Output:  
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_MarkVirtualCells(pRegion,c2cCSR,c2cCSRInfo,cellStatus)                                   

    USE ModSortSearch
    
    USE RFLU_ModCellFaceEdgeInfo

    IMPLICIT NONE
      
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, DIMENSION(:) :: cellStatus,c2cCSR
    INTEGER, DIMENSION(:,:) :: c2cCSRInfo           
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: c1,c1k,c2,c2k,i,icg,icl,ifg
    INTEGER, DIMENSION(1) :: iMaxLoc
    INTEGER, DIMENSION(6) :: iRegList,iRegCntr        
    TYPE(t_global), POINTER :: global     
    TYPE(t_grid), POINTER :: pGrid       
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PART_MarkVirtualCells',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   
! ******************************************************************************

	DO icl = 1,pGrid%nCellsTot
	  cellStatus(icl) = CELL_STATUS_NOT_ADDED
	END DO ! icl
				  	  
	DO ifg = pGrid%nFaces+1,pGrid%nFacesTot	
	  c1 = pGrid%f2c(1,ifg)
	  c2 = pGrid%f2c(2,ifg)
				  
	  c1k = RFLU_GetGlobalCellKind(global,pGrid,c1)
	  c2k = RFLU_GetGlobalCellKind(global,pGrid,c2)            
	
	  SELECT CASE ( RFLU_GetFaceKind(global,c1k,c2k) ) 
		CASE ( FACE_KIND_VV )	
		  IF ( (pGrid%sc2r(c2) /= CRAZY_VALUE_INT) .AND. &
			   (pGrid%sc2r(c1) == CRAZY_VALUE_INT) ) THEN 				 
			icl = c1 - pGrid%nCells

			i = c2cCSRInfo(1,icl) + c2cCSRInfo(3,icl)
			c2cCSR(i) = pGrid%sc2r(c2)
			c2cCSRInfo(3,icl) = c2cCSRInfo(3,icl) + 1
			
			cellStatus(c1) = CELL_STATUS_TO_BE_ADDED
		  END IF ! pGrid%sc2r			  
			 
		  IF ( (pGrid%sc2r(c1) /= CRAZY_VALUE_INT) .AND. & 
			   (pGrid%sc2r(c2) == CRAZY_VALUE_INT) ) THEN 
			icl = c2 - pGrid%nCells

			i = c2cCSRInfo(1,icl) + c2cCSRInfo(3,icl)
			c2cCSR(i) = pGrid%sc2r(c1)
			c2cCSRInfo(3,icl) = c2cCSRInfo(3,icl) + 1
			
			cellStatus(c2) = CELL_STATUS_TO_BE_ADDED
		  END IF ! pGrid%sc2r
	  END SELECT ! RFLU_GetFaceKind	            
	END DO ! ifg		  

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_MarkVirtualCells






! ******************************************************************************
!
! Purpose: Set region for virtual cells in serial region (arising from SYPE
!   boundaries).
!
! Description: None.
!
! Input:   
!   pRegion		Pointer to region
!   c2cCSR		CSR data structure containing AV and VV face neighbours of 
!               each cell
!   c2cCSRInfo	Access data structure for c2cCSR
!   cellStatus	Cell status array
!
! Output:  
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_PART_SetVirtualCellRegionIndex(pRegion,c2cCSR,c2cCSRInfo, &
                                                 cellStatus)                                   

    USE ModSortSearch

    IMPLICIT NONE
      
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, DIMENSION(:) :: cellStatus,c2cCSR
    INTEGER, DIMENSION(:,:) :: c2cCSRInfo           
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: i,iBeg,icg,icl,iEnd,iLoc,iReg,nIns
    INTEGER, DIMENSION(1) :: iMaxLoc
    INTEGER, DIMENSION(6) :: iRegList,iRegCntr        
    TYPE(t_global), POINTER :: global     
    TYPE(t_grid), POINTER :: pGrid      
      
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_PART_SetVirtualCellRegionIndex', &
         __FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   
! ******************************************************************************

	DO icg = pGrid%nCells+1,pGrid%nCellsTot ! Loop over virtual cells
	  IF ( cellStatus(icg) == CELL_STATUS_TO_BE_ADDED ) THEN 
		icl = icg - pGrid%nCells
		
		iBeg = c2cCSRInfo(1,icl)
		iEnd = c2cCSRInfo(1,icl) + c2cCSRInfo(3,icl) - 1
		
		nIns = 0
		iRegList(:) = 0  
		iRegCntr(:) = 0 ! Must be non-positive
		
		DO i = iBeg,iEnd ! Loop over cells in CSR list
		  iReg = c2cCSR(i) ! region of neighboring cell

		  IF ( nIns > 1 ) THEN 
			CALL BinarySearchInteger(iRegList,6,iReg,iLoc)
							
			IF ( iLoc == ELEMENT_NOT_FOUND ) THEN 
			  nIns = nIns + 1
			  iRegList(nIns) = iReg
			  CALL QuickSortIntegerInteger(iRegList,iRegCntr,nIns)
			ELSE 
			  iRegCntr(iLoc) = iRegCntr(iLoc) + 1
			END IF ! iloc
		  ELSE 
			nIns = 1
			iRegList(1) = iReg
			iRegCntr(1) = 1
		  END IF ! i
		END DO ! i
		
		iMaxLoc = MAXLOC(iRegCntr) ! Find most frequently appearing region
		pGrid%sc2r(icg) = iRegList(iMaxLoc(1))
	  END IF ! cellStatus
	END DO ! icg
                          
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PART_SetVirtualCellRegionIndex





! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModPartitionRegionUtils


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModPartitionRegionUtils.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.1.1.1  2012/09/26 06:49:00  haselbac
! Initial revision at ETH
!
! Revision 1.4  2008/12/06 08:43:56  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:10  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/04/11 16:13:25  mtcampbe
! Updating for Rocflu split from RocfluidMP (filename too long, needed 2nd line)
!
! Revision 1.1  2007/04/09 18:56:52  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.6  2007/03/06 18:10:09  haselbac
! Adapted to building 1d stencils for higher-order scheme
!
! Revision 1.5  2006/01/06 22:17:20  haselbac
! Adapted to name changes
!
! Revision 1.4  2005/10/27 19:21:26  haselbac
! Adapted to changes in stencil routine names
!
! Revision 1.3  2005/10/05 14:25:11  haselbac
! Adapted to changes in stencil modules
!
! Revision 1.2  2005/06/14 17:48:14  haselbac
! Adapted RFLU_AddVirtualCellsInv2 to adaptive memory allocation
!
! Revision 1.1  2005/04/15 15:09:13  haselbac
! Initial revision
!
! Revision 1.1  2005/01/17 19:47:47  haselbac
! Initial revision
!
! ******************************************************************************

