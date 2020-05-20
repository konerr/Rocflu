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
! Purpose: Clone grid.
!
! Description: None.
!
! Input:
!   pRegionSerial    	Pointer to serial region
!   pRegion     	Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_CloneGrid.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_CloneGrid(pRegionSerial,pRegion)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI
  
  USE RFLU_ModCellMapping, ONLY: RFLU_CreateCellMapping

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================
   
  TYPE(t_region), POINTER :: pRegion,pRegionSerial 
   
! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,icl,ifl,iPatch,iPatchSerial,iPatchSerialPeriodic1, &
             iPatchSerialPeriodic2, &
             ivg,nPatchSerialPeriodic,shiftDir
  INTEGER, DIMENSION(:), ALLOCATABLE :: ip2ipsMap
  REAL(RFREAL) :: shiftDist
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid,pGridSerial
  TYPE(t_patch), POINTER :: pPatch,pPatchSerial

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_CloneGrid.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global
  
  CALL RegisterFunction(global,'RFLU_CloneGrid',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Cloning grid...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                     pRegion%iRegionGlobal       
  END IF ! global%verbLevel

  pGrid       => pRegion%grid  
  pGridSerial => pRegionSerial%grid  

! ******************************************************************************
! Clone volume grid
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Cloning volume grid...'
  END IF ! global%verbLevel
    
! ==============================================================================
! Clone dimensions
! ==============================================================================

  pGrid%nVert    = pGridSerial%nVert
  pGrid%nVertTot = pGridSerial%nVertTot
  pGrid%nVertMax = pGridSerial%nVertMax  
  
  pGrid%nCells    = pGridSerial%nCells   
  pGrid%nCellsTot = pGridSerial%nCellsTot
  pGrid%nCellsMax = pGridSerial%nCellsMax  
  
  pGrid%nTets    = pGridSerial%nTets   
  pGrid%nTetsTot = pGridSerial%nTetsTot 
  pGrid%nTetsMax = pGridSerial%nTetsMax 
  
  pGrid%nHexs    = pGridSerial%nHexs   
  pGrid%nHexsTot = pGridSerial%nHexsTot   
  pGrid%nHexsMax = pGridSerial%nHexsMax          

  pGrid%nPris    = pGridSerial%nPris   
  pGrid%nPrisTot = pGridSerial%nPrisTot 
  pGrid%nPrisMax = pGridSerial%nPrisMax   

  pGrid%nPyrs    = pGridSerial%nPyrs   
  pGrid%nPyrsTot = pGridSerial%nPyrsTot
  pGrid%nPyrsMax = pGridSerial%nPyrsMax  

! ==============================================================================
! Allocate memory
! ==============================================================================

  ALLOCATE(pGrid%xyz(XCOORD:ZCOORD,pGrid%nVertMax),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%xyz')
  END IF ! global%error        

  IF ( pGrid%nTetsMax > 0 ) THEN 
    ALLOCATE(pGrid%tet2v(4,pGrid%nTetsMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%tet2v')
    END IF ! global%error      
  ELSE 
    NULLIFY(pGrid%tet2v)
  END IF ! pGrid%nTetsMax

  IF ( pGrid%nHexsMax > 0 ) THEN 
    ALLOCATE(pGrid%hex2v(8,pGrid%nHexsMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%hex2v')
    END IF ! global%error
  ELSE 
    NULLIFY(pGrid%hex2v)
  END IF ! pGrid%nHexsMax

  IF ( pGrid%nPrisMax > 0 ) THEN 
    ALLOCATE(pGrid%pri2v(6,pGrid%nPrisMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pri2v')
    END IF ! global%error
  ELSE 
    NULLIFY(pGrid%pri2v)
  END IF ! pGrid%nPrisMax    

  IF ( pGrid%nPyrsMax > 0 ) THEN 
    ALLOCATE(pGrid%pyr2v(5,pGrid%nPyrsMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pyr2v')
    END IF ! global%error
  ELSE 
    NULLIFY(pGrid%pyr2v)
  END IF ! pGrid%nPyrsMax 
  
! ==============================================================================
! Clone coordinates. Frist determine shifting direction and distance. 
! NOTE the shifting direction is given by periodic boundary conditions.
! ==============================================================================

  DO iPatchSerial = 1,pGridSerial%nPatches
    pPatchSerial => pRegionSerial%patches(iPatchSerial)
  
    IF ( pPatchSerial%bcType == BC_PERIODIC ) THEN 
      shiftDir = pPatchSerial%axisRelated
      
      EXIT
    END IF ! pPatchSerial%bcType
  END DO ! iPatchSerial

  shiftDist = MAXVAL(pGridSerial%xyz(shiftDir,1:pGridSerial%nVert)) & 
            - MINVAL(pGridSerial%xyz(shiftDir,1:pGridSerial%nVert))

  DO ivg = 1,pGrid%nVertTot
    pGrid%xyz(XCOORD,ivg) = pGridSerial%xyz(XCOORD,ivg)
    pGrid%xyz(YCOORD,ivg) = pGridSerial%xyz(YCOORD,ivg)
    pGrid%xyz(ZCOORD,ivg) = pGridSerial%xyz(ZCOORD,ivg)        

    pGrid%xyz(shiftDir,ivg) = pGrid%xyz(shiftDir,ivg) & 
                            + (pRegion%iRegionGlobal-1)*shiftDist
  END DO ! ivg

! ==============================================================================
! Clone connectivity
! ==============================================================================

  DO icl = 1,pGrid%nTetsTot
    pGrid%tet2v(1,icl) = pGridSerial%tet2v(1,icl)
    pGrid%tet2v(2,icl) = pGridSerial%tet2v(2,icl)
    pGrid%tet2v(3,icl) = pGridSerial%tet2v(3,icl)
    pGrid%tet2v(4,icl) = pGridSerial%tet2v(4,icl)            
  END DO ! icl
  
  DO icl = 1,pGrid%nHexsTot
    pGrid%hex2v(1,icl) = pGridSerial%hex2v(1,icl)
    pGrid%hex2v(2,icl) = pGridSerial%hex2v(2,icl)
    pGrid%hex2v(3,icl) = pGridSerial%hex2v(3,icl)
    pGrid%hex2v(4,icl) = pGridSerial%hex2v(4,icl)  
    pGrid%hex2v(5,icl) = pGridSerial%hex2v(5,icl)
    pGrid%hex2v(6,icl) = pGridSerial%hex2v(6,icl)
    pGrid%hex2v(7,icl) = pGridSerial%hex2v(7,icl)
    pGrid%hex2v(8,icl) = pGridSerial%hex2v(8,icl)                          
  END DO ! icl  
  
  DO icl = 1,pGrid%nPrisTot
    pGrid%pri2v(1,icl) = pGridSerial%pri2v(1,icl)
    pGrid%pri2v(2,icl) = pGridSerial%pri2v(2,icl)
    pGrid%pri2v(3,icl) = pGridSerial%pri2v(3,icl)
    pGrid%pri2v(4,icl) = pGridSerial%pri2v(4,icl)  
    pGrid%pri2v(5,icl) = pGridSerial%pri2v(5,icl) 
    pGrid%pri2v(6,icl) = pGridSerial%pri2v(6,icl)                            
  END DO ! icl
  
  DO icl = 1,pGrid%nPyrsTot
    pGrid%pyr2v(1,icl) = pGridSerial%pyr2v(1,icl)
    pGrid%pyr2v(2,icl) = pGridSerial%pyr2v(2,icl)
    pGrid%pyr2v(3,icl) = pGridSerial%pyr2v(3,icl)
    pGrid%pyr2v(4,icl) = pGridSerial%pyr2v(4,icl)  
    pGrid%pyr2v(5,icl) = pGridSerial%pyr2v(5,icl)                        
  END DO ! icl      

! ==============================================================================
! Clone cell mapping
! ==============================================================================

  CALL RFLU_CreateCellMapping(pRegion)
  
  DO icl = 1,pGrid%nTetsTot
    pGrid%tet2CellGlob(icl) = pGridSerial%tet2CellGlob(icl)
    pGrid%tet2CellGlob(icl) = pGridSerial%tet2CellGlob(icl)
    pGrid%tet2CellGlob(icl) = pGridSerial%tet2CellGlob(icl)
    pGrid%tet2CellGlob(icl) = pGridSerial%tet2CellGlob(icl)            
  END DO ! icl
  
  DO icl = 1,pGrid%nHexsTot
    pGrid%hex2CellGlob(icl) = pGridSerial%hex2CellGlob(icl)
    pGrid%hex2CellGlob(icl) = pGridSerial%hex2CellGlob(icl)
    pGrid%hex2CellGlob(icl) = pGridSerial%hex2CellGlob(icl)
    pGrid%hex2CellGlob(icl) = pGridSerial%hex2CellGlob(icl)  
    pGrid%hex2CellGlob(icl) = pGridSerial%hex2CellGlob(icl)
    pGrid%hex2CellGlob(icl) = pGridSerial%hex2CellGlob(icl)
    pGrid%hex2CellGlob(icl) = pGridSerial%hex2CellGlob(icl)
    pGrid%hex2CellGlob(icl) = pGridSerial%hex2CellGlob(icl)                          
  END DO ! icl  
  
  DO icl = 1,pGrid%nPrisTot
    pGrid%pri2CellGlob(icl) = pGridSerial%pri2CellGlob(icl)
    pGrid%pri2CellGlob(icl) = pGridSerial%pri2CellGlob(icl)
    pGrid%pri2CellGlob(icl) = pGridSerial%pri2CellGlob(icl)
    pGrid%pri2CellGlob(icl) = pGridSerial%pri2CellGlob(icl)  
    pGrid%pri2CellGlob(icl) = pGridSerial%pri2CellGlob(icl) 
    pGrid%pri2CellGlob(icl) = pGridSerial%pri2CellGlob(icl)                            
  END DO ! icl
  
  DO icl = 1,pGrid%nPyrsTot
    pGrid%pyr2CellGlob(icl) = pGridSerial%pyr2CellGlob(icl)
    pGrid%pyr2CellGlob(icl) = pGridSerial%pyr2CellGlob(icl)
    pGrid%pyr2CellGlob(icl) = pGridSerial%pyr2CellGlob(icl)
    pGrid%pyr2CellGlob(icl) = pGridSerial%pyr2CellGlob(icl)  
    pGrid%pyr2CellGlob(icl) = pGridSerial%pyr2CellGlob(icl)                        
  END DO ! icl   

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Cloning volume grid done.'
  END IF ! global%verbLevel

! ******************************************************************************
! Clone patches
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Cloning patches...'
  END IF ! global%verbLevel
  
! ==============================================================================
! Determine number of patches. NOTE first and last regions have one periodic 
! patch less, in-between regions do not have either of the periodic patches 
! ==============================================================================

  IF ( (pRegion%iRegionGlobal == 1) .OR. & 
       (pRegion%iRegionGlobal == global%nRegionsLocal) ) THEN
    pGrid%nPatches = pGridSerial%nPatches - 1
  ELSE 
    pGrid%nPatches = pGridSerial%nPatches - 2  
  END IF ! pRegion%iRegionGlobal 

! ==============================================================================
! Build mapping from serial patches to cloned patches
! ==============================================================================

  ALLOCATE(pRegion%patches(pGrid%nPatches),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%patches')
  END IF ! global%error  

  nPatchSerialPeriodic = 0

  DO iPatchSerial = 1,pGridSerial%nPatches
    pPatchSerial => pRegionSerial%patches(iPatchSerial)
    
    IF ( pPatchSerial%bcType == BC_PERIODIC ) THEN 
      nPatchSerialPeriodic = nPatchSerialPeriodic + 1
      
      IF ( nPatchSerialPeriodic == 1 ) THEN 
        iPatchSerialPeriodic1 = iPatchSerial
      ELSE IF ( nPatchSerialPeriodic == 2 ) THEN 
        iPatchSerialPeriodic2 = iPatchSerial      
      ELSE 
        CALL ErrorStop(global,ERR_CLONABILITY,__LINE__)    
      END IF ! nPatchSerialPeriodic                               
    END IF ! pPatchSerial%bcType                    
  END DO ! iPatchSerial

  ALLOCATE(ip2ipsMap(pGrid%nPatches),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'ip2ipsMap')
  END IF ! global%error 

  iPatch = 0

  IF ( pRegion%iRegionGlobal == 1 ) THEN 
    DO iPatchSerial = 1,pGridSerial%nPatches
      IF ( iPatchSerial /= iPatchSerialPeriodic2 ) THEN 
        iPatch = iPatch + 1
      
        ip2ipsMap(iPatch) = iPatchSerial                        
      END IF ! iPatchSerial
    END DO ! iPatchSerial
  ELSE IF ( pRegion%iRegionGlobal == global%nRegionsLocal ) THEN 
    DO iPatchSerial = 1,pGridSerial%nPatches
      IF ( iPatchSerial /= iPatchSerialPeriodic1 ) THEN 
        iPatch = iPatch + 1
      
        ip2ipsMap(iPatch) = iPatchSerial 
      END IF ! iPatchSerial
    END DO ! iPatchSerial   
  ELSE 
    DO iPatchSerial = 1,pGridSerial%nPatches
      IF ( (iPatchSerial /= iPatchSerialPeriodic1) .AND. & 
           (iPatchSerial /= iPatchSerialPeriodic2) ) THEN
        iPatch = iPatch + 1
      
        ip2ipsMap(iPatch) = iPatchSerial 
      END IF ! iPatchSerial
    END DO ! iPatchSerial    
  END IF ! pRegion%iRegionGlobal

! ==============================================================================
! Copy patch data
! ==============================================================================

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
    
    iPatchSerial = ip2ipsMap(iPatch)
    
    pPatchSerial => pRegionSerial%patches(iPatchSerial)      
    
    pPatch%iPatchGlobal = iPatchSerial
    
    pPatch%nBTris    = pPatchSerial%nBTris
    pPatch%nBTrisTot = pPatchSerial%nBTrisTot
    pPatch%nBTrisMax = pPatchSerial%nBTrisMax  
    
    pPatch%nBQuads    = pPatchSerial%nBQuads
    pPatch%nBQuadsTot = pPatchSerial%nBQuadsTot
    pPatch%nBQuadsMax = pPatchSerial%nBQuadsMax            
    
    pPatch%nBQuads    = pPatchSerial%nBQuads
    pPatch%nBQuadsTot = pPatchSerial%nBQuadsTot
    pPatch%nBQuadsMax = pPatchSerial%nBQuadsMax
    
    pPatch%nBCellsVirt = pPatchSerial%nBCellsVirt
    
    pPatch%bcCoupled    = pPatchSerial%bcCoupled
    pPatch%movePatchDir = pPatchSerial%movePatchDir
    pPatch%flatFlag     = pPatchSerial%flatFlag    
    
    IF ( pPatch%nBTrisTot > 0 ) THEN 
      ALLOCATE(pPatch%bTri2v(3,pPatch%nBTrisTot),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bTri2v')
      END IF ! global%error   
    ELSE 
      NULLIFY(pPatch%bTri2v)
    END IF ! pPatch%nBTrisTot

    IF ( pPatch%nBQuadsTot > 0 ) THEN 
      ALLOCATE(pPatch%bQuad2v(4,pPatch%nBQuadsTot),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bQuad2v')
      END IF ! global%error   
    ELSE 
      NULLIFY(pPatch%bQuad2v)
    END IF ! pPatch%nBQuadTot
    
    IF ( pPatch%nBCellsVirt > 0 ) THEN 
      ALLOCATE(pPatch%bvc(pPatch%nBCellsVirt),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bvc')
      END IF ! global%error   
    ELSE 
      NULLIFY(pPatch%bvc)
    END IF ! pPatch%nBCellsVirt   
    
    DO ifl = 1,pPatch%nBTrisTot
      pPatch%bTri2v(1,ifl) = pPatchSerial%bTri2v(1,ifl)
      pPatch%bTri2v(2,ifl) = pPatchSerial%bTri2v(2,ifl)
      pPatch%bTri2v(3,ifl) = pPatchSerial%bTri2v(3,ifl)            
    END DO ! ifl
    
    DO ifl = 1,pPatch%nBQuadsTot
      pPatch%bQuad2v(1,ifl) = pPatchSerial%bQuad2v(1,ifl)
      pPatch%bQuad2v(2,ifl) = pPatchSerial%bQuad2v(2,ifl)
      pPatch%bQuad2v(3,ifl) = pPatchSerial%bQuad2v(3,ifl)   
      pPatch%bQuad2v(4,ifl) = pPatchSerial%bQuad2v(4,ifl)                     
    END DO ! ifl 
    
    DO icl = 1,pPatch%nBCellsVirt
      pPatch%bvc(icl) = pPatchSerial%bvc(icl)                    
    END DO ! icl        
  END DO ! iPatchSerial

! ==============================================================================
! Destroy mapping 
! ==============================================================================

  DEALLOCATE(ip2ipsMap,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'ip2ipsMap')
  END IF ! global%error 
  
  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Cloning patches done.'
  END IF ! global%verbLevel

! ******************************************************************************
! End
! ******************************************************************************
 
  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Cloning grid done.'
  END IF ! global%verbLevel 
 
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CloneGrid


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CloneGrid.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:54  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:06  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/08/07 17:13:59  haselbac
! Initial revision
!
! ******************************************************************************

