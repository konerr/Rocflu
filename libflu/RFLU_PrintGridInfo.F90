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
! Purpose: Display information on grid.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
! 
! Output: N/A.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_PrintGridInfo.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2000-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PrintGridInfo(pRegion)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError

  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region   

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Parameters
! ==============================================================================   

  TYPE(t_region), POINTER :: pRegion 

! ==============================================================================
! Local variables
! ==============================================================================
 
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iPatch,nBQuads,nBQuadsTot,nBTris,nBTrisTot
  INTEGER :: dummy(1)  
  INTEGER :: loc(XCOORD:ZCOORD,MIN_VAL:MAX_VAL)
  TYPE(t_grid), POINTER :: pGrid  
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global
  
  RCSIdentString = '$RCSfile: RFLU_PrintGridInfo.F90,v $ $Revision: 1.1.1.1 $'

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_PrintGridInfo',__FILE__)

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing grid information...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                     pRegion%iRegionGlobal 
    IF ( global%flowType == FLOW_UNSTEADY ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                          global%currentTime
    END IF ! global%flowType    
  END IF ! global%verbLevel

! ==============================================================================
! Set pointers
! ==============================================================================

  pGrid => pRegion%grid

! ==============================================================================
! Compute patch dimensions here to be independent
! ==============================================================================

  nBTris     = 0
  nBTrisTot  = 0
  nBQuads    = 0
  nBQuadsTot = 0
  
  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
    
    nBTris     = nBTris     + pPatch%nBTris 
    nBTrisTot  = nBTrisTot  + pPatch%nBTrisTot     
    nBQuads    = nBQuads    + pPatch%nBQuads 
    nBQuadsTot = nBQuadsTot + pPatch%nBQuadsTot            
  END DO ! iPatch  

! ==============================================================================
! Grid dimensions
! ==============================================================================

  WRITE(STDOUT,'(A,3X,A)')                SOLVER_NAME,'Grid statistics:'
  WRITE(STDOUT,'(A,5X,A,3X,2(1X,I9))')    SOLVER_NAME,'Vertices:      ', &
        pGrid%nVert,pGrid%nVertTot-pGrid%nVert
  WRITE(STDOUT,'(A,5X,A,3X,2(1X,I9))')    SOLVER_NAME,'Cells:         ', &
        pGrid%nCells,pGrid%nCellsTot-pGrid%nCells
  WRITE(STDOUT,'(A,7X,A,1X,2(1X,I9))')    SOLVER_NAME,'Tetrahedra:    ', &
        pGrid%nTets,pGrid%nTetsTot-pGrid%nTets
  WRITE(STDOUT,'(A,7X,A,1X,2(1X,I9))')    SOLVER_NAME,'Hexahedra:     ', &
        pGrid%nHexs,pGrid%nHexsTot-pGrid%nHexs
  WRITE(STDOUT,'(A,7X,A,1X,2(1X,I9))')    SOLVER_NAME,'Prisms:        ', & 
        pGrid%nPris,pGrid%nPrisTot-pGrid%nPris
  WRITE(STDOUT,'(A,7X,A,1X,2(1X,I9))')    SOLVER_NAME,'Pyramids:      ', & 
        pGrid%nPyrs,pGrid%nPyrsTot-pGrid%nPyrs
  WRITE(STDOUT,'(A,5X,A,4X,I9)')          SOLVER_NAME,'Patches:       ', & 
        pGrid%nPatches
  WRITE(STDOUT,'(A,5X,A,3X,2(1X,I9))')    SOLVER_NAME,'Patch faces:   ', & 
        nBTris+nBQuads,nBTrisTot+nBQuadsTot-nBTris-nBQuads
  WRITE(STDOUT,'(A,7X,A,1X,2(1X,I9))')    SOLVER_NAME,'Triangles:     ', & 
        nBTris,nBTrisTot-nBTris
  WRITE(STDOUT,'(A,7X,A,1X,2(1X,I9))')    SOLVER_NAME,'Quadrilaterals:', &
        nBQuads,nBQuadsTot-nBQuads
  
! ==============================================================================
! Patch dimensions
! ==============================================================================  

  IF ( pGrid%nPatches > 0 ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Patch statistics:'

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      WRITE(STDOUT,'(A,5X,A,2X,I4)')       SOLVER_NAME,'Patch:',iPatch   
      WRITE(STDOUT,'(A,7X,A,6X,2(1X,I9))') SOLVER_NAME,'Triangles:', & 
            pPatch%nBTris,pPatch%nBTrisTot-pPatch%nBTris
      WRITE(STDOUT,'(A,7X,A,1X,2(1X,I9))') SOLVER_NAME,'Quadrilaterals:', &
            pPatch%nBQuads,pPatch%nBQuadsTot-pPatch%nBQuads    
    END DO ! iPatch
  END IF ! pGrid%nPatches

! ==============================================================================
! Find locations of extrema and print information on extrema: NOTE Asinine 
! coding needed because of poor FORTRAN interface for MINLOC and MAXLOC
! functions... 
! ==============================================================================

  dummy               = MINLOC(pGrid%xyz(XCOORD,1:pGrid%nVert))
  loc(XCOORD,MIN_VAL) = dummy(1)

  dummy               = MINLOC(pGrid%xyz(YCOORD,1:pGrid%nVert))
  loc(YCOORD,MIN_VAL) = dummy(1)

  dummy               = MINLOC(pGrid%xyz(ZCOORD,1:pGrid%nVert))
  loc(ZCOORD,MIN_VAL) = dummy(1)


  dummy               = MAXLOC(pGrid%xyz(XCOORD,1:pGrid%nVert))
  loc(XCOORD,MAX_VAL) = dummy(1)

  dummy               = MAXLOC(pGrid%xyz(YCOORD,1:pGrid%nVert))
  loc(YCOORD,MAX_VAL) = dummy(1)

  dummy               = MAXLOC(pGrid%xyz(ZCOORD,1:pGrid%nVert))
  loc(ZCOORD,MAX_VAL) = dummy(1)


  WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Coordinate extrema:'
  WRITE(STDOUT,'(A,5X,A,2(1X,E23.16),2(1X,I9))') & 
                SOLVER_NAME,'X-coordinate:', & 
                MINVAL(pGrid%xyz(XCOORD,1:pGrid%nVert)), & 
                MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVert)), & 
                loc(XCOORD,MIN_VAL),loc(XCOORD,MAX_VAL)
  WRITE(STDOUT,'(A,5X,A,2(1X,E23.16),2(1X,I9))') & 
                SOLVER_NAME,'Y-coordinate:', & 
                MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVert)), & 
                MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVert)), & 
                loc(YCOORD,MIN_VAL),loc(YCOORD,MAX_VAL)
  WRITE(STDOUT,'(A,5X,A,2(1X,E23.16),2(1X,I9))') & 
                SOLVER_NAME,'Z-coordinate:', & 
                MINVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert)), & 
                MAXVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert)), & 
                loc(ZCOORD,MIN_VAL),loc(ZCOORD,MAX_VAL)
                  
! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing grid information done.'
  END IF ! global%verbLevel    
  
  CALL DeregisterFunction(global)  

END SUBROUTINE RFLU_PrintGridInfo

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PrintGridInfo.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:35  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:53  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.13  2004/10/19 19:24:45  haselbac
! Removed printing of derived dims, cosmetics
!                                    
! Revision 1.12  2004/01/22 16:04:48  haselbac                        
! Changed declaration to eliminate warning on ALC                     
!
! Revision 1.11  2003/03/15 17:03:38  haselbac                        
! Made changes for parallel calcs, added bVert info                   
!
! Revision 1.10  2003/02/01 00:27:44  haselbac                        
! Increased precision for printing coordinate extrema                 
!
! Revision 1.9  2003/01/28 15:39:16  haselbac                         
! Compute total no of boundary faces (deleted from pGrid), cosmetics  
!
! Revision 1.8  2002/10/27 18:52:56  haselbac                         
! Cosmetic changes, added iRegionGlobal                               
!
! Revision 1.7  2002/10/17 19:57:52  haselbac                         
! Cosmetic changes to output                                          
!
! Revision 1.6  2002/09/09 14:15:01  haselbac                         
! global now under regions                                            
!
! Revision 1.5  2002/06/27 15:38:21  haselbac                         
! Added dummy cells and vertex output                                 
!
! Revision 1.4  2002/06/17 13:31:22  haselbac                         
! Prefixed SOLVER_NAME to all screen output                           
!
! Revision 1.3  2002/06/14 20:09:28  haselbac                         
! Added writing out of iPatchGlobal                                   
!
! Revision 1.2  2002/06/11 22:27:43  haselbac                         
! Added writing of patch statistics                                   
!
! Revision 1.1  2002/06/05 18:34:12  haselbac                         
! Initial revision                                                    
!
! ******************************************************************************

