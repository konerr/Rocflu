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
! Purpose: Build stencil of four cells so can use bilinear interpolation to 
!   get solution at intersection point.
!
! Description: None.
!
! Input: 
!   pRegion	Pointer to region
!   ifg		Index of face which was intersected
!   xLoc	x-location of intersection
!   yLoc	y-location of intersection
!   zLoc	z-location of intersection
!
! Output: 
!   cs		Cell stencil
!
! Notes: 
!   1. Only applicable to structured-like quadrilateral grids in 2d.
!   2. Only applicable to 2d grids in x-y plane.
!
! ******************************************************************************
!
! $Id: RFLU_BuildStencilInterp.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_BuildStencilInterp(pRegion,ifg,xLoc,yLoc,zLoc,cs)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  USE ModSortSearch
      
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Parameters
! ==============================================================================

  INTEGER, INTENT(IN) :: ifg
  INTEGER, DIMENSION(4), INTENT(OUT) :: cs
  REAL(RFREAL), INTENT(IN) :: xLoc,yLoc,zLoc
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: c1,c2,errorFlag,icl,icl2,icg,icg2,ict,ifg2,ifl,iLoc,iPatch, & 
             nCsTempMax,nCsTempMerged,nCsTempMergedMax,nFaces,nFacesPool
  INTEGER, DIMENSION(2) :: iFacePool,nCsTemp
  INTEGER, DIMENSION(4) :: key
  INTEGER, DIMENSION(:), ALLOCATABLE :: csTempMerged
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: csTemp
  INTEGER, DIMENSION(:,:), POINTER :: pC2f
  REAL(RFREAL) :: fcx,fcy,xMax,xMin,yMax,yMin
  REAL(RFREAL), DIMENSION(2) :: dist
  REAL(RFREAL), DIMENSION(4) :: sc,xc,yc
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: RFLU_BuildStencilInterp.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BuildStencilInterp', &
                        __FILE__)

! ******************************************************************************
! Set pointers 
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  nCsTempMax       = 8
  nCsTempMergedMax = 8

  nCsTemp(1:2) = 0

  ALLOCATE(csTemp(2,nCsTempMax),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'csTemp')
  END IF ! global%error
  
  ALLOCATE(csTempMerged(nCsTempMergedMax),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'csTempMerged')
  END IF ! global%error

! ******************************************************************************
! Loop over cells adjacent to intersected face and store face-neighbors of each
! cell in sorted list. NOTE this part can be generalized, but at the moment only
! applicable to structured-like quadrilateral grids
! ******************************************************************************

  DO icl = 1,2
    icg = pGrid%f2c(icl,ifg)

    ict  = pGrid%cellGlob2Loc(1,icg) ! cell type
    icl2 = pGrid%cellGlob2Loc(2,icg) ! local cell index

    SELECT CASE ( ict )
      CASE ( CELL_TYPE_HEX )
        pC2f => pGrid%hex2f(:,:,icl2)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! ict

    nFaces = SIZE(pC2f,2)

    DO ifl = 1,nFaces
      iPatch = pC2f(1,ifl)
      ifg2   = pC2f(2,ifl)

      IF ( iPatch == 0 ) THEN 
        c1 = pGrid%f2c(1,ifg2)
        c2 = pGrid%f2c(2,ifg2)

        IF ( c1 == icg ) THEN 
          icg2 = c2
        ELSE 
          icg2 = c1
        END IF ! c1

        IF ( icg2 /= pGrid%f2c(3-icl,ifg) ) THEN
          nCsTemp(icl) = nCsTemp(icl) + 1

          csTemp(icl,nCsTemp(icl)) = icg2
        END IF ! icg2
      END IF ! iPatch
    END DO ! ifl

    CALL QuickSortInteger(csTemp(icl,1:nCsTemp(icl)),nCsTemp(icl))
  END DO ! icl 

! ******************************************************************************
! Merge sorted list of neighboring cells and deallocate temporary memory
! ******************************************************************************

  CALL MergeSortedIntegers(global,nCsTemp(1),nCsTemp(2), & 
                           csTemp(1,1:nCsTemp(1)),csTemp(2,1:nCsTemp(2)), & 
                           nCsTempMergedMax,nCsTempMerged,csTempMerged)

  DEALLOCATE(csTemp,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'csTemp')
  END IF ! global%error

! ******************************************************************************
! Loop over cells in merged list and extract all faces in grid which have both
! of its adjacent cells in merged list, they are put in a pool of faces. NOTE 
! with the restriction to structured-like quadrilateral cells, can have only 
! two faces in pool.
! ******************************************************************************

  nFacesPool = 0

  iFacePool(1:2) = 0

  DO ifl = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifl)
    c2 = pGrid%f2c(2,ifl)

    IF ( (MIN(c1,c2) >= csTempMerged(1)) .AND. & 
         (MAX(c1,c2) <= csTempMerged(nCsTempMerged)) ) THEN 
      CALL BinarySearchInteger(csTempMerged(1:nCsTempMerged), &
                               nCsTempMerged,c1,iLoc)  
      
      IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
        CALL BinarySearchInteger(csTempMerged(1:nCsTempMerged), &
                                 nCsTempMerged,c2,iLoc)

        IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
          IF ( nFacesPool < 2 ) THEN
            nFacesPool = nFacesPool + 1
          ELSE 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END IF ! nFacesPool

          iFacePool(nFacesPool) = ifl
        END IF ! iLoc
      END IF ! iLoc
    END IF ! MIN 
  END DO ! ifl
  
  DEALLOCATE(csTempMerged,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'csTempMerged')
  END IF ! global%error

! ******************************************************************************
! Sort faces in pool according to their distance from the intersection point,
! then pick face with smaller distance
! ******************************************************************************

  dist(1:2) = HUGE(1.0_RFREAL)

  DO ifl = 1,nFacesPool
    ifg2 = iFacePool(ifl)

    fcx = pGrid%fc(XCOORD,ifg2)
    fcy = pGrid%fc(YCOORD,ifg2)

    dist(ifl) = (fcx-xLoc)**2 + (fcy-yLoc)**2
  END DO ! ifl 

  IF ( dist(1) < dist(2) ) THEN 
    ifg2 = iFacePool(1) 
  ELSE 
    ifg2 = iFacePool(2)
  END IF ! dist

! ******************************************************************************
! Build list of cells which belong to intersected face and other close-by face.
! Sort according to function which allows to order cells into increasing order
! of x and y so can be addressed with (i,j) notation
! ******************************************************************************

  cs(1) = pGrid%f2c(1,ifg)
  cs(2) = pGrid%f2c(2,ifg)
  cs(3) = pGrid%f2c(1,ifg2)
  cs(4) = pGrid%f2c(2,ifg2)

  DO icl = 1,4
    xc(icl)  = pGrid%cofg(XCOORD,cs(icl))
    yc(icl)  = pGrid%cofg(YCOORD,cs(icl))
    key(icl) = icl
  END DO ! icl

  xMin = MINVAL(xc(1:4))
  xMax = MAXVAL(xc(1:4))
  yMin = MINVAL(yc(1:4))
  yMax = MAXVAL(yc(1:4))

  DO icl = 1,4
    xc(icl)  = (xc(icl)-xMin)/(xMax-xMin)
    yc(icl)  = (yc(icl)-yMin)/(yMax-yMin)
    sc(icl)  = xc(icl) + 10.0_RFREAL*yc(icl)
  END DO ! icl

  CALL QuickSortRFREALInteger(sc,cs,4)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BuildStencilInterp

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_BuildStencilInterp.F90,v $
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
! ******************************************************************************

