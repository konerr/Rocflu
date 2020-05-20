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
! Purpose: Initialize particle solution in serial region for 2d cases.
!
! Description: None.
!
! Input: 
!   pRegion	Pointer to serial region
!
! Output: None.
!
! Notes: 
!   1. Assume grid spacing is constant!
!   2. Assume cells are numbered in ascending order in direction of increasing
!      x-coordinate, then in direction of increasing y-coordinate.
!
! ******************************************************************************
!
! $Id: PLAG_RFLU_InitSolSerial_2D_Subbu_01Mar2014.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_InitSolSerial_2D(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY: t_plag
  USE ModParameters

  USE RFLU_ModFaceList, ONLY: RFLU_CreateCell2FaceList, & 
                              RFLU_BuildCell2FaceList, & 
                              RFLU_DestroyCell2FaceList    
  USE RFLU_ModGeometryTools, ONLY: RFLU_GetGridSpacingCartesian
  USE RFLU_ModInCellTest, ONLY: RFLU_ICT_TestInCell
  USE RFLU_ModTopologyUtils, ONLY: RFLU_GetNCellsCartesian

  USE PLAG_ModParameters    
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  ! Subbu - Vars Pcl Init
  INTEGER :: iRow,iDep,iCol,layer,nCellsCharge,nCellsSecI,Sgn
  REAL(RFREAL) :: nTol,rad,radLoc,the,theLoc,theLocActCell,Tol,z
  ! Subbu - End Vars Pcl Init

  LOGICAL :: foundFlag
  CHARACTER(CHRLEN) :: errorString,RCSIdentString
  INTEGER :: errorFlag,icg,icl,iPcl,ix,ix0,iy,iy0,nCellsX
  INTEGER, DIMENSION(2,8) :: dIcl 
  REAL(RFREAL) :: dx,dy,idx,idy,xMax,xMin,xPcl,yMin,yPcl,zPcl
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag

! TEMPORARY: Manoj: Guiding to find correct cell for particles  
  INTEGER :: icgOffset
! END TEMPORARY

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: PLAG_RFLU_InitSolSerial_2D_Subbu_01Mar2014.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_InitSolSerial_2D', &
                        __FILE__)
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing particle solution '// & 
                             'for serial region in 2d...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                     pRegion%iRegionGlobal
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and values
! ******************************************************************************
  
  pGrid => pRegion%grid
  pPlag => pRegion%plag  

  dIcl(1,1) = -1 ! Set local cell increments for alternative search
  dIcl(2,1) = -1
  dIcl(1,2) =  0
  dIcl(2,2) = -1
  dIcl(1,3) =  1
  dIcl(2,3) = -1
  dIcl(1,4) = -1
  dIcl(2,4) =  0
  dIcl(1,5) =  1
  dIcl(2,5) =  0
  dIcl(1,6) = -1
  dIcl(2,6) =  1
  dIcl(1,7) =  0
  dIcl(2,7) =  1
  dIcl(1,8) =  1
  dIcl(2,8) =  1

! ******************************************************************************
! Build cell-to-face list (needed for in-cell test and getting number of cells)
! ******************************************************************************

  CALL RFLU_CreateCell2FaceList(pRegion)
  CALL RFLU_BuildCell2FaceList(pRegion)

! ******************************************************************************
! Compute (inverse) grid spacings and number of cells in x-direction
! ******************************************************************************

  xMin = MINVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
  yMin = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))

  CALL RFLU_GetGridSpacingCartesian(pRegion,XCOORD,dx)
  CALL RFLU_GetGridSpacingCartesian(pRegion,YCOORD,dy)

  idx = 1.0_RFREAL/dx
  idy = 1.0_RFREAL/dy

 Tol = 1.0E-10_RFREAL
 ! Note: pRegion is still pRegionSerial here
 !pGrid%Jmax = pGrid%nCells/pRegion%patches(3)%nBQuads
 !pGrid%Kmax = pGrid%nCells/pRegion%patches(1)%nBQuads
 !pGrid%Imax = pGrid%nCells/(pGrid%Jmax * pGrid%Kmax)

 pGrid%Jmax = pGrid%nCells/pRegion%patches(2)%nBQuads
 pGrid%Imax = pRegion%patches(1)%nBQuads/pGrid%Jmax
 nCellsCharge = (pGrid%Imax/4)*(pGrid%Imax/4) ! No of cells within charge for every z-layer
 pGrid%Kmax = (pRegion%patches(2)%nBQuads - nCellsCharge)/ pGrid%Imax


 icg = 0
 DO layer = 1,pGrid%Kmax
   icg = (layer-1)*(pGrid%Imax*pGrid%Jmax) + 1
   rad = SQRT(pGrid%cofg(XCOORD,icg)**2.0_RFREAL + pGrid%cofg(YCOORD,icg)**2.0_RFREAL)
   IF (ABS(rad - pRegion%plagInput%iniRandYMin) .LE. Tol) THEN
     pGrid%radMin = rad 
     pGrid%theMin = ATAN2(pGrid%cofg(YCOORD,icg),pGrid%cofg(XCOORD,icg))
     pGrid%SgntheMin = INT(DSIGN(1.0_RFREAL,pGrid%theMin))
     pGrid%zMin = pGrid%cofg(ZCOORD,icg)
     pGrid%dz = 2.0_RFREAL*pGrid%zMin
     icg = icg + 1 
     the = ATAN2(pGrid%cofg(YCOORD,icg),pGrid%cofg(XCOORD,icg))
     pGrid%dthe = the - pGrid%theMin
     icg = icg-1
     icg = icg + (pGrid%Imax*pGrid%Jmax)
     rad = SQRT(pGrid%cofg(XCOORD,icg)**2.0_RFREAL + pGrid%cofg(YCOORD,icg)**2.0_RFREAL)
     pGrid%drad = rad - pGrid%radMin
     !== Finding number of cells in Sec I ==!
     !nCellsSecI = (layer-1)*pGrid%Imax*pGrid%Jmax
   EXIT
   END IF
 END DO

! Subbu - Initialize particle in a given cell
! ******************************************************************************
! Loop over particles in serial region
! ******************************************************************************
  nTol = 1.0E-14_RFREAL
  DO iPcl = 1,pPlag%nPcls    

! ==============================================================================
!   Get particle location and determine cell which should contain this particle
! ==============================================================================

    xPcl = pPlag%cv(CV_PLAG_XPOS,iPcl)
    yPcl = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zPcl = pPlag%cv(CV_PLAG_ZPOS,iPcl)

    radLoc = DSQRT(xPcl**2.0_RFREAL + yPcl**2.0_RFREAL)
    theLoc = DATAN2(yPcl,xPcl)

    iDep = FLOOR(DABS((zPcl - (pGrid%zMin-0.5_RFREAL*pGrid%dz))/(pGrid%dz)) + 1)
    iRow = FLOOR(DABS((radLoc -(pGrid%radMin-0.5_RFREAL*pGrid%drad))/(pGrid%drad)) + 1) 

    !WRITE(*,*) 'iPcl,iRow,the=',iPcl,iRow,theLoc*180/global%pi, &
    !                            (pGrid%theMin-0.5_RFREAL*pGrid%dthe)*180/global%pi,pGrid%dthe*180/global%pi
    Sgn = INT(DSIGN(1.0_RFREAL,theLoc))
    IF (Sgn /= pGrid%SgntheMin) & !THEN
        theLoc = theLoc + pGrid%SgntheMin*2.0_RFREAL*global%pi
    IF ((DABS(theLoc) - DABS(pGrid%theMin-0.5_RFREAL*pGrid%dthe)) .LE. nTol) &
        theLoc = theLoc + pGrid%SgntheMin*2.0_RFREAL*global%pi
    iCol = FLOOR(DABS((theLoc -(pGrid%theMin-0.5_RFREAL*pGrid%dthe))/(pGrid%dthe)) + 1)
    icg  = (iDep-1)*(pGrid%Imax*pGrid%Kmax + nCellsCharge) + (iRow-1)*pGrid%Imax + iCol
    !icg = icg + nCellsSecI
    pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg
    pPlag%aiv(AIV_PLAG_REGINI,iPcl) = pRegion%iRegionGlobal
    !WRITE(*,*) 'iPcl,iRow,iCol,the,icg=',iPcl,iRow,iCol,theLoc,pGrid%theMin-0.5_RFREAL*pGrid%dthe,pGrid%dthe
  END DO ! iPcl
   
! Subbu - End Initialize particle in a given cell

! ******************************************************************************
! Destroy cell-to-face list
! ******************************************************************************

  CALL RFLU_DestroyCell2FaceList(pRegion)        

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Initializing particle solution for serial region in 2d done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_InitSolSerial_2D

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_InitSolSerial_2D_Subbu_01Mar2014.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.5  2008/12/06 08:43:51  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2008/04/05 18:26:30  haselbac
! Made determination of cell containing particle more general
!
! Revision 1.2  2007/04/26 20:23:44  haselbac
! Added search of nearby cells, fixed bug in RegisterFunction call
!
! Revision 1.1  2007/04/15 02:33:23  haselbac
! Initial revision
!
! ******************************************************************************

