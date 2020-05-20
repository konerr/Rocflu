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
! Purpose: Compute number of cells containing particles
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
!
! $Id: PLAG_RFLU_ComputeCellsContainingPcls.F90,v 1.2 2015/07/27 04:45:42 brollin Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_ComputeCellsContainingPcls(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModParameters
  USE ModMPI
  
  USE ModRandom, ONLY: Rand1Uniform
  
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

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,icl,rMinicl(1),rMaxicl(1)
  REAL(RFREAL) :: rad,radMin,radMax,rMinReg,rMaxReg,theMin,theMax,tol,zMin,zMax
  TYPE(t_global), POINTER :: global
  TYPE(t_grid),   POINTER :: pGrid
!******************************************************************************
!Saptarshi Locals
!******************************************************************************
 INTEGER ::xMinicl(1),xMaxicl(1),yMinicl(1),yMaxicl(1),zMinicl(1),zMaxicl(1)
 REAL(RFREAL) :: xMin,xMax,yMin,yMax,X,Y,dx,xMinReg,xMaxReg,yMinReg,yMaxReg
 REAL(RFREAL) :: Z,zMinReg,zMaxReg
!******************************************************************************
!Saptarshi Ends
!******************************************************************************

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: PLAG_RFLU_ComputeCellsContainingPcls.F90,v $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_ComputeCellsContainingPcls', &
                        __FILE__)
 
  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
       WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal          
       WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                          'Computing number of cells that can be seeded with &
                           particles ...'
  END IF ! global%verbLevel
  
! ******************************************************************************
! Set pointers and values
! ******************************************************************************

  pGrid => pRegion%grid

  SELECT CASE ( global%casename )

  CASE("cyldet")

! ******************************************************************************  
! Build bounding box. NOTE not from grid, but from user input
! ******************************************************************************

  theMin = pRegion%plagInput%iniRandXMin
  theMax = pRegion%plagInput%iniRandXMax
  radMin = pRegion%plagInput%iniRandYMin
  radMax = pRegion%plagInput%iniRandYMax
  zMin = pRegion%plagInput%iniRandZMin
  zMax = pRegion%plagInput%iniRandZMax

  IF (pRegion%mixtInput%dimens == 2) THEN
    icl = 1 ! Could be any value between 1 and pGrid%nCells
    zMin = pGrid%cofg(ZCOORD,icl)
    zMax = pGrid%cofg(ZCOORD,icl)
  END IF

! ******************************************************************************
! Compute how many cells can be seeded with particles given radMin,radMax etc.
! from input file.
! Based on the above computation set the logical var initPclPresent to true or
! false
! ******************************************************************************

  ! Push the following two lines into routine where vars are initialized
  pGrid%initPclPresent = .FALSE.
  pGrid%nCellsSeedPcls = 0

  tol = 1.0E-08
 
! ******************************************************************************
! Allocate memory to store a random number in each of the cells 
! If the min and max radii of the entire region fall within the range specified 
! in the input file generate a random number for all cells in region
! Else if only a section of cells fall within the range, seed those cells with a
! a random number and seed the remaining cells with 0.
! If there are no pcls in a given region deallocate memory
! ******************************************************************************

  rMinicl = MINLOC(ABS(pGrid%cofg(YCOORD,:)))
  rMaxicl = MAXLOC(ABS(pGrid%cofg(YCOORD,:)))
  rMinReg = SQRT(pGrid%cofg(XCOORD,rMinicl(1))**2.0_RFREAL + &
                 pGrid%cofg(YCOORD,rMinicl(1))**2.0_RFREAL) 
  rMaxReg = SQRT(pGrid%cofg(XCOORD,rMaxicl(1))**2.0_RFREAL + &
                 pGrid%cofg(YCOORD,rMaxicl(1))**2.0_RFREAL)

! ------------------------ Allocating memory for all regions ------------------- 

  ALLOCATE(pGrid%nPclsPerCell(pGrid%nCells),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'nPclsPerCell')
  END IF ! global%error
 
  IF ( (radMin-rMinReg .LE. tol) .AND. &
       (rMaxReg-radMax .LE. tol) ) THEN 
     pGrid%nCellsSeedPcls = pGrid%nCells 
     DO icl = 1,pGrid%nCells ! random number will be [0,1)
       pGrid%nPclsPerCell(icl) = Rand1Uniform(pRegion%randData) 
     END DO
  ELSEIF ( (radMax-rMinReg .LE. tol) .OR. &
         (rMaxReg-radMin .LE. tol) ) THEN 
         pGrid%initPclPresent = .FALSE.
  ELSE
     DO icl = 1,pGrid%nCells
       rad = SQRT(pGrid%cofg(XCOORD,icl)**2.0_RFREAL + &
                 pGrid%cofg(YCOORD,icl)**2.0_RFREAL)
       IF ( (radMin-rad .LE. tol) .AND. &
            (rad-radMax .LE. tol) ) THEN 
          pGrid%nCellsSeedPcls = pGrid%nCellsSeedPcls + 1
          pGrid%nPclsPerCell(icl) = Rand1Uniform(pRegion%randData) 
       ELSE
          pGrid%nPclsPerCell(icl) = -1E6 ! Some large -ve number
       END IF
     END DO  
  END IF
  
  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A,1X,I16)') SOLVER_NAME, &
          'Number of cells that can be seeded with particles is', &
           pGrid%nCellsSeedPcls
  END IF ! global%verbLevel

  IF (pGrid%nCellsSeedPcls .GE. 1) THEN
    pGrid%initPclPresent = .TRUE.
  ELSE
! -------------------- Deallocating memory for regions w/o pcls ----------------
    DEALLOCATE(pGrid%nPclsPerCell,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'nPclsPerCell')
    END IF ! global%error
  END IF

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Computing number of cells that can be seeded with particles done.'
  END IF ! global%verbLevel

!*******************************************************************************
!Saptarshi 30th January 2015
!******************************************************************************
 CASE("shktb")

  xMin = pRegion%plagInput%iniRandXMin
  xMax = pRegion%plagInput%iniRandXMax
  yMin = pRegion%plagInput%iniRandYMin
  yMax = pRegion%plagInput%iniRandYMax
  zMin = pRegion%plagInput%iniRandZMin
  zMax = pRegion%plagInput%iniRandZMax

  IF (pRegion%mixtInput%dimens == 2) THEN
    icl = 1 ! Could be any value between 1 and pGrid%nCells
    zMin = pGrid%cofg(ZCOORD,icl)
    zMax = pGrid%cofg(ZCOORD,icl)
  END IF

! ******************************************************************************
! Compute how many cells can be seeded with particles given radMin,radMax etc.
! from input file.
! Based on the above computation set the logical var initPclPresent to true or
! false
! ******************************************************************************

  ! Push the following two lines into routine where vars are initialized
  pGrid%initPclPresent = .FALSE.
  pGrid%nCellsSeedPcls = 0

  tol = 1.0E-08

! ******************************************************************************
! Allocate memory to store a random number in each of the cells
! If the min and max length of the entire region fall within the range specified
! in the input file generate a random number for all cells in region
! Else if only a section of cells fall within the range, seed those cells with a
! a random number and seed the remaining cells with 0.
! If there are no pcls in a given region deallocate memory
! ******************************************************************************

  xMinicl = MINLOC(ABS(pGrid%cofg(XCOORD,:)))
  xMaxicl = MAXLOC(ABS(pGrid%cofg(XCOORD,:)))
  xMinReg = pGrid%cofg(XCOORD,xMinicl(1))
  xMaxReg = pGrid%cofg(XCOORD,xMaxicl(1))
  yMinicl = MINLOC(ABS(pGrid%cofg(YCOORD,:)))
  yMaxicl = MAXLOC(ABS(pGrid%cofg(YCOORD,:)))
  yMinReg = pGrid%cofg(YCOORD,yMinicl(1))
  yMaxReg = pGrid%cofg(YCOORD,yMaxicl(1))
  zMinicl = MINLOC(ABS(pGrid%cofg(ZCOORD,:)))
  zMaxicl = MAXLOC(ABS(pGrid%cofg(ZCOORD,:)))
  zMinReg = pGrid%cofg(ZCOORD,zMinicl(1))
  zMaxReg = pGrid%cofg(ZCOORD,zMaxicl(1))

! ------------------------ Allocating memory for all regions -------------------

! ===== Florian =====> Begins

  ALLOCATE(pGrid%nPclsPerCell(pGrid%nCells),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'nPclsPerCell')
  END IF ! global%error

! Addition of condition on Y and Z

  IF (pRegion%mixtInput%dimens == 2) THEN
    IF ((xMin-xMinReg .LE. tol) .AND. (xMaxReg-xMax .LE. tol) .AND. &
       (yMin-yMinReg .LE. tol) .AND. (yMaxReg-yMax .LE. tol)) THEN
       pGrid%nCellsSeedPcls = pGrid%nCells
       DO icl = 1,pGrid%nCells ! random number will be [0,1)
         pGrid%nPclsPerCell(icl) = Rand1Uniform(pRegion%randData)
       END DO
    ELSEIF ((xMax-xMinReg .LE. tol) .OR. (xMaxReg-xMin .LE. tol) .OR. &
           (yMax-yMinReg .LE. tol) .OR. (yMaxReg-yMin .LE. tol)) THEN
       pGrid%initPclPresent = .FALSE.
    ELSE
       DO icl = 1,pGrid%nCells
          X = pGrid%cofg(XCOORD,icl)
          Y = pGrid%cofg(YCOORD,icl)
          Z = pGrid%cofg(ZCOORD,icl)
          IF (((xMin-X < tol) .AND. (X-xMax < tol)) .AND. &
             ((yMin-Y < tol) .AND. (Y-yMax < tol))) THEN
             pGrid%nCellsSeedPcls = pGrid%nCellsSeedPcls + 1
             pGrid%nPclsPerCell(icl) = Rand1Uniform(pRegion%randData)
          ELSE
             pGrid%nPclsPerCell(icl) = -1E6 ! Some large -ve number
          END IF
       END DO
    END IF 

    ELSEIF (pRegion%mixtInput%dimens == 3) THEN
      IF ((xMin-xMinReg .LE. tol) .AND. (xMaxReg-xMax .LE. tol) .AND. &
         (yMin-yMinReg .LE. tol) .AND. (yMaxReg-yMax .LE. tol) .AND. &
         (zMin-zMinReg .LE. tol) .AND. (zMaxReg-zMax .LE. tol)) THEN
         pGrid%nCellsSeedPcls = pGrid%nCells
         DO icl = 1,pGrid%nCells ! random number will be [0,1)
            pGrid%nPclsPerCell(icl) = Rand1Uniform(pRegion%randData)
         END DO
      ELSEIF ((xMax-xMinReg .LE. tol) .OR. (xMaxReg-xMin .LE. tol) .OR. &
             (yMax-yMinReg .LE. tol) .OR. (yMaxReg-yMin .LE. tol) .OR. &
             (zMax-zMinReg .LE. tol) .OR. (zMaxReg-zMin .LE. tol)) THEN
         pGrid%initPclPresent = .FALSE.
      ELSE
         DO icl = 1,pGrid%nCells
           X = pGrid%cofg(XCOORD,icl)
           Y = pGrid%cofg(YCOORD,icl)
           Z = pGrid%cofg(ZCOORD,icl)
           IF (((xMin-X < tol) .AND. (X-xMax < tol)) .AND. &
              ((yMin-Y < tol) .AND. (Y-yMax < tol)) .AND. &
              ((zMin-Z <= tol) .AND. (Z-zMax <= tol))) THEN
              pGrid%nCellsSeedPcls = pGrid%nCellsSeedPcls + 1
              pGrid%nPclsPerCell(icl) = Rand1Uniform(pRegion%randData)
           ELSE
              pGrid%nPclsPerCell(icl) = -1E6 ! Some large -ve number
           END IF
         END DO
       END IF
     END IF

! ===== Florian =====> Ends

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A,1X,I16)') SOLVER_NAME, &
          'Number of cells that can be seeded with particles is', &
           pGrid%nCellsSeedPcls
  END IF ! global%verbLevel

  IF (pGrid%nCellsSeedPcls .GE. 1) THEN
    pGrid%initPclPresent = .TRUE.
  ELSE

! -------------------- Deallocating memory for regions w/o pcls ----------------
    DEALLOCATE(pGrid%nPclsPerCell,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'nPclsPerCell')
    END IF ! global%error
  END IF

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Computing number of cells that can be seeded with particles done.'
  END IF ! global%verbLevel

! ------------------------------------------------------------------------------
!       Default - must be due to input error
! ------------------------------------------------------------------------------

        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! global%casename
!-------------------------------------------------------------------------------

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_ComputeCellsContainingPcls

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_ComputeCellsContainingPcls.F90,v $
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
!
! ******************************************************************************

