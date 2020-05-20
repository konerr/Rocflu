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
! Purpose: Pick regions by coordinate range.
!
! Description: Check whether a given region is inside a user-specified 
!   bounding box; if yes, flag region as active, otherwise as inactive.
!
! Input: 
!   pRegion		Pointer to region data
!
! Output: None.
!
! Notes: 
!   1. If a given region is not inside the bounding box, then postActiveFlag
!      is set to FALSE.
!
! ******************************************************************************
!
! $Id: RFLU_PickRegionsCoord.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PickRegionsCoord(pRegion)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModMPI
  USE ModParameters
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  REAL(RFREAL) :: xMax,xMin,yMax,yMin,zMax,zMin
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PickRegionsCoord.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_PickRegionsCoord',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Picking region by coordinates...'
  END IF ! global%verbLevel

! ******************************************************************************
! Pick by bounding box
! ******************************************************************************

  pGrid => pRegion%grid

  xMin = MINVAL(pGrid%xyz(XCOORD,1:pGrid%nVertTot))
  xMax = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVertTot))
  yMin = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVertTot))
  yMax = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVertTot))
  zMin = MINVAL(pGrid%xyz(ZCOORD,1:pGrid%nVertTot))
  zMax = MAXVAL(pGrid%xyz(ZCOORD,1:pGrid%nVertTot))

  IF ( xMin > global%pickXCoordLow .AND. xMax < global%pickXCoordUpp .AND. & 
       yMin > global%pickYCoordLow .AND. yMax < global%pickYCoordUpp .AND. & 
       zMin > global%pickZCoordLow .AND. zMax < global%pickZCoordUpp ) THEN 
    pRegion%postActiveFlag = .TRUE.

    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME,'Picked region:', &
                                     pRegion%iRegionGlobal
    END IF ! global%verbLevel
  ELSE 
    pRegion%postActiveFlag = .FALSE.  
    
    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I6,1X,A)') SOLVER_NAME,'Region:', &
                                          pRegion%iRegionGlobal,'not picked.'
    END IF ! global%verbLevel    
  END IF ! xMin    
  
! ******************************************************************************
! End
! ******************************************************************************

  IF ( (global%myProcid == MASTERPROC) .AND. & 
       (global%verbLevel > VERBOSE_NONE) ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Picking region by coordinates done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PickRegionsCoord

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PickRegionsCoord.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:57  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:11  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:57:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2005/12/10 23:29:34  haselbac
! Renamed geom post variables, cosmetics
!
! Revision 1.2  2004/10/19 19:30:08  haselbac
! Added output statement, cosmetics
!
! Revision 1.1  2003/08/07 15:14:40  haselbac
! Initial revision
!
! ******************************************************************************

