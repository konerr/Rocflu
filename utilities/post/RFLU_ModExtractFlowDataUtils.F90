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
! Purpose: Collect utility routines for extraction of data from flow solution.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModExtractFlowDataUtils.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2007-2008 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModExtractFlowDataUtils

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region

#ifdef PLAG
  USE PLAG_ModParameters
#endif

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PARAMETER, PRIVATE :: &
    RCSIdentString = '$RCSfile: RFLU_ModExtractFlowDataUtils.F90,v $ $Revision: 1.1.1.1 $'

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: RFLU_ExtractDiscontLocation1D

! ==============================================================================
! Private functions
! ==============================================================================



! ******************************************************************************
! Subroutines and functions
! ******************************************************************************

  CONTAINS




! ******************************************************************************
!
! Purpose: Extract discontinuity location for 1D cases.
!
! Description: Look for location of largest density gradient and then, 
!   starting from this position, determine position of zero second derivative, 
!   which is then taken to be shock location.
!
! Input:
!   pRegion		Pointer to region
!   icgBeg		Beginning cell index
!   icgEnd		Ending cell index
!   nCellsX		Number of cells in x-direction
!   pXv			Pointer to solution data
!   iXv			Index of variable in pXv
!   xDiscMin		Lower value of x
!   xDiscMax		Upper value of x
!
! Output: 
!   xs			Discontinuity location
!
! Notes: 
!   1. This routine can be used for 2d cases, but then the array numbering must
!      be such that the indices follow the longest direction.
!   2. When no suitable position could be found, CRAZY_VALUE_INT is returned.
!
! ******************************************************************************

SUBROUTINE RFLU_ExtractDiscontLocation1D(pRegion,icgBeg,icgEnd,nCellsX, &
                                         pXv,iXv,xs,xDiscMin,xDiscMax)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: icgBeg,icgEnd,iXv,nCellsX
  REAL(RFREAL), INTENT(IN), OPTIONAL :: xDiscMax,xDiscMin
  REAL(RFREAL), INTENT(OUT) :: xs  
  REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: pXv
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: foundFlag
  INTEGER :: errorFlag,icg,icl,iclDiscont,iclMax,iclMin,iclOffs
  INTEGER :: dummy(1)
  REAL(RFREAL) :: dist,distMax,distMin,idx,r,rm1,rp1,rxx,rxxp1,x,xp1
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: gradx,gradxx
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ExtractDiscontLocation1D', &
                        __FILE__)

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Set variables
! ******************************************************************************

  iclOffs = 4
 
  IF ( (PRESENT(xDiscMax) .EQV. .TRUE.) .AND. & 
       (PRESENT(xDiscMin) .EQV. .TRUE.) ) THEN 
    distMax = HUGE(1.0_RFREAL)
    distMin = HUGE(1.0_RFREAL)

    DO icl = 1,nCellsX-2
      icg = icgBeg + icl

      dist = ABS(pGrid%cofg(XCOORD,icg)-xDiscMin)

      IF ( dist < distMin ) THEN
        distMin = dist
        iclMin  = icl
      END IF ! dist

      dist = ABS(pGrid%cofg(XCOORD,icg)-xDiscMax)

      IF ( dist < distMax ) THEN
        distMax = dist
        iclMax  = icl
      END IF ! dist
    END DO ! icl 
  ELSE 
    iclMin = 1
    iclMax = nCellsX-2 
  END IF ! PRESENT

  xs = REAL(CRAZY_VALUE_INT,KIND=RFREAL)

! ******************************************************************************
! Allocate temporary memory
! ******************************************************************************

  ALLOCATE(gradx(nCellsX-2),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gradx')
  END IF ! global%error

  ALLOCATE(gradxx(nCellsX-2),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gradxx')
  END IF ! global%error

! ******************************************************************************
! Compute first and second density derivatives and find location of maximum 
! (absolute) value of first derivative
! ******************************************************************************

  DO icl = iclMin,iclMax
    icg = icgBeg + icl 

    rp1 = pXv(iXv,icg+1)       
    r   = pXv(iXv,icg  )
    rm1 = pXv(iXv,icg-1)       

    gradx(icl)  = 0.5_RFREAL*(rp1-rm1)
    gradxx(icl) = rp1-2.0_RFREAL*r+rm1
  END DO ! icl

  dummy = MAXLOC(ABS(gradx(iclMin:iclMax)))
  iclDiscont = dummy(1) 

! ******************************************************************************
! Starting from this location, search for zero crossing of second derivative 
! and then use linear approximation to compute location of zero crossing. NOTE 
! this search is carried out for iclOffs points on either side of location of 
! maximum (absolute) value of first derivative, so need to make sure do not step 
! out of bounds. If that is the case, simply pick shock location to be that cell
! with maximum (absolute) value of first derivative.
! ******************************************************************************

  foundFlag = .FALSE.

  IF ( (iclDiscont <= (nCellsX-iclOffs-2)) .AND. & 
       (iclDiscont >= (iclOffs+1)) ) THEN  
    iclLoop: DO icl = iclDiscont-iclOffs,iclDiscont+iclOffs-1
      icg = icgBeg + icl

! TEMPORARY - NOTE the following criterion may need to be evaluated in more 
! detail. First criterion does not work at particle contact because of bias
! introduced by magnitude of 1 (and not -1).  
!      IF ( SIGN(1.0_RFREAL,gradxx(icl)) /= SIGN(1.0_RFREAL,gradxx(icl+1)) ) THEN 
      IF ( gradxx(icl)*gradxx(icl+1) < 0.0_RFREAL ) THEN 
! END TEMPORARY
        rxxp1 = gradxx(icl+1)
        rxx   = gradxx(icl)
  
        xp1 = pGrid%cofg(XCOORD,icg+1)
        x   = pGrid%cofg(XCOORD,icg)
  
        xs = (x*rxxp1-xp1*rxx)/(rxxp1-rxx)

        foundFlag = .TRUE.

        EXIT iclLoop
      END IF ! SIGN
    END DO iclLoop
  ELSE 
    foundFlag = .TRUE.

    xs = pGrid%cofg(XCOORD,icgBeg+iclDiscont)
  END IF ! icl

  IF ( foundFlag .EQV. .FALSE. ) THEN 
    global%warnCounter = global%warnCounter + 1

    WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'*** WARNING ***', & 
                                  'Discontinuity could not be located!'
  END IF ! foundFlag 

! ******************************************************************************
! Deallocate temporary memory
! ******************************************************************************

  DEALLOCATE(gradx,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gradx')
  END IF ! global%error

  DEALLOCATE(gradxx,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gradxx')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ExtractDiscontLocation1D





END MODULE RFLU_ModExtractFlowDataUtils

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModExtractFlowDataUtils.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.7  2008/12/06 08:43:58  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:12  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2008/01/20 22:55:40  haselbac
! Fix problem when disc cannot be located (with imposed extrema)
!
! Revision 1.4  2008/01/17 13:13:19  haselbac
! Improved shock extraction routine
!
! Revision 1.3  2008/01/10 18:47:54  haselbac
! Converted shock extraction to discontinuity extraction
!
! Revision 1.2  2007/04/12 12:11:57  haselbac
! Changed shock extraction
!
! Revision 1.1  2007/04/05 01:38:07  haselbac
! Initial revision
!
! ******************************************************************************

