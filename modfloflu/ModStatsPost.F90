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
! Purpose: Suite for statistics postprocessing routines.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ModStatsPost.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE ModStatsPost

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY: t_grid
  USE ModMPI
 
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: STAT_RFLU_AllocMemPost, &
            STAT_RFLU_AllocMemPostVert, &
            STAT_RFLU_DeallocMemPost, &
            STAT_RFLU_DeallocMemPostVert
        
! ******************************************************************************
! Declarations and definitions
! ****************************************************************************** 
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: ModStatsPost.F90,v $ $Revision: 1.1.1.1 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS


! ******************************************************************************
!
! Purpose: Allocate memory of statistics post quantities.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE STAT_RFLU_AllocMemPost(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: nStat, errorFlag
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'STAT_RFLU_AllocMemPost',__FILE__)

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid  => pRegion%grid

! ******************************************************************************
! Allocate memory
! ******************************************************************************

! ==============================================================================
! Mixture statistics
! ==============================================================================

  nStat = global%mixtNStat
  IF (nStat > 0) THEN
    ALLOCATE(pRegion%mixt%tav(nStat,pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%mixt%tav')
    END IF ! global%error
  ENDIF

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE STAT_RFLU_AllocMemPost



! ******************************************************************************
!
! Purpose: Deallocate memory of statistics post quantities.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE STAT_RFLU_DeallocMemPost(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: nStat, errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'STAT_RFLU_DeallocMemPost',__FILE__)

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

! ==============================================================================
! Mixture statistics
! ==============================================================================

  nStat = global%mixtNStat
  IF (nStat > 0) THEN
    DEALLOCATE(pRegion%mixt%tav,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%tav')
    END IF ! global%error
  ENDIF

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE STAT_RFLU_AllocMemPostVert



! ******************************************************************************
!
! Purpose: Deallocate memory of statistics arrays at vertices.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE STAT_RFLU_DeallocMemPostVert(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: nStat, errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'STAT_RFLU_DeallocMemPostVert',__FILE__)

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

! ==============================================================================
! Mixture statistics
! ==============================================================================

  nStat = global%mixtNStat
  IF (nStat > 0) THEN
    DEALLOCATE(pRegion%mixt%tavVert,STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%tavVert')
    END IF ! global%error
  ENDIF

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE STAT_RFLU_DeallocMemPostVert


! ******************************************************************************
! End
! ******************************************************************************
      
END MODULE ModStatsPost

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModStatsPost.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:38  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:11  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:17  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2006/01/12 09:51:49  wasistho
! initial import
!
! ******************************************************************************

