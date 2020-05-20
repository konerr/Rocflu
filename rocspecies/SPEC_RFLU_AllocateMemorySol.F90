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
! Purpose: Allocate memory for species solution.
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
!
! $Id: SPEC_RFLU_AllocateMemorySol.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_RFLU_AllocateMemorySol(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

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
  INTEGER :: errorFlag,iSpec,nSpecies
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_RFLU_AllocateMemorySol.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_AllocateMemorySol',__FILE__)

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid  => pRegion%grid

  nSpecies = pRegion%specInput%nSpecies

! ******************************************************************************
! Allocate memory
! ******************************************************************************

! ==============================================================================
! State vectors
! ==============================================================================

  ALLOCATE(pRegion%spec%cv(nSpecies,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%spec%cv')
  END IF ! global%error

  ALLOCATE(pRegion%spec%cvInfo(nSpecies),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%spec%cvInfo')
  END IF ! global%error

! ==============================================================================
! Transport variables
! ==============================================================================

  IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN
    ALLOCATE(pRegion%spec%tv(nSpecies,pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%spec%tv')
    END IF ! global%error
  ELSE
    NULLIFY(pRegion%spec%tv)
  END IF ! pMixtInput%flowModel

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_RFLU_AllocateMemorySol

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_RFLU_AllocateMemorySol.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:53  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:05  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:51:23  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:50  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2004/11/02 02:33:52  haselbac
! Removed init of cvInfo, cosmetics
!
! Revision 1.3  2004/03/03 23:54:41  jferry
! Corrected typo
!
! Revision 1.2  2004/01/29 22:59:37  haselbac
! Added allocation of cvInfo and tv
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
! ******************************************************************************

