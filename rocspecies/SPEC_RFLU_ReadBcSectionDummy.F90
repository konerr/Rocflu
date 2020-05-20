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
! Purpose: Read in user input related to boundary condition for species for 
!   patches which do not actually require any data. 
!
! Description: None.
!
! Input: 
!   pRegion     Region pointer
!
! Output: None.
!
! Notes: This routine is necessary because the allocation of arrays related to
!   time-dependent boundary conditions requires that all patches have had their
!   values of nData initialized. 
!
! ******************************************************************************
!
! $Id: SPEC_RFLU_ReadBcSectionDummy.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_RFLU_ReadBcSectionDummy(pRegion)

  USE ModDataTypes
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModParameters
  
  USE ModInterfaces, ONLY: MakeNumberedKeys, & 
                           ReadPatchSection 

  IMPLICIT NONE

! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region), POINTER :: pRegion

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: bcName,RCSIdentString
  CHARACTER(10), DIMENSION(:), ALLOCATABLE :: keys
  CHARACTER(256) :: fileName
  LOGICAL, DIMENSION(:), ALLOCATABLE :: defined
  INTEGER :: checkSum,distrib,errorFlag,iKey,iPatch,iPatchBeg,iPatchEnd, &
             iReg,iVal,nKeys
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: vals
  TYPE(t_grid) :: grid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: SPEC_RFLU_ReadBcSectionDummy.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_ReadBcSectionDummy',__FILE__)

! *****************************************************************************
! Allocate memory 
! *****************************************************************************

  nKeys = pRegion%specInput%nSpecies + 1 ! Add one because of default key

  ALLOCATE(keys(nKeys),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'keys')
  END IF ! global%error

  ALLOCATE(vals(nKeys),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vals')
  END IF ! global%error  

  ALLOCATE(defined(nKeys),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'defined')
  END IF ! global%error 

! *****************************************************************************
! Generate keys. NOTE first key is default key.
! *****************************************************************************

  keys(1) = 'SPEC_'
  
  CALL MakeNumberedKeys(keys,2,'SPEC',1,nKeys,1)

! *****************************************************************************
! Read section
! *****************************************************************************

  CALL ReadPatchSection(global,IF_INPUT,nKeys,keys,vals,iPatchBeg,iPatchEnd, &
                        distrib,fileName,bcName,defined)

! *****************************************************************************
! Check if specified number of patches exceeds available ones
! *****************************************************************************

  IF ( iPatchEnd > global%nPatches ) THEN 
    CALL ErrorStop(global,ERR_PATCH_RANGE,__LINE__)
  END IF ! iPatchEnd

! *****************************************************************************
! Set options and check if all necessary values defined
! *****************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

! =============================================================================
!   Check whether this global patch exists in this region
! =============================================================================

    IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. & 
         pPatch%iPatchGlobal <= iPatchEnd ) THEN
      pPatch%spec%nData     = 0
      pPatch%spec%nSwitches = 0
      pPatch%spec%distrib   = BCDAT_CONSTANT
    END IF ! pPatch%iPatchGlobal
  END DO ! iPatch

! *****************************************************************************
! Deallocate memory 
! *****************************************************************************

  DEALLOCATE(keys,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'keys')
  END IF ! global%error

  DEALLOCATE(vals,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vals')
  END IF ! global%error  

  DEALLOCATE(defined,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'defined')
  END IF ! global%error 
  
! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_RFLU_ReadBcSectionDummy

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_RFLU_ReadBcSectionDummy.F90,v $
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
! Revision 1.3  2006/08/19 15:40:34  mparmar
! Renamed patch variables
!
! Revision 1.2  2006/04/07 15:19:25  haselbac
! Removed tabs
!
! Revision 1.1  2005/03/31 17:22:25  haselbac
! Initial revision
!
! ******************************************************************************

