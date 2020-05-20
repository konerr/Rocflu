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
!******************************************************************************
!
! Purpose: Allocate memory for species related to time stepping.
!
! Description: None.
!
! Input: 
!   pRegion     Region pointer
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: SPEC_RFLU_AllocateMemoryTStep.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE SPEC_RFLU_AllocateMemoryTStep(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModParameters  
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch  
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input  
  USE ModMPI

  USE ModInterfaces, ONLY: RFLU_DecideNeedBGradFace  
 
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
  INTEGER :: errorFlag,iPatch,nBFaces,nBFacesTot,nSpecies
  TYPE(t_global), POINTER :: global  
  TYPE(t_grid), POINTER :: pGrid 
  TYPE(t_mixt_input), POINTER :: pMixtInput  
  TYPE(t_patch), POINTER :: pPatch 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_RFLU_AllocateMemoryTStep.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_AllocateMemoryTStep',__FILE__)

! ******************************************************************************
! Set pointers and variables 
! ******************************************************************************

  pGrid      => pRegion%grid
  pMixtInput => pRegion%mixtInput

  nSpecies = pRegion%specInput%nSpecies

! ******************************************************************************
! Allocate memory
! ******************************************************************************

! ==============================================================================
! Compute number of boundary faces for later use
! ==============================================================================

  nBFaces    = 0 
  nBFacesTot = 0 
     
  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)

    nBFaces    = nBFaces    + pPatch%nBTris    + pPatch%nBQuads
    nBFacesTot = nBFacesTot + pPatch%nBTrisTot + pPatch%nBQuadsTot
  END DO ! iPatch
  
! ==============================================================================     
! Old state vector
! ==============================================================================

  ALLOCATE(pRegion%spec%cvOld(nSpecies,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%spec%cvOld')
  END IF ! global%error  

! ==============================================================================     
! Residuals 
! ==============================================================================
  
  ALLOCATE(pRegion%spec%rhs(nSpecies,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%spec%rhs')
  END IF ! global%error  
  
  ALLOCATE(pRegion%spec%diss(nSpecies,pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%spec%diss')
  END IF ! global%error  

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    ALLOCATE(pRegion%spec%rhsSum(nSpecies,pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%spec%rhsSum')
    END IF ! global%error              
  ELSE
    NULLIFY(pRegion%spec%rhsSum)
  END IF ! global%flowType   

! ============================================================================== 
! Gradients 
! ==============================================================================
  
! ------------------------------------------------------------------------------
! Cell gradients
! ------------------------------------------------------------------------------  
  
  IF ( pMixtInput%spaceOrder > 1 ) THEN 
    ALLOCATE(pRegion%spec%gradCell(XCOORD:ZCOORD,nSpecies,pGrid%nCellsTot), &
             STAT=errorFlag)
    global%error = errorFlag         
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%spec%gradCell')
    END IF ! global%error     
  ELSE 
    NULLIFY(pRegion%spec%gradCell)
  END IF ! pMixtInput%spaceOrder  
  
! ------------------------------------------------------------------------------
! Face gradients
! ------------------------------------------------------------------------------  
  
  IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN 
    ALLOCATE(pRegion%spec%gradFace(XCOORD:ZCOORD,nSpecies,pGrid%nFaces), &
             STAT=errorFlag)
    global%error = errorFlag         
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%spec%gradFace')
    END IF ! global%error
  ELSE 
    NULLIFY(pRegion%spec%gradFace)
  END IF ! pMixtInput%flowModel  
   
  IF ( pGrid%nPatches > 0 ) THEN         
    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

      IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
        ALLOCATE(pPatch%spec%gradFace(XCOORD:ZCOORD,nSpecies,pPatch%nBFaces), &
                                        STAT=errorFlag)
        global%error = errorFlag         
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%spec%gradFace')
        END IF ! global%error
      ELSE
        NULLIFY(pPatch%spec%gradFace)         
      END IF ! RFLU_DecideNeedBGradFace
    END DO ! iPatch
  END IF ! pGrid%nPatches
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_RFLU_AllocateMemoryTStep

!******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_RFLU_AllocateMemoryTStep.F90,v $
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
! Revision 1.5  2006/10/20 21:34:33  mparmar
! Bug fix in allocation of pPatch%spec%gradFace
!
! Revision 1.4  2006/08/19 15:40:27  mparmar
! Moved region%spec%bGradFace to patch%spec%gradFace
!
! Revision 1.3  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.2  2004/01/29 22:59:38  haselbac
! Added allocation of cell and face gradients
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
!******************************************************************************

