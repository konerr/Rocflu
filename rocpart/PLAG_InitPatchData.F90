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
! Purpose: Initialize patch data for Rocpart.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Only treat injection boundaries right now.  
!   2. Avoid error if randUnif is 0, and set timefactor such that 
!      EXP(-50) = 1.9E-22
!
!******************************************************************************
!
! $Id: PLAG_InitPatchData.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_InitPatchData(pRegion)

  USE ModDataTypes
  USE ModPartLag, ONLY: t_tile_plag
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModRandom, ONLY: Rand1Uniform
  USE ModError
  USE ModParameters
  USE ModMPI

  USE PLAG_ModParameters
  
  USE PLAG_ModInflow, ONLY: PLAG_RFLU_InflowMakeParticle

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
  INTEGER :: inflowDiamDist,iPatch,iTile,nPatches,nTiles
  REAL(RFREAL) :: randUnif
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_tile_plag), POINTER :: pTilePlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InitPatchData.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_InitPatchData',__FILE__)

  nPatches = pRegion%grid%nPatches

! ******************************************************************************
! Loop over patches
! ******************************************************************************

  DO iPatch=1,nPatches
    pPatch => pRegion%patches(iPatch)

! ==============================================================================
!   Select patch type
! ==============================================================================

    SELECT CASE ( pPatch%bcType )
    
! ------------------------------------------------------------------------------    
!     Injection boundary
! ------------------------------------------------------------------------------    
    
      CASE ( BC_INJECTION,                &
             BC_INFLOW, BC_INFLOW_VELTEMP ) 

! ----- Get dimensions ---------------------------------------------------------

        nTiles = pPatch%nBFaces

        pTilePlag => pPatch%tilePlag


! ----- Set distribution type --------------------------------------------------

        inflowDiamDist = pPatch%plag%inflowDiamDist

! ----- Loop over tiles --------------------------------------------------------

        DO iTile = 1,nTiles
          CALL PLAG_RFLU_InflowMakeParticle( pRegion,pPatch,inflowDiamDist,     &
                                             pTilePlag%dv(DV_TILE_DIAM,iTile),  &
                                             pTilePlag%dv(DV_TILE_SPLOAD,iTile) )

          randUnif = Rand1Uniform(pRegion%randData) 

          IF ( randUnif <= 0.0_RFREAL ) THEN          
            pTilePlag%dv(DV_TILE_COUNTDOWN,iTile) = 50.0_RFREAL 
          ELSE
            pTilePlag%dv(DV_TILE_COUNTDOWN,iTile) = -LOG(randUnif)
          END IF ! randUnif
        END DO ! iTile
    END SELECT ! pPatch%bcType
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_InitPatchData

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InitPatchData.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.6  2008/12/06 08:43:50  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2007/08/07 21:54:01  fnajjar
! Removed statements with BC_RANGE since obsolete for RocfluMP
!
! Revision 1.3  2007/05/16 22:28:59  fnajjar
! Modified calls for new bc datastructure
!
! Revision 1.2  2007/04/16 23:20:36  fnajjar
! Deleted all RFLO-pertinent calls and statements
!
! Revision 1.1  2007/04/09 18:50:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/09/18 20:29:07  fnajjar
! Activated tile infrastructure for inflow bc
!
! Revision 1.3  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.2  2004/06/16 22:58:13  fnajjar
! Renamed injcModel to injcDiamDist and dv(TIMEFCTR) to dv(COUNTDOWN) for CRE kernel
!
! Revision 1.1  2004/03/05 23:15:54  haselbac
! Initial revision
!
!******************************************************************************

