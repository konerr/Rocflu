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
! Purpose: Update solution and sums residuals for tile data infrastructure.
!
! Description: none.
!
! Input: region = current region.
!        iStage = RK stage
!
! Output: region%levels%tilePlag%cv
!         region%levels%tilePlag%rhsSum
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_InjcTileRKUpdate.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_InjcTileRKUpdate( region, iStage )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag_input, t_tile_plag
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMixture, ONLY    : t_mixt_input
  USE ModError  
  USE ModParameters

  USE PLAG_ModParameters
  
  USE ModInterfaces, ONLY : RkUpdateGeneric
  
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  
  INTEGER        :: iStage

! ... loop variables
  INTEGER :: iCont, iPatch, iTile

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: bcType, nCont, nPatches, nTiles
  INTEGER :: ivTileBeg, ivTileEnd

  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv, pCvOld, pRhs, pRhsSum
  
  TYPE(t_patch),     POINTER :: pPatch
  TYPE(t_tile_plag), POINTER :: pTilePlag    
  TYPE(t_global),    POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InjcTileRKUpdate.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_InjcTileRKUpdate',__FILE__ )
 
! Get dimensions --------------------------------------------------------------

  nPatches = region%grid%nPatches
  nCont    = region%plagInput%nCont 
 
! Loop over patches -----------------------------------------------------------

  DO iPatch=1,nPatches

    pPatch  => region%patches(iPatch)

    bcType = pPatch%bcType

! - Select injection boundary condition ---------------------------------------

    IF ( bcType == BC_INJECTION     .OR. &
         bcType == BC_INFLOW        .OR. &
	 bcType == BC_INFLOW_VELTEMP     ) THEN 
    
! -- Get tile dimensions and pointers -----------------------------------------

      nTiles  = pPatch%nBFaces
      
      pTilePlag   => pPatch%tilePlag

      pCv     => pTilePlag%cv
      pCvOld  => pTilePlag%cvOld
      pRhs    => pTilePlag%rhs
      pRhsSum => pTilePlag%rhsSum

! -- Use Generic RKUpdate routine ---------------------------------------------
      
      ivTileBeg = CV_TILE_MOMNRM
      ivTileEnd = CV_TILE_LAST+nCont
      
      CALL RkUpdateGeneric( region,VAR_TYPE_POINT,iStage,1,nTiles, &
                            ivTileBeg,ivTileEnd,pCv,pCvOld,pRhs,pRhsSum )
      
    END IF !bcType
    
  END DO ! iPatch
  
! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_InjcTileRKUpdate

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InjcTileRKUpdate.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.5  2008/12/06 08:43:50  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/08/07 21:54:01  fnajjar
! Removed statements with BC_RANGE since obsolete for RocfluMP
!
! Revision 1.2  2007/04/16 23:20:36  fnajjar
! Deleted all RFLO-pertinent calls and statements
!
! Revision 1.1  2007/04/09 18:50:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/09/18 20:30:18  fnajjar
! Activated tile datastructure for inflow bc
!
! Revision 1.1  2004/12/01 20:57:47  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2004/11/17 16:43:17  haselbac
! Replaced RKUpdate call with call to RkUpdateGeneric
!
! Revision 1.4  2004/03/08 22:23:14  fnajjar
! Modified routine to be RFLU-aware and added call to PLAG_rkUpdateGeneric
!
! Revision 1.3  2004/02/25 21:56:15  fnajjar
! Moved tile pointers outside do-loop
!
! Revision 1.2  2003/01/16 20:15:11  f-najjar
! Removed iRegionGlobal
!
! Revision 1.1  2002/10/25 14:16:32  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************

