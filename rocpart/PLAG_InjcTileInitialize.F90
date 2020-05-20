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
! Purpose: initialize the tiles for the multiphase injection algorithm.
!
! Description: none.
!
! Input: region    = current region
!
! Output: regions(iReg)%levels%patch%tile = initial tile values
!                                           of current region.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_InjcTileInitialize.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_InjcTileInitialize( region )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag_input, t_tile_plag
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModRandom, ONLY     : Rand1Uniform

  USE ModError
  USE ModParameters
  USE ModMPI
  USE PLAG_ModParameters
  USE PLAG_ModInflow, ONLY: PLAG_RFLU_InflowMakeParticle
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iPatch, iTile

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: bcType, nPatches, nTiles
  INTEGER :: inflowDiamDist
  INTEGER :: nCv, nDv
  INTEGER :: inflowModel,iCont,nCont
  
  REAL(RFREAL) :: randUnif
  REAL(RFREAL) :: poolVolumeInit, ratioPhiDensInv   
  REAL(RFREAL), POINTER, DIMENSION(:) :: inflowMassFluxRatio, densPlag

  TYPE(t_patch),     POINTER :: pPatch
  TYPE(t_tile_plag), POINTER :: pTilePlag
  TYPE(t_global),    POINTER :: global
    
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InjcTileInitialize.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_InjcTileInitialize',__FILE__ )

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Initializing tile memory for PLAG...'
  END IF ! global%verbLevel
  
! Get dimensions --------------------------------------------------------------

  nPatches = region%grid%nPatches

  nCont        = region%plagInput%nCont  
  densPlag          => region%plagInput%dens

! Loop over patches -----------------------------------------------------------

  DO iPatch=1,nPatches

    pPatch  => region%patches(iPatch)

    bcType = pPatch%bcType

! - Select injection boundary condition ---------------------------------------

    IF (bcType>=BC_INJECTION) THEN

! - Set pointers -------------------------------------------------------------

      pTilePlag => pPatch%tilePlag

! - Get variables ------------------------------------------------------------

      nTiles = pPatch%nBFaces

      nCv = pTilePlag%nCv
      nDv = pTilePlag%nDv

      inflowDiamDist = region%patch%plag%inflowDiamDist
      inflowModel    = region%patch%plag%inflowModel                

      inflowMassFluxRatio => region%patch%plag%inflowMassFluxRatio

! - Initialize variables -------------------------------------------------------

      pTilePlag%nPclsInjc(:nTiles)   = 0
      
      pTilePlag%cv(:nCv,:nTiles)     = 0.0_RFREAL
      pTilePlag%dv(:nDv,:nTiles)     = 0.0_RFREAL
        
      pTilePlag%rhs(:nCv,:nTiles)    = 0.0_RFREAL
        
      pTilePlag%cvOld(:nCv,:nTiles)  = 0.0_RFREAL
      pTilePlag%rhsSum(:nCv,:nTiles) = 0.0_RFREAL
                        
! --- Loop over tile patch ----------------------------------------------------

      DO iTile = 1, nTiles

        CALL PLAG_RFLU_InflowMakeParticle( region, injcDiamDist,         &
                                           pTilePlag%dv(DV_TILE_DIAM,iTile), &
                                           pTilePlag%dv(DV_TILE_SPLOAD,iTile) )
 
        randUnif = Rand1Uniform(region%randData) 

! ---- Avoid error if randUnif is 0 -------------------------------------------
! ---- and set timefactor such that EXP(-50) = 1.9E-22 ------------------------

        IF ( randUnif <= 0.0_RFREAL) THEN          
          pTilePlag%dv(DV_TILE_COUNTDOWN,iTile) = 50.0_RFREAL 
        ELSE
          pTilePlag%dv(DV_TILE_COUNTDOWN,iTile) = -LOG(randUnif)
        END IF ! randUnif

! ---- Set initial pool volume

       IF ( inflowModel == PLAG_EJEC_CRE ) THEN
         poolVolumeInit  = 0.0_RFREAL
         ratioPhiDensInv = poolVolumeInit  &
                         *1.0_RFREAL/SUM(inflowMassFluxRatio(:)/densPlag(:))

         DO iCont = 1, nCont
           pTilePlag%cv(pTilePlag%cvTileMass(iCont),iTile) = poolVolumeInit    &
                                                           * inflowMassFluxRatio(iCont)
         END DO ! iCont
       ENDIF ! ejecModel       
      END DO ! iTile
                                       
    ENDIF  ! bcType

  ENDDO    ! iPatch
  
! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_injcTileInitialize

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InjcTileInitialize.F90,v $
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
! Revision 1.3  2007/05/16 22:29:59  fnajjar
! Modified calls for new bc datastructure
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
! Revision 1.2  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 20:57:45  fnajjar
! Initial revision after changing case
!
! Revision 1.8  2004/06/30 15:43:27  fnajjar
! Added initialization step for CRE model
!
! Revision 1.7  2004/06/16 23:06:33  fnajjar
! Renamed variabled for CRE kernel
!
! Revision 1.6  2004/03/03 00:30:52  fnajjar
! Incorrect ifdef construct for RFLU pPatch pointer
!
! Revision 1.5  2003/11/26 22:00:28  fnajjar
! Removed improper comment symbol after USE ModRandom
!
! Revision 1.4  2003/11/21 22:35:51  fnajjar
! Update Random Number Generator
!
! Revision 1.3  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.2  2003/01/16 20:15:11  f-najjar
! Removed iRegionGlobal
!
! Revision 1.1  2002/10/25 14:16:31  f-najjar
! Initial Import of Rocpart
!
!******************************************************************************

