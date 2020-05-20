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
! Purpose: computes the RHS for the multiphase injection algorithm
!          specific to RFLU.
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
! $Id: PLAG_RFLU_InjcTileCalcRhs.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_InjcTileCalcRhs( pRegion )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag_input, t_tile_plag
  USE ModBndPatch, ONLY   : t_patch, t_bcvalues_plag
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global  
  USE ModError
  USE ModParameters

  USE PLAG_ModParameters

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
 
  INTEGER           :: bcType,c1,distrib,errorFlag,i2dVals,iCont,iCvMass,ifc, &
                       iPatch,iTile,nCont,nPatches
                       
  INTEGER, POINTER, DIMENSION(:) :: bf2c,pCvTileMass

  REAL(RFREAL)  :: area,cp,heatCapSum,heatCapSumR,inflowVelRatio,&
                   massFluxLimit,massFluxSum,massFluxSumR,minj, &
                   nm,nx,ny,nz,rhoMixtbCond,tinj,tileTemp,tileVelNrm               
  REAL(RFREAL), POINTER, DIMENSION(:)   :: inflowMassFluxRatio, inflowTemp, &
                                           specHeat
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: cv, fn, pRhs, vals
    
  TYPE(t_global), POINTER    :: global
  TYPE(t_patch), POINTER     :: pPatch
  TYPE(t_tile_plag), POINTER :: pTilePlag  
  TYPE(t_bcvalues_plag), POINTER :: pPlagBc 

  REAL(RFREAL)  :: heatCapRatio,massFluxRatio,rhoPlag

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_InjcTileCalcRhs.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_InjcTileCalcRhs',__FILE__)

! ******************************************************************************
! Set pointers and variables 
! ******************************************************************************

  nCont    = pRegion%plagInput%nCont
  nPatches = pRegion%grid%nPatches

  specHeat => pRegion%plagInput%spht

  cv => pRegion%mixt%cv

  massFluxLimit = 1.0E-10_RFREAL

! ******************************************************************************  
! Loop over patches 
! ******************************************************************************

  DO iPatch=1,nPatches

    pPatch => pRegion%patches(iPatch)
    vals   => pPatch%mixt%vals
    
    bf2c   => pPatch%bf2c
    fn     => pPatch%fn
    
    bcType = pPatch%bcType
    distrib = pPatch%mixt%distrib

    pPlagBc => pPatch%plag

!===============================================================================
!   Select injection boundary condition 
!===============================================================================

    IF ( bcType == BC_INJECTION      .OR. &
         bcType == BC_INFLOW         .OR. &
	 bcType == BC_INFLOW_VELTEMP      ) THEN 

!-------------------------------------------------------------------------------
!    Set bc data
!-------------------------------------------------------------------------------

     inflowVelRatio = pPlagBc%inflowVelRatio  
     inflowMassFluxRatio => pPlagBc%inflowMassFluxRatio
     inflowTemp          => pPlagBc%inflowTemp

!-------------------------------------------------------------------------------
!     Loop over face patches
!-------------------------------------------------------------------------------

      DO ifc = 1, pPatch%nBFaces
        iTile = ifc
        c1 = bf2c(ifc)

        nx = fn(XCOORD,ifc)
        ny = fn(YCOORD,ifc)
        nz = fn(ZCOORD,ifc)
        nm = fn(XYZMAG,ifc)

        i2dVals = distrib * ifc
        mInj       = vals(BCDAT_INJECT_MFRATE,i2dVals)
        tInj       = vals(BCDAT_INJECT_TEMP  ,i2dVals)

!-------------------------------------------------------------------------------
!       Compute mixture density on boundary 
!        Using a low-order projection based on cell value
!-------------------------------------------------------------------------------
        
        rhoMixtbCond = cv(CV_MIXT_DENS,c1)

!-------------------------------------------------------------------------------
!       Update tile Rhs datastructure 
!-------------------------------------------------------------------------------

        pTilePlag => pPatch%tilePlag
 
        pRhs        => pTilePlag%rhs
        pCvTileMass => pTilePlag%cvTileMass

        SELECT CASE(bcType)

!-------------------------------------------------------------------------------
!         Injection 
!-------------------------------------------------------------------------------

          CASE(BC_INJECTION) 
            massFluxSum = SUM ( mInj *inflowMassFluxRatio(:))
            heatCapSum  = DOT_PRODUCT( specHeat, mInj *inflowMassFluxRatio(:) )

            tileTemp    = tInj
            tileVelNrm  = inflowVelRatio *mInj/rhoMixtbCond
            area        = nm

            DO iCont = 1, nCont
              iCvMass = pCvTileMass(iCont)
              pRhs(iCvMass,iTile) = -area *mInj *inflowMassFluxRatio(iCont) 
            END DO ! iCont            

            pRhs(CV_TILE_MOMNRM,iTile) = -area *massFluxSum *tileVelNrm

            pRhs(CV_TILE_ENER  ,iTile) = -area *                             &
                                      ( 0.5_RFREAL*massFluxSum*tileVelNrm**2 &
                                      +            heatCapSum *tileTemp      )

!-------------------------------------------------------------------------------
!         Inflow based on velocity and temperature 
!-------------------------------------------------------------------------------

          CASE(BC_INFLOW_VELTEMP) 
            area        = nm

            massFluxSum = 0.0_RFREAL
            heatCapSum  = 0.0_RFREAL
            tileVelNrm  = 0.0_RFREAL
            tileTemp    = 0.0_RFREAL

            DO iCont = 1, nCont
              rhoPlag =  pRegion%plagInput%dens(iCont)
              massFluxSum = massFluxSum + inflowMassFluxRatio(iCont) *rhoPlag
              heatCapSum  = heatCapSum  + specHeat(iCont) &
	                                * inflowMassFluxRatio(iCont) *rhoPlag 
            END DO ! iCont

! ==============================================================================
!           Set inverse of massFluxSum to avoid division by zero 
! ==============================================================================
  
            IF ( massFluxSum > massFluxLimit ) THEN
              massFluxSumR = 1.0_RFREAL/massFluxSum
              heatCapSumR  = 1.0_RFREAL/heatCapSum
            ELSE
              massFluxSumR = 1.0_RFREAL
              heatCapSumR  = 1.0_RFREAL
            END IF ! massFluxSum

            DO iCont = 1, nCont
              iCvMass = pCvTileMass(iCont)
              rhoPlag = pRegion%plagInput%dens(iCont)

              massFluxRatio = inflowMassFluxRatio(iCont) *rhoPlag *massFluxSumR
              heatCapRatio  = inflowMassFluxRatio(iCont) *rhoPlag &
                            * specHeat(iCont) *heatCapSumR

              tileVelNrm    = tileVelNrm +massFluxRatio *inflowMassFluxRatio(iCont)
              tileTemp      = tileTemp   +heatCapRatio  *inflowTemp(iCont)

              pRhs(iCvMass,iTile) = -area * inflowMassFluxRatio(iCont) *rhoPlag
            END DO ! iCont

            pRhs(CV_TILE_MOMNRM,iTile) = -area *massFluxSum *tileVelNrm

            pRhs(CV_TILE_ENER  ,iTile) = -area *                              &
                                       ( 0.5_RFREAL*massFluxSum*tileVelNrm**2 &
                                       +            heatCapSum *tileTemp      )

!-------------------------------------------------------------------------------
!         Trap error for default 
!-------------------------------------------------------------------------------

          CASE DEFAULT 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! bcType

      END DO ! ifc                                                        
    ENDIF    ! bcType 
  ENDDO      ! iPatch

! ******************************************************************************
! finalize
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLU_InjcTileCalcRhs

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_InjcTileCalcRhs.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.5  2008/12/06 08:43:52  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/08/07 21:54:01  fnajjar
! Removed statements with BC_RANGE since obsolete for RocfluMP
!
! Revision 1.2  2007/05/16 22:37:16  fnajjar
! Modified routines to be aligned with new bc datastructure
!
! Revision 1.1  2007/04/09 18:50:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:34  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.6  2007/03/09 15:16:31  fnajjar
! Fixed bug for massFluxSum being zero and avoiding division by zero
!
! Revision 1.5  2006/09/18 20:35:30  fnajjar
! Activated tile datastructure for inflow bc and created proper particle bc
!
! Revision 1.4  2006/08/19 15:40:07  mparmar
! Renamed patch variables
!
! Revision 1.3  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.2  2004/12/07 22:56:24  fnajjar
! Removed rhoVrel being obsolete
!
! Revision 1.1  2004/03/08 22:36:08  fnajjar
! Initial import of RFLU-specific injection routines
!
!******************************************************************************

