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
! Purpose: Collection of routines for program burn. 
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: SPEC_RFLU_ModPBA.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2006-2007 by the University of Illinois
!
! ******************************************************************************

MODULE SPEC_RFLU_ModPBA

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModError

#ifdef SPEC
  USE ModSpecies, ONLY: t_spec_input, &
                        t_spec_type
  USE ModInterfacesSpecies, ONLY: SPEC_GetSpeciesIndex
  USE SPEC_RFLU_ModPBAUtils, ONLY: SPEC_RFLU_PBA_ComputeY, &
                                   SPEC_RFLU_PBA_GetSolution
#endif

  USE RFLU_ModConvertCv, ONLY: RFLU_ScalarConvertCvCons2Prim, &
                               RFLU_ScalarConvertCvPrim2Cons

  USE ModInterfaces, ONLY: MixtPerf_G_CpR, &
                           MixtPerf_R_M

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: SPEC_RFLU_ModPBA.F90,v $'
  
! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: SPEC_RFLU_PBA_ProgramBurn, &
            SPEC_RFLU_PBA_ReactionZone, &
            SPEC_RFLU_PBA_ZeroResiduals

! ==============================================================================
! Private functions
! ==============================================================================
  
!  PRIVATE ::

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  




! ******************************************************************************
!
! Purpose: Program burn algorithm.
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

  SUBROUTINE SPEC_RFLU_PBA_ProgramBurn(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
   
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================

    CHARACTER(CHRLEN) :: errorString
    INTEGER :: icg,iCvSpecAir,iCvSpecExplosive,iCvSpecProducts,indCp,indMol, &
               iVarScal
    REAL(RFREAL) :: a,cpProducts,detonVel,e,Eo,flameLoc,flameWidth,gcProducts, &
                    gProducts,mwProducts,nTol,p,r,ro,T,u,Vm2,x,xFinal,xFront, &
                    xHalfW,xInit,xTail,xWidth, &
                    y,YProducts
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvMixt,pCvOldSpec,pCvSpec, &
                                             pDvMixt,pGvMixt,pRhsSpec, &
                                             pRhsSumSpec
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_mixt_input), POINTER :: pMixtInput
    TYPE(t_spec_input), POINTER :: pSpecInput

! ******************************************************************************
!   Start, if pbaBurnFlag is FALSE then no need to go further
! ******************************************************************************

    global => pRegion%global

    IF ( global%pbaBurnFlag .EQV. .FALSE. ) THEN
      RETURN
    END IF ! flameLoc

    CALL RegisterFunction(global,'SPEC_RFLU_PBA_ProgramBurn',__FILE__)

! ******************************************************************************
!   Set pointers and variables 
! ******************************************************************************

    pGrid      => pRegion%grid
    pCvMixt    => pRegion%mixt%cv
    pCvOldSpec => pRegion%spec%cvOld
    pCvSpec    => pRegion%spec%cv
    pDvMixt    => pRegion%mixt%dv
    pGvMixt    => pRegion%mixt%gv
    pMixtInput => pRegion%mixtInput
    pSpecInput => pRegion%specInput

    pRhsSpec    => pRegion%spec%rhs
    pRhsSumSpec => pRegion%spec%rhsSum

    indCp  = pRegion%mixtInput%indCp
    indMol = pRegion%mixtInput%indMol

    nTol = 1.0E-14_RFREAL

! ******************************************************************************
!   Need species in primitive form
! ******************************************************************************

    CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv, &
                                       pRegion%spec%cvState)

! ******************************************************************************
!   Use program burn algorithm based on case name
! ******************************************************************************

    SELECT CASE ( global%casename )

! ==============================================================================
!     Jenkins 2D planar case
! ==============================================================================

      CASE ( "jenkins2d" )
        iCvSpecAir       = SPEC_GetSpeciesIndex(global,pSpecInput,'AIR')
        iCvSpecExplosive = SPEC_GetSpeciesIndex(global,pSpecInput,'EXPLOSIVE')
        iCvSpecProducts  = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')

        mwProducts = pSpecInput%specType(iCvSpecProducts)%pMaterial%molw
        cpProducts = pSpecInput%specType(iCvSpecProducts)%pMaterial%spht

        gcProducts = MixtPerf_R_M(mwProducts)
        gProducts  = MixtPerf_G_CpR(cpProducts,gcProducts)

        ro = pMixtInput%prepRealVal6

        detonVel   = pSpecInput%specType(iCvSpecExplosive)%pMaterial%detonVel
        flameWidth = pMixtInput%prepRealVal8
        xInit      = pMixtInput%prepRealVal1
        xInit      = xInit + flameWidth ! Initialize whole reaction zone inside explosive region
        xFinal     = pMixtInput%prepRealVal2
        flameLoc   = xInit + detonVel*(global%currentTimeRK)

        xFront = xInit + detonVel*(global%currentTimeRK)
        xWidth = pMixtInput%prepRealVal8
        xTail  = xFront - xWidth
!        xHalfW = 0.0_RFREAL
        xHalfW = 5.0E-05_RFREAL

        DO icg = 1,pGrid%nCellsTot
          x = pGrid%cofg(XCOORD,icg)
          y = pGrid%cofg(YCOORD,icg)

          IF ( x >= pMixtInput%prepRealVal1 .AND. &
               x <= pMixtInput%prepRealVal2 .AND. &
               y >= pMixtInput%prepRealVal3 .AND. &
               y <= pMixtInput%prepRealVal4 ) THEN             
! ------------------------------------------        
! TEMPORARY: Manoj-PBA, Using volume averaged Y for jenkins2d case
!            Need to improve logic for general case
!            For now no air is assumed in initial explosive location

!            pCvSpec(iCvSpecProducts,icg)  = SPEC_RFLU_PBA_ComputeY(x,xFront,xHalfW,xWidth,xTail)
!            pCvSpec(iCvSpecExplosive,icg) = 1.0_RFREAL - pCvSpec(iCvSpecProducts,icg)

! --------- Transition zone ----------------------------------------------------
            IF ( (flameLoc > pGrid%cofg(XCOORD,icg)) .AND. &
                 (flameLoc - flameWidth < pGrid%cofg(XCOORD,icg)) ) THEN
              pCvSpec(iCvSpecProducts,icg)  = &
                                    (flameLoc-pGrid%cofg(XCOORD,icg))/flameWidth
              pCvSpec(iCvSpecExplosive,icg) = 1.0_RFREAL &
                                              - pCvSpec(iCvSpecProducts,icg)
! --------- Just formed products -----------------------------------------------
            ELSEIF ( (flameLoc - flameWidth > pGrid%cofg(XCOORD,icg)) .AND. &
                     (ABS(pCvSpec(iCvSpecExplosive,icg)) > nTol) ) THEN
              pCvSpec(iCvSpecProducts,icg)  = 1.0_RFREAL &
                                              - pCvSpec(iCvSpecAir,icg)
              pCvSpec(iCvSpecExplosive,icg) = 0.0_RFREAL
              
! ----------- Delete all history for this cell species -------------------------              
              DO iVarScal = 1,pSpecInput%nSpecies
                pCvOldSpec(iVarScal,icg)  = pCvSpec(iVarScal,icg)
                pRhsSpec(iVarScal,icg)    = 0.0_RFREAL
                pRhsSumSpec(iVarScal,icg) = 0.0_RFREAL
              END DO ! iVarScal

            END IF ! flameLoc
! END TEMPORARY
! ------------------------------------------        

! ------- Outside initial explosive region -------------------------------------
          ELSE
            IF ( ABS(pCvSpec(iCvSpecExplosive,icg)) > nTol ) THEN
              pCvSpec(iCvSpecProducts,icg)  = 1.0_RFREAL - pCvSpec(iCvSpecAir,icg)
              pCvSpec(iCvSpecExplosive,icg) = 0.0_RFREAL
            END IF ! pCvSpec
          END IF
        END DO ! icg

        IF ( flameLoc - flameWidth > xFinal ) THEN
          global%pbaBurnFlag = .FALSE.
        END IF ! flameLoc

! ==============================================================================
!     PBA 1D planar case
! ==============================================================================

      CASE ( "pba1d" )
        iCvSpecExplosive = SPEC_GetSpeciesIndex(global,pSpecInput,'EXPLOSIVE')
        iCvSpecProducts  = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')

        mwProducts = pSpecInput%specType(iCvSpecProducts)%pMaterial%molw
        cpProducts = pSpecInput%specType(iCvSpecProducts)%pMaterial%spht

        gcProducts = MixtPerf_R_M(mwProducts)
        gProducts  = MixtPerf_G_CpR(cpProducts,gcProducts)

        ro = pMixtInput%prepRealVal6

! DEBUG: Manoj-PBA1D, Notes: To match Jianghui's implementation
!                         1. Remove volume averaged Y
!               Adding use of volume averaged Y
! END DEBUG

        detonVel   = pSpecInput%specType(iCvSpecExplosive)%pMaterial%detonVel
        xInit      = pMixtInput%prepRealVal1
        flameLoc   = xInit + detonVel*(global%currentTimeRK)
        flameWidth = pMixtInput%prepRealVal8

        xFront = xInit + detonVel*(global%currentTimeRK)
        xWidth = pMixtInput%prepRealVal8
        xTail  = xFront - xWidth
!        xHalfW = 0.0_RFREAL
        xHalfW = 5.0E-05_RFREAL

        DO icg = 1,pGrid%nCellsTot
          x = pGrid%cofg(XCOORD,icg)

          pCvSpec(iCvSpecProducts,icg)  = SPEC_RFLU_PBA_ComputeY(x,xFront,xHalfW,xWidth,xTail)
          pCvSpec(iCvSpecExplosive,icg) = 1.0_RFREAL - pCvSpec(iCvSpecProducts,icg)
        END DO ! icg

! ==============================================================================
!     Default
! ==============================================================================  
  
      CASE DEFAULT     
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! global%casename        

! ******************************************************************************
!   Convert species back to conservative form
! ******************************************************************************

    CALL RFLU_ScalarConvertCvPrim2cons(pRegion,pRegion%spec%cv, &
                                       pRegion%spec%cvState)

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE SPEC_RFLU_PBA_ProgramBurn








! ******************************************************************************
!
! Purpose: Set solution in reaction zone.
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

  SUBROUTINE SPEC_RFLU_PBA_ReactionZone(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
   
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================

    CHARACTER(CHRLEN) :: errorString
    INTEGER :: icg,iCvSpecAir,iCvSpecExplosive,iCvSpecProducts,indCp,indMol
    REAL(RFREAL) :: a,cpProducts,d,detonVel,e,Eo,flameLoc,flameWidth, &
                    gcProducts,gProducts,mwProducts,nTol,p,po,r,ro,T,u,v,Vm2, &
                    w,x,xFinal,xInit,y,YProducts
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvMixt,pCvSpec,pDvMixt,pGvMixt
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_mixt_input), POINTER :: pMixtInput
    TYPE(t_spec_input), POINTER :: pSpecInput

! ******************************************************************************
!   Start, if pbaBurnFlag is FALSE then no need to go further
! ******************************************************************************

    global => pRegion%global

    IF ( global%pbaBurnFlag .EQV. .FALSE. ) THEN
      RETURN
    END IF ! flameLoc

    CALL RegisterFunction(global,'SPEC_RFLU_PBA_ReactionZone',__FILE__)

! ******************************************************************************
!   Set pointers and variables 
! ******************************************************************************

    pGrid      => pRegion%grid
    pCvMixt    => pRegion%mixt%cv
    pCvSpec    => pRegion%spec%cv
    pDvMixt    => pRegion%mixt%dv
    pGvMixt    => pRegion%mixt%gv
    pMixtInput => pRegion%mixtInput
    pSpecInput => pRegion%specInput

    indCp  = pRegion%mixtInput%indCp
    indMol = pRegion%mixtInput%indMol

    nTol = 1.0E-14_RFREAL

! ******************************************************************************
!   Need species in primitive form
! ******************************************************************************

    CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv, &
                                       pRegion%spec%cvState)

! ******************************************************************************
!   Use program burn algorithm based on case name
! ******************************************************************************

    SELECT CASE ( global%casename )

! ==============================================================================
!     Jenkins 2D planar case
! ==============================================================================

      CASE ( "jenkins2d" )
        iCvSpecAir       = SPEC_GetSpeciesIndex(global,pSpecInput,'AIR')
        iCvSpecExplosive = SPEC_GetSpeciesIndex(global,pSpecInput,'EXPLOSIVE')
        iCvSpecProducts  = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')

        mwProducts = pSpecInput%specType(iCvSpecProducts)%pMaterial%molw
        cpProducts = pSpecInput%specType(iCvSpecProducts)%pMaterial%spht

        gcProducts = MixtPerf_R_M(mwProducts)
        gProducts  = MixtPerf_G_CpR(cpProducts,gcProducts)

        ro = pMixtInput%prepRealVal6

        detonVel   = pSpecInput%specType(iCvSpecExplosive)%pMaterial%detonVel
        xInit      = pMixtInput%prepRealVal1
        xFinal     = pMixtInput%prepRealVal2
        flameLoc   = xInit + detonVel*(global%currentTime)
        flameWidth = pMixtInput%prepRealVal8

        DO icg = 1,pGrid%nCellsTot
          x = pGrid%cofg(XCOORD,icg)
          y = pGrid%cofg(YCOORD,icg)

          IF ( x >= pMixtInput%prepRealVal1 .AND. &
               x <= pMixtInput%prepRealVal2 .AND. &
               y >= pMixtInput%prepRealVal3 .AND. &
               y <= pMixtInput%prepRealVal4 ) THEN             
            IF ( (flameLoc > pGrid%cofg(XCOORD,icg)) .AND. &
                 (flameLoc - flameWidth < pGrid%cofg(XCOORD,icg)) ) THEN
!             Do nothing for now.
            END IF ! flameLoc
          END IF
        END DO ! icg

! ==============================================================================
!     PBA 1D planar case
! ==============================================================================

      CASE ( "pba1d" )
        iCvSpecExplosive = SPEC_GetSpeciesIndex(global,pSpecInput,'EXPLOSIVE')
        iCvSpecProducts  = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')

        mwProducts = pSpecInput%specType(iCvSpecProducts)%pMaterial%molw
        cpProducts = pSpecInput%specType(iCvSpecProducts)%pMaterial%spht

        gcProducts = MixtPerf_R_M(mwProducts)
        gProducts  = MixtPerf_G_CpR(cpProducts,gcProducts)

        ro = pMixtInput%prepRealVal6
        po = pMixtInput%prepRealVal7

        detonVel   = pSpecInput%specType(iCvSpecExplosive)%pMaterial%detonVel
        xInit      = pMixtInput%prepRealVal1
        flameLoc   = xInit + detonVel*(global%currentTime)
        flameWidth = pMixtInput%prepRealVal8

        DO icg = 1,pGrid%nCellsTot
          x = pGrid%cofg(XCOORD,icg)

          IF ( x <= flameLoc-flameWidth ) THEN
!           Region before reaction zone. Do nothing                  
          ELSEIF ( x >= flameLoc ) THEN
!           Region after reaction zone. Do nothing                  
          ELSE
            v = 0.0_RFREAL
            w = 0.0_RFREAL

            YProducts = pCvSpec(iCvSpecProducts,icg)
 
            CALL SPEC_RFLU_PBA_GetSolution(gProducts,gcProducts,ro,po,detonVel, &
                                           YProducts,e,p,d,T,u)

            pCvMixt(CV_MIXT_DENS,icg) = d
            pCvMixt(CV_MIXT_XMOM,icg) = d*u
            pCvMixt(CV_MIXT_YMOM,icg) = d*v
            pCvMixt(CV_MIXT_ZMOM,icg) = d*w
            pCvMixt(CV_MIXT_ENER,icg) = d*e+0.5_RFREAL*d*(u*u+v*v+w*w)
          END IF
        END DO ! icg

! ==============================================================================
!     Default
! ==============================================================================  
  
      CASE DEFAULT     
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! global%casename        

! ******************************************************************************
!   Convert species back to conservative form
! ******************************************************************************

    CALL RFLU_ScalarConvertCvPrim2cons(pRegion,pRegion%spec%cv, &
                                       pRegion%spec%cvState)

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE SPEC_RFLU_PBA_ReactionZone








! ******************************************************************************
!
! Purpose: Set rhs for explosive to zero.
!
! Description: None.
!
! Input:
!   region       Data of current region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE SPEC_RFLU_PBA_ZeroResiduals(region)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
   
    TYPE(t_region), TARGET :: region
    
! ==============================================================================
!   Locals
! ==============================================================================

    CHARACTER(CHRLEN) :: errorString
    INTEGER :: icg,iCvSpecExplosive,iVarScal
    REAL(RFREAL) :: nTol,x,y,YExplosive
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvMixt,pCvSpec,pRhsMixt, &
                                             pRhsSpec,pRhsSumSpec
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_mixt_input), POINTER :: pMixtInput
    TYPE(t_spec_input), POINTER :: pSpecInput
    TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
!   Start, if pbaBurnFlag is FALSE then no need to go further
! ******************************************************************************

    pRegion => region%pRegion
    global  => pRegion%global

    IF ( global%pbaBurnFlag .EQV. .FALSE. ) THEN
      RETURN
    END IF ! flameLoc

    CALL RegisterFunction(global,'SPEC_RFLU_PBA_ZeroResiduals',__FILE__)

! ******************************************************************************
!   Set pointers and variables 
! ******************************************************************************

    pGrid       => pRegion%grid
    pCvMixt     => pRegion%mixt%cv
    pCvSpec     => pRegion%spec%cv
    pMixtInput  => pRegion%mixtInput
    pSpecInput  => pRegion%specInput
    pRhsMixt    => pRegion%mixt%rhs
    pRhsSpec    => pRegion%spec%rhs
    pRhsSumSpec => pRegion%spec%rhsSum

    nTol = 1.0E-14_RFREAL

! ******************************************************************************
!   Use program burn algorithm based on case name
! ******************************************************************************

    SELECT CASE ( global%casename )

! ==============================================================================
!     Jenkins 2D planar case
! ==============================================================================

      CASE ( "jenkins2d" )
        iCvSpecExplosive = SPEC_GetSpeciesIndex(global,pSpecInput,'EXPLOSIVE')

        DO icg = 1,pGrid%nCellsTot
          x = pGrid%cofg(XCOORD,icg)
          y = pGrid%cofg(YCOORD,icg)

          IF ( x >= pMixtInput%prepRealVal1 .AND. &
               x <= pMixtInput%prepRealVal2 .AND. &
               y >= pMixtInput%prepRealVal3 .AND. &
               y <= pMixtInput%prepRealVal4 ) THEN             
            YExplosive = pCvSpec(iCvSpecExplosive,icg)/pCvMixt(CV_MIXT_DENS,icg)

! --------- Mixture Rhs in pure explosive should be zero -----------------------
            IF ( ABS(1.0_RFREAL - YExplosive) < nTol ) THEN
              pRhsMixt(CV_MIXT_DENS,icg) = 0.0_RFREAL 
              pRhsMixt(CV_MIXT_XMOM,icg) = 0.0_RFREAL
              pRhsMixt(CV_MIXT_YMOM,icg) = 0.0_RFREAL
              pRhsMixt(CV_MIXT_ZMOM,icg) = 0.0_RFREAL
              pRhsMixt(CV_MIXT_ENER,icg) = 0.0_RFREAL
            END IF ! YExplosive

! --------- Speceies Rhs in reaction zone and pure explosive should be zero ----
            IF ( ABS(YExplosive) > nTol ) THEN
              DO iVarScal = 1,pSpecInput%nSpecies
                pRhsSpec(iVarScal,icg)    = 0.0_RFREAL
                pRhsSumSpec(iVarScal,icg) = 0.0_RFREAL
              END DO ! iVarScal
            END IF ! YExplosive
          END IF ! x,y
        END DO ! icg

! ==============================================================================
!     PBA 1D planar case
! ==============================================================================

      CASE ( "pba1d" )
        iCvSpecExplosive = SPEC_GetSpeciesIndex(global,pSpecInput,'EXPLOSIVE')

        DO icg = 1,pGrid%nCellsTot
          x = pGrid%cofg(XCOORD,icg)

          YExplosive = pCvSpec(iCvSpecExplosive,icg)/pCvMixt(CV_MIXT_DENS,icg)

! ------- Mixture Rhs in pure explosive should be zero -------------------------
          IF ( ABS(1.0_RFREAL - YExplosive) < nTol ) THEN
            pRhsMixt(CV_MIXT_DENS,icg) = 0.0_RFREAL 
            pRhsMixt(CV_MIXT_XMOM,icg) = 0.0_RFREAL
            pRhsMixt(CV_MIXT_YMOM,icg) = 0.0_RFREAL
            pRhsMixt(CV_MIXT_ZMOM,icg) = 0.0_RFREAL
            pRhsMixt(CV_MIXT_ENER,icg) = 0.0_RFREAL
          END IF ! YExplosive

! ------- Species Rhs in reaction zone and pure explosive should be zero -------
          IF ( ABS(YExplosive) > nTol ) THEN
            DO iVarScal = 1,pSpecInput%nSpecies
              pRhsSpec(iVarScal,icg) = 0.0_RFREAL
            END DO ! iVarScal
          END IF ! YExplosive
        END DO ! icg

! ==============================================================================
!     Default
! ==============================================================================  
  
      CASE DEFAULT     
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! global%casename        

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE SPEC_RFLU_PBA_ZeroResiduals








END MODULE SPEC_RFLU_ModPBA

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_RFLU_ModPBA.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
! ******************************************************************************

