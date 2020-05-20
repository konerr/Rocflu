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
! Purpose: Correct the mixture properties using higher-order
!          interpolation schemed onto particle locations
!
! Description: none.
!
! Input: pRegion = current region.
!
! Output: region%levels%plag%dv 
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_RFLU_CorrectMixtProperties.F90,v 1.3 2016/02/21 17:18:16 fred Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_CorrectMixtProperties(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModMixture, ONLY    : t_mixt
  USE ModError  
  USE ModParameters
  USE PLAG_ModParameters

  USE ModInterfaces, ONLY: MixtPerf_R_M, & 
                           MixtPerf_T_DPR, &
                           MixtPerf_G_CpR
  USE RFLU_ModConvertCv, ONLY: RFLU_ScalarConvertCvCons2Prim, &
                               RFLU_ScalarConvertCvPrim2Cons
  USE RFLU_ModJWL

#ifdef SPEC
  USE ModInterfacesSpecies, ONLY: SPEC_GetSpeciesIndex
  USE ModSpecies,           ONLY: t_spec_input
#endif
!Fred - added JWL EOS capabilities - 9/22/15

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: icg,indMol,intrplMixtModel,iPcl,nPcls
  INTEGER, POINTER, DIMENSION(:,:)    :: pAiv
  
  REAL(RFREAL) :: dx,dy,dz,gc,mm,p,r,T,a,e,g,gcAir,ir,cpAir
  REAL(RFREAL), POINTER, DIMENSION(:,:)    :: pCv,pDv,pDvMixt,pGvMixt
  REAL(RFREAL), POINTER, DIMENSION(:,:,:) :: pGrad
  
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag),   POINTER :: pPlag 
  TYPE(t_mixt),   POINTER :: pMixt 

#ifdef SPEC
   LOGICAL :: scalarConvFlag
   INTEGER :: iCvSpecProducts,iCvSpecAir
   REAL(RFREAL) :: Yproducts
   REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
   TYPE(t_spec_input), POINTER :: pSpecInput
#endif

 
! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_CorrectMixtProperties.F90,v $ $Revision: 1.3 $'

  global => pRegion%global
  
  CALL RegisterFunction( global, 'PLAG_RFLU_CorrectMixtProperties',__FILE__ )

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pMixt => pRegion%mixt 
  pGrid => pRegion%grid
  pPlag => pRegion%plag

  pDvMixt => pMixt%dv
  pGvMixt => pMixt%gv

  nPcls  = pRegion%plag%nPcls
  indMol = pRegion%mixtInput%indMol

  intrplMixtModel = pRegion%plagInput%intrplMixtModel

! ******************************************************************************
! Set pointers for mixture 
! ******************************************************************************

  pDvMixt => pMixt%dv
  pGvMixt => pMixt%gv

  pGrad => pRegion%mixt%gradCell

#ifdef SPEC
   pSpecInput => pRegion%specInput
   pCvSpec => pRegion%spec%cv
#endif


! ******************************************************************************
! Set pointers for discrete particles 
! ******************************************************************************

  pAiv    => pPlag%aiv
  pCv     => pPlag%cv  
  pDv     => pPlag%dv

! ******************************************************************************
! Check that have primitive state vector
! ******************************************************************************

  IF ( pMixt%cvState /= CV_MIXT_STATE_DUVWP ) THEN 
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! pMixt%cvState

! ******************************************************************************
! Trap error for inconsistent input deck
! ******************************************************************************

  IF ( intrplMixtModel /= ZEROTH_ORDER .AND. &
       pRegion%mixtInput%spaceOrder == 1     ) THEN 
    CALL ErrorStop(global,ERR_PLAG_INTRPLMODEL,__LINE__)
  END IF ! intrplMixtModel

! ******************************************************************************        
! Correct discrete particle dv
! ******************************************************************************
    
  SELECT CASE ( intrplMixtModel )    
    CASE ( ZEROTH_ORDER )
    CASE ( FIRST_ORDER  )        
      DO iPcl = 1, nPcls
        icg = pAiv(AIV_PLAG_ICELLS,iPcl)

        dx = pCv(CV_PLAG_XPOS,iPcl) - pGrid%cofg(XCOORD,icg)
        dy = pCv(CV_PLAG_YPOS,iPcl) - pGrid%cofg(YCOORD,icg)
        dz = pCv(CV_PLAG_ZPOS,iPcl) - pGrid%cofg(ZCOORD,icg)                

        pDv(DV_PLAG_UVELMIXT,iPcl) = pDv(DV_PLAG_UVELMIXT,iPcl) & 
                                   + pGrad(XCOORD,GRC_MIXT_XVEL,icg)*dx &
                                   + pGrad(YCOORD,GRC_MIXT_XVEL,icg)*dy &
                                   + pGrad(ZCOORD,GRC_MIXT_XVEL,icg)*dz
        pDv(DV_PLAG_VVELMIXT,iPcl) = pDv(DV_PLAG_VVELMIXT,iPcl) & 
                                   + pGrad(XCOORD,GRC_MIXT_YVEL,icg)*dx &
                                   + pGrad(YCOORD,GRC_MIXT_YVEL,icg)*dy &
                                   + pGrad(ZCOORD,GRC_MIXT_YVEL,icg)*dz
        pDv(DV_PLAG_WVELMIXT,iPcl) = pDv(DV_PLAG_WVELMIXT,iPcl) & 
                                   + pGrad(XCOORD,GRC_MIXT_ZVEL,icg)*dx &
                                   + pGrad(YCOORD,GRC_MIXT_ZVEL,icg)*dy &
                                   + pGrad(ZCOORD,GRC_MIXT_ZVEL,icg)*dz
        pDv(DV_PLAG_PRESMIXT,iPcl) = pDv(DV_PLAG_PRESMIXT,iPcl) & 
                                   + pGrad(XCOORD,GRC_MIXT_PRES,icg)*dx &
                                   + pGrad(YCOORD,GRC_MIXT_PRES,icg)*dy &
                                   + pGrad(ZCOORD,GRC_MIXT_PRES,icg)*dz
        pDv(DV_PLAG_DENSMIXT,iPcl) = pDv(DV_PLAG_DENSMIXT,iPcl) & 
                                   + pGrad(XCOORD,GRC_MIXT_DENS,icg)*dx &
                                   + pGrad(YCOORD,GRC_MIXT_DENS,icg)*dy &
                                   + pGrad(ZCOORD,GRC_MIXT_DENS,icg)*dz 
      END DO ! iPcl
      
      SELECT CASE ( pRegion%mixtInput%gasModel ) 
        CASE ( GAS_MODEL_TCPERF,GAS_MODEL_MIXT_TCPERF )
          DO iPcl = 1,nPcls      
            icg = pAiv(AIV_PLAG_ICELLS,iPcl)
          
            mm = pGvMixt(GV_MIXT_MOL,indMol*icg)
            gc = MixtPerf_R_M(mm)
                  
            r = pDv(DV_PLAG_DENSMIXT,iPcl)
            p = pDv(DV_PLAG_PRESMIXT,iPcl)        
                  
            pDv(DV_PLAG_TEMPMIXT,iPcl) = MixtPerf_T_DPR(r,p,gc)
          END DO ! iPcl
#ifdef SPEC
        CASE (GAS_MODEL_MIXT_JWL)

         IF (global%specUsed .EQV. .TRUE.) THEN
         IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
             scalarConvFlag = .FALSE.
         ELSE
             scalarConvFlag = .TRUE.
    CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
         END IF

          DO iPcl = 1,nPcls
            icg = pAiv(AIV_PLAG_ICELLS,iPcl)

            mm = pGvMixt(GV_MIXT_MOL,indMol*icg)
            gc = MixtPerf_R_M(mm)

            r = pDv(DV_PLAG_DENSMIXT,iPcl)
            p = pDv(DV_PLAG_PRESMIXT,iPcl)

            ir = 1.0_RFREAL/r
            iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')
            iCvSpecAir = SPEC_GetSpeciesIndex(global,pSpecInput,'AIR')
            cpAir = pSpecInput%specType(iCvSpecAir)%pMaterial%spht
            Yproducts = pCvSpec(iCvSpecProducts,icg)

            g = MixtPerf_G_CpR(cpAir,gc)

            CALL RFLU_JWL_ComputeEnergyMixt(pRegion,icg,g,gc,p,r,Yproducts,a,e,T)
            pDv(DV_PLAG_TEMPMIXT,ipcl) = T

          END DO ! iPcl

         IF (scalarConvFlag .EQV. .TRUE.) THEN
  CALL RFLU_ScalarConvertcvPrim2Cons(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
        END IF
      END IF
#endif

        CASE DEFAULT
          CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )
      END SELECT ! pRegion%mixtInput%gasModel
      
    CASE ( SECOND_ORDER )
      CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )
    CASE DEFAULT
      CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )     
  END SELECT  ! intrplMixtModel

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLU_CorrectMixtProperties

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_CorrectMixtProperties.F90,v $
! Revision 1.3  2016/02/21 17:18:16  fred
! Fixing an incorrect precompiler flag
!
! Revision 1.2  2016/02/04 20:00:57  fred
! Adding iterative JWL EOS capabilities for the cylindrical detonation case
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:51  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:50:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.1  2005/11/30 22:26:54  fnajjar
! Initial version
!
!******************************************************************************

