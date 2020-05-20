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
! Purpose: Suite of routines to implement axisymmetry.
!
! Description: None.
!
! Notes: 
!
! ******************************************************************************
!
! $Id: RFLU_ModAxisymmetry.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2007-2008 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModAxisymmetry

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  IMPLICIT NONE
   
  PRIVATE
  PUBLIC :: RFLU_AXI_DescaleGeometry, &
            RFLU_AXI_SourceTermsNSCBC, &
            RFLU_AXI_ScaleGeometry, & 
            RFLU_AXI_SourceTerms 
                                   
  SAVE    
     
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
     
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModAxisymmetry.F90,v $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
 





! ******************************************************************************
!
! Purpose: Descale faces and volume of mesh by ycoord of centroid. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Scaling/Descaling of geometry must be done with ABS(y) coord of centroid. 
!
! ******************************************************************************
    
  SUBROUTINE RFLU_AXI_DescaleGeometry(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    INTEGER :: errorFlag,ic,ifc,iPatch
    REAL(RFREAL) :: factor
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_AXI_DescaleGeometry',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Descaling geometry ...' 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel

    pGrid => pRegion%grid  

! ==============================================================================
!   Descale cell volumes
! ==============================================================================

    DO ic = 1,pGrid%nCellsTot
      factor        = 1.0_RFREAL/(ABS(pGrid%cofg(YCOORD,ic))+1.0E-15_RFREAL) 
      pGrid%vol(ic) = factor*pGrid%vol(ic)
    END DO ! ic

! ==============================================================================
!   Descale interior face areas
! ==============================================================================

    DO ifc = 1,pGrid%nFacesTot 
      factor             = 1.0_RFREAL/(ABS(pGrid%fc(YCOORD,ifc))+1.0E-15_RFREAL)
      pGrid%fn(XYZMAG,ifc) = factor*pGrid%fn(XYZMAG,ifc) 
    END DO ! ifc      

! ==============================================================================
!   Descale patch face areas
! ==============================================================================

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)  

      DO ifc = 1,pPatch%nBFacesTot 
        factor          = 1.0_RFREAL/(ABS(pPatch%fc(YCOORD,ifc))+1.0E-15_RFREAL)
        pPatch%fn(XYZMAG,ifc) = factor*pPatch%fn(XYZMAG,ifc)                    
      END DO ! ifc 
    END DO ! iPatch                  

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Descaling geometry done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_AXI_DescaleGeometry








! ******************************************************************************
!
! Purpose: Computing the source terms for flow equations solved at boundary in 
!          NSCBC implementation in axisymmetric framework. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Source terms must be compted with y coord of centroid and not with ABS(y) 
!   2. Equations are solved with finite difference discretization and scaling of
!      geometry does not affect source terms.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_AXI_SourceTermsNSCBC(pRegion,pPatch)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    INTEGER :: errorFlag,ifl
    REAL(RFREAL) :: divTerm,dudx,dvdy,dwdz,p,r,rE,ru,rv,rw,S_Rho,S_RhoU, &
                    S_RhoV,S_RhoW,S_RhoE,v,visc,y
    REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv,pDv,pRhs
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_AXI_SourceTermsNSCBC',__FILE__)

    pCv  => pPatch%mixt%cv
    pDv  => pPatch%mixt%dv
    pRhs => pPatch%mixt%rhs

! ==============================================================================
!   Loop over patch faces and compute inviscid source terms 
! ==============================================================================

    DO ifl = 1,pPatch%nBFaces
      y  = pPatch%fc(YCOORD,ifl) 

      r  = pCv(CV_MIXT_DENS,ifl)
      ru = pCv(CV_MIXT_XMOM,ifl)
      rv = pCv(CV_MIXT_YMOM,ifl)
      rw = pCv(CV_MIXT_ZMOM,ifl)
      rE = pCv(CV_MIXT_ENER,ifl)
      p  = pDv(DV_MIXT_PRES,ifl)
      v  = rv/r

      S_rho  = -r*v/y
      S_rhoU = -ru*v/y
      S_rhoV = -rv*v/y
      S_rhoW = -rw*v/y
      S_rhoE = -(rE+p)*v/y 

      pRhs(CV_MIXT_DENS,ifl) = pRhs(CV_MIXT_DENS,ifl) + S_Rho 
      pRhs(CV_MIXT_XMOM,ifl) = pRhs(CV_MIXT_XMOM,ifl) + S_RhoU 
      pRhs(CV_MIXT_YMOM,ifl) = pRhs(CV_MIXT_YMOM,ifl) + S_RhoV 
      pRhs(CV_MIXT_ZMOM,ifl) = pRhs(CV_MIXT_ZMOM,ifl) + S_RhoW 
      pRhs(CV_MIXT_ENER,ifl) = pRhs(CV_MIXT_ENER,ifl) + S_RhoE 
    END DO ! ifl

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_AXI_SourceTermsNSCBC








! ******************************************************************************
!
! Purpose: Scale faces and volume of mesh by ycoord of centroid. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Scaling/Descaling of geometry must be done with ABS(y) coord of centroid. 
!
! ******************************************************************************
    
  SUBROUTINE RFLU_AXI_ScaleGeometry(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    INTEGER :: errorFlag,ic,ifc,iPatch
    REAL(RFREAL) :: factor
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_AXI_ScaleGeometry',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Scaling geometry ...' 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel

    pGrid => pRegion%grid  

! ==============================================================================
!   Scale cell volumes
! ==============================================================================

    IF ( (global%solverType == SOLV_IMPLICIT_HM) .AND. &
         (pRegion%mixtInput%axiFlag .EQV. .TRUE.) ) THEN
      DO ic = 1,pGrid%nCellsTot
        pGrid%volus(ic) = pGrid%vol(ic)
      END DO ! ic
    END IF ! global%solverType

    DO ic = 1,pGrid%nCellsTot
      factor          = ABS(pGrid%cofg(YCOORD,ic))
      pGrid%vol(ic)   = factor*pGrid%vol(ic)
    END DO ! ic

! ==============================================================================
!   Scale interior face areas
! ==============================================================================

    IF ( (global%solverType == SOLV_IMPLICIT_HM) .AND. &
         (pRegion%mixtInput%axiFlag .EQV. .TRUE.) ) THEN
      DO ifc = 1,pGrid%nFacesTot 
        pGrid%fnmus(ifc) = pGrid%fn(XYZMAG,ifc) 
      END DO ! ifc      
    END IF ! global%solverType

    DO ifc = 1,pGrid%nFacesTot 
      factor               = ABS(pGrid%fc(YCOORD,ifc))
      pGrid%fn(XYZMAG,ifc) = factor*pGrid%fn(XYZMAG,ifc) 
    END DO ! ifc      

! ==============================================================================
!   Scale patch face areas
! ==============================================================================

    IF ( (global%solverType == SOLV_IMPLICIT_HM) .AND. &
         (pRegion%mixtInput%axiFlag .EQV. .TRUE.) ) THEN
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)  

        DO ifc = 1,pPatch%nBFacesTot 
          pPatch%fnmus(ifc)     = pPatch%fn(XYZMAG,ifc)                    
        END DO ! ifc 
      END DO ! iPatch                  
    END IF ! global%solverType

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)  

      DO ifc = 1,pPatch%nBFacesTot 
        factor                = ABS(pPatch%fc(YCOORD,ifc))
        pPatch%fn(XYZMAG,ifc) = factor*pPatch%fn(XYZMAG,ifc)                    
      END DO ! ifc 
    END DO ! iPatch                  

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Scaling geometry done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_AXI_ScaleGeometry






! ******************************************************************************
!
! Purpose: Computing the source terms due to solving axisymmetric equations. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Source terms must be compted with y coord of centroid and not with ABS(y) 
!
! ******************************************************************************
    
  SUBROUTINE RFLU_AXI_SourceTerms(region)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_region), TARGET :: region

! ==============================================================================
!   Locals      
! ==============================================================================

    INTEGER :: errorFlag,ic,ifc,iPatch
    REAL(RFREAL), PARAMETER :: TWO_THIRDS = 2.0_RFREAL/3.0_RFREAL
    REAL(RFREAL) :: divTerm,dudx,dvdy,dwdz,Sy,v,visc,y
    TYPE(t_global), POINTER :: global
    TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
!   Start
! ******************************************************************************

    pRegion => region
    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_AXI_SourceTerms',__FILE__)

! ==============================================================================
!   Loop over cells and compute inviscid source terms 
! ==============================================================================

    DO ic = 1,pRegion%grid%nCellsTot
      y  = pRegion%grid%cofg(YCOORD,ic)

      Sy = pRegion%mixt%dv(DV_MIXT_PRES,ic)/y

      pRegion%mixt%rhs(CV_MIXT_YMOM,ic) = pRegion%mixt%rhs(CV_MIXT_YMOM,ic) &
                                        - Sy*pRegion%grid%vol(ic)
    END DO ! ic

! ==============================================================================
!   Loop over cells and compute viscous source terms 
! ==============================================================================

    IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN
      DO ic = 1,pRegion%grid%nCellsTot
        visc = pRegion%mixt%tv(TV_MIXT_MUEL,ic)

        y    = pRegion%grid%cofg(YCOORD,ic)
        v    = pRegion%mixt%cv(CV_MIXT_YMOM,ic)/pRegion%mixt%cv(CV_MIXT_DENS,ic)

        dudx = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_XVEL,ic)
        dvdy = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_YVEL,ic)
        dwdz = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_ZVEL,ic)

        divTerm = TWO_THIRDS*(dudx + dvdy + dwdz + v/y)

        Sy = -visc*(2.0_RFREAL*v/y - divTerm)/y

        pRegion%mixt%rhs(CV_MIXT_YMOM,ic) = pRegion%mixt%rhs(CV_MIXT_YMOM,ic) &
                                          - Sy*pRegion%grid%vol(ic)
      END DO ! ic
    END IF ! flowModel

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_AXI_SourceTerms



! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModAxisymmetry


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModAxisymmetry.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.6  2009/09/28 14:21:34  mparmar
! Added initialization of unscaled geometry
!
! Revision 1.5  2009/07/08 19:11:44  mparmar
! Added RFLU_AXI_SourceTermsNSCBC
!
! Revision 1.4  2008/12/06 08:43:39  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:16:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2008/05/29 01:35:26  mparmar
! Adapted axi-symmetric computations for Navier Stokes equations
!
! Revision 1.1  2008/03/27 12:03:34  haselbac
! Initial revision
!
!
! ******************************************************************************

