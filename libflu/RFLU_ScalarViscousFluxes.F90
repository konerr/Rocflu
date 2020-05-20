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
! Purpose: Compute viscous fluxes of scalar for actual faces.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!  nVarScal     Number of scalar variables
!  tvScal       Scalar transport variables
!  gradScal     Scalar face gradients
!  resScal      Scalar residuals
!
! Output: None.
!
! Notes: 
!   1. The viscosity must be a dynamic viscosity!
!
!******************************************************************************
!
! $Id: RFLU_ScalarViscousFluxes.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_ScalarViscousFluxes(pRegion,nVarScal,tvScal,gradScal,resScal)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
   
  IMPLICIT NONE

! *****************************************************************************
! Declarations
! *****************************************************************************

! =============================================================================  
! Arguments 
! =============================================================================  

  INTEGER, INTENT(IN) :: nVarScal
  REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: tvScal
  REAL(RFREAL), DIMENSION(:,:), INTENT(INOUT) :: resScal 
  REAL(RFREAL), DIMENSION(:,:,:), INTENT(IN) :: gradScal   
  TYPE(t_region), POINTER :: pRegion

! =============================================================================  
! Locals
! =============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: c1,c2,ifg,iVarScal
  REAL(RFREAL) :: beta,flx,mu,nm,nx,ny,nz
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ScalarViscousFluxes.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ScalarViscousFluxes',__FILE__)

! *****************************************************************************
! Set variables and pointers
! *****************************************************************************

  pGrid => pRegion%grid
    
  beta = pRegion%mixtInput%betrk(pRegion%irkStep)  
  
! *****************************************************************************  
! Loop over faces and compute viscous fluxes
! *****************************************************************************

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

! =============================================================================   
!   Get face geometry
! =============================================================================    
    
    nx = pGrid%fn(XCOORD,ifg)
    ny = pGrid%fn(YCOORD,ifg)
    nz = pGrid%fn(ZCOORD,ifg) 
    nm = pGrid%fn(XYZMAG,ifg)        

! =============================================================================   
!   Get face state. NOTE this simple average is not accurate for non-uniform 
!   grids, so this will have to be replaced by a proper interpolation.
! =============================================================================

    DO iVarScal = 1,nVarScal
      mu = 0.5_RFREAL*(tvScal(iVarScal,c1) + tvScal(iVarScal,c2)) 
     
! -----------------------------------------------------------------------------
!     Compute fluxes
! -----------------------------------------------------------------------------
    
      flx = mu*(  gradScal(XCOORD,iVarScal,ifg)*nx     & 
                + gradScal(YCOORD,iVarScal,ifg)*ny     & 
                + gradScal(ZCOORD,iVarScal,ifg)*nz)*nm
  
! -----------------------------------------------------------------------------
!     Accumulate into residual     
! -----------------------------------------------------------------------------

      resScal(iVarScal,c1) = resScal(iVarScal,c1) + beta*flx    
      resScal(iVarScal,c2) = resScal(iVarScal,c2) - beta*flx
    END DO ! iVarScal
  END DO ! ifg

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ScalarViscousFluxes

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ScalarViscousFluxes.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:36  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:53  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.1  2004/01/29 22:56:17  haselbac
! Initial revision
!
!******************************************************************************

