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
! Purpose: Compute first-order accurate discretization of scalar inviscid flux.
!
! Description: None.
!
! Input: 
!  pRegion              Pointer to region data
!  nVarScal             Number of scalar variables
!  cvScal               Conserved scalar variables
!  resScal              Residual due to central scalar fluxes
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_ScalarFirst.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_ScalarFirst(pRegion,nVarScal,cvScal,resScal)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModParameters
    
  IMPLICIT NONE

! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: nVarScal
  REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cvScal
  REAL(RFREAL), DIMENSION(:,:), INTENT(INOUT) :: resScal
  TYPE(t_region), POINTER :: pRegion

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: c1,c2,ifc,iVarScal
  REAL(RFREAL) :: flx,mf,mfn,mfp,sl,sr
  REAL(RFREAL), DIMENSION(:), POINTER :: pMf  
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ScalarFirst.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ScalarFirst',__FILE__)
  
! *****************************************************************************
! Checks: Defensive coding, should never occur
! *****************************************************************************
 
  IF ( pRegion%mixtInput%indMfMixt /= 1 ) THEN 
    CALL ErrorStop(global,ERR_INDMFMIXT_INVALID,__LINE__)
  END IF ! pRegion%mixtInput%indMfMixt

! *****************************************************************************
! Set dimensions and pointers 
! *****************************************************************************

  pGrid => pRegion%grid
  pMf   => pRegion%mixt%mfMixt

! *****************************************************************************
! Compute fluxes
! *****************************************************************************

  DO ifc = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifc)
    c2 = pGrid%f2c(2,ifc) 

! =============================================================================  
!   Get mass flux 
! =============================================================================  

    mf  = pMf(ifc)
    mfp = MAX(mf,0.0_RFREAL)
    mfn = MIN(mf,0.0_RFREAL)
         
! =============================================================================  
!   Compute flux and accumulate into residual
! =============================================================================  
     
    DO iVarScal = 1,nVarScal
      sl = cvScal(iVarScal,c1) 
      sr = cvScal(iVarScal,c2) 
     
      flx = mfp*sl + mfn*sr

      resScal(iVarScal,c1) = resScal(iVarScal,c1) + flx
      resScal(iVarScal,c2) = resScal(iVarScal,c2) - flx
    END DO ! iVarScal
  END DO  ! ifc
  
! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ScalarFirst

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ScalarFirst.F90,v $
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
! Revision 1.1  2004/01/29 22:56:07  haselbac
! Initial revision
!
!******************************************************************************

