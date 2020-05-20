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
! Purpose: Compute first-order accurate discretization of scalar inviscid flux 
!   through boundary faces.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to data of current region
!   pPatch      Pointer to data of current patch
!   nVarScal    Number of scalars
!   cvScal      Vector of conserved scalar variables
!   valScal     Boundary values of scalar variables     
!
! Output: 
!   resScal     Residual of scalar variables
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ScalarFirstPatch.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ScalarFirstPatch(pRegion,pPatch,nVarScal,cvScal,valScal,resScal)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModBndPatch, ONLY: t_bcvalues,t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModParameters
    
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: nVarScal
  REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cvScal
  REAL(RFREAL), DIMENSION(:,:), INTENT(INOUT) :: resScal
  TYPE(t_bcvalues) :: valScal
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: c1,bcType,distScal,ifc,iVarScal
  REAL(RFREAL) :: flx,mf
  REAL(RFREAL), DIMENSION(:), POINTER :: pMfMixt
  REAL(RFREAL), DIMENSION(:,:), POINTER :: rhs
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  
! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ScalarFirstPatch.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ScalarFirstPatch',__FILE__)

! ******************************************************************************
! Checks: Defensive coding, should never occur
! ******************************************************************************
 
  IF ( pRegion%mixtInput%indMfMixt /= 1 ) THEN 
    CALL ErrorStop(global,ERR_INDMFMIXT_INVALID,__LINE__)
  END IF ! pRegion%mixtInput%indMfMixt
  
! ******************************************************************************
! Set pointers and variables
! ******************************************************************************
      
  pGrid   => pRegion%grid    
  pMfMixt => pPatch%mfMixt
  
  bcType   = pPatch%bcType  
  distScal = valScal%distrib  

! ******************************************************************************
! Select boundary type
! ******************************************************************************

  SELECT CASE ( bcType )

! ==============================================================================  
!   Inflow
! ==============================================================================

    CASE ( BC_INFLOW_TOTANG,BC_INFLOW_VELTEMP ) 
      DO ifc = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifc)

        mf = pMfMixt(ifc) 
                                                                
        DO iVarScal = 1,nVarScal                                      
          flx = mf*valScal%vals(iVarScal,distScal*ifc)
          
          resScal(iVarScal,c1) = resScal(iVarScal,c1) + flx
        END DO ! iVarScal
      END DO ! ifc
            
! ==============================================================================  
!   Outflow
! ==============================================================================

    CASE ( BC_OUTFLOW )        
      DO ifc = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifc)
        
        mf = pMfMixt(ifc)                   
                                          
        DO iVarScal = 1,nVarScal                
          flx = mf*cvScal(iVarScal,c1)
                                  
          resScal(iVarScal,c1) = resScal(iVarScal,c1) + flx
        END DO ! iVarScal 
      END DO ! ifc

! ==============================================================================  
!   Slip wall 
! ==============================================================================
    
    CASE ( BC_SLIPWALL ) 

! ==============================================================================  
!   No-slip wall
! ==============================================================================
 
    CASE ( BC_NOSLIPWALL_HFLUX,BC_NOSLIPWALL_TEMP ) 

! ==============================================================================  
!   Farfield
! ==============================================================================

    CASE ( BC_FARFIELD )
      DO ifc = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifc)                                              

        mf = pMfMixt(ifc)
                   
        IF ( mf > 0.0_RFREAL ) THEN ! Outflow
          DO iVarScal = 1,nVarScal        
            flx = mf*cvScal(iVarScal,c1)
                        
            resScal(iVarScal,c1) = resScal(iVarScal,c1) + flx
          END DO ! iVarScal                      
        ELSE ! Inflow
          DO iVarScal = 1,nVarScal 
            flx = mf*valScal%vals(iVarScal,distScal*ifc)
                      
            resScal(iVarScal,c1) = resScal(iVarScal,c1) + flx
          END DO ! iVarScal
        END IF ! mf   
      END DO ! ifc      
    
! ==============================================================================  
!   Injection
! ==============================================================================

    CASE ( BC_INJECTION )
      DO ifc = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifc)

        mf = pMfMixt(ifc)
                                                
        DO iVarScal = 1,nVarScal                
          flx = mf*valScal%vals(iVarScal,distScal*ifc)
                                  
          resScal(iVarScal,c1) = resScal(iVarScal,c1) + flx
        END DO ! iVarScal                           
      END DO ! ifc

! ==============================================================================
!   Boundaries for which fluxes must not or need not be computed
! ==============================================================================

    CASE ( BC_PERIODIC, &
           BC_SYMMETRY, & 
           BC_VIRTUAL )
           
! ==============================================================================  
!   Default
! ==============================================================================

    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! bcType

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ScalarFirstPatch

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ScalarFirstPatch.F90,v $
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
! Revision 1.5  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.4  2006/03/25 21:43:47  haselbac
! Added CASEs for sype patches
!
! Revision 1.3  2005/11/10 02:03:16  haselbac
! Added virtual boundary, cleaned up CASE statements
!
! Revision 1.2  2005/04/20 14:40:15  haselbac
! Removed CHECK_UNIFLOW code section, cosmetics
!
! Revision 1.1  2004/01/29 22:56:10  haselbac
! Initial revision
!
! ******************************************************************************

