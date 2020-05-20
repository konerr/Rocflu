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
! Purpose: compute convective fluxes, add them to dissipation in order
!          to obtain the residual.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%mixt%rhs = convective fluxes + dissipation.
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: ConvectiveFluxes.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ConvectiveFluxes( region )

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE ModParameters

  USE RFLU_ModNSCBC

  USE ModInterfaces, ONLY: RFLU_ComputeFluxInv

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), TARGET :: region

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: icg,spaceDiscr,spaceOrder
  REAL(RFREAL), POINTER :: diss(:,:), rhs(:,:)
  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
! Start
! ******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'ConvectiveFluxes',__FILE__ )

! ******************************************************************************
! Get dimensions and pointers
! ******************************************************************************

  pRegion => region
  
  diss => region%mixt%diss
  rhs  => region%mixt%rhs  

  spaceDiscr = region%mixtInput%spaceDiscr
  spaceOrder = region%mixtInput%spaceOrder

! ******************************************************************************
! Initialize residual (set = dissipation). NOTE only need to do this for 
! unsteady flows because for steady flows this is now done by a dedicated 
! routine (RFLU_EMS_SetRhs).
! ******************************************************************************

  IF ( global%flowType == FLOW_UNSTEADY ) THEN ! Unsteady
    DO icg = 1,region%grid%nCellsTot
      rhs(CV_MIXT_DENS,icg) = -diss(CV_MIXT_DENS,icg)
      rhs(CV_MIXT_XMOM,icg) = -diss(CV_MIXT_XMOM,icg)
      rhs(CV_MIXT_YMOM,icg) = -diss(CV_MIXT_YMOM,icg)
      rhs(CV_MIXT_ZMOM,icg) = -diss(CV_MIXT_ZMOM,icg)
      rhs(CV_MIXT_ENER,icg) = -diss(CV_MIXT_ENER,icg)
    END DO ! icg
  END IF ! global%flowType

! ******************************************************************************
! Compute fluxes
! ******************************************************************************

! ==============================================================================
! Check state of conserved variable vector
! ==============================================================================

  IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN 
    CALL ErrorStop(pRegion%global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! pRegion%mixt%cvState

! ==============================================================================
! Call flux functions 
! ==============================================================================

  CALL RFLU_ComputeFluxInv(pRegion,FLUX_PART_BOTH)

! ******************************************************************************
! Compute Boundary Patch Rhs 
! ******************************************************************************

  IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
    CALL RFLU_NSCBC_CompRhs(pRegion)
  END IF ! RFLU_NSCBC_DecideHaveNSCBC(pRegion)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( region%global )

END SUBROUTINE ConvectiveFluxes

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ConvectiveFluxes.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:24  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.8  2006/08/19 15:38:26  mparmar
! Added computing of fluxes for boundary variables
!
! Revision 1.7  2006/05/01 20:57:43  haselbac
! Converted to new flux routine
!
! Revision 1.6  2006/03/26 20:21:12  haselbac
! Rewrite bcos of of GL model
!
! Revision 1.5  2005/11/14 16:53:37  haselbac
! Added AUSM flux computation for pseudo-gas
!
! Revision 1.4  2005/11/10 01:55:16  haselbac
! Added new logic for AUSM scheme
!
! Revision 1.3  2005/07/14 21:39:28  haselbac
! Added AUSM flux function
!
! Revision 1.2  2005/05/16 20:37:56  haselbac
! Changed args to pRegion, rmvd OLES call, changed init of rhs, cosmetics
!
! Revision 1.1  2004/12/01 16:48:25  haselbac
! Initial revision after changing case
!
! Revision 1.17  2004/02/13 02:56:23  haselbac
! Added calls (commented out for now) to fast Roe routines
!
! Revision 1.16  2004/01/29 22:52:38  haselbac
! Rewrote for positive species update
!
! Revision 1.15  2003/12/04 03:22:57  haselbac
! Added second-order schemes
!
! Revision 1.14  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.11  2003/05/29 17:28:42  jblazek
! Implemented Roe scheme.
!
! Revision 1.10  2003/05/16 22:05:40  haselbac
! Added HLLC flux
!
! Revision 1.9  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.8  2002/09/09 13:56:51  haselbac
! mixtInput under region, bug fix for iec, added state check
!
! Revision 1.7  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.6  2002/07/29 17:13:08  jblazek
! Clean up after RFLU and TURB.
!
! Revision 1.5  2002/07/25 14:46:49  haselbac
! Added optimal LES flux routine
!
! Revision 1.4  2002/05/04 16:20:24  haselbac
! Added RFLU statements
!
! Revision 1.3  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.2  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.1  2002/01/23 03:51:24  jblazek
! Added low-level time-stepping routines.
!
! ******************************************************************************

