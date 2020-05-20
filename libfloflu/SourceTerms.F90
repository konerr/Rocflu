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
! Purpose: Add source terms to the residual.
!
! Description: None.
!
! Input: 
!   region      data of current region.
!
! Output: None.
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: SourceTerms.F90,v 1.6 2016/05/06 00:36:42 rahul Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SourceTerms(region)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region,t_grid,t_mixt,t_mixt_input
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region) :: region
  TYPE(t_grid)   :: grid
  TYPE(t_mixt)   :: mixt
  TYPE(t_mixt_input) :: mixtInput
! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: ic,errorflag
  REAL(RFREAL) :: rhoVol
                  
  REAL(RFREAL), POINTER :: cv(:,:),rhs(:,:),vol(:)
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'SourceTerms',__FILE__ )

! ******************************************************************************
! Get dimensions and pointers
! ******************************************************************************

  cv  => region%mixt%cv
  rhs => region%mixt%rhs
  vol => region%grid%vol 

! ******************************************************************************
! Source term due to acceleration
! ******************************************************************************

!  IF ( global%accelOn .EQV. .TRUE. ) THEN
    DO ic = 1,region%grid%nCells
      rhoVol = vol(ic)*cv(CV_MIXT_DENS,ic)
      rhs(CV_MIXT_XMOM,ic) = rhs(CV_MIXT_XMOM,ic) - global%accelX*rhoVol
      rhs(CV_MIXT_YMOM,ic) = rhs(CV_MIXT_YMOM,ic) - global%accelY*rhoVol
      rhs(CV_MIXT_ZMOM,ic) = rhs(CV_MIXT_ZMOM,ic) - global%accelZ*rhoVol
      rhs(CV_MIXT_ENER,ic) = rhs(CV_MIXT_ENER,ic) - &
                             vol(ic)*(global%accelX*cv(CV_MIXT_XMOM,ic)+ &
                                      global%accelY*cv(CV_MIXT_YMOM,ic)+ &
                                      global%accelZ*cv(CV_MIXT_ZMOM,ic))
    END DO ! ic
!  END IF ! global%accelOn

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SourceTerms

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SourceTerms.F90,v $
! Revision 1.6  2016/05/06 00:36:42  rahul
! 1. Moved non-conservative terms from E-L AUSM+up out of this subroutine to
!    SourceTermsPint.F90.
! 2. Commented out accelOn flag for gravity. Added this check in
!    SourceTermsMP.F90.
!
! Revision 1.5  2016/02/08 22:17:51  rahul
! Fixed a sign error in rhs of energy term.
!
! Revision 1.4  2016/02/04 23:12:15  rahul
! Fixed a minor bug in ALLOCATE statement.
!
! Revision 1.3  2016/01/31 04:53:54  rahul
! Fixed a bug from earlier check-in.
!
! Revision 1.2  2015/12/18 23:08:57  rahul
! Added computation of non-conservative source terms in multiphase AUSM+up
! formulation.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.3  2005/11/10 01:58:37  haselbac
! Added Rocflu support for accel terms, clean-up
!
! Revision 1.2  2004/12/04 06:13:31  wasistho
! cvOld to cv in updating rhs by DualTST source terms
!
! Revision 1.1  2004/12/01 16:51:24  haselbac
! Initial revision after changing case
!
! Revision 1.9  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/08/28 20:05:38  jblazek
! Added acceleration terms.
!
! Revision 1.5  2003/07/03 21:48:44  jblazek
! Implemented dual-time stepping.
!
! Revision 1.4  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.3  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/01/23 03:51:24  jblazek
! Added low-level time-stepping routines.
!
! *****************************************************************************

