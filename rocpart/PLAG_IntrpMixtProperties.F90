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
! Purpose: Interpolate mixture properties onto particle locations
!
! Description: none.
!
! Input: region = current region.
!
! Output: region%levels%plag%dv 
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_IntrpMixtProperties.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_IntrpMixtProperties( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModMixture, ONLY    : t_mixt
  USE ModError  
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region
  
! ... loop variables
  INTEGER :: iPcls

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iCell, nPcls
  INTEGER, POINTER, DIMENSION(:,:)    :: pAiv

  REAL(RFREAL), POINTER, DIMENSION(:,:)    :: pCv, pDv, pTv
  REAL(RFREAL), POINTER, DIMENSION(:,:)    :: pCvMixt, pDvMixt, pTvMixt
  
  TYPE(t_global), POINTER :: global
  TYPE(t_plag),   POINTER :: pPlag 
  TYPE(t_mixt),   POINTER :: pMixt 
 
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_IntrpMixtProperties.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_IntrpMixtProperties',__FILE__ )

! Get dimensions --------------------------------------------------------------

  nPcls = region%plag%nPcls

! Set pointers ----------------------------------------------------------------

  pMixt => region%mixt
  pPlag => region%plag

! Set pointers for mixture ----------------------------------------------------

  pCvMixt => pMixt%cv
  pDvMixt => pMixt%dv
  pTvMixt => pMixt%tv

! Set pointers for discrete particles -----------------------------------------

  pAiv    => pPlag%aiv
  pCv     => pPlag%cv  
  pDv     => pPlag%dv
  pTv     => pPlag%tv

! Check that have primitive state vector --------------------------------------

  IF ( pMixt%cvState /= CV_MIXT_STATE_DUVWP ) THEN 
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! region%mixt%cvState
     
! Assume a piece-wise constant distribution in a cell -------------------------
!  Note: Distribution has been properly corrected for RFLU
!        with the routine PLAG_RFLU_CorrectMixtProperties.
!        This routine is called after the cell gradients are computed
    
  DO iPcls = 1, nPcls

    iCell = pAiv(AIV_PLAG_ICELLS,iPcls)
    pDv(DV_PLAG_UVELMIXT,iPcls) = pCvMixt(CV_MIXT_XVEL,iCell)
    pDv(DV_PLAG_VVELMIXT,iPcls) = pCvMixt(CV_MIXT_YVEL,iCell)  
    pDv(DV_PLAG_WVELMIXT,iPcls) = pCvMixt(CV_MIXT_ZVEL,iCell)
    pDv(DV_PLAG_PRESMIXT,iPcls) = pDvMixt(DV_MIXT_PRES,iCell)
    pDv(DV_PLAG_TEMPMIXT,iPcls) = pDvMixt(DV_MIXT_TEMP,iCell)
    pDv(DV_PLAG_DENSMIXT,iPcls) = pCvMixt(CV_MIXT_DENS,iCell)

    pTv(TV_PLAG_MUELMIXT,iPcls) = pTvMixt(TV_MIXT_MUEL,iCell)
    pTv(TV_PLAG_TCOLMIXT,iPcls) = pTvMixt(TV_MIXT_TCOL,iCell)
  END DO ! iPcls
   
! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_IntrpMixtProperties

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_IntrpMixtProperties.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:51  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:03  mtcampbe
! Added Illinois Open Source License/Copyright
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
! Revision 1.2  2005/11/30 22:18:57  fnajjar
! Changed bcos of corr routine, now only do piecewise const here
!
! Revision 1.1  2004/12/01 20:57:51  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2004/02/26 21:02:15  haselbac
! Fixed RFLU bug, removed superfluous statements
!
! Revision 1.4  2004/02/06 21:18:20  fnajjar
! Initial Integration of Rocpart with Rocflu
!
! Revision 1.3  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.2  2003/01/16 20:18:53  f-najjar
! Removed iRegionGlobal
!
! Revision 1.1  2002/10/25 14:16:32  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************

