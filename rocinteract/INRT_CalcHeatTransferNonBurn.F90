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
! Purpose: compute interaction source for thermal forces on Lagrangian particles
!
! Description: none.
!
! Input: region  = current region.
!
! Output: region%levels(iLev)%plag%inrtSources
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_CalcHeatTransferNonBurn.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_CalcHeatTransferNonBurn( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartLag,    ONLY : t_plag
  USE ModError
  USE INRT_ModParameters

#ifdef PLAG
  USE PLAG_ModParameters
#endif
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region

! ... loop variables
  INTEGER :: iPcls

! ... local variables
  CHARACTER(CHRLEN)  :: RCSIdentString

  INTEGER :: thermalModel, iCont, nCont, nPcls
  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv

  REAL(RFREAL) :: accelL, diamL, heatCapL, mixtVolR, oneThird, pi, &
                  prLam, psiL, relVelMagL, relTemp, reyL, tauLR

  REAL(RFREAL),          DIMENSION(3)   :: relVel
  REAL(RFREAL), POINTER, DIMENSION(:)   :: pSpcHeat
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv, pDv, pTv

  TYPE(t_plag)  , POINTER :: pPlag
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: INRT_CalcHeatTransferNonBurn.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_CalcHeatTransferNonBurn',__FILE__ )

#ifdef PLAG
! begin -----------------------------------------------------------------------

! Check if there are any particles

  nPcls = 0

  IF (global%plagUsed) nPcls = region%plag%nPcls

  IF (nPcls < 1) GO TO 999

! Get dimensions --------------------------------------------------------------

  oneThird = 1.0_RFREAL/3.0_RFREAL
  pi    = global%pi
  prLam = global%prLam

  nCont        = region%plagInput%nCont
  thermalModel = region%inrtInput%inrts(INRT_TYPE_HTRANSNB)%switches( &
    INRT_SWI_HTRANSNB_MODEL)

! Set pointers ----------------------------------------------------------------

  pPlag     => region%plag
  pAiv      => pPlag%aiv
  pCv       => pPlag%cv
  pDv       => pPlag%dv
  pTv       => pPlag%tv

  pCvPlagMass => pPlag%cvPlagMass
  pSpcHeat    => region%plagInput%spht

! Loop over all the particles -------------------------------------------------

  DO iPcls = 1,nPcls

    IF (pAiv(AIV_PLAG_BURNSTAT,iPcls) /= INRT_BURNSTAT_OFF) CYCLE

    diamL  = pDv(DV_PLAG_DIAM,iPcls)

    heatCapL = SUM( pCv(pCvPlagMass(:),iPcls) * pSpcHeat(:) )

! - Compute thermal timescale for nonburning particles ------------------------

    tauLR = 2.0_RFREAL*pi*pTv(TV_PLAG_TCOLMIXT,iPcls) * diamL/heatCapL

    relTemp = pDv(DV_PLAG_TEMPMIXT,iPcls)-pDv(DV_PLAG_TEMP,iPcls)

    SELECT CASE (thermalModel)

      CASE (INRT_HTRANSNB_MODEL_STOKES)

        psiL = 1.0_RFREAL

      CASE (INRT_HTRANSNB_MODEL_RM)

        relVel(1)   = pDv(DV_PLAG_UVELMIXT,iPcls)-pDv(DV_PLAG_UVEL,iPcls)
        relVel(2)   = pDv(DV_PLAG_VVELMIXT,iPcls)-pDv(DV_PLAG_VVEL,iPcls)
        relVel(3)   = pDv(DV_PLAG_WVELMIXT,iPcls)-pDv(DV_PLAG_WVEL,iPcls)

        relVelMagL  = SQRT( relVel(1)*relVel(1)+ &
                            relVel(2)*relVel(2)+ &
                            relVel(3)*relVel(3)  )

        reyL = diamL * relVelMagL * pDv(DV_PLAG_DENSMIXT,iPcls) / &
                                    pTv(TV_PLAG_MUELMIXT,iPcls)

        psiL = 1.0_RFREAL + 0.3_RFREAL* (reyL**0.5_RFREAL) *(prLam**oneThird)

      CASE DEFAULT
        CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

    END SELECT ! thermalModel

    accelL = psiL*relTemp*tauLR

! - fill the effect of the particles (L) on the gas (G).
! - the source terms will be added to the particle Rhs.

    pPlag%inrtSources(INRT_HTRANSNB_L_ENER_G,iPcls) = - heatCapL*accelL
  END DO ! iPcls

! finalize --------------------------------------------------------------------

999  CONTINUE
#endif
  CALL DeregisterFunction( global )

END SUBROUTINE INRT_CalcHeatTransferNonBurn

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_CalcHeatTransferNonBurn.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:49  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:01  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:50:11  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:14  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 21:56:15  fnajjar
! Initial revision after changing case
!
! Revision 1.7  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.6  2004/03/02 21:47:29  jferry
! Added After Update interactions
!
! Revision 1.5  2004/01/31 03:59:22  haselbac
! Initial integration for Rocflu and Rocpart
!
! Revision 1.4  2003/04/03 21:10:18  jferry
! implemented additional safety checks for rocinteract
!
! Revision 1.3  2003/04/02 22:32:03  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.2  2003/03/11 16:09:39  jferry
! Added comments
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!******************************************************************************

