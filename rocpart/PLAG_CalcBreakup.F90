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
! Purpose: invoke breakup model case for Lagrangian particles.
!
! Description: none.
!
! Input: region  = current region.
!
! Output: region%levels(iLev)%plag%cv
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_CalcBreakup.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CalcBreakup( region, iReg )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartLag,    ONLY : t_plag, t_plag_input
  USE ModError
  USE PLAG_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region

  INTEGER, INTENT(IN) :: iReg

! ... loop variables
  INTEGER :: iCont, iPcls

! ... local variables
  CHARACTER(CHRLEN)  :: RCSIdentString

  INTEGER :: breakupModel, breakupWebSwi, nCont, nPcls
  INTEGER, POINTER, DIMENSION(:) :: pCvPlagMass

  REAL(RFREAL) :: breakupFac, breakupFacR, densG, diamL, diamLSplit, &
                  oneThird, pi, presG, relVelMagL, surfTensL,        &
                  surfTensSum, voluSum, voluSumR, weberL, weberCrit

  REAL(RFREAL),          DIMENSION(3)   :: relVel
  REAL(RFREAL), POINTER, DIMENSION(:)   :: pDens,pSurfTens
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv,pCv,pDv

  TYPE(t_plag),   POINTER :: pPlag
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CalcBreakup.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction( global,'PLAG_CalcBreakup',__FILE__ )

! begin =======================================================================

! Check if there are any particles

  nPcls = 0
  IF (global%plagUsed) nPcls = region%plag%nPcls
  IF (nPcls < 1) GO TO 999

! Set pointers ----------------------------------------------------------------

  pPlag => region%plag

  pArv  => pPlag%arv
  pCv   => pPlag%cv
  pDv   => pPlag%dv

  pCvPlagMass => pPlag%cvPlagMass
  pSurfTens   => region%plagInput%surftens
  pDens       => region%plagInput%dens

! Get dimensions --------------------------------------------------------------

  pi = global%pi
  oneThird = 1.0_RFREAL/3.0_RFREAL

  nCont  = region%plagInput%nCont

  breakupModel  = region%plagInput%breakupModel
  breakupWebSwi = region%plagInput%breakupWebSwi
  breakupFac    = region%plagInput%breakupFac
  breakupFacR   = 1.0_RFREAL/breakupFac

! Set appropriate coefficients pertinent to breakup Model ---------------------

  SELECT CASE (breakupModel)

    CASE (PLAG_BREAKUP_MODEL1)
      weberCrit = 10.0_RFREAL

    CASE DEFAULT
      CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

  END SELECT ! breakupModel

! Loop over all the particles -------------------------------------------------

  DO iPcls = 1,nPcls

! - Extract gas properties ----------------------------------------------------

    densG = pDv(DV_PLAG_DENSMIXT,iPcls)

!- Compute particle volume for each constituent -------------------------------

    voluSum = 0.0_RFREAL

    DO iCont = 1, nCont
      voluSum = voluSum + pCv(pCvPlagMass(iCont),iPcls)/pDens(iCont)
    END DO ! iCont

    voluSumR = 1.0_RFREAL/voluSum

! - Compute Particle surface tension ------------------------------------------

    surfTensSum = 0.0_RFREAL

    DO iCont = 1, nCont
      surfTensSum = surfTensSum + pCv(pCvPlagMass(iCont),iPcls)/pDens(iCont) &
                  * pSurfTens(iCont)
    END DO ! iCont

    surfTensL =  surfTensSum * voluSumR

    diamL = pDv(DV_PLAG_DIAM,iPcls)

! - Compute relative velocities and its magnitude ----------------------------

    relVel(1) = pDv(DV_PLAG_UVELMIXT,iPcls)-pDv(DV_PLAG_UVEL,iPcls)
    relVel(2) = pDv(DV_PLAG_VVELMIXT,iPcls)-pDv(DV_PLAG_VVEL,iPcls)
    relVel(3) = pDv(DV_PLAG_WVELMIXT,iPcls)-pDv(DV_PLAG_WVEL,iPcls)

    relVelMagL  = relVel(1)*relVel(1) &
                + relVel(2)*relVel(2) &
                + relVel(3)*relVel(3)

! - Compute weber number ------------------------------------------------------

    weberL = densG * diamL * relVelMagL / surfTensL

! - Check if critical weber number is met -------------------------------------

    IF ( weberL >= weberCrit ) THEN

#ifdef PLAG_DEBUG
      WRITE(*,'(A,3X,I3,3X,I4,3X,1PE12.5)') &
      'PLAG_CalcBreakup-Critical We reached: iReg, iPcls, We = ',&
       iReg, iPcls, weberL
#endif

! -- Redefine breakupFactor based on critical Weber, if needed ----------------

      IF ( breakupWebSwi == PLAG_BREAKUP_WEBSWI1 ) THEN
        breakupFac  = ( densG * diamL * relVelMagL /( surfTensL *weberCrit ) ) **3
        breakupFacR = 1.0_RFREAL/breakupFac

#ifdef PLAG_DEBUG
      WRITE(*,'(A,3X,I3,3X,I4,3X,1PE12.5)') &
      'PLAG_CalcBreakup-Breakup Switch Active: iReg, iPcls, breakupFac = ',&
       iReg, iPcls, breakupFac
#endif

      END IF ! breakupSwi

! -- Update cv and arv values --------------------------------------------------

      DO iCont = 1, nCont
       pCv(pCvPlagMass(iCont),iPcls) = pCv(pCvPlagMass(iCont),iPcls)*breakupFacR
      END DO ! iCont

      pCv(CV_PLAG_XMOM,iPcls) = pCv(CV_PLAG_XMOM,iPcls) * breakupFacR
      pCv(CV_PLAG_YMOM,iPcls) = pCv(CV_PLAG_YMOM,iPcls) * breakupFacR
      pCv(CV_PLAG_ZMOM,iPcls) = pCv(CV_PLAG_ZMOM,iPcls) * breakupFacR
      pCv(CV_PLAG_ENER,iPcls) = pCv(CV_PLAG_ENER,iPcls) * breakupFacR
      pCv(CV_PLAG_ENERVAPOR,iPcls) = pCv(CV_PLAG_ENERVAPOR,iPcls) * breakupFacR

      pArv(ARV_PLAG_SPLOAD,iPcls) = pArv(ARV_PLAG_SPLOAD,iPcls)* breakupFac
    END IF ! weberL

  END DO ! iPcls

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CalcBreakup

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CalcBreakup.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.5  2008/12/06 08:43:50  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:02  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/04/26 20:45:34  fnajjar
! Removed dv-based calculations that are irrelevant and modified the needed ones
!
! Revision 1.2  2007/04/16 23:19:35  fnajjar
! Deleted all RFLO-pertinent calls and statements
!
! Revision 1.1  2007/04/09 18:50:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 20:57:00  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2004/03/25 21:16:43  jferry
! fixed Vapor Energy bug
!
! Revision 1.4  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.3  2004/02/06 21:18:20  fnajjar
! Initial Integration of Rocpart with Rocflu
!
! Revision 1.2  2003/09/15 20:26:54  fnajjar
! Corrected breakupFac and removed cubeRootFac
!
! Revision 1.1  2003/09/13 20:15:13  fnajjar
! Initialimport of breakup model
!
!******************************************************************************

