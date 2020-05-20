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
! Purpose: compute interaction source for drag forces on Lagrangian particles.
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
! $Id: INRT_CalcDrag.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_CalcDrag( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartLag,    ONLY : t_plag
  USE ModError
  USE ModParameters
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

  INTEGER :: icg, dragModel, iCont, nCont, nPcls
  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass

  REAL(RFREAL) :: C1,C2,C3,C4,CdFinal,CdSTD,CdM1,CdM2,CdMcr,DeltaCd,diamL,f1M, &
                  f1M1,f1M2,f2M,f2M1,f2M2,f3M,f3M1,f3M2,f4M,factor,gamma,lre, &
                  Mach,Mach1,Mach2,massL,mixtVolR,pi,psiL,relVelMagL,reyL, &
                  tauLR,vFrac,vFracCorr
  REAL(RFREAL),          DIMENSION(3)   :: relVel, accelL
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv, pDv, pTv

  TYPE(t_plag)  , POINTER :: pPlag
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_CalcDrag.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_CalcDrag',__FILE__ )

#ifdef PLAG
! begin -----------------------------------------------------------------------

! Check if there are any particles

  nPcls = 0

  IF (global%plagUsed) nPcls = region%plag%nPcls

  IF (nPcls < 1) GO TO 999

! Get dimensions --------------------------------------------------------------

  pi = global%pi

  nCont     = region%plagInput%nCont
  dragModel = region%inrtInput%inrts(INRT_TYPE_DRAG)%switches( &
    INRT_SWI_DRAG_MODEL)

! Set pointers ----------------------------------------------------------------

  pPlag     => region%plag

  pCv       => pPlag%cv
  pDv       => pPlag%dv
  pTv       => pPlag%tv

  pCvPlagMass => pPlag%cvPlagMass

! Loop over all cells and initialize force total ------------------------------
  DO icg = 1,region%grid%nCellsTot
    pPlag%forceTotal(XCOORD,icg) = 0.0_RFREAL
    pPlag%forceTotal(YCOORD,icg) = 0.0_RFREAL
    pPlag%forceTotal(ZCOORD,icg) = 0.0_RFREAL
  END DO ! icg

! Loop over all the particles -------------------------------------------------

  DO iPcls = 1,nPcls
    icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcls)

    vFrac     = pPlag%vFracL(1,iPcls)
    vFracCorr = (1.0_RFREAL + 2.0_RFREAL*vFrac)/(1.0_RFREAL-vFrac)**3.0_RFREAL

    diamL = pDv(DV_PLAG_DIAM,iPcls)

    massL = SUM( pCv(pCvPlagMass(:),iPcls) )

    tauLR = 3.0_RFREAL*pi*pTv(TV_PLAG_MUELMIXT,iPcls)* &
            diamL/massL

    relVel(1)   = pDv(DV_PLAG_UVELMIXT,iPcls)-pDv(DV_PLAG_UVEL,iPcls)
    relVel(2)   = pDv(DV_PLAG_VVELMIXT,iPcls)-pDv(DV_PLAG_VVEL,iPcls)
    relVel(3)   = pDv(DV_PLAG_WVELMIXT,iPcls)-pDv(DV_PLAG_WVEL,iPcls)

    relVelMagL  = SQRT( relVel(1)*relVel(1)+ &
                        relVel(2)*relVel(2)+ &
                        relVel(3)*relVel(3)  )

    reyL = diamL * relVelMagL * pDv(DV_PLAG_DENSMIXT,iPcls) / &
                                (pTv(TV_PLAG_MUELMIXT,iPcls)*(1.0_RFREAL-vFrac))

! TEMPORARY: Need proper method to interpolate mixture gamma at particle location. 
        gamma = 1.667_RFREAL
        gamma = 1.4_RFREAL   ! For shock-particle-cloud, rectangle to triangle.
! END TEMPORARY

    Mach  = relVelMagL/SQRT(gamma*pDv(DV_PLAG_PRESMIXT,iPcls) &
            *(1.0_RFREAL-vFrac)/pDv(DV_PLAG_DENSMIXT,iPcls))


    SELECT CASE (dragModel)

      CASE (INRT_DRAG_MODEL_STOKES)

        psiL = 1.0_RFREAL

      CASE (INRT_DRAG_MODEL_SN)

        psiL = 1.0_RFREAL + 0.15_RFREAL* (reyL**0.687_RFREAL)

      CASE (INRT_DRAG_MODEL_SMRFLD)

        psiL = 112.0_RFREAL/24.0_RFREAL *reyL**(+0.02_RFREAL)

      CASE (INRT_DRAG_MODEL_PARMAR)

        lre = LOG(reyL)

        IF ( reyL < 1.0E-14_RFREAL ) THEN
          psiL = 1.0_RFREAL
        ELSE
          IF (Mach <= 0.6_RFREAL) THEN
            CdSTD = (24.0_RFREAL/reyL)*(1.0_RFREAL+0.15_RFREAL*reyL**0.687_RFREAL) &
                    +0.42_RFREAL/(1.0_RFREAL+42500.0_RFREAL/reyL**1.16_RFREAL)
            CdMcr = (24.0_RFREAL/reyL)*(1.0_RFREAL+0.15_RFREAL*reyL**0.684_RFREAL) &
                    +0.513_RFREAL/(1.0_RFREAL+483.0_RFREAL/reyL**.669_RFREAL)
            factor  = Mach/0.6_RFREAL
            CdFinal = CdSTD + factor*(CdMcr-CdSTD)
          ELSE IF (Mach <= 1.0_RFREAL) THEN
            CdMcr = (24.0_RFREAL/reyL)*(1.0_RFREAL+0.15_RFREAL*reyL**0.684_RFREAL) &
                    +0.513_RFREAL/(1.0_RFREAL+483.0_RFREAL/reyL**.669_RFREAL)
            CdM1  = (24.0_RFREAL/reyL)*(1.0_RFREAL+0.118_RFREAL*reyL**0.813_RFREAL) &
                    +0.69_RFREAL/(1.0_RFREAL+3550.0_RFREAL/reyL**.793_RFREAL)

            C1 =  6.48_RFREAL
            C2 =  9.28_RFREAL
            C3 = 12.21_RFREAL

            f1M = -1.884_RFREAL +8.422_RFREAL*Mach -13.70_RFREAL*Mach*Mach +8.162_RFREAL*Mach*Mach*Mach
            f2M = -2.228_RFREAL +10.35_RFREAL*Mach -16.96_RFREAL*Mach*Mach +9.840_RFREAL*Mach*Mach*Mach
            f3M =  4.362_RFREAL -16.91_RFREAL*Mach +19.84_RFREAL*Mach*Mach -6.296_RFREAL*Mach*Mach*Mach

            factor = f1M*(lre-C2)*(lre-C3)/((C1-C2)*(C1-C3)) &
                    +f2M*(lre-C1)*(lre-C3)/((C2-C1)*(C2-C3)) &
                    +f3M*(lre-C1)*(lre-C2)/((C3-C1)*(C3-C2)) 
            CdFinal = CdMcr + factor*(CdM1-CdMcr)
          ELSE IF (Mach < 1.75_RFREAL) THEN
            CdM1 = (24.0_RFREAL/reyL)*(1.0_RFREAL+0.118_RFREAL*reyL**0.813_RFREAL) &
                   +0.69_RFREAL/(1.0_RFREAL+3550.0_RFREAL/reyL**.793_RFREAL)
            CdM2 = (24.0_RFREAL/reyL)*(1.0_RFREAL+0.107_RFREAL*reyL**0.867_RFREAL) &
                   +0.646_RFREAL/(1.0_RFREAL+861.0_RFREAL/reyL**.634_RFREAL)

            C1 =  6.48_RFREAL
            C2 =  8.93_RFREAL
            C3 = 12.21_RFREAL

            f1M = -2.963_RFREAL +4.392_RFREAL*Mach -1.169_RFREAL*Mach*Mach -0.027_RFREAL*Mach*Mach*Mach &
                  -0.233_RFREAL*exp((1.0_RFREAL-Mach)/0.011_RFREAL)
            f2M = -6.617_RFREAL +12.11_RFREAL*Mach -6.501_RFREAL*Mach*Mach +1.182_RFREAL*Mach*Mach*Mach &
                  -0.174_RFREAL*exp((1.0_RFREAL-Mach)/0.010_RFREAL)
            f3M = -5.866_RFREAL +11.57_RFREAL*Mach -6.665_RFREAL*Mach*Mach +1.312_RFREAL*Mach*Mach*Mach &
                  -0.350_RFREAL*exp((1.0_RFREAL-Mach)/0.012_RFREAL)

            factor = f1M*(lre-C2)*(lre-C3)/((C1-C2)*(C1-C3)) &
                    +f2M*(lre-C1)*(lre-C3)/((C2-C1)*(C2-C3)) &
                    +f3M*(lre-C1)*(lre-C2)/((C3-C1)*(C3-C2)) 
            CdFinal = CdM1 + factor*(CdM2-CdM1)
          ELSE
            CdFinal = (24.0_RFREAL/reyL)*(1.0_RFREAL+0.107_RFREAL*reyL**0.867_RFREAL) &
                      +0.646_RFREAL/(1.0_RFREAL+861.0_RFREAL/reyL**.634_RFREAL)
          END IF ! Mach

          psiL = reyL*CdFinal/24.0_RFREAL 
        END IF ! reyL

      CASE DEFAULT
        CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

    END SELECT ! dragModel

    accelL(1) = psiL*relVel(1)*tauLR
    accelL(2) = psiL*relVel(2)*tauLR
    accelL(3) = psiL*relVel(3)*tauLR

! - Correction due to volume fraction
    accelL(1) = accelL(1)*vFracCorr
    accelL(2) = accelL(2)*vFracCorr
    accelL(3) = accelL(3)*vFracCorr

! - fill the effect of the particles (L) on the gas (G).
! - the source terms will be added to the particle Rhs.

    pPlag%inrtSources(INRT_DRAG_L_XMOM_G,iPcls) = -massL*accelL(1)
    pPlag%inrtSources(INRT_DRAG_L_YMOM_G,iPcls) = -massL*accelL(2)
    pPlag%inrtSources(INRT_DRAG_L_ZMOM_G,iPcls) = -massL*accelL(3)
! HARDCODE: Manoj
!IF ( region%iRegionGlobal == 8 .AND. icg==12867 ) THEN
!   WRITE(*,'(A,I4,1X,5(E12.6,2X))') "QSDrag: ",iPcls,reyL,Mach,CdFinal,massL*accelL(1),massL*accelL(2)
!END IF
!IF ( iPcls == 1) THEN
!   WRITE(101,'(101(E12.6,2X))') global%currentTime,region%plag%dv(DV_PLAG_UVEL,iPcls),region%plag%dv(DV_PLAG_UVELMIXT,iPcls),CdFinal,reyL,Mach,massL*accelL(1),pTv(TV_PLAG_MUELMIXT,1),pDv(DV_PLAG_DENSMIXT,1),diamL
!END IF ! iPcls
!   Mach,reyL,CdFinal,diamL,pDv(DV_PLAG_DENSMIXT,1),pTv(TV_PLAG_MUELMIXT,1),region%mixt%tv(TV_MIXT_MUEL,icg)
! END HARDCODE           

! - Update force total --------------------------------------------------------

    pPlag%forceTotal(XCOORD,icg) = pPlag%forceTotal(XCOORD,icg) + massL*accelL(1)
    pPlag%forceTotal(YCOORD,icg) = pPlag%forceTotal(YCOORD,icg) + massL*accelL(2)
    pPlag%forceTotal(ZCOORD,icg) = pPlag%forceTotal(ZCOORD,icg) + massL*accelL(3)

  END DO ! iPcls

! finalize --------------------------------------------------------------------

999  CONTINUE
#endif
  CALL DeregisterFunction( global )

END SUBROUTINE INRT_CalcDrag

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_CalcDrag.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2009/07/09 20:44:19  mparmar
! Added Parmar drag law
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
! Revision 1.3  2007/03/08 14:59:36  fnajjar
! Fixed bug in Sommerfeld drag law for psiL
!
! Revision 1.2  2007/03/07 22:16:47  fnajjar
! Added Sommerfeld drag law
!
! Revision 1.1  2004/12/01 21:56:14  fnajjar
! Initial revision after changing case
!
! Revision 1.6  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.5  2004/01/31 03:59:22  haselbac
! Initial integration for Rocflu and Rocpart
!
! Revision 1.4  2003/04/03 21:10:17  jferry
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

