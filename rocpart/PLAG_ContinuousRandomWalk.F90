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
! Purpose: compute particle velocity fluctuation using CRW.
!
! Description: none.
!
! Input: region  = data of current region.
!
! Output: region%levels%plag%cv
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_ContinuousRandomWalk.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_ContinuousRandomWalk( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters

  USE ModRandom, ONLY: Rand1Uniform

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region

! ... loop variables
  INTEGER :: iPcls

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: nPcls
  REAL(RFREAL) :: delU,delV,delW,diameter,factor,heatCapSum,ke,keNew, &
                  massSum,massSumR, &
                  phi,rn,temp,TxScale,TyScale,TzScale,u, &
                  uNew,uOld,uPrime,uPrimeOld,uScale,v,vNew,vOld,vPrime, &
                  vPrimeOld,vScale,w,wNew,wOld,wPrime,wPrimeOld,wScale
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv,pCv,pDv,pDvOld
  
  TYPE(t_plag) ,  POINTER :: pPlag  
  TYPE(t_global), POINTER :: global 
   
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ContinuousRandomWalk.F90,v $'

  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_ContinuousRandomWalk',__FILE__ )

! Get dimensions --------------------------------------------------------------

  nPcls = region%plag%nPcls

! Set pointers ----------------------------------------------------------------

  pPlag  => region%plag
  pArv   => pPlag%arv
  pCv    => pPlag%cv
  pDv    => pPlag%dv 
  pDvOld => pPlag%dvOld

  factor = 1.0_RFREAL

! Calculate particle velocity fluctuation -------------------------------------

  DO iPcls = 1, nPcls
    massSum  = SUM( pCv(pPlag%cvPlagMass(:),iPcls) )    
    massSumR = 1.0_RFREAL/massSum
    diameter = pDv(DV_PLAG_DIAM,iPcls)

    u = pCv(CV_PLAG_XMOM,iPcls)*massSumR
    v = pCv(CV_PLAG_YMOM,iPcls)*massSumR
    w = pCv(CV_PLAG_ZMOM,iPcls)*massSumR

    uOld = pDvOld(DV_PLAG_UVEL,iPcls)
    vOld = pDvOld(DV_PLAG_VVEL,iPcls)
    wOld = pDvOld(DV_PLAG_WVEL,iPcls)

! TEMPORARY: Manoj: July 22 2012: Trying different definition for velocity scale
!    uScale = ABS(u-uOld)
!    vScale = ABS(v-vOld)
!    wScale = ABS(w-wOld)
    uScale = ABS(u)
    vScale = ABS(v)
    wScale = ABS(w)

    phi = pPlag%vFracL(1,iPcls)

    TxScale = 0.5_RFREAL*(11.0_RFREAL**0.5_RFREAL) &
              *(diameter/(uScale+1.0E-15_RFREAL))*(phi**(-2.0_RFREAL/3.0_RFREAL))
    TyScale = 0.5_RFREAL*(11.0_RFREAL**0.5_RFREAL) &
              *(diameter/(vScale+1.0E-15_RFREAL))*(phi**(-2.0_RFREAL/3.0_RFREAL))
    TzScale = 0.5_RFREAL*(11.0_RFREAL**0.5_RFREAL) &
              *(diameter/(wScale+1.0E-15_RFREAL))*(phi**(-2.0_RFREAL/3.0_RFREAL))

    factor = ((11.0_RFREAL**0.5_RFREAL)*(1-phi**0.125_RFREAL) + phi**0.125_RFREAL) &
             *(phi**(1.0_RFREAL/3.0_RFREAL))

    rn   = 2.0_RFREAL*Rand1Uniform(region%randData) - 1.0_RFREAL
    delU = factor*uScale*SIGN(1.0_RFREAL,rn)
    rn   = 2.0_RFREAL*Rand1Uniform(region%randData) - 1.0_RFREAL
    delV = factor*vScale*SIGN(1.0_RFREAL,rn)
    rn   = 2.0_RFREAL*Rand1Uniform(region%randData) - 1.0_RFREAL
    delW = factor*wScale*SIGN(1.0_RFREAL,rn)

    uPrimeOld = pArv(ARV_PLAG_UPRIME,iPcls)
    vPrimeOld = pArv(ARV_PLAG_VPRIME,iPcls)
    wPrimeOld = pArv(ARV_PLAG_WPRIME,iPcls)

    uPrime = EXP(-global%dtMin/TxScale)*uPrimeOld &
             + SQRT(1 - EXP(-2.0_RFREAL*global%dtMin/TxScale))*delU
    vPrime = EXP(-global%dtMin/TyScale)*vPrimeOld &
             + SQRT(1 - EXP(-2.0_RFREAL*global%dtMin/TyScale))*delV
    wPrime = EXP(-global%dtMin/TzScale)*wPrimeOld &
             + SQRT(1 - EXP(-2.0_RFREAL*global%dtMin/TzScale))*delW

    pArv(ARV_PLAG_UPRIME,iPcls) = uPrime
    pArv(ARV_PLAG_VPRIME,iPcls) = vPrime
    pArv(ARV_PLAG_WPRIME,iPcls) = wPrime

    uNew = u + uPrime - uPrimeOld
    vNew = v + vPrime - vPrimeOld
    wNew = w + wPrime - wPrimeOld

    ke    = 0.5_RFREAL*massSum*(u*u+v*v+w*w)
    keNew = 0.5_RFREAL*massSum*(uNew*uNew+vNew*vNew+wNew*wNew)

    pCv(CV_PLAG_XMOM,iPcls) = uNew*massSum
    pCv(CV_PLAG_YMOM,iPcls) = vNew*massSum
    pCv(CV_PLAG_ZMOM,iPcls) = wNew*massSum

! TEMPORARY: Manoj: 2012-07-27: Fixing error with -ve energy
    IF ( pCv(CV_PLAG_ENER,iPcls) <= ke ) THEN
      temp        = pDv(DV_PLAG_TEMP,iPcls)
      IF ( temp < 0.0_RFREAL ) THEN
        WRITE(*,*) "CRW: Temperature is -ve for particle #",iPcls,temp
        STOP
      END IF ! temp

      heatCapSum  = SUM(pCv(pPlag%cvPlagMass(:),iPcls) * region%plagInput%spht(:) )

      pCv(CV_PLAG_ENER,iPcls) = temp*heatCapSum + ke 

      WRITE(*,'(A,I8,A)') "  CRW:----  Fixed -ve energy in particle",iPcls,"  --------  "
    END IF ! ener
! END TEMPORARY

    pCv(CV_PLAG_ENER,iPcls) = pCv(CV_PLAG_ENER,iPcls) - ke + keNew
  END DO  ! iPcls

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

! DEBUG: Manoj: 2012-06-25    
!WRITE(*,*) "Stopping here in CRW...."
!STOP
! END DEBUG

END SUBROUTINE PLAG_ContinuousRandomWalk

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ContinuousRandomWalk.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
!
!******************************************************************************

