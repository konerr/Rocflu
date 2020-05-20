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
! Purpose: update PLAG derived variables (dv) arrays.
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
! $Id: PLAG_CalcDerivedVariables.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CalcDerivedVariables( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region

! ... loop variables
  INTEGER :: iCont, iPcls

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: nCont, nPcls
  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv

  REAL(RFREAL) :: heatCapSum, massSum, massSumR, &
                  oneThird, pi, piR, sphtL, voluSum
  REAL(RFREAL), POINTER, DIMENSION(:)   :: pDens, pMixtVol, pSpcHeat
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv, pDv
  
  TYPE(t_plag)  , POINTER :: pPlag
  TYPE(t_global), POINTER :: global  
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CalcDerivedVariables.F90,v $ $Revision: 1.1.1.1 $'
  
  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_CalcDerivedVariables',__FILE__ )

! Get dimensions --------------------------------------------------------------

  oneThird = 1.0_RFREAL/3.0_RFREAL  
  pi       = global%pi
  piR      = 1.0_RFREAL/pi

  nPcls = region%plag%nPcls
  nCont = region%plagInput%nCont
  
! Set pointers ----------------------------------------------------------------

  pMixtVol  => region%grid%vol
  pPlag     => region%plag

  pAiv   => pPlag%aiv  
  pCv    => pPlag%cv
  pDv    => pPlag%dv
  
  pCvPlagMass => pPlag%cvPlagMass
  
  pDens     => region%plagInput%dens
  pSpcHeat  => region%plagInput%spht
  
! Calculate derived variables for Lagrangian field ----------------------------

  DO iPcls = 1, nPcls

!- Compute particle mass ------------------------------------------------------

    massSum  = SUM( pCv(pPlag%cvPlagMass(:),iPcls) )    
    massSumR = 1.0_RFREAL/massSum
    
!- Compute particle volume ----------------------------------------------------

    voluSum = 0.0_RFREAL 

    DO iCont = 1, nCont
      voluSum = voluSum + pCv(pCvPlagMass(iCont),iPcls)/pDens(iCont)
    END DO ! iCont   
      
!- Compute particle heat capacity ---------------------------------------------
     
    heatCapSum  = SUM(pCv(pCvPlagMass(:),iPcls) * pSpcHeat(:) )

!- Compute particle specific heat ---------------------------------------------

    sphtL = heatCapSum * massSumR

!- Extract derived variables -------------------------------------------------- 
          
    pDv(DV_PLAG_UVEL,iPcls) = pCv(CV_PLAG_XMOM,iPcls) * massSumR
    pDv(DV_PLAG_VVEL,iPcls) = pCv(CV_PLAG_YMOM,iPcls) * massSumR
    pDv(DV_PLAG_WVEL,iPcls) = pCv(CV_PLAG_ZMOM,iPcls) * massSumR
    
    pDv(DV_PLAG_DIAM,iPcls) = (6.0_RFREAL*piR*voluSum)**oneThird
    
    pDv(DV_PLAG_TEMP,iPcls) = (           pCv(CV_PLAG_ENER,iPcls)*massSumR &
                            - 0.5_RFREAL*(pDv(DV_PLAG_UVEL,iPcls)**2       &
                                         +pDv(DV_PLAG_VVEL,iPcls)**2       &
                                         +pDv(DV_PLAG_WVEL,iPcls)**2 ) )/  &
                               sphtL
  END DO  ! iPcls 

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CalcDerivedVariables

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CalcDerivedVariables.F90,v $
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
! Revision 1.1  2004/12/01 20:57:01  fnajjar
! Initial revision after changing case
!
! Revision 1.6  2004/03/26 21:27:03  fnajjar
! Set particle volume in dv array
!
! Revision 1.5  2004/03/02 21:50:03  jferry
! Changed name of DV_PLAG_HTCP to DV_PLAG_SPHT
!
! Revision 1.4  2004/02/06 21:18:20  fnajjar
! Initial Integration of Rocpart with Rocflu
!
! Revision 1.3  2003/01/16 20:15:11  f-najjar
! Removed iRegionGlobal
!
! Revision 1.2  2002/12/03 19:55:37  f-najjar
! Fixed subroutine string name in RegisterFunction
!
! Revision 1.1  2002/10/25 14:14:44  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************

