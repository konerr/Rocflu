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
! Purpose: Collect relations for static and total temperature for perfect 
!   gases.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: MixtPerf_T.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

FUNCTION MixtPerf_T_CGR(C,G,R)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: C,G,R
  REAL(RFREAL) :: MixtPerf_T_CGR
   
  MixtPerf_T_CGR = C*C/(G*R)

END FUNCTION MixtPerf_T_CGR

!------------------------------------------------------------------------------

FUNCTION MixtPerf_T_CpHoVm2(Cp,Ho,Vm2)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: Cp,Ho,Vm2
  REAL(RFREAL) :: MixtPerf_T_CpHoVm2
   
  MixtPerf_T_CpHoVm2 = (Ho-0.5_RFREAL*Vm2)/Cp

END FUNCTION MixtPerf_T_CpHoVm2

!------------------------------------------------------------------------------

FUNCTION MixtPerf_T_CvEoVm2(Cv,Eo,Vm2)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: Cv,Eo,Vm2
  REAL(RFREAL) :: MixtPerf_T_CvEoVm2
   
  MixtPerf_T_CvEoVm2 = (Eo-0.5_RFREAL*Vm2)/Cv

END FUNCTION MixtPerf_T_CvEoVm2

! ------------------------------------------------------------------------------

FUNCTION MixtPerf_T_DPR(D,P,R)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: D,P,R
  REAL(RFREAL) :: MixtPerf_T_DPR
   
  MixtPerf_T_DPR = P/(D*R)

END FUNCTION MixtPerf_T_DPR

! ------------------------------------------------------------------------------

FUNCTION MixtPerf_HM_T_DGMP(D,G,M,P)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: D,G,M,P
  REAL(RFREAL) :: MixtPerf_HM_T_DGMP
   
  MixtPerf_HM_T_DGMP = (G*M*M*P + 1.0_RFREAL)/D 

END FUNCTION MixtPerf_HM_T_DGMP

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_T_GMaTo(G,Ma,To)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: G,Ma,To
  REAL(RFREAL) :: MixtPerf_T_GMaTo
   
  MixtPerf_T_GMaTo = To/(1.0_RFREAL + 0.5_RFREAL*(G - 1.0_RFREAL)*Ma*Ma)

END FUNCTION MixtPerf_T_GMaTo

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_To_CpTUVW(Cp,T,U,V,W)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: Cp,T,U,V,W
  REAL(RFREAL) :: MixtPerf_To_CpTUVW
   
  MixtPerf_To_CpTUVW = T + 0.5_RFREAL*(U*U + V*V + W*W)/Cp

END FUNCTION MixtPerf_To_CpTUVW

!******************************************************************************
!
! RCS Revision history:
!
! $Log: MixtPerf_T.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:31  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:16:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/11/28 23:17:42  mparmar
! Added MixtPerf_HM_T_DGMP function
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/05/01 20:58:15  haselbac
! Added MixtPerf_T_CpHoVm2 function
!
! Revision 1.3  2006/03/26 20:21:17  haselbac
! Added fuction
!
! Revision 1.2  2005/07/14 21:39:56  haselbac
! Added MixtPerf_To_CpTUVW
!
! Revision 1.1  2004/12/01 16:49:30  haselbac
! Initial revision after changing case
!
! Revision 1.3  2002/06/10 21:16:27  haselbac
! Added MixtPerf_T_GMaTo function
!
! Revision 1.2  2002/06/05 18:31:14  haselbac
! Added function mixtPerf_T_DPR
!
! Revision 1.1  2002/05/04 16:16:52  haselbac
! Initial revision
!
! ******************************************************************************

