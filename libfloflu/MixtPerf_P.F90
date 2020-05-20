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
! Purpose: Collect relations for static and total pressure for perfect gases.
!
! Description: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: MixtPerf_P.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

FUNCTION MixtPerf_P_DEoGVm2(D,Eo,G,Vm2)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: D,Eo,G,Vm2
  REAL(RFREAL) :: MixtPerf_P_DEoGVm2
   
  MixtPerf_P_DEoGVm2 = (G - 1.0_RFREAL)*D*(Eo - 0.5_RFREAL*Vm2)

END FUNCTION MixtPerf_P_DEoGVm2

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_P_DRT(D,R,T)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: D,R,T
  REAL(RFREAL) :: MixtPerf_P_DRT
   
  MixtPerf_P_DRT = D*R*T

END FUNCTION MixtPerf_P_DRT

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_P_GMaPo(G,Ma,Po)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: G,Ma,Po
  REAL(RFREAL) :: MixtPerf_P_GMaPo
   
  MixtPerf_P_GMaPo = Po/ & 
    ((1.0_RFREAL + 0.5_RFREAL*(G - 1.0_RFREAL)*Ma*Ma)**(G/(G - 1.0_RFREAL)))

END FUNCTION MixtPerf_P_GMaPo

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_P_DDoGPo(D,Do,G,Po)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: D,Do,G,Po
  REAL(RFREAL) :: MixtPerf_P_DDoGPo
   
  MixtPerf_P_DDoGPo = Po*(D/Do)**G

END FUNCTION MixtPerf_P_DDoGPo

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_P_GPoTTo(G,Po,T,To)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: G,Po,T,To
  REAL(RFREAL) :: MixtPerf_P_GPoTTo
   
  MixtPerf_P_GPoTTo = Po*(T/To)**(G/(G - 1.0_RFREAL))

END FUNCTION MixtPerf_P_GPoTTo

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_Po_GPTTo(G,P,T,To)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: G,P,T,To
  REAL(RFREAL) :: MixtPerf_Po_GPTTo
   
  MixtPerf_Po_GPTTo = P*(To/T)**(G/(G - 1.0_RFREAL))

END FUNCTION MixtPerf_Po_GPTTo

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_Po_CGPUVW(C,G,P,U,V,W)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: C,G,P,U,V,W
  REAL(RFREAL) :: MixtPerf_Po_CGPUVW
   
  MixtPerf_Po_CGPUVW = & 
    P*(1.0_RFREAL + 0.5_RFREAL*(G - 1.0_RFREAL)*(U*U+V*V+W*W)/(C*C))**(G/(G - 1.0_RFREAL))

END FUNCTION MixtPerf_Po_CGPUVW

!******************************************************************************
!
! RCS Revision history:
!
! $Log: MixtPerf_P.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 16:49:18  haselbac
! Initial revision after changing case
!
! Revision 1.7  2003/09/16 01:24:08  haselbac
! Added function MixtPerf_P_DDoGPo
!
! Revision 1.6  2003/06/04 21:54:21  haselbac
! Added MixtPerf_P_DRT
!
! Revision 1.5  2002/07/05 23:20:46  jblazek
! Corrected bug in perfgasDependentVars.F90; did some cosmetics.
!
! Revision 1.4  2002/06/10 21:15:53  haselbac
! Added MixtPerf_P_GMaPo function
!
! Revision 1.3  2002/06/05 18:31:38  haselbac
! Added function mixtPerf_P_DEoGVm2
!
! Revision 1.2  2002/05/28 13:45:09  haselbac
! Added new functions
!
! Revision 1.1  2002/05/04 16:16:52  haselbac
! Initial revision
!
!******************************************************************************

