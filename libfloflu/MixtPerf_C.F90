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
! Purpose: Collect relations for static and total speed of sound for perfect
!   gases.
!
! Description: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: MixtPerf_C.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

FUNCTION MixtPerf_C_Co2GUVW(Co2,G,U,V,W)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: Co2,G,U,V,W
  REAL(RFREAL) :: MixtPerf_C_Co2GUVW
   
  MixtPerf_C_Co2GUVW = & 
    SQRT(Co2 - 0.5_RFREAL*(G - 1.0_RFREAL)*(U*U + V*V + W*W))

END FUNCTION MixtPerf_C_Co2GUVW

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_C_DGP(D,G,P)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: D,G,P
  REAL(RFREAL) :: MixtPerf_C_DGP
   
  MixtPerf_C_DGP = SQRT(G*P/D)

END FUNCTION MixtPerf_C_DGP

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_C_GHoVm2(G,Ho,Vm2)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: G,Ho,Vm2
  REAL(RFREAL) :: MixtPerf_C_GHoVm2
   
  MixtPerf_C_GHoVm2 = SQRT((G - 1.0_RFREAL)*(Ho - 0.5_RFREAL*Vm2))

END FUNCTION MixtPerf_C_GHoVm2

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_C_GRT(G,R,T)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: G,R,T
  REAL(RFREAL) :: MixtPerf_C_GRT
   
  MixtPerf_C_GRT = SQRT(G*R*T)

END FUNCTION MixtPerf_C_GRT

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_C2_GRT(G,R,T)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: G,R,T
  REAL(RFREAL) :: MixtPerf_C2_GRT
   
  MixtPerf_C2_GRT = G*R*T

END FUNCTION MixtPerf_C2_GRT

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_Co2_CGUVW(C,G,U,V,W)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: C,G,U,V,W
  REAL(RFREAL) :: MixtPerf_Co2_CGUVW
   
  MixtPerf_Co2_CGUVW = C*C + 0.5_RFREAL*(G - 1.0_RFREAL)*(U*U + V*V + W*W)

END FUNCTION MixtPerf_Co2_CGUVW

!******************************************************************************
!
! RCS Revision history:
!
! $Log: MixtPerf_C.F90,v $
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
! Revision 1.1  2004/12/01 16:48:47  haselbac
! Initial revision after changing case
!
! Revision 1.2  2002/05/28 13:44:44  haselbac
! Added new functions
!
! Revision 1.1  2002/05/04 16:16:52  haselbac
! Initial revision
!
!******************************************************************************

