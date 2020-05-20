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
! Purpose: Collect relations for mixture speed of sound for gas-liquid model.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: MixtGasLiq_C.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

FUNCTION MixtGasLiq_C(Cvm,D,P,Dl,Dv,Dg,VFl,VFv,VFg,Cl2,Cv2,Cg2,Bl2,Bv2,Bg2)

   USE ModDataTypes

   IMPLICIT NONE

   REAL(RFREAL), INTENT(IN) :: Bg2,Bl2,Bv2,Cl2,Cv2,Cg2,Cvm,D,Dg,Dl,Dv,P,VFg, & 
                               VFl,VFv
   REAL(RFREAL) :: MixtGasLiq_C
   
   REAL(RFREAL) ::  denom,numer,term1,term2,term3

   term1 = Bl2*VFl*Dv*Cv2*Dg*Cg2
   term2 = Bv2*VFv*Dl*Cl2*Dg*Cg2
   term3 = Bg2*VFg*Dl*Cl2*Dv*Cv2

   numer  = D*Cvm*Dl*Cl2*Dv*Cv2*Dg*Cg2 + P*(term1 + term2 + term3)
   denom  = D*D*Cvm*(VFl*Dv*Cv2*Dg*Cg2 + VFv*Dl*Cl2*Dg*Cg2 + VFg*Dl*Cl2*Dv*Cv2)

   MixtGasLiq_C = SQRT(numer/denom)

END FUNCTION MixtGasLiq_C

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: MixtGasLiq_C.F90,v $
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
! Revision 1.1  2006/03/26 20:20:46  haselbac
! Initial revision
!
! ******************************************************************************

