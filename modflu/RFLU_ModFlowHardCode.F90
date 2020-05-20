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
! Purpose: Suite of routines to support hard-coded flow solutions.
!
! Description: None.
!
! Notes: 
!   1. Collected routines in single module because setting of parameters is
!      needed in at least three places: initialization of solution, setting 
!      of boundary profiles, and computation of errors. 
!
! ******************************************************************************
!
! $Id: RFLU_ModFlowHardCode.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModFlowHardCode

  USE ModParameters
  USE ModDataTypes
      
  IMPLICIT NONE
     
  PRIVATE
  PUBLIC :: RFLU_GetParamsHardCodePAcoust, & 
            RFLU_GetParamsHardCodeProudman, &
            RFLU_GetParamsHardCodeRingleb, & 
            RFLU_GetParamsHardCodeSsVortex

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModFlowHardCode.F90,v $ $Revision: 1.1.1.1 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  


! ******************************************************************************
!
! Purpose: Get parameters for pipe acoustics.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   pTot        Total pressure
!   aTot        Total speed of sound
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_GetParamsHardCodePAcoust(pTot,aTot)

    REAL(RFREAL), INTENT(OUT) :: aTot,pTot
    
    aTot = 340.0_RFREAL
    pTot = 1.0E+5_RFREAL
    
  END SUBROUTINE RFLU_GetParamsHardCodePAcoust
  
  
  
! ******************************************************************************
!
! Purpose: Get parameters for Proudman-Culick flow.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   dInc        Density
!   mInj        Injection mass flux
!   vInj        Injection velocity
!   pTot        Total pressure
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_GetParamsHardCodeProudman(dInc,mInj,vInj,pTot)

    REAL(RFREAL), INTENT(OUT) :: dInc,mInj,pTot,vInj
    
    dInc = 1.0_RFREAL
    mInj = 2.42_RFREAL
    vInj = mInj/dInc
    pTot = 1.0E+5_RFREAL
    
  END SUBROUTINE RFLU_GetParamsHardCodeProudman



! ******************************************************************************
!
! Purpose: Get parameters for Ringleb flow.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   pTot        Total pressure
!   tTot        Total temperature
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_GetParamsHardCodeRingleb(pTot,tTot)

    REAL(RFREAL), INTENT(OUT) :: pTot,tTot
    
    pTot = 1.0E+5_RFREAL
    tTot = 288.15_RFREAL
    
  END SUBROUTINE RFLU_GetParamsHardCodeRingleb



! ******************************************************************************
!
! Purpose: Get parameters for supersonic vortex flow.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   ri          Inner radius
!   Mi          Mach number at inner radius    
!   pTot        Total pressure
!   tTot        Total temperature
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_GetParamsHardCodeSsVortex(ri,Mi,pTot,tTot)

    REAL(RFREAL), INTENT(OUT) :: Mi,pTot,ri,tTot
    
    ri   = 1.0_RFREAL 
    Mi   = 2.25_RFREAL
    pTot = 1.0E+5_RFREAL
    tTot = 288.15_RFREAL
    
  END SUBROUTINE RFLU_GetParamsHardCodeSsVortex




! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModFlowHardCode


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModFlowHardCode.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:40  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:24  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:40  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.3  2005/03/15 20:45:06  haselbac
! Added routine to get parameters for pipe acoustics
!
! Revision 1.2  2004/07/06 15:14:39  haselbac
! Cosmetics only
!
! Revision 1.1  2004/02/23 23:01:49  haselbac
! Initial revision
!
! ******************************************************************************

