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
! Purpose: set particle velocity (dependent variables) at earlier time.
!
! Description: none.
!
! Input: region  = data of current region.
!
! Output: region%levels%plag%dvOld
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_SetDependentVarsOld.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_SetDependentVarsOld( region )

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

  INTEGER :: iVars,nPcls
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pDv,pDvOld
  
  TYPE(t_plag) ,  POINTER :: pPlag  
  TYPE(t_global), POINTER :: global 
   
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_SetDependentVarsOld.F90,v $'

  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_SetDependentVarsOld',__FILE__ )

! Get dimensions --------------------------------------------------------------

  nPcls = region%plag%nPcls

! Set pointers ----------------------------------------------------------------

  pPlag  => region%plag
  pDv    => pPlag%dv 
  pDvOld => pPlag%dvOld

! Set old depedent variables for particles ------------------------------------

  DO iPcls = 1, nPcls
    DO iVars = 1, DV_PLAG_LAST 
      pDvOld(iVars,iPcls) = pDv(iVars,iPcls)
    END DO ! iVars
  END DO  ! iPcls

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_SetDependentVarsOld

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_SetDependentVarsOld.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
!******************************************************************************

