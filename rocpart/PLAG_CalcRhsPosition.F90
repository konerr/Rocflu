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
! Purpose: compute RHS for Lagrangian particle position field.
!
! Description: none.
!
! Input: region  = data of current region.
!
! Output: region%levels%plag%rhs
!
! Notes: Use negative values of rhs for consistent RK Updating step.
!
!******************************************************************************
!
! $Id: PLAG_CalcRhsPosition.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CalcRhsPosition( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region

! ... loop variables
  INTEGER :: iPcls

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: nPcls
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pDv, pRhs
  
  TYPE(t_plag) ,  POINTER :: pPlag  
  TYPE(t_global), POINTER :: global 
   
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CalcRhsPosition.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_CalcRhsPosition',__FILE__ )

! Get dimensions --------------------------------------------------------------

  nPcls = region%plag%nPcls

! Set pointers ----------------------------------------------------------------

  pPlag => region%plag
  pDv   => pPlag%dv 
  pRhs  => pPlag%rhs

! Calculate RHS for position vector field -------------------------------------

  DO iPcls = 1, nPcls      
    pRhs(CV_PLAG_XPOS,iPcls) = -pDv(DV_PLAG_UVEL,iPcls)    
    pRhs(CV_PLAG_YPOS,iPcls) = -pDv(DV_PLAG_VVEL,iPcls)       
    pRhs(CV_PLAG_ZPOS,iPcls) = -pDv(DV_PLAG_WVEL,iPcls)         
  END DO  ! iPcls

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CalcRhsPosition

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CalcRhsPosition.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:50  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:02  mtcampbe
! Added Illinois Open Source License/Copyright
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
! Revision 1.1  2004/12/01 20:57:03  fnajjar
! Initial revision after changing case
!
! Revision 1.2  2004/02/06 21:18:20  fnajjar
! Initial Integration of Rocpart with Rocflu
!
! Revision 1.1  2002/10/25 14:14:44  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************

