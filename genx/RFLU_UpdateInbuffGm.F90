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
! Purpose: Dummy routine to fill the geometry buffer.
!
! Description: None.
!
! Input: 
!   globalGenx	Pointer to global data
!   dAlpha    	Relative time step
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_UpdateInbuffGm.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_UpdateInbuffGm(globalGenx,dAlpha)

  USE ModDataTypes
  USE ModGenx, ONLY: t_globalGenx
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  
  IMPLICIT NONE

! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  DOUBLE PRECISION, INTENT(IN) :: dAlpha
  TYPE(t_globalGenx), POINTER :: globalGenx

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iPatch,iReg,iv
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_UpdateInbuffGm.F90,v $ $Revision: 1.1.1.1 $'

  global => globalGenx%global

  CALL RegisterFunction(global,'RFLU_UpdateInbuffGm',__FILE__)

! *****************************************************************************
! Fill displacement buffer
! *****************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => globalGenx%levels(1)%regions(iReg)    
  
    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)
    
      DO iv = 1,pPatch%nBVert
        pPatch%duAlp(XCOORD,iv) = 0.0_RFREAL
        pPatch%duAlp(YCOORD,iv) = 0.0_RFREAL
        pPatch%duAlp(ZCOORD,iv) = 0.0_RFREAL
      END DO ! iv

      DO iv = pPatch%nBVert+1,pPatch%nBVertTot
        pPatch%duAlp(XCOORD,iv) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%duAlp(YCOORD,iv) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%duAlp(ZCOORD,iv) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      END DO ! iv
    END DO ! iPatch
  END DO ! iReg

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_UpdateInbuffGm

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_UpdateInbuffGm.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:30  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:47:51  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:58:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2003/05/13 23:45:02  haselbac
! Added init of dummy vertices
!
! Revision 1.3  2003/04/12 21:35:01  haselbac
! Cosmetics only
!
! Revision 1.2  2002/12/03 22:58:27  haselbac
! Added proper code
!
! Revision 1.1  2002/10/05 18:28:18  haselbac
! Initial revision
!
!******************************************************************************

