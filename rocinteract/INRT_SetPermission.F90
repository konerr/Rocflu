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
! Purpose: converts a real number input to a Permission level
!
! Description: none.
!
! Input: val = real number from a ReadSection call
!
! Output: perm = integer used to designate Permission level
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_SetPermission.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************
SUBROUTINE INRT_SetPermission( global,val,perm )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  USE INRT_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER     :: global
  REAL(RFREAL),   INTENT(IN)  :: val
  INTEGER,        INTENT(OUT) :: perm

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_SetPermission.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction( global,'INRT_SetPermission',__FILE__ )

! begin -----------------------------------------------------------------------

  SELECT CASE (NINT(val))

  CASE (0)                 ! value used in input deck
    perm = INRT_PERM_BLOCK ! value used inside the code

  CASE (1)                 ! value used in input deck
    perm = INRT_PERM_PMASS ! value used inside the code

  CASE (2)                 ! value used in input deck
    perm = INRT_PERM_PMOME ! value used inside the code

  CASE (3)                 ! value used in input deck
    perm = INRT_PERM_PALL  ! value used inside the code

  CASE DEFAULT
    CALL ErrorStop( global,ERR_INRT_BADPERM,__LINE__ )

  END SELECT ! val

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_SetPermission

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_SetPermission.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:50  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:02  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:50:12  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:15  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 21:56:45  fnajjar
! Initial revision after changing case
!
! Revision 1.1  2003/04/02 22:32:04  jferry
! codified Activeness and Permission structures for rocinteract
!
!******************************************************************************

