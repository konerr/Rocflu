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
! Purpose: searches list of materials for one with the given name
!
! Description: none.
!
! Input:  name = input name of material
!
! Output: material points to the correct element of global%materials(:)
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_SetMaterial.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_SetMaterial(global,material,name)

  USE ModDataTypes
  USE ModGlobal,    ONLY : t_global
  USE ModMaterials, ONLY : t_material
  USE ModError
  IMPLICIT NONE

! ... parameters
  TYPE(t_global),   POINTER    :: global
  TYPE(t_material), POINTER    :: material
  CHARACTER(*),     INTENT(in) :: name

! ... loop variables
  INTEGER :: iMat

! ... local variables
  LOGICAL :: found

  CHARACTER(CHRLEN) :: RCSIdentString

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_SetMaterial.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction( global,'INRT_SetMaterial',__FILE__ )

! begin -----------------------------------------------------------------------

  found = .FALSE.

  DO iMat = 1,global%nMaterials

    IF (TRIM(name) == TRIM(global%materials(iMat)%name)) THEN
      material => global%materials(iMat)
      found = .TRUE.
      EXIT
    END IF ! name

  END DO ! iMat

  IF (.NOT.found) CALL ErrorStop( global,ERR_INRT_BADMAT,__LINE__,TRIM(name))

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_SetMaterial

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_SetMaterial.F90,v $
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
! Revision 1.1  2004/12/01 21:56:42  fnajjar
! Initial revision after changing case
!
! Revision 1.1  2003/03/24 23:23:25  jferry
! converted from libfloflu routine to rocinteract routine
!
!******************************************************************************

