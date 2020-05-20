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
! Purpose: write out user input for materials
!
! Description: none.
!
! Input: global = user input.
!
! Output: to standard output.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_PrintMaterialInput.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_PrintMaterialInput( global ) ! PUBLIC

  USE ModDataTypes
  USE ModGlobal,    ONLY : t_global
  USE ModError
  USE ModParameters
  USE ModMaterials, ONLY : t_material
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER :: imat

! ... local variables
  TYPE(t_material), POINTER :: material

  CHARACTER(CHRLEN) :: RCSIdentString

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_PrintMaterialInput.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction( global,'INRT_PrintMaterialInput',__FILE__ )

! begin -----------------------------------------------------------------------

  WRITE(STDOUT,*)
  WRITE(STDOUT,1030) SOLVER_NAME//' Materials:'
  WRITE(STDOUT,1010) SOLVER_NAME//'   Number of materials defined', &
    global%nMaterials
  WRITE(STDOUT,*)

  DO iMat = 1,global%nMaterials

    material => global%materials(iMat)

    WRITE(STDOUT,1010) SOLVER_NAME//'   *** Material',iMat
    WRITE(STDOUT,1030) SOLVER_NAME//'     name    =  '//TRIM(material%name)
    SELECT CASE (material%phase)
    CASE (1)
      WRITE(STDOUT,1030) SOLVER_NAME//'     phase   =  GAS'
    CASE (2)
      WRITE(STDOUT,1030) SOLVER_NAME//'     phase   =  LIQUID'
    CASE (3)
      WRITE(STDOUT,1030) SOLVER_NAME//'     phase   =  SOLID'
    CASE DEFAULT
      WRITE(STDOUT,1030) SOLVER_NAME//'***WARNING: no phase specified.'// &
                                      '  Setting phase = GAS'
      material%phase = 1
    END SELECT
    WRITE(STDOUT,1010) SOLVER_NAME//'     index  ', material%index
    WRITE(STDOUT,1021) SOLVER_NAME//'     mol.wt.', material%molw
    WRITE(STDOUT,1021) SOLVER_NAME//'     density', material%dens
    WRITE(STDOUT,1021) SOLVER_NAME//'     sp.heat', material%spht
    WRITE(STDOUT,1021) SOLVER_NAME//'     surf.tension', material%surftens
    WRITE(STDOUT,1021) SOLVER_NAME//'     boiling pt.', material%Tboil
    WRITE(STDOUT,1021) SOLVER_NAME//'     melting pt.', material%Tmelt
    WRITE(STDOUT,1021) SOLVER_NAME//'     ref. visc.', material%refVisc
    WRITE(STDOUT,1021) SOLVER_NAME//'     suth. temp.', material%suthTemp
    WRITE(STDOUT,1021) SOLVER_NAME//'     suth. coef.', material%suthCoef
    WRITE(STDOUT,1021) SOLVER_NAME//'     Prandtl no.', material%pr
    WRITE(STDOUT,*)

  END DO ! iMat

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

1010 FORMAT(A,' =',I3)
1020 FORMAT(A,' =',ES15.5)
1021 FORMAT(A,' =',EN15.5)
1030 FORMAT(A)

END SUBROUTINE INRT_PrintMaterialInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_PrintMaterialInput.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2009/07/09 20:44:21  mparmar
! Added printing for viscosity and Prandtl number for material
!
! Revision 1.3  2008/12/06 08:43:49  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:02  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:50:11  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:15  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 21:56:29  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2004/07/23 22:43:16  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.4  2004/04/15 16:04:21  jferry
! minor formatting (removed trailing spaces)
!
! Revision 1.3  2004/03/02 21:49:46  jferry
! Added melting and boiling point to material definitions
!
! Revision 1.2  2003/09/13 20:17:31  fnajjar
! Added surface tension to Materials datastructure
!
! Revision 1.1  2003/03/24 23:23:25  jferry
! converted from libfloflu routine to rocinteract routine
!
!******************************************************************************

