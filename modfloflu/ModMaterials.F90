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
! Purpose: define the properties of materials used throughout the code
!
! Description: none.
!
! Notes:
!
!   this module will grow as the modeling for materials becomes more
!   sophisticated
!
!******************************************************************************
!
! $Id: ModMaterials.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

MODULE ModMaterials

  USE ModDataTypes
  IMPLICIT NONE

! data types ------------------------------------------------------------------

  TYPE t_material
    CHARACTER(CHRLEN) :: name     ! Name of material
    INTEGER           :: phase    ! 1 = Gas, 2 = Liquid, 3 = Solid
    INTEGER           :: index    ! Index of material in global%materials(:)
    REAL(RFREAL)      :: molw     ! Molecular weight
    REAL(RFREAL)      :: dens     ! Density
    REAL(RFREAL)      :: spht     ! Specific heat
    REAL(RFREAL)      :: surftens ! Surface tension
    REAL(RFREAL)      :: Tboil    ! Boiling (or condensation) temperature
    REAL(RFREAL)      :: Tmelt    ! Melting (or freezing) temperature
    REAL(RFREAL)      :: refVisc  ! Reference viscosity used in Sutherland's law
    REAL(RFREAL)      :: suthTemp ! Reference temperature used in Sutherland's law
    REAL(RFREAL)      :: suthCoef ! Sutherland's coefficient
    REAL(RFREAL)      :: pr       ! Sutherland's coefficient
    REAL(RFREAL)      :: cond     ! Sutherland's coefficient
    REAL(RFREAL)      :: detonVel ! Detonation velocity
  END TYPE t_material

END MODULE ModMaterials

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModMaterials.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2009/07/09 20:44:03  mparmar
! Added viscosity, conductivity and Prandtl number to material
!
! Revision 1.3  2008/12/06 08:43:37  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:10  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:17  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.5  2004/07/23 22:43:16  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.4  2004/03/02 21:49:45  jferry
! Added melting and boiling point to material definitions
!
! Revision 1.3  2003/09/13 20:15:55  fnajjar
! Added surface tension to Materials datastructure
!
! Revision 1.2  2003/03/12 18:06:13  jferry
! added index component to t_material
!
! Revision 1.1  2003/03/11 15:53:25  jferry
! Created data type for material properties
!
!******************************************************************************

