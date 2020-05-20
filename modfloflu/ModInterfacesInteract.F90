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
! Purpose: set explicit interfaces to subroutines and functions.
!
! Description: none
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ModInterfacesInteract.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

MODULE ModInterfacesInteract

  IMPLICIT NONE

  INTERFACE

! =============================================================================
! Multiphase Interactions
! =============================================================================

  SUBROUTINE INRT_BuildVersionString( versionString )
    CHARACTER(*) :: versionString
  END SUBROUTINE INRT_BuildVersionString

  SUBROUTINE INRT_BurnStatusUpdate( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT), TARGET :: region
  END SUBROUTINE INRT_BurnStatusUpdate

  SUBROUTINE INRT_PrintMaterialInput( global )
    USE ModGlobal,    ONLY : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE INRT_PrintMaterialInput

  SUBROUTINE INRT_PrintUserInput( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(IN) :: region
  END SUBROUTINE INRT_PrintUserInput

  SUBROUTINE INRT_ReadMaterialInput( global )
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE INRT_ReadMaterialInput

  SUBROUTINE INRT_SetMaterial(global,material,name)
    USE ModGlobal,    ONLY : t_global
    USE ModMaterials, ONLY : t_material
    TYPE(t_global),   POINTER    :: global
    TYPE(t_material), POINTER    :: material
    CHARACTER(*),     INTENT(in) :: name
  END SUBROUTINE INRT_SetMaterial

  SUBROUTINE INRT_SetParticleTemp( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE INRT_SetParticleTemp

  SUBROUTINE INRT_SourceTerms( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE INRT_SourceTerms

  SUBROUTINE INRT_UserInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE INRT_UserInput

  SUBROUTINE INRT_VaporEnergyConversion( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT), TARGET :: region
  END SUBROUTINE INRT_VaporEnergyConversion

  END INTERFACE

END MODULE ModInterfacesInteract

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModInterfacesInteract.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
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
! Revision 1.6  2004/12/01 00:09:25  wasistho
! added BuildVersionString
!
! Revision 1.5  2004/07/27 21:27:14  jferry
! removed rocinteract allocation routines (moved to rocpart)
!
! Revision 1.4  2004/03/02 21:47:28  jferry
! Added After Update interactions
!
! Revision 1.3  2003/03/24 23:25:48  jferry
! moved some routines from libfloflu to rocinteract
!
! Revision 1.2  2003/03/11 15:59:16  jferry
! Added routines to rocinteract
!
! Revision 1.1  2003/03/04 22:12:34  jferry
! Initial import of Rocinteract
!
!******************************************************************************

