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
! Description: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: SPEC_ModInterfaces.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001-2003 by the University of Illinois
!
!******************************************************************************

MODULE SPEC_ModInterfaces

  IMPLICIT NONE

  INTERFACE

    SUBROUTINE SPEC_CheckUserInput(regions)
      USE ModDataStruct, ONLY: t_region
      TYPE(t_region), DIMENSION(:), POINTER :: regions
    END SUBROUTINE SPEC_CheckUserInput

    SUBROUTINE SPEC_DerivedInputValues(regions)
      USE ModDataStruct, ONLY: t_region
      TYPE(t_region), DIMENSION(:), POINTER :: regions
    END SUBROUTINE SPEC_DerivedInputValues

    SUBROUTINE SPEC_EqEulCorrPatch(pRegion,pPatch,iSpec)
      USE ModDataStruct, ONLY : t_region
      USE ModBndPatch, ONLY: t_patch
      TYPE(t_region), POINTER :: pRegion
      TYPE(t_patch), POINTER :: pPatch
      INTEGER, INTENT(IN) :: iSpec
    END SUBROUTINE SPEC_EqEulCorrPatch

    SUBROUTINE SPEC_InitInputValues(regions)
      USE ModDataStruct, ONLY: t_region
      TYPE(t_region), DIMENSION(:), POINTER :: regions
    END SUBROUTINE SPEC_InitInputValues

    SUBROUTINE SPEC_InitInputValuesSpecType(region,iSpecType)
      USE ModDataStruct, ONLY: t_region
      INTEGER, INTENT(IN) :: iSpecType
      TYPE(t_region) :: region
    END SUBROUTINE SPEC_InitInputValuesSpecType

    SUBROUTINE SPEC_ReadInputFile(regions)
      USE ModDataStruct, ONLY: t_region
      TYPE(t_region), DIMENSION(:), POINTER :: regions
    END SUBROUTINE SPEC_ReadInputFile

    SUBROUTINE SPEC_ReadSpecSection(regions)
      USE ModDataStruct, ONLY: t_region
      TYPE(t_region), DIMENSION(:), POINTER :: regions
    END SUBROUTINE SPEC_ReadSpecSection

    SUBROUTINE SPEC_ReadSpecTypeSection(regions,iSpecType)
      USE ModDataStruct, ONLY: t_region
      INTEGER :: iSpecType
      TYPE(t_region), DIMENSION(:), POINTER :: regions
    END SUBROUTINE SPEC_ReadSpecTypeSection

  END INTERFACE

END MODULE SPEC_ModInterfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_ModInterfaces.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:53  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:51:23  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2004/07/30 22:47:37  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.2  2003/11/25 21:08:33  haselbac
! Added interfaces
!
! Revision 1.1.1.1  2001/12/03 21:44:04  jblazek
! Import of RocfluidMP
!
!******************************************************************************

