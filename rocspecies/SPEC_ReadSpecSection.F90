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
! Purpose: Read user input applicable to all species.
!
! Description: None.
!
! Input:
!   regions                Data associated with regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: SPEC_ReadSpecSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_ReadSpecSection(regions)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters

#ifdef RFLU
  USE ModInterfaces, ONLY: ReadSection
#endif

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER, PARAMETER :: NKEYS = 2

  CHARACTER(CHRLEN) :: keys(NKEYS)
  CHARACTER(CHRLEN) :: RCSIdentString
  LOGICAL :: defined(NKEYS),usedFlag
  INTEGER :: errorFlag
  REAL(RFREAL) :: vals(NKEYS)
  TYPE(t_global), POINTER :: global

#ifdef RFLU
  INTEGER :: iReg
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_ReadSpecSection.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'SPEC_ReadSpecSection',__FILE__)

! ******************************************************************************
! Read user input for species
! ******************************************************************************

  keys(1) = 'USED'
  keys(2) = 'NSPECIES'

#ifdef RFLU
  CALL ReadSection(regions(1)%global,IF_INPUT,NKEYS,keys,vals,defined)

  IF ( defined(1) .EQV. .TRUE. ) THEN
    usedFlag = (NINT(vals(1)) /= 0)   ! if keyword USED appears, evaluate
  ELSE
    usedFlag = .TRUE.                 ! if it does not appear, default = true
  END IF ! defined(1)

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
    regions(iReg)%specInput%usedFlag = usedFlag
  END DO ! iReg

  IF ( (defined(2) .EQV. .TRUE.) .AND. (usedFlag .EQV. .TRUE.) ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%specInput%nSpecies = NINT(vals(2))
      regions(iReg)%global%nSpecies    = regions(iReg)%specInput%nSpecies
    END DO ! iReg
  END IF ! defined
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_ReadSpecSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_ReadSpecSection.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:53  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:05  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:51:23  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:50  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2005/10/05 20:10:38  haselbac
! Added setting of global%nSpecies
!
! Revision 1.2  2004/07/28 15:31:34  jferry
! added USED field to SPECIES input section
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
! ******************************************************************************

