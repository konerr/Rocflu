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
! Purpose: Returns species index given the material name
!
! Description: none.
!
! Input:  name = input name of material
!
! Output: Index number of species in spec%cv(:,icg)
!
! Notes: It is assumed that same material is not used in 2 different species.
!        A check for this has already been added to SPEC_ReadInputFile().
!
!******************************************************************************
!
! $Id: SPEC_GetSpeciesIndex.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2011 by the University of Illinois
!
!******************************************************************************

FUNCTION SPEC_GetSpeciesIndex(global,pSpecInput,name)

  USE ModGlobal,    ONLY : t_global
  USE ModDataTypes
  USE ModError
  USE ModParameters

  USE ModSpecies, ONLY: t_spec_input

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global
  TYPE(t_spec_input), POINTER :: pSpecInput
  CHARACTER(*), INTENT(in) :: name
  INTEGER :: SPEC_GetSpeciesIndex

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: found
  INTEGER :: iSpec
  CHARACTER(CHRLEN) :: RCSIdentString

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_GetSpeciesIndex.F90,v $'

  CALL RegisterFunction( global,'SPEC_GetSpeciesIndex',__FILE__ )

  found = .FALSE.

  DO iSpec = 1,pSpecInput%nSpecies
    IF (TRIM(name) == TRIM(pSpecInput%specType(iSpec)%pMaterial%name)) THEN
      SPEC_GetSpeciesIndex = iSpec
      found = .TRUE.
      EXIT
    END IF ! name
  END DO ! iSpec

  IF (.NOT.found) CALL ErrorStop( global,ERR_SPEC_BADMAT,__LINE__,TRIM(name))

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)  

END FUNCTION SPEC_GetSpeciesIndex

!******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_GetSpeciesIndex.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
!******************************************************************************
