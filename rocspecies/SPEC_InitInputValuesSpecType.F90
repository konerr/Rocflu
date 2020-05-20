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
! Purpose: Initialize user input for species to default values.
!
! Description: None.
!
! Input:
!   regions        Region data
!   iSpecType        Species type
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: SPEC_InitInputValuesSpecType.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE SPEC_InitInputValuesSpecType(region,iSpecType)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModSpecies, ONLY: t_spec_input
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters
  USE ModMPI

  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: iSpecType
  TYPE(t_region) :: region

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = &
    '$RCSfile: SPEC_InitInputValuesSpecType.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction(global,'SPEC_InitInputValuesSpecType',__FILE__)

! ******************************************************************************
! Initialize
! ******************************************************************************

! ==============================================================================
! Common to all species
! ==============================================================================

  region%specInput%specType(iSpecType)%frozenFlag     = .FALSE.
  region%specInput%specType(iSpecType)%initVal        = 1.0_RFREAL
  region%specInput%specType(iSpecType)%sourceType     = SPEC_SOURCE_TYPE_NONE
  region%specInput%specType(iSpecType)%schmidtNumber  = 1.0_RFREAL

! ==============================================================================
! Specific to species representing particles
! ==============================================================================
  
  region%specInput%specType(iSpecType)%iSpec2iSpecEEv = CRAZY_VALUE_INT
  region%specInput%specType(iSpecType)%iSpecEEv2iSpec = CRAZY_VALUE_INT  
  region%specInput%specType(iSpecType)%diameter       = 0.0_RFREAL
  region%specInput%specType(iSpecType)%puffFactor     = 1.0_RFREAL
  region%specInput%specType(iSpecType)%velocityMethod = SPEC_METHV_FLUIDVEL
  region%specInput%specType(iSpecType)%settlingFlag   = .FALSE.

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_InitInputValuesSpecType

!******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_InitInputValuesSpecType.F90,v $
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
! Revision 1.1  2007/04/09 18:51:22  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.6  2005/11/27 01:54:46  haselbac
! Added init for EEv variables, cosmetics
!
! Revision 1.5  2005/11/10 02:35:13  haselbac
! Added init for settlingFlag
!
! Revision 1.4  2004/07/30 22:47:37  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.3  2004/07/23 22:43:17  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.2  2004/01/29 22:59:34  haselbac
! Added Schmidt number
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
!******************************************************************************

