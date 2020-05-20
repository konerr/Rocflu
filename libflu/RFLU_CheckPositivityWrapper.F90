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
! Purpose: Wrapper for check for posivity of variables
!
! Description: None.
!
! Input:
!   pRegion        Region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_CheckPositivityWrapper.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_CheckPositivityWrapper(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModParameters
  USE ModMPI

  USE ModInterfaces, ONLY: RFLU_CheckPositivity, &
                           RFLU_CheckPositivity_GL, &
                           RFLU_CheckValidity_GL 

#ifdef SPEC
  USE ModInterfaces, ONLY: RFLU_ScalarCheckPositivity
#endif

#ifdef PLAG
  USE PLAG_ModCheckVars, ONLY: PLAG_CheckPositivity
#endif

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: RFLU_CheckPositivityWrapper.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_CheckPositivityWrapper',__FILE__)

! ******************************************************************************
! Mixture
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%fluidModel )

! ==============================================================================
!   Incompressible fluid model
! ==============================================================================

     CASE ( FLUID_MODEL_INCOMP )

! ==============================================================================
!   Compressible fluid model. NOTE check state of solution vector.
! ==============================================================================

     CASE ( FLUID_MODEL_COMP )
       SELECT CASE ( pRegion%mixtInput%gasModel ) 
         CASE ( GAS_MODEL_TCPERF, & 
                GAS_MODEL_MIXT_TCPERF, & 
                GAS_MODEL_MIXT_TPERF, & 
                GAS_MODEL_MIXT_PSEUDO, &
                GAS_MODEL_MIXT_JWL )
           CALL RFLU_CheckPositivity(pRegion)
         CASE ( GAS_MODEL_MIXT_GASLIQ )
           CALL RFLU_CheckPositivity_GL(pRegion)
         CASE DEFAULT
           CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
       END SELECT ! pRegion%mixtInput%gasModel

! ==============================================================================
!   Default
! ==============================================================================

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixtInput%fluidModel


! ******************************************************************************
! Multiphysics modules
! ******************************************************************************

#ifdef SPEC
! ==============================================================================
! Species
! ==============================================================================

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    CALL RFLU_ScalarCheckPositivity(pRegion,FTYPE_SPEC, &
                                    pRegion%specInput%nSpecies, &
                                    pRegion%spec%cv)
  END IF ! global%specUsed
#endif

#ifdef PLAG
! ==============================================================================
! Particles
! ==============================================================================

  IF ( global%plagUsed .EQV. .TRUE. ) THEN
    CALL PLAG_CheckPositivity(pRegion)
  END IF ! global%plagUsed
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CheckPositivityWrapper

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CheckPositivityWrapper.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:34  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:52  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:48  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/03/26 20:21:27  haselbac
! Rewrite bcos of GL model, cosmetics
!
! Revision 1.3  2005/12/01 21:54:43  fnajjar
! Added call to PLAG_CheckPositivity
!
! Revision 1.2  2004/07/28 15:29:19  jferry
! created global variable for spec use
!
! Revision 1.1  2004/01/29 22:55:50  haselbac
! Initial revision
!
! ******************************************************************************

