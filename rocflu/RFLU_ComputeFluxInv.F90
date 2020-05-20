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
! Purpose: Wrapper function for computing inviscid fluxes of mixture.
!
! Description: None.
!
! Input: 
!   pRegion	Pointer to region data
!   fluxPart	Part of flux (central or dissipation or both)
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ComputeFluxInv.F90,v 1.2 2015/12/18 23:02:15 rahul Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************


SUBROUTINE RFLU_ComputeFluxInv(pRegion,fluxPart)

  USE ModDataTypes
  USE ModParameters  
  USE ModError  
  USE ModGlobal, ONLY: t_global  
  USE ModDataStruct, ONLY: t_region
  
  USE RFLU_ModAUSMFlux
  USE RFLU_ModAUSMPlusUpFlux
  USE RFLU_ModHLLCFlux
  USE RFLU_ModRoeFlux

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: fluxPart
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ComputeFluxInv.F90,v $ $Revision: 1.2 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeFluxInv',__FILE__)
          
! ******************************************************************************
! Compute fluxes
! ******************************************************************************
         
  SELECT CASE ( pRegion%mixtInput%spaceDiscr )
    CASE ( DISCR_UPW_ROE )        
      CALL RFLU_ROE_ComputeFlux(pRegion,fluxPart)
    CASE ( DISCR_UPW_HLLC )
      IF ( fluxPart == FLUX_PART_CENTRAL .OR. fluxPart == FLUX_PART_BOTH ) THEN
        CALL RFLU_HLLC_ComputeFlux(pRegion)      
      END IF ! fluxPart
    CASE ( DISCR_UPW_AUSMPLUS )     
      IF ( fluxPart == FLUX_PART_CENTRAL .OR. fluxPart == FLUX_PART_BOTH ) THEN     
        CALL RFLU_AUSM_ComputeFlux(pRegion)
      END IF ! fluxPart
    CASE ( DISCR_UPW_AUSMPLUSUP )     
      IF ( fluxPart == FLUX_PART_CENTRAL .OR. fluxPart == FLUX_PART_BOTH ) THEN     
        CALL RFLU_AUSMPlusUp_ComputeFlux(pRegion)
      END IF ! fluxPart
  CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
  END SELECT ! pMixtInput%spaceDiscr

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputeFluxInv

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeFluxInv.F90,v $
! Revision 1.2  2015/12/18 23:02:15  rahul
! Added call to AUSM+up scheme.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:47  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:59  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:56  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:01  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2006/05/01 21:02:58  haselbac
! Initial revision
!
! ******************************************************************************

