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
! Purpose: Compute EE approximation to velocity and temperature of particulate
!   species.
!
! Description: None.
!
! Input:
!   pRegion        Pointer to data of current region
!   iSpec          Index of species
!
! Output: None.
!
! Notes: 
!   1. Use nomenclature (gas instead of mixture) of document describing
!      governing equations for compressible multiphase flow.
!
! ******************************************************************************
!
! $Id: SPEC_RFLU_SetEEv.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_RFLU_SetEEv(pRegion,iSpec)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModParameters

  USE ModSpecies, ONLY: t_spec_type
  
  USE SPEC_ModParameters

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: iSpec
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: icg,iSpecEEv
  REAL(RFREAL) :: gx,gy,gz,irg,tau,tauCoef,Tee,Tg,ug,uee,vg,vee,wg,wee            
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvG,pDvG,pSd,pTvG
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pEEv  
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_spec_type), POINTER :: pSpecType

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_RFLU_SetEEv.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_SetEEv',__FILE__)

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pSd   => pRegion%mixt%sd
  pCvG  => pRegion%mixt%cv
  pDvG  => pRegion%mixt%dv    
  pTvG  => pRegion%mixt%tv
  pEEv  => pRegion%spec%eev

  pSpecType => pRegion%specInput%specType(iSpec)  

  tauCoef  = pSpecType%tauCoefficient
  iSpecEEv = pSpecType%iSpec2iSpecEEv

! ******************************************************************************
! Checks: Defensive coding, should never occur
! ******************************************************************************

  IF ( pSpecType%discreteFlag .EQV. .FALSE. ) THEN
    CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__)
  END IF ! pSpecType%discreteFlag

  IF ( pSpecType%velocityMethod /= SPEC_METHV_EQEUL ) THEN
    CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__)
  END IF ! pSpecType%velocityMethod

  IF ( tauCoef <= 0.0_RFREAL ) THEN
    CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__)
  END IF ! tauCoef

  IF ( pRegion%mixtInput%indSd /= 1 ) THEN
    CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__)
  END IF ! pRegion%mixtInput%indSd

  IF ( pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! pRegion%mixt%cvState

! ******************************************************************************
! Set gravity vector for settling velocity contribution
! ******************************************************************************

  IF ( pSpecType%settlingFlag .EQV. .TRUE. ) THEN 
    IF ( global%accelOn .EQV. .TRUE. ) THEN 
      gx = global%accelX
      gy = global%accelY
      gz = global%accelZ
    ELSE 
      gx = 0.0_RFREAL
      gy = 0.0_RFREAL
      gz = 0.0_RFREAL
    END IF ! global%accelOn
  ELSE 
    gx = 0.0_RFREAL
    gy = 0.0_RFREAL
    gz = 0.0_RFREAL              
  END IF ! pSpecType%settlingFlag

! ******************************************************************************
! Compute EE variables
! ******************************************************************************

  DO icg = 1,pGrid%nCellsTot
    irg = 1.0_RFREAL/pCvG(CV_MIXT_DENS,icg)
    ug  = irg*pCvG(CV_MIXT_XMOM,icg)
    vg  = irg*pCvG(CV_MIXT_YMOM,icg)
    wg  = irg*pCvG(CV_MIXT_ZMOM,icg)
    Tg  =     pDvG(DV_MIXT_TEMP,icg)

    tau = tauCoef/pTvG(TV_MIXT_MUEL,icg)

    uee = ug - tau*(pSd(SD_XMOM,icg) - gx)
    vee = vg - tau*(pSd(SD_YMOM,icg) - gy)
    wee = wg - tau*(pSd(SD_ZMOM,icg) - gz)

    Tee = Tg ! TEMPORARY               
    
    pEEv(EEV_SPEC_XVEL,iSpecEEv,icg) = uee
    pEEv(EEV_SPEC_YVEL,iSpecEEv,icg) = vee
    pEEv(EEV_SPEC_ZVEL,iSpecEEv,icg) = wee
    pEEv(EEV_SPEC_TEMP,iSpecEEv,icg) = Tee 
  END DO ! icg

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_RFLU_SetEEv

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_RFLU_SetEEv.F90,v $
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
! Revision 1.3  2006/04/07 15:19:25  haselbac
! Removed tabs
!
! Revision 1.2  2006/01/06 22:43:14  haselbac
! Bug fix: Incorrect declaration of iSpecEEv
!
! Revision 1.1  2005/11/27 01:47:26  haselbac
! Initial revision
!
! ******************************************************************************

