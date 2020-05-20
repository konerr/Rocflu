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
! Purpose: Initialize particle solution in a region from scratch.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: PLAG_RFLU_InitSolutionScratch.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_InitSolutionScratch(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModParameters
   
  USE PLAG_ModParameters    

  USE PLAG_RFLU_ModFindCells
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iCont,iPcl
  INTEGER, POINTER, DIMENSION(:) :: pCvPlagMass
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv
  REAL(RFREAL) :: heatCapSum,massSum,massRatio,massFluxRatioLimit, &
                  massFluxRatioSum,massFluxRatioSumR
  REAL(RFREAL), POINTER, DIMENSION(:)   :: pDens,pIniComp,pSpcHeat
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv, pCv
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: PLAG_RFLU_InitSolutionScratch.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_InitSolutionScratch', &
                        __FILE__)
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Initializing particle solution from scratch...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and values
! ******************************************************************************

  pRegion%plag%nextIdNumber = pRegion%plag%nPcls
  
  pDens       => pRegion%plagInput%dens
  pSpcHeat    => pRegion%plagInput%spht
  pIniComp    => pRegion%plagInput%iniComp
  
  pCvPlagMass => pRegion%plag%cvPlagMass
  pAiv        => pRegion%plag%aiv
  pArv        => pRegion%plag%arv
  pCv         => pRegion%plag%cv
  
  massFluxRatioLimit = 1.0E-10_RFREAL
  massFluxRatioSum = SUM(pIniComp)

! ==============================================================================
! Set inverse of massFluxRatioSum to avoid division by zero
! ==============================================================================
  
  IF ( massFluxRatioSum > massFluxRatioLimit ) THEN
    massFluxRatioSumR = 1.0_RFREAL/massFluxRatioSum
  ELSE
    massFluxRatioSumR = 1.0_RFREAL
  END IF ! massFluxRatioSum 

! ******************************************************************************
! Define initial solution
! ******************************************************************************

  DO iPcl = 1,pRegion%plag%nPcls
    pCv(CV_PLAG_XPOS,iPcl) = pRegion%plagInput%iniPosX(iPcl)
    pCv(CV_PLAG_YPOS,iPcl) = pRegion%plagInput%iniPosY(iPcl)
    pCv(CV_PLAG_ZPOS,iPcl) = pRegion%plagInput%iniPosZ(iPcl)

    DO iCont = 1,pRegion%plagInput%nCont
      massRatio = pIniComp(iCont)*massFluxRatioSumR
      pCv(pCvPlagMass(iCont),iPcl) = pDens(iCont)*massRatio*global%pi/ & 
        6.0_RFREAL*pRegion%plagInput%iniDiam(iPcl)**3
    END DO ! iCont

    heatCapSum = SUM(pCv(pCvPlagMass(:),iPcl)*pSpcHeat(:))
    massSum    = SUM(pCv(pCvPlagMass(:),iPcl)) 

    pCv(CV_PLAG_XMOM,iPcl) = massSum*pRegion%plagInput%iniVelX(iPcl)
    pCv(CV_PLAG_YMOM,iPcl) = massSum*pRegion%plagInput%iniVelY(iPcl)
    pCv(CV_PLAG_ZMOM,iPcl) = massSum*pRegion%plagInput%iniVelZ(iPcl)
    pCv(CV_PLAG_ENER,iPcl) = heatCapSum*pRegion%plagInput%iniTemp(iPcl)
    
    pArv(ARV_PLAG_SPLOAD,iPcl) = pRegion%plagInput%iniSpLoad(iPcl)
    pAiv(AIV_PLAG_PIDINI,iPcl) = iPcl
    pAiv(AIV_PLAG_ICELLS,iPcl) = CRAZY_VALUE_INT ! Value used in initialization
    pAiv(AIV_PLAG_REGINI,iPcl) = CRAZY_VALUE_INT     
  END DO ! iPcl

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Initializing particle solution from scratch done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_InitSolutionScratch

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_InitSolutionScratch.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:52  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/05/16 22:33:46  fnajjar
! Modified to align with new bc datastructure
!
! Revision 1.1  2007/04/09 18:50:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:34  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.14  2007/03/08 15:04:19  fnajjar
! Fixed bug for massFluxRatioSum being zero and avoiding division by zero
!
! Revision 1.13  2006/05/05 18:05:26  haselbac
! Added init of ICELLS and REGINI
!
! Revision 1.12  2005/04/27 14:57:58  fnajjar
! Included module call to PLAG_RFLU_ModFindCells
!
! Revision 1.11  2005/03/31 20:25:59  fnajjar
! Added initial particle velocities for scratch solution
!
! Revision 1.10  2005/03/11 02:27:29  haselbac
! Simplified locating particles
!
! Revision 1.9  2004/11/04 16:29:19  fnajjar
! Deleted nPcls definition
!
! Revision 1.8  2004/10/11 22:11:27  haselbac
! Bug fix and renamed procedures
!
! Revision 1.7  2004/10/10 17:05:54  fnajjar
! Cleanup of duplicate variables
!
! Revision 1.6  2004/10/09 21:22:07  haselbac
! Bug fix: Brute force approach also requires c2f data structure
!
! Revision 1.5  2004/10/08 22:10:30  haselbac
! Added brute-force search option for finding particle positions
!
! Revision 1.4  2004/08/20 23:28:46  fnajjar
! Aligned with Plag prep tool
!
! Revision 1.3  2004/03/15 21:08:48  haselbac
! Changed name of particle location routine
!
! Revision 1.2  2004/03/08 22:18:30  fnajjar
! Commented off section specific to CFVEL test case
!
! Revision 1.1  2004/02/26 21:00:47  haselbac
! Initial revision
!
! ******************************************************************************

