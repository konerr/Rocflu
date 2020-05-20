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
! Purpose: Initialize particle solution in a region from random state
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
! $Id: PLAG_RFLU_InitSolutionRandom.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_InitSolutionRandom(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModParameters
  
  USE ModRandom, ONLY: Rand1Uniform
  
  USE PLAG_ModParameters    
     
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
  INTEGER :: errorFlag,iCont,infLoopCntr,infLoopCntrMax,iPcl,nPcls
  INTEGER, POINTER, DIMENSION(:) :: pCvPlagMass
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv
  REAL(RFREAL) :: delFrac,diamDel,diamMax,diamMin,diam,heatCapSum,massRatio, &
                  massFluxRatioSum,massFluxRatioSumR,massFluxRatioLimit, &
                  massSum,nPclsFrac,rn,spLoadDel,spLoadMax,spLoadMin,spLoad, &
                  tempDel,tempMax,tempMin,temp,uDel,uMax,uMin,u,vDel,vMin, &
                  vMax,v,wDel,wMin,wMax,w,xDel,xLoc,xMax,xMin,yDel,yLoc,yMax, &
                  yMin,zDel,zLoc,zMax,zMin
  REAL(RFREAL), POINTER, DIMENSION(:) :: pDens,pIniComp,pSpcHeat
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv,pCv
  TYPE(t_global), POINTER :: global
  TYPE(t_grid),   POINTER :: pGrid

! Subbu - Vars for scaling superloading
  REAL(RFREAL) :: rad,scal,the
! Subbu - End Vars for scaling superloading

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: PLAG_RFLU_InitSolutionRandom.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_InitSolutionRandom', &
                        __FILE__)
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                          'Initializing particle solution from random state...'
  END IF ! global%verbLevel
  
! ******************************************************************************
! Set pointers and values
! ******************************************************************************

  pGrid => pRegion%grid

  delFrac = 0.01_RFREAL
  
  pRegion%plag%nextIdNumber = pRegion%plag%nPcls
  
  nPcls       = pRegion%plag%nPcls
  nPclsFrac   = 0.1_RFREAL    
  iPcl        = 0
  infLoopCntr = 0 

  IF ( nPcls > HUGE(1)/100 ) THEN
    infLoopCntrMax = HUGE(1)
  ELSE
    infLoopCntrMax = 100*nPcls
  END IF ! nPcls
  
  massFluxRatioLimit = 1.0E-10_RFREAL

! ******************************************************************************  
! Build bounding box. NOTE not from grid, but from user input
! ******************************************************************************

  xMin = pRegion%plagInput%iniRandXMin
  xMax = pRegion%plagInput%iniRandXMax
  yMin = pRegion%plagInput%iniRandYMin
  yMax = pRegion%plagInput%iniRandYMax
  zMin = pRegion%plagInput%iniRandZMin
  zMax = pRegion%plagInput%iniRandZMax

  xDel = xMax - xMin 
  yDel = yMax - yMin
  zDel = zMax - zMin

! ******************************************************************************
! Set pointers and values for PLAG infrastructure
! ******************************************************************************

  pDens    => pRegion%plagInput%dens
  pSpcHeat => pRegion%plagInput%spht
  pIniComp => pRegion%plagInput%iniComp
  
  pCvPlagMass => pRegion%plag%cvPlagMass
  pAiv => pRegion%plag%aiv
  pArv => pRegion%plag%arv
  pCv  => pRegion%plag%cv
  
  massFluxRatioSum = SUM(pIniComp)

! ==============================================================================  
! Set inverse of massFluxRatioSum to avoid division by zero 
! ==============================================================================  

  IF ( massFluxRatioSum > massFluxRatioLimit ) THEN 
    massFluxRatioSumR = 1.0_RFREAL/massFluxRatioSum
  ELSE
    massFluxRatioSumR = 1.0_RFREAL
  END IF ! massFluxRatioSum 

  diamMin   = pRegion%plagInput%iniRandDiamMin
  diamMax   = pRegion%plagInput%iniRandDiamMax
  tempMin   = pRegion%plagInput%iniRandTempMin
  tempMax   = pRegion%plagInput%iniRandTempMax 
  spLoadMin = pRegion%plagInput%iniRandSpLoadMin
  spLoadMax = pRegion%plagInput%iniRandSpLoadMax

  uMin      = pRegion%plagInput%iniRandUMin
  uMax      = pRegion%plagInput%iniRandUMax
  vMin      = pRegion%plagInput%iniRandVMin
  vMax      = pRegion%plagInput%iniRandVMax
  wMin      = pRegion%plagInput%iniRandWMin
  wMax      = pRegion%plagInput%iniRandWMax
  
  diamDel   = diamMax   - diamMin
  tempDel   = tempMax   - tempMin
  spLoadDel = spLoadMax - spLoadMin

  uDel      = uMax      - uMin
  vDel      = vMax      - vMin
  wDel      = wMax      - wMin

! ******************************************************************************
! Define initial solution through an infinite loop. Exit loop once nPcls is 
! reached
! ******************************************************************************

  infLoop: DO 
    infLoopCntr = infLoopCntr + 1

    IF ( iPcl == nPcls ) THEN 
      EXIT infLoop
    END IF ! iPcl
   
    iPcl = iPcl + 1

! ==============================================================================  
!   Set random particle data
! ==============================================================================  
     
    rn   = Rand1Uniform(pRegion%randData) 
    xLoc = xMin + xDel*rn

    rn   = Rand1Uniform(pRegion%randData) 
! DEBUG: Manoj 2012-03-21: more probability near centerline
!    rn   = rn-0.5_RFREAL
!    rn   = SIGN(1.0_RFREAL,rn)*rn*rn
!    rn   = rn+0.5_RFREAL
!    rn   = rn + 0.1_RFREAL*SIN(rn*2.0_RFREAL*global%pi)
! END DEBUG
    yloc = yMin + yDel*rn

    rn   = Rand1Uniform(pRegion%randData) 
! DEBUG: Manoj 2012-03-21: more probability near centerline
!    rn   = rn + 0.1_RFREAL*SIN(rn*2.0_RFREAL*global%pi)
! END DEBUG
    zLoc = zMin + zDel*rn

! TEMPORARY
!   Allows to create a cylindrical-configuration
!
!    IF ( yLoc*yLoc + zLoc*zLoc > 0.250_RFREAL ) THEN 
!      iPcl = iPcl - 1
!      infLoopCntr = infLoopCntr - 1
!      
!      CYCLE
!     END IF ! yLoc
! END TEMPORARY

    rn     = Rand1Uniform(pRegion%randData) 
    diam   = diamMin + diamDel*rn

    rn     = Rand1Uniform(pRegion%randData) 
    temp   = tempMin  + tempDel*rn

    rn     = Rand1Uniform(pRegion%randData) 
    spLoad = spLoadMin + spLoadDel*rn

    rn     = Rand1Uniform(pRegion%randData) 
    u      = uMin + uDel*rn

    rn     = Rand1Uniform(pRegion%randData) 
    v      = vMin + vDel*rn

    rn     = Rand1Uniform(pRegion%randData) 
    w      = wMin + wDel*rn

! TEMPORARY: Manoj
! ==============================================================================  
! Modifying due to axi-symmetry
! ==============================================================================  

  IF ( pRegion%mixtInput%axiFlag .EQV. .TRUE. ) THEN
    spLoad = spLoad*ABS(yLoc)
  END IF ! pRegion%mixtInput%axiFlag
! END TEMPORARY

! Subbu - Scale superloading for cyldet case
! Redefine xLoc and yLoc to correspond to a cylindrical geometry - Blastwave
! problem -->  xMin(rMin) <=> rMin(rMax); yMin(yMax) <=> thetaMin(thetaMax)
! xLoc identical to rLoc and yLoc indential to thetaLoc
! Define XPOS = rLocCOS(thetaLoc) and YPOS = rLocSIN(thetaLoc) 

  rad  = yLoc
  the  = xLoc
  xLoc = rad * SIN(the)
  yLoc = rad * COS(the)
  scal = 0.5_RFREAL*(yMin-yMax)
  spLoad = spLoad*ABS(rad/scal)
! Subbu - End Scale superloading for cyldet case


! ==============================================================================  
!   Store data
! ==============================================================================

    pCv(CV_PLAG_XPOS,iPcl) = xLoc
    pCv(CV_PLAG_YPOS,iPcl) = yLoc
    pCv(CV_PLAG_ZPOS,iPcl) = zLoc

    DO iCont = 1,pRegion%plagInput%nCont
      massRatio = pIniComp(iCont)*massFluxRatioSumR
      pCv(pCvPlagMass(iCont),iPcl) = pDens(iCont)*massRatio*global%pi/6.0_RFREAL &
                                   * diam**3
    END DO ! iCont

    heatCapSum = SUM(pCv(pCvPlagMass(:),iPcl)*pSpcHeat(:))
    massSum    = SUM(pCv(pCvPlagMass(:),iPcl))
    
    pCv(CV_PLAG_XMOM,iPcl) = massSum*u
    pCv(CV_PLAG_YMOM,iPcl) = massSum*v 
    pCv(CV_PLAG_ZMOM,iPcl) = massSum*w 
    pCv(CV_PLAG_ENER,iPcl) = heatCapSum*temp
    
    pArv(ARV_PLAG_SPLOAD,iPcl) = spLoad
    pAiv(AIV_PLAG_PIDINI,iPcl) = iPcl   
    pAiv(AIV_PLAG_ICELLS,iPcl) = CRAZY_VALUE_INT ! Value used in initialization
    pAiv(AIV_PLAG_REGINI,iPcl) = CRAZY_VALUE_INT ! Value used in initialization

! ==============================================================================  
!   Write some info
! ==============================================================================  

    IF ( iPcl/REAL(nPcls,KIND=RFREAL) > nPclsFrac ) THEN 
      nPclsFrac = nPclsFrac + 0.1_RFREAL 

      IF ( global%verbLevel > VERBOSE_LOW ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I10,1X,A)') SOLVER_NAME,'Generated',iPcl, & 
                                             'particles.'
      END IF ! global%verbLevel 
    END IF ! iPcl

! ------------------------------------------------------------------------------
!   Check for infinite loop
! ------------------------------------------------------------------------------        

    IF ( infLoopCntr >= infLoopCntrMax ) THEN 
      CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
    END IF ! infLoopCntr
  END DO infLoop

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Initializing particle solution from random state done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_InitSolutionRandom

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_InitSolutionRandom.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.5  2008/12/06 08:43:52  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/08/15 13:36:00  haselbac
! Added check for infinite loop counter to avoid wrap-around
!
! Revision 1.2  2007/05/16 22:33:45  fnajjar
! Modified to align with new bc datastructure
!
! Revision 1.1  2007/04/09 18:50:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.11  2007/03/27 00:22:31  haselbac
! Some simplifications
!
! Revision 1.10  2007/03/08 15:04:19  fnajjar
! Fixed bug for massFluxRatioSum being zero and avoiding division by zero
!
! Revision 1.9  2006/10/26 15:00:14  fnajjar
! Added computation of particle momentum based on random velocity field
!
! Revision 1.8  2006/09/18 20:34:26  fnajjar
! Increased FORMAT output to I10 for particle size
!
! Revision 1.7  2006/05/05 18:39:46  haselbac
! Changed logic so serial grid no longer needed
!
! Revision 1.4  2004/11/04 16:29:09  fnajjar
! Deleted call to PLAG_SetMaxDimensions and nPcls definition
!
! Revision 1.3  2004/10/11 22:10:54  haselbac
! Bug fix and cosmetics
!
! Revision 1.2  2004/10/11 19:38:48  fnajjar
! Renamed ininPlag to nPclsIni to follow naming convention
!
! Revision 1.1  2004/10/10 20:06:35  fnajjar
! Initial import
!
! ******************************************************************************

