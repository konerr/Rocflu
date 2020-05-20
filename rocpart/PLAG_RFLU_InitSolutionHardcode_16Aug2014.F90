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
! $Id: PLAG_RFLU_InitSolutionHardcode_16Aug2014.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_InitSolutionHardcode(pRegion)

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
                  massSum,massSumR,nPclsFrac,rn,spLoadDel,spLoadMax,spLoadMin,spLoad, &
                  tempDel,tempMax,tempMin,temp,uDel,uMax,uMin,u,vDel,vMin, &
                  vMax,v,wDel,wMin,wMax,w,xDel,xLoc,xMax,xMin,yDel,yLoc,yMax, &
                  yMin,zDel,zLoc,zMax,zMin
  REAL(RFREAL), POINTER, DIMENSION(:) :: pDens,pIniComp,pSpcHeat
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv,pCv
  TYPE(t_global), POINTER :: global
  TYPE(t_grid),   POINTER :: pGrid

! TEMPORARY: Manoj: Hardcoded particle distribution
  INTEGER :: iIndex,jIndex,nPclsX,nPclsY,nPclsUniform
! END TEMPORARY

! Subbu - Vars Pcl Init
  REAL(RFREAL) :: BaseVol,delr,rad,scal,the
  INTEGER :: icg,kIndex,layer,nPclsZ,nCellsCharge,nCellsSecI,rMinLayer,rMaxLayer
! Subbu - End Vars Pcl Init
! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: PLAG_RFLU_InitSolutionHardcode_16Aug2014.F90,v $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_InitSolutionHardcode', &
                        __FILE__)
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                          'Initializing particle solution from hardcode...'
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

  ! Subbu - If 2d run then set zMin and zMax to z-centroid location of the cells
    IF (pRegion%mixtInput%dimens == 2) THEN
      icg = 1 ! Could be any value between 1 and pGrid%nCells
      zMin = pGrid%cofg(ZCOORD,icg)
      zMax = pGrid%cofg(ZCOORD,icg)
    END IF
  ! Subbu - End set zMin and zMax for 2d

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

! TEMPORARY: Manoj: Hardcoded particle distribution

  WRITE(*,*) "========================================="
  WRITE(*,*) "nPcls        =",nPcls

  ! TEMPORARY: Manoj
  !WRITE(*,*) "Enter nPclsX:"
  !READ(*,*) nPclsX
  !WRITE(*,*) "Enter nPclsY:"
  !READ(*,*) nPclsY
  !WRITE(*,*) "Enter nPclsZ:"
  !READ(*,*) nPclsZ

  nPclsX = pRegion%mixtInput%prepIntVal14 
  nPclsY = pRegion%mixtInput%prepIntVal15 
  nPclsZ = pRegion%mixtInput%prepIntVal16 
  nPclsUniform = nPclsX*nPclsY*nPclsZ

  WRITE(*,*) "nPclsX       =",nPclsX
  WRITE(*,*) "nPclsY       =",nPclsY
  WRITE(*,*) "nPclsZ       =",nPclsZ
  WRITE(*,*) "nPclsX.nPclsY.nPclsZ=",nPclsUniform
  WRITE(*,*) "nPcls        =",nPcls
 
  IF (nPcls /= nPclsUniform) THEN
    WRITE(*,*)'ERROR-Pcl Initialization : nPcls do not match nPclsUniform'
    WRITE(*,*)'Error in file at line number:',__FILE__,__LINE__
    STOP
  END IF
! END TEMPORARY: Manoj

  !pGrid%Jmax = pGrid%nCells/pRegion%patches(3)%nBQuads
  !pGrid%Kmax = pGrid%nCells/pRegion%patches(1)%nBQuads
  !pGrid%Imax = pGrid%nCells/(pGrid%Jmax * pGrid%Kmax)
  
  pGrid%Jmax = pGrid%nCells/pRegion%patches(2)%nBQuads
  pGrid%Imax = pRegion%patches(1)%nBQuads/pGrid%Jmax
  nCellsCharge = (pGrid%Imax/4)*(pGrid%Imax/4) ! No of cells within charge for every z-layer
  pGrid%Kmax = (pRegion%patches(2)%nBQuads - nCellsCharge)/ pGrid%Imax
  WRITE(*,*) 'Imax,Jmax,Kmax=',pGrid%Imax,pGrid%Jmax,pGrid%Kmax

! Subbu - Recast pcl initialization
  
  icg = 0
  DO layer = 1,pGrid%Kmax
    !icg = (layer-1)*(pGrid%Imax*pGrid%Jmax) + 1
    icg = (layer-1)*pGrid%Imax + 1
    rad = SQRT(pGrid%cofg(XCOORD,icg)**2.0_RFREAL + pGrid%cofg(YCOORD,icg)**2.0_RFREAL)
    IF (ABS(rad - yMin) .LE. 1.0E-14) THEN
      BaseVol = pGrid%vol(icg)
      !nCellsSecI = (layer-1)*(pGrid%Imax*pGrid%Jmax) 
      nCellsSecI = (layer-1)*pGrid%Imax 
      rMinLayer = layer
    END IF
    IF (ABS(rad - yMax) .LE. 1.0E-14) THEN
      rMaxLayer = layer
      EXIT
    END IF
  END DO
  delr = (yMax-yMin)/(rMaxLayer-rMinLayer)
  WRITE(*,*) 'rMinLayer,rMaxLayer,nCellsI=',rMinLayer,rMaxLayer,nCellsSecI
 
  DO iIndex = 1,nPclsX
    DO jIndex = 1,nPclsY
      DO kIndex = 1,nPclsZ
        iPcl = iPcl + 1
        xLoc = xMin + (xDel/(nPclsX-0)) * (iIndex-1)
        yLoc = yMin + (yDel/(nPclsY-1)) * (jIndex-1)
        IF (nPclsZ > 1) THEN
          zLoc = zMin + (zDel/(nPclsZ-1)) * (kIndex-1)
        ELSE
          zLoc = zMin
        END IF

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

! Subbu - Redefine x,y as the,rad
! Redefine xLoc and yLoc to correspond to a cylindrical geometry - Blastwave
! problem -->  xMin(rMin) <=> rMin(rMax); yMin(yMax) <=> thetaMin(thetaMax)
! xLoc identical to rLoc and yLoc indential to thetaLoc
! Define XPOS = rLocCOS(thetaLoc) and YPOS = rLocSIN(thetaLoc) 
  
       rad  = yLoc
       the  = xLoc
       xLoc = rad * DCOS(the)
       yLoc = rad * DSIN(the)
       layer = CEILING((rad-yMin+0.5*delr)/delr)
       !icg = nCellsSecI + (layer-1)*(pGrid%Imax*pGrid%Jmax) + 1
       icg = nCellsSecI + (layer-1)*pGrid%Imax + 1
       scal = pGrid%vol(icg)/BaseVol
       spLoad = spLoad*scal
       !WRITE(*,*) 'iPcl,r,t=',iPcl,rad,the*180/global%pi
! Subbu - End redefine x,y as the,rad

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
       massSumR = 1.0_RFREAL/massSum
     
       pCv(CV_PLAG_XMOM,iPcl) = massSum*u
       pCv(CV_PLAG_YMOM,iPcl) = massSum*v 
       pCv(CV_PLAG_ZMOM,iPcl) = massSum*w 
       pCv(CV_PLAG_ENER,iPcl) = heatCapSum*temp + & 
                                massSum*0.5_RFREAL*((massSumR*pCv(CV_PLAG_XMOM,iPcl))**2.0_RFREAL + &
                                            (massSumR*pCv(CV_PLAG_YMOM,iPcl))**2.0_RFREAL + &
                                            (massSumR*pCv(CV_PLAG_ZMOM,iPcl))**2.0_RFREAL)

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

      END DO
    END DO
  END DO
! Subbu - End Recast pcl initialization

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Initializing particle solution from hardcode done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_InitSolutionHardcode

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_InitSolutionHardcode_16Aug2014.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
! ******************************************************************************

