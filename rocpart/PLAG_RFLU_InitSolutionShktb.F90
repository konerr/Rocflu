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
! Saptarshi
!
! Purpose: Initialize particle solution for Shock Tube cases
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
!   nPclsSumReg Sum of particles in regions 1 to n-1, n being the current region
!
! Output: None.
!
! ******************************************************************************
!
! $Id: PLAG_RFLU_InitSolutionShktb.F90,v 1.1 2015/08/12 03:55:40 brollin Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_InitSolutionShktb(pRegion,nPclsSumReg)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt_input
  USE ModMPI
  USE ModPartLag, ONLY: t_plag
  USE ModParameters  
  USE ModRandom, ONLY: Rand1Uniform,Rand1Normal
  
  USE PLAG_ModParameters    

  USE RFLU_ModFaceList, ONLY: RFLU_CreateCell2FaceList, &
                              RFLU_BuildCell2FaceList, &
                              RFLU_DestroyCell2FaceList
  USE RFLU_ModInCellTest, ONLY: RFLU_ICT_TestInCell
     
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

  INTERFACE
    SUBROUTINE CurtainProfile(xs,xmid,phi)
      INTEGER, PARAMETER :: RFREAL = selected_real_kind(15,307)
      REAL(RFREAL), INTENT(IN) :: xs,xmid
      REAL(RFREAL), INTENT(OUT) :: phi
    END SUBROUTINE CurtainProfile
  END INTERFACE

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion
  INTEGER :: nPclsSumReg

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,icg,icl,iCont,iPcl,loopCounter,m,n,nPclsBeg,nPclsEnd
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv
  INTEGER, POINTER, DIMENSION(:) :: pCvPlagMass
  LOGICAL :: foundFlag
  REAL(RFREAL) :: dia,gaussAmp,heatCapSum,massRatio,massFluxRatioSum, &
                  massFluxRatioSumR,massFluxRatioLimit,massSum,massSumR,meanDia, &
                  meanVfrac,perturb,rad,rand,rMinCell,rMaxCell,radtmp,the,spLoad, &
                  tMinCell,thetmp,tMaxCell,tol,T,u,v,vFrac, &
                  volPcl,volPclsSum,w,xLoc,yLoc,z,zLoc,zMinCell,&
                  zMaxCell 
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv,pCv
  REAL(RFREAL), POINTER, DIMENSION(:) :: pDens,pIniComp,pSpcHeat
  TYPE(t_global), POINTER :: global
  TYPE(t_grid),   POINTER :: pGrid
  TYPE(t_plag),   POINTER :: pPlag
  REAL(RFREAL) :: x,y,xMinCell,xMaxCell,yMinCell,yMaxCell,dev
  REAL(RFREAL) :: xMin,xMax,xmid,phi
!=============================================================================  
! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: PLAG_RFLU_InitSolutionShktb.F90,v $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_InitSolutionShktb', &
                        __FILE__)
 
  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
       WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal          
       WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                          'Initializing particle solution...' 
                           
  END IF ! global%verbLevel
  
! ******************************************************************************
! Set pointers and values
! ******************************************************************************

  pAiv  => pRegion%plag%aiv
  pArv  => pRegion%plag%arv
  pCv   => pRegion%plag%cv
  pCvPlagMass => pRegion%plag%cvPlagMass
  pDens => pRegion%plagInput%dens
  pGrid => pRegion%grid
  pIniComp => pRegion%plagInput%iniComp
  pPlag => pRegion%plag
  pSpcHeat => pRegion%plagInput%spht

! ******************************************************************************
! Build cell-to-face list (needed for in-cell test and getting number of cells)
! ******************************************************************************

  CALL RFLU_CreateCell2FaceList(pRegion)
  CALL RFLU_BuildCell2FaceList(pRegion)

! ==============================================================================  
! Set inverse of massFluxRatioSum to avoid division by zero 
! ==============================================================================  

  massFluxRatioLimit = 1.0E-10_RFREAL
  massFluxRatioSum = SUM(pIniComp)
  IF ( massFluxRatioSum > massFluxRatioLimit ) THEN
    massFluxRatioSumR = 1.0_RFREAL/massFluxRatioSum
  ELSE
    massFluxRatioSumR = 1.0_RFREAL
  END IF ! massFluxRatioSum 

! ******************************************************************************  
! Read initial pcl velocities,temperature,diameter and vFrac from user input
! ******************************************************************************

  u = pRegion%plagInput%iniRandUMin
  v = pRegion%plagInput%iniRandVMin
  w = pRegion%plagInput%iniRandWMin
  T = pRegion%plagInput%iniRandTempMin
  xMin = pRegion%plagInput%iniRandXMin
  xMax = pRegion%plagInput%iniRandXMax
  
  meanDia = pRegion%plagInput%iniRandDiamMin
  meanVfrac = pRegion%plagInput%iniRandSpLoadMin

  ! Temporarily set pPlag%nPcls = 0 
  ! It will be recomputed as we place particles one by one in this routine

  pPlag%nPcls = 0
  tol = 1.0E-14_RFREAL
  xmid = 0.5_RFREAL*(xMin+xMax)
 
! ******************************************************************************
! Sum up the number of particles in the all other regions except the current one
! to assign proper global initial PCL_ID.
! ******************************************************************************

  DO icl = 1, pGrid%nCells
    IF (INT(pGrid%nPclsPerCell(icl)) .GT. 0) THEN

      x = pGrid%cofg(XCOORD,icl)
      y = pGrid%cofg(YCOORD,icl)
      z = pGrid%cofg(ZCOORD,icl)

! ******************************************************************************
! Compute the min and max r,the,z coordiantes for each cell
! ******************************************************************************
      ! pGrid%dx,pGrid%dy,pGrid%dz computed in 
      ! libflu/RFLU_ComputeGridSpacingShktb.F90       

      xMinCell = x - 0.5_RFREAL*pGrid%dx
      xMaxCell = x + 0.5_RFREAL*pGrid%dx
      yMinCell = y - 0.5_RFREAL*pGrid%dy
      yMaxCell = y + 0.5_RFREAL*pGrid%dy
      zMinCell = z   - 0.5_RFREAL*pGrid%dz
      zMaxCell = z   + 0.5_RFREAL*pGrid%dz

      volPclsSum = 0.0_RFREAL
      nPclsBeg = pPlag%nPcls + 1
      nPclsEnd = nPclsBeg + INT(pGrid%nPclsPerCell(icl)) - 1
! ******************************************************************************
! Use random number generator to generate a random location for each particle
! within a cell
! ******************************************************************************

      loopCounter = 0
      iPcl = nPclsBeg
      DO WHILE (iPcl .LE. nPclsEnd)
        xLoc = xMinCell &
             + Rand1Uniform(pRegion%randData) * (xMaxCell - xMinCell)
        yLoc = yMinCell &
             + Rand1Uniform(pRegion%randData) * (yMaxCell - yMinCell)
        IF (pRegion%mixtInput%dimens == 2) THEN
          zLoc = z
        ELSE
          zLoc = zMinCell &
               + Rand1Uniform(pRegion%randData) * (zMaxCell - zMinCell)
        END IF

        pCv(CV_PLAG_XPOS,iPcl) = xLoc
        pCv(CV_PLAG_YPOS,iPcl) = yLoc
        pCv(CV_PLAG_ZPOS,iPcl) = zLoc

! ******************************************************************************
! Check if the particle is placed within the cell of interest, if not generate
! a random number until particle (x,y,z) falls inside cell
! The nth particle's initial cell location and regID is also initialized
! ******************************************************************************

        icg = pGrid%hex2CellGlob(icl)
        foundFlag = .FALSE.
        IF ( RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,icg) .EQV. .TRUE. ) THEN
          pPlag%aiv(AIV_PLAG_PIDINI,iPcl) = iPcl + nPclsSumReg
          pAiv(AIV_PLAG_REGINI,iPcl) = pRegion%iRegionGlobal
          pAiv(AIV_PLAG_ICELLS,iPcl) = icg
          foundFlag = .TRUE.

          dia = meanDia !+ Rand1Normal(meanDia,stdDev) ! Add stdDev later
          volPcl = global%pi*dia**3/6.0_RFREAL
! ******************************************************************************
! Compute mass of each particle and initialize each particle with x,y and z
! momentum in addition to total energy
! ******************************************************************************

          DO iCont = 1,pRegion%plagInput%nCont
            massRatio = pIniComp(iCont)*massFluxRatioSumR
            pCv(pCvPlagMass(iCont),iPcl) = pDens(iCont)*massRatio*volPcl
          END DO ! iCont
          heatCapSum = SUM(pCv(pCvPlagMass(:),iPcl)*pSpcHeat(:))
          massSum    = SUM(pCv(pCvPlagMass(:),iPcl))
          massSumR = 1.0_RFREAL/massSum

          pCv(CV_PLAG_XMOM,iPcl) = massSum*u
          pCv(CV_PLAG_YMOM,iPcl) = massSum*v
          pCv(CV_PLAG_ZMOM,iPcl) = massSum*w
          pCv(CV_PLAG_ENER,iPcl) = heatCapSum*T + &
                                   massSum*0.5_RFREAL* &
                                      ((massSumR*pCv(CV_PLAG_XMOM,iPcl))**2+ &
                                       (massSumR*pCv(CV_PLAG_YMOM,iPcl))**2+ &
                                       (massSumR*pCv(CV_PLAG_ZMOM,iPcl))**2)
          volPclsSum = volPclsSum + volPcl

          iPcl = iPcl + 1
          pPlag%nPcls = pPlag%nPcls + 1

        ELSE   ! Increment loop counter if pcl falls outside cell

          loopCounter = loopCounter + 1
          IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN ! Prevent infinite loop
            CALL ErrorStop(global,ERR_INFINITE_LOOP ,__LINE__)
          END IF ! loopCounter

        END IF ! RFLU_ICT_TestInCell
      END DO ! iPcl

! ******************************************************************************
! Compute superparticle loading : spLoad for all particles within a given cell
! is taken to be the same
! ******************************************************************************
     IF (pRegion%mixtInput%prepIntVal2 == 1) THEN
      dev = pRegion%mixtInput%prepRealVal12 ! Wave no in theta
      vFrac = Rand1Normal(meanVFrac,dev,pRegion%randData) ! Add stdDev late
     ELSE IF (pRegion%mixtInput%prepIntVal2 == 2) THEN
      IF ((x .GE. xMin) .AND. (x .LE. xMax)) THEN
!       write(198,*) meanVfrac
       CALL CurtainProfile(x*1000.0_RFREAL,xmid*1000.0_RFREAL,phi)
       vFrac = meanVfrac*phi
!       write(199,*) vFrac
      END IF ! x .GE. xMin
     ELSE
!      vFrac = 0.25_RFREAL*( (DTANH(k*( (x-b)+1.0_RFREAL))+1.0_RFREAL)* &
!                            (DTANH(k*(-(x-b)+1.0_RFREAL))+1.0_RFREAL) )
      vFrac = meanVFrac !*vFrac !Rand1Normal(meanVFrac,dev,pRegion%randData) ! Add stdDev late
     END IF
!      gaussAmp = pRegion%mixtInput%prepRealVal5
!      m = pRegion%mixtInput%prepIntVal1 ! Wave no in theta
!      n = pRegion%mixtInput%prepIntVal2 ! Wave number in z
!      perturb = gaussAmp*DCOS(m*the)
!      vFrac = meanVfrac*(1.0_RFREAL+perturb)

      spLoad = vFrac * pGrid%vol(icl)/volPclsSum
      DO iPcl = nPclsBeg,nPclsEnd
        pArv(ARV_PLAG_SPLOAD,iPcl) = spLoad
      END DO
    END IF ! nPclsPerCell > 0
  END DO ! icl

! ******************************************************************************
! Destroy cell-to-face list
! ******************************************************************************

  CALL RFLU_DestroyCell2FaceList(pRegion)

! ******************************************************************************
! Destroy memory for number of particles per cell (only if a region contains
! pcls)
! NOTE memory deallocated in PLAG_RFLU_ComputeCellsContiningPcls if no particles
! present
! ******************************************************************************

  IF (pGrid%initPclPresent .EQV. .TRUE.) THEN
    DEALLOCATE(pGrid%nPclsPerCell,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'nPclsPerCell')
    END IF ! global%error
  END IF

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN
       WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
            'Initializing particle solution done.'
  END IF ! global%verbLevel

!*******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_InitSolutionShktb 

SUBROUTINE CurtainProfile(xs,xmid,phi)

  IMPLICIT NONE

  INTEGER, PARAMETER :: RFREAL = selected_real_kind(15,307)

  REAL(RFREAL), INTENT(IN) :: xs,xmid
  REAL(RFREAL), INTENT(OUT) :: phi

  INTEGER :: i
  REAL(RFREAL) :: centers(14),theta(14),betas(14),sqrdDists(14),phis(14), &
                  iden(14)
  REAL(RFREAL) :: sigma,beta

  DATA centers / -2.0, -1.7, -1.3, -1.1, -0.5, -0.3, -0.05, 0.05, 0.3, 0.5, &
                  1.1, 1.3, 1.7, 2.0 /
  DATA theta / -695.00, 2440.96, -8530.24, 10307.40, -26162.70, 45596.39, &
               -22973.34, -22973.34, 45596.39, -26162.70, 10307.40, -8530.24, &
                2440.96, -695.00 /
  DO i = 1,14
    sqrdDists(i) = 0.0_RFREAL
    phis(i) = 0.0_RFREAL
    iden(i) = 1.0_RFREAL
  END DO

  sigma = 0.65_RFREAL
  beta = 1.0_RFREAL/(2.0_RFREAL* sigma*sigma)

  DO i = 1,14
    !sqrdDists(i) = (centers(i)-(xs-xmid)*iden)**2.0
    sqrdDists(i) = (centers(i)-(xs-xmid))**2 !.0_RFREAL
    !sqrdDists(i) = (centers(i)-xs)**2.0_RFREAL
    phis(i) = EXP(-beta * sqrdDists(i))
  END DO

  phis = phis/SUM(phis)
  phi = DOT_PRODUCT(theta,phis)

END SUBROUTINE CurtainProfile

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_InitSolutionShktb.F90,v $
! Revision 1.1  2015/08/12 03:55:40  brollin
! New Subroutine to Initialize Shock Tube related problems.
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
!******************************************************************************

