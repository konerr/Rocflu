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
! Purpose: Initialize particle solution in a region read from an ASCII file
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
! $Id: PLAG_RFLU_InitSolutionFile.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_InitSolutionFile(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModParameters
  
  !USE ModRandom, ONLY: Rand1Uniform
  
  USE PLAG_ModParameters    
  USE ModBuildFileNames, ONLY: BuildFileNamePlain
   
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

  CHARACTER(CHRLEN) :: RCSIdentString,iFileName
  INTEGER :: errorFlag,icg,iCont,iFile,iPcl,nPcls
  INTEGER, POINTER, DIMENSION(:) :: pCvPlagMass
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv
  REAL(RFREAL) :: heatCapSum,massRatio,massFluxRatioSum,massFluxRatioSumR, &
                  massFluxRatioLimit,massSum,massSumR,nPclsFrac
  REAL(RFREAL), POINTER, DIMENSION(:) :: pDens,pIniComp,pSpcHeat
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv,pCv
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: diam
  TYPE(t_global), POINTER :: global
  TYPE(t_grid),   POINTER :: pGrid


! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: PLAG_RFLU_InitSolutionFile.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_InitSolutionFile', &
                        __FILE__)
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                          'Initializing particle solution - Reading from ASCII &
                           file...'
  END IF ! global%verbLevel
 
! ******************************************************************************
!   Start, open file and read title
! ******************************************************************************

  CALL BuildFileNamePlain(global,FILEDEST_INDIR,'_ptcl.dat',iFileName)

  iFile = 10 
  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
  END IF

  READ(iFile,'(I9.1)') nPcls

! == Check if nPcls is consistent with pRegion%plag%nPcls read earlier == !
  IF (nPcls /= pRegion%plag%nPcls) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
  END IF
 
! Allocate memory for dia
  ALLOCATE(diam(nPcls),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'diam')
  END IF 

! ******************************************************************************
! Set pointers and values
! ******************************************************************************

  pGrid => pRegion%grid

  pRegion%plag%nextIdNumber = pRegion%plag%nPcls
  
  nPclsFrac   = 0.1_RFREAL    
  iPcl        = 0
  
  massFluxRatioLimit = 1.0E-10_RFREAL

! ******************************************************************************
! Set pointers and values for PLAG infrastructure
! ******************************************************************************

  pDens    => pRegion%plagInput%dens
  pSpcHeat => pRegion%plagInput%spht
  pIniComp => pRegion%plagInput%iniComp
  
  pAiv => pRegion%plag%aiv
  pArv => pRegion%plag%arv
  pCvPlagMass => pRegion%plag%cvPlagMass
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


! ******************************************************************************
! Read the initial particle location,dia,temp,spl and velocities from file 
! ******************************************************************************

  DO iPcl = 1,nPcls
      READ(iFile,'(I9.1,9(1X,E16.9))') pAiv(AIV_PLAG_PIDINI,iPcl), & 
                                       pCv(CV_PLAG_XPOS:CV_PLAG_ZPOS,iPcl),diam(iPcl), &
                                       pCv(CV_PLAG_ENER,iPcl),pArv(ARV_PLAG_SPLOAD,iPcl), &
                                       pCv(CV_PLAG_XMOM:CV_PLAG_ZMOM,iPcl)
  END DO

  CLOSE(iFile)

  ! If 2d run then set pcl z-loc to z-centroid location of the cells
    IF (pRegion%mixtInput%dimens == 2) THEN
      icg = 1 ! Could be any value between 1 and pGrid%nCells
      pCv(CV_PLAG_ZPOS,1:nPcls) = pGrid%cofg(ZCOORD,icg)
    END IF
  ! End set pcl z-loc for 2d

  DO iPcl = 1,nPcls

! ==============================================================================  
! Modifying due to axi-symmetry
! ==============================================================================  

    IF ( pRegion%mixtInput%axiFlag .EQV. .TRUE. ) THEN
      pArv(ARV_PLAG_SPLOAD,iPcl) = pArv(ARV_PLAG_SPLOAD,iPcl)*ABS(pCv(CV_PLAG_YPOS,iPcl))
    END IF ! pRegion%mixtInput%axiFlag

! ==============================================================================  
!   Store data
! ==============================================================================

    DO iCont = 1,pRegion%plagInput%nCont
      massRatio = pIniComp(iCont)*massFluxRatioSumR
      pCv(pCvPlagMass(iCont),iPcl) = pDens(iCont)*massRatio*global%pi/6.0_RFREAL &
                                   * diam(iPcl)**3
    END DO ! iCont

    heatCapSum = SUM(pCv(pCvPlagMass(:),iPcl)*pSpcHeat(:))
    massSum    = SUM(pCv(pCvPlagMass(:),iPcl))
    massSumR = 1.0_RFREAL/massSum

    pCv(CV_PLAG_XMOM,iPcl) = massSum*pCv(CV_PLAG_XMOM,iPcl)
    pCv(CV_PLAG_YMOM,iPcl) = massSum*pCv(CV_PLAG_YMOM,iPcl)
    pCv(CV_PLAG_ZMOM,iPcl) = massSum*pCv(CV_PLAG_ZMOM,iPcl) 
    pCv(CV_PLAG_ENER,iPcl) = heatCapSum*pCv(CV_PLAG_ENER,iPcl) + &
                             massSum*0.5_RFREAL*( &
                             (massSumR*pCv(CV_PLAG_XMOM,iPcl))**2.0_RFREAL + &
                             (massSumR*pCv(CV_PLAG_YMOM,iPcl))**2.0_RFREAL + &
                             (massSumR*pCv(CV_PLAG_ZMOM,iPcl))**2.0_RFREAL)
    
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

! Deallocate memory for dia
  DEALLOCATE(diam,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'diam')
  END IF

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Initializing particle solution - Reading from ASCII file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_InitSolutionFile

! ******************************************************************************
!
! RCS Revision history:

! Subbu
! Read the particle coordinates, dia, superparticle loading, temperature 
! and initial velocity from an ASCII file
!
! ******************************************************************************

