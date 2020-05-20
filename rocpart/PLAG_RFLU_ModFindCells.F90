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
! Purpose: Suite of routines for particle tracking on Eulerian grid.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: PLAG_RFLU_ModFindCells.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE PLAG_RFLU_ModFindCells
  
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModDataTypes
  USE ModGlobal,     ONLY: t_global
  USE ModGrid,       ONLY: t_grid
  USE ModPartLag,    ONLY: t_plag, t_plag_input
  USE ModBndPatch,   ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModBorder,     ONLY: t_border
  USE ModError
  USE ModMPI
  USE PLAG_ModParameters
    
  USE ModTools, ONLY: FloatEqual
  
  USE RFLU_ModCellFaceEdgeInfo, ONLY: RFLU_GetGlobalCellKind, & 
                                      RFLU_GetFaceKind  
  USE RFLU_ModGeometryTools
  USE RFLU_ModInCellTest, ONLY: RFLU_ICT_TestInCell, & 
                                RFLU_ICT_TestInCellFancy, & 
                                RFLU_ICT_TestInCellLohner
  USE RFLU_ModOctree
  USE RFLU_ModRelatedPatches, ONLY: RFLU_RELP_TransformVector

  USE PLAG_ModInterfaces, ONLY: PLAG_ReflectParticleData 
  USE PLAG_ModSurfStats, ONLY: PLAG_GatherSurfStats

  USE RFLU_ModMPI, ONLY: RFLU_MPI_RecreateBufferIPclSend  

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: PLAG_RFLU_ComputeDistTot,    &
            PLAG_RFLU_FindCellsBrute,    &
            PLAG_RFLU_FindCellsBruteMod, &
            PLAG_RFLU_FindCellsHardCode, & ! Subbu - Add Hardcoded Pcl Tracking
            PLAG_RFLU_FindCellsLohner,   & 
            PLAG_RFLU_FindCellsOct,      &
            PLAG_RFLU_FindCellsOctMod,   &
            PLAG_RFLU_FindCellsTrajFast, &
            PLAG_RFLU_FindCellsTrajSafe
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: PLAG_RFLU_ModFindCells.F90,v $ $Revision: 1.1.1.1 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS











! ******************************************************************************
!
! Purpose: Compute total distance travelled.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_ComputeDistTot(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: errorString
  INTEGER ::iPcl
  REAL(RFREAL) :: xLocNew,xLocOld,xTraj,yLocNew,yLocOld,yTraj,zLocNew, &
                  zLocOld,zTraj
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_ComputeDistTot',__FILE__)

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

! ******************************************************************************
! Loop over particles
! ******************************************************************************

  DO iPcl = 1,pPlag%nPcls
 
! ==============================================================================  
!   Compute distance travelled and trajectory
! ==============================================================================  
  
    xLocNew = pPlag%cv(CV_PLAG_XPOS,iPcl)
    yLocNew = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLocNew = pPlag%cv(CV_PLAG_ZPOS,iPcl)  

    xLocOld = pPlag%cvOld(CV_PLAG_XPOS,iPcl)
    ylocOld = pPlag%cvOld(CV_PLAG_YPOS,iPcl)
    zLocOld = pPlag%cvOld(CV_PLAG_ZPOS,iPcl)
        
    xTraj = xLocNew - xLocOld
    yTraj = yLocNew - yLocOld      
    zTraj = zLocNew - zLocOld

    pPlag%arv(ARV_PLAG_DISTOT,iPcl)  = SQRT( xTraj*xTraj + yTraj*yTraj &
                                           + zTraj*zTraj )
          
  END DO ! iPcl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_ComputeDistTot









! ******************************************************************************
!
! Purpose: Determine cells which contain particles using brute force approach.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsBrute(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  LOGICAL :: foundFlag
  CHARACTER(CHRLEN) :: errorString
  INTEGER :: errorFlag,icg,iPcl
  REAL(RFREAL) :: xLoc,yLoc,zLoc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsBrute',__FILE__)

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

! ******************************************************************************
! Loop over particles
! ******************************************************************************

  DO iPcl = 1,pPlag%nPcls  

! ==============================================================================  
!   Get particle position
! ==============================================================================  

    xLoc = pPlag%cv(CV_PLAG_XPOS,iPcl)
    yloc = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLoc = pPlag%cv(CV_PLAG_ZPOS,iPcl)

! ==============================================================================  
!   Loop over actual cells
! ==============================================================================  
  
    foundFlag = .FALSE.
  
    cellLoop: DO icg = 1,pGrid%nCells         
      IF ( RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,icg) .EQV. .TRUE. ) THEN
        pPlag%aiv(AIV_PLAG_REGINI,iPcl) = pRegion%iRegionGlobal
        pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg

        foundFlag = .TRUE.

        EXIT cellLoop
      END IF ! RFLU_ICT_TestInCell        
    END DO cellLoop
      
! ==============================================================================  
!   Check whether particle was found
! ==============================================================================  
            
    IF ( foundFlag .EQV. .FALSE. ) THEN
      WRITE(errorString,'(I6)') iPcl 
      CALL ErrorStop(global,ERR_PLAG_PCL_NOT_FOUND,__LINE__,TRIM(errorString))
    END IF ! foundFlag
  END DO ! iPcl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsBrute







! ******************************************************************************
!
! Purpose: Kernel for brute force search.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to region
!   xLoc,yLoc,zLoc      Location
!
! Output:
!   icgOut              Cell containing location
!
! Notes: 
!   1. icgOut is set to CRAZY_VALUE_INT if no cell contained location.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsBruteKernel(pRegion,xLoc,yLoc,zLoc,icgOut)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(OUT) :: icgOut
  REAL(RFREAL), INTENT(IN) :: xLoc,yLoc,zLoc
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  INTEGER :: errorFlag,icg
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsBruteKernel',__FILE__)

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Loop over cells
! ******************************************************************************

  icgOut = CRAZY_VALUE_INT

  cellLoop: DO icg = 1,pGrid%nCells         
    IF ( RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,icg) .EQV. .TRUE. ) THEN
      icgOut = icg

      EXIT cellLoop
    END IF ! RFLU_ICT_TestInCell        
  END DO cellLoop
      
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsBruteKernel







! ******************************************************************************
!
! Purpose: Determine cells which contain particles using modified brute 
!  force approach.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine is suitable to determine the cell containing a particle 
!      after the latters position has been updated. 
!   2. This routine cannot be used in parallel runs and for particles which 
!      leave the domain or bounce from solid walls.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsBruteMod(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: errorString
  INTEGER :: errorFlag,icg,iPcl
  REAL(RFREAL) :: xLoc,yLoc,zLoc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsBruteMod',__FILE__)

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

! ******************************************************************************
! Loop over particles
! ******************************************************************************

  DO iPcl = 1,pPlag%nPcls  

! ==============================================================================  
!   Get particle cell and position
! ==============================================================================  

    icg  = pPlag%aivOld(AIV_PLAG_ICELLS,iPcl) ! NOTE OLD cell index

    xLoc = pPlag%cv(CV_PLAG_XPOS,iPcl) ! NOTE already NEW particle position
    yloc = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLoc = pPlag%cv(CV_PLAG_ZPOS,iPcl)

! ==============================================================================  
!   Check if particle still in same cell 
! ==============================================================================  

    IF ( RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,icg) .EQV. .FALSE. ) THEN    

! ------------------------------------------------------------------------------  
!     Find cell containing particle
! ------------------------------------------------------------------------------  
  
      CALL PLAG_RFLU_FindCellsBruteKernel(pRegion,xLoc,yLoc,zLoc,icg)
  
      IF ( icg /= CRAZY_VALUE_INT ) THEN 
        pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg
      ELSE 

! TEMPORARY
        WRITE(*,*) 'timeCurrent,iPcl,xLocOld,yLocOld,zLocOld,xLoc,yLoc,zLoc = ',&
        global%currentTime,iPcl,pPlag%cvOld(CV_PLAG_XPOS:CV_PLAG_ZPOS,iPcl),&
        xLoc,yLoc,zLoc 
! END TEMPORARY

        WRITE(errorString,'(I6)') iPcl 
        CALL ErrorStop(global,ERR_PLAG_PCL_NOT_FOUND,__LINE__,TRIM(errorString))     
      END IF ! icg
            
! ==============================================================================  
!   Particle still in same cell, copy cell index
! ==============================================================================        
      
    ELSE   
      pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg
    END IF ! RFLU_ICT_TestInCell
  END DO ! iPcl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsBruteMod






! Subbu - Add hardcoded pcl tracking subroutine
! ******************************************************************************
!
! Purpose: Determine cells which contain particles using hardcoded routines
!
! Description: Locate the cell containing particle based on delrad, deltheta
!   and delz info 
!
! Input: 
!  pRegion      Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Currently works only for cyldet case in 2D/3D provided the points
!      are equisapced
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsHardCode(pRegion,iPclBeg,iPclEnd)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  
  INTEGER :: iPclBeg,iPclEnd
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  TYPE(t_border), POINTER :: pBorder
  CHARACTER(CHRLEN) :: errorString
  INTEGER :: flagVirtCol,flagVirtDep,flagVirtRow,iBorder,icg,icl, &
             icgVirt,iCol,iDep,ifl,ifgOut,iLayer,iPatch,iPcl,iRow,iReg, & 
             nRows,nRegs,Sgn,SgnVirt
  REAL(RFREAL) :: nTol,nr,nx,ny,nz,radLoc,theLoc,theLocActCell,theta,xLoc, &
                  xLocOld,yLoc,yLocOld,zLoc,zLocOld
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsHardCode',__FILE__)

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag
  nTol = 1.0E-14_RFREAL

   !radMax = 1.997E-02_RFREAL
! ******************************************************************************
! Loop over particles
! ******************************************************************************
 DO iPcl = iPclBeg,iPclEnd

! ==============================================================================  
!   Get particle position
! ==============================================================================  

    xLocOld = pPlag%cvOld(CV_PLAG_XPOS,iPcl)
    ylocOld = pPlag%cvOld(CV_PLAG_YPOS,iPcl)
    zLocOld = pPlag%cvOld(CV_PLAG_ZPOS,iPcl)

    xLoc = pPlag%cv(CV_PLAG_XPOS,iPcl)
    yloc = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLoc = pPlag%cv(CV_PLAG_ZPOS,iPcl)
     
    radLoc = SQRT(xLoc**2.0_RFREAL + yLoc**2.0_RFREAL)
    theLoc = ATAN2(yLoc,xLoc)

    IF (pGrid%radMaxflag == 1 .AND. radLoc-pGrid%radMaxGlob .GT. nTol ) THEN ! Pcls outside domain
    !IF ( (pGrid%zMinflag == 1 .AND. pGrid%zMinGlob-zLoc .GT. nTol) .OR. &
    !     (pGrid%zMaxflag == 1 .AND. zLoc-pGrid%zMaxGlob .GT. nTol) ) THEN ! Pcls outside domain

        !pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = pPlag%aivOld(AIV_PLAG_ICELLS,iPcl)
        icg = pPlag%aivOld(AIV_PLAG_ICELLS,iPcl)
        icl = pGrid%cellGlob2Loc(2,icg)
        DO ifl = 1,6
          iPatch = pGrid%hex2f(1,ifl,icl)
          IF (iPatch .NE. 0) THEN
            pPatch => pRegion%patches(iPatch)
            nx = pPatch%fn(XCOORD,ifl)
            ny = pPatch%fn(YCOORD,ifl)
            nz = DABS(pPatch%fn(ZCOORD,ifl))
            nr = DSQRT(nx**2.0_RFREAL + ny**2.0_RFREAL)
            IF (DABS(nr - 1.0_RFREAL) .LE. nTol) THEN ! nTol can be reduced and condition for radLoc>radMax
            !IF (DABS(nz - 1.0_RFREAL) .LE. nTol) THEN ! nTol can be reduced
              ifgOut = pGrid%hex2f(2,ifl,icl)
              EXIT
            END IF
          END IF
        END DO
        CALL PLAG_ReflectParticleData(pPatch,pPlag,ifgOut,iPcl,xLocOld,yLocOld, &
                                      zLocOld,xLoc,yLoc,zLoc)
        IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
          CALL PLAG_GatherSurfStats(pRegion,pPlag,pPatch%statsPlag,&
                                    ifgOut,iPcl,theta)
        END IF ! pPatch%plotStatsFlag
        GOTO 200
    END IF

    radLoc = SQRT(xLoc**2.0_RFREAL + yLoc**2.0_RFREAL)
    theLoc = ATAN2(yLoc,xLoc)
    theLocActCell = theLoc

    IF (theLoc .LT. nTol) &
      theLocActCell = theLoc + 2.0_RFREAL*global%pi

    iDep = FLOOR(ABS((zLoc - (pGrid%zMin-0.5_RFREAL*pGrid%dz))/(pGrid%dz)) + 1)
    iRow = FLOOR(ABS((radLoc -(pGrid%radMin-0.5_RFREAL*pGrid%drad))/(pGrid%drad)) + 1)
   
    IF ( iDep .GE. 1 .AND. iDep .LE. pGrid%Jmax ) THEN
      flagVirtDep = 0
      IF ( iRow .GE. 1 .AND. iRow .LE. pGrid%Kmax ) THEN
        flagVirtRow = 0
        !IF ((ABS(theLocActCell) - ABS(pGrid%theLowLim)) .GT. nTol .AND.  &
        !    (ABS(pGrid%theUppLim) - ABS(theLocActCell)) .GT. nTol) THEN ! Pcl lies in actual cell
        IF (theLocActCell .GT. pGrid%theLowLim-1.0E-08 .AND.  &
            theLocActCell .LT. pGrid%theUppLim+1.0E-08) THEN ! Pcl lies in actual cell
            flagVirtCol = 0 
            Sgn = INT(DSIGN(1.0_RFREAL,theLoc))
            IF (Sgn /= pGrid%SgntheMin) & !THEN
              theLoc = theLoc + pGrid%SgntheMin*2.0_RFREAL*global%pi
            iCol = FLOOR(ABS((theLoc -(pGrid%theMin-0.5_RFREAL*pGrid%dthe))/(pGrid%dthe)) + 1)
            icg  = (iDep-1)*(pGrid%Imax*pGrid%Kmax) + (iRow-1)*pGrid%Imax + iCol
            GOTO 100
         ELSE
            flagVirtCol = 1 ! Pcl lies in virtual cell
            !WRITE(*,*) 'iPcl,r,t,z=',iPcl,radLoc,theLoc,zLoc
            !WRITE(*,*) 'iPcl,tAct,tLow,tUpp',iPcl,theLocActCell,pGrid%theLowLim,pGrid%theUppLim
         END IF
      ELSE
        flagVirtRow = 1
      END IF
    ELSE
      flagVirtDep = 1
    END IF

    IF (flagVirtDep == 1) THEN
      iDep = FLOOR(ABS((zLoc - (pGrid%zMinVirtExt-0.5_RFREAL*pGrid%dz))/(pGrid%dz)) + 1)
      iRow = FLOOR(ABS((radLoc -(pGrid%radMinVirtExt-0.5_RFREAL*pGrid%drad))/(pGrid%drad)) + 1)
      iCol = FLOOR(ABS((theLoc -(pGrid%theMinVirtExt-0.5_RFREAL*pGrid%dthe))/(pGrid%dthe)) + 1)
      Sgn = INT(DSIGN(1.0_RFREAL,theLoc))
      IF (Sgn /= pGrid%SgntheMinVirtExt) THEN
        theLoc = theLoc + pGrid%SgntheMinVirtExt*2.0_RFREAL*global%pi
        iCol = FLOOR(ABS((theLoc -(pGrid%theMinVirtExt-0.5_RFREAL*pGrid%dthe))/(pGrid%dthe)) + 1)
        IF (iCol .GT. pGrid%ImaxVirt) & 
          iCol = iCol - global%nthe + pGrid%ImaxVirt
      END IF
      icg  = (iDep-1)*(pGrid%ImaxVirt*pGrid%KmaxVirt) + (iRow-1)*pGrid%ImaxVirt + iCol
    END IF

    ! *******************  Pcl in iRow < 1 or iRow > Kmax ****************************
    IF (flagVirtRow == 1) THEN
      IF (iRow .GT. pGrid%Kmax) THEN
      iCol = FLOOR(ABS((theLoc -(pGrid%theMinVirtExt-0.5_RFREAL*pGrid%dthe))/(pGrid%dthe)) + 1)
      Sgn = INT(DSIGN(1.0_RFREAL,theLoc))
      IF (Sgn /= pGrid%SgntheMinVirtExt) THEN
        theLoc = theLoc + pGrid%SgntheMinVirtExt*2.0_RFREAL*global%pi
        iCol = FLOOR(ABS((theLoc -(pGrid%theMinVirtExt-0.5_RFREAL*pGrid%dthe))/(pGrid%dthe)) + 1)
        IF (iCol .GT. pGrid%ImaxVirt) & 
          iCol = iCol - global%nthe + pGrid%ImaxVirt
      END IF
      icg = (iDep-1)*(pGrid%ImaxVirt*pGrid%KmaxVirt - pGrid%Imax*pGrid%Kmax) + &
            pGrid%nCellsShare + pGrid%Kmax*(pGrid%ImaxVirt-pGrid%Imax) + &
            (iRow-pGrid%Kmax-1)*pGrid%ImaxVirt + iCol
      icg = icg + pGrid%nCells + pGrid%nCellsVirt
      ELSE
      iRow = iRow + 4
      iCol = FLOOR(ABS((theLoc -(pGrid%theMinVirtExt-0.5_RFREAL*pGrid%dthe))/(pGrid%dthe)) + 1)
      Sgn = INT(DSIGN(1.0_RFREAL,theLoc))
      IF (Sgn /= pGrid%SgntheMinVirtExt) THEN
        theLoc = theLoc + pGrid%SgntheMinVirtExt*2.0_RFREAL*global%pi
        iCol = FLOOR(ABS((theLoc -(pGrid%theMinVirtExt-0.5_RFREAL*pGrid%dthe))/(pGrid%dthe)) + 1)
        IF (iCol .GT. pGrid%ImaxVirt) & 
          iCol = iCol - global%nthe + pGrid%ImaxVirt
      END IF
      icg = (iDep-1)*(pGrid%ImaxVirt*pGrid%KmaxVirt - pGrid%Imax*pGrid%Kmax) + &
            (iRow-1)*pGrid%ImaxVirt + iCol
      icg = icg + pGrid%nCells + pGrid%nCellsVirt
      END IF
    END IF

    ! *******************  Pcl in iCol < 1 or iRow > Imax ****************************
    IF (flagVirtCol == 1) THEN
      iCol = FLOOR(ABS((theLoc -(pGrid%theMinVirtInt-0.5_RFREAL*pGrid%dthe))/(pGrid%dthe)) + 1)
      Sgn = INT(DSIGN(1.0_RFREAL,theLoc))
      IF (Sgn /= pGrid%SgntheMinVirtInt) THEN
        theLoc = theLoc + pGrid%SgntheMinVirtInt*2.0_RFREAL*global%pi
        iCol = FLOOR(ABS((theLoc -(pGrid%theMinVirtInt-0.5_RFREAL*pGrid%dthe))/(pGrid%dthe)) + 1)
      END IF
      IF (iCol .GT. 3) THEN
        iCol = iCol - pGrid%Imax
        IF (iCol .GT. 6) & !THEN
          iCol = (iCol+pGrid%Imax) - (global%nthe-pGrid%Imax) + 6
      END IF
      icg = (iDep-1)*(pGrid%ImaxVirt*pGrid%KmaxVirt - pGrid%Imax*pGrid%Kmax) + &
            pGrid%nCellsShare + (iRow-1)*(pGrid%ImaxVirt-pGrid%Imax) + iCol
      icg = icg + pGrid%nCells + pGrid%nCellsVirt
    END IF

    IF (icg .GT. pGrid%nCells) THEN
      !WRITE(*,*) '==After=='
      !WRITE(*,*) 'iPcl,Virtflag-R,C,D',iPcl,flagVirtRow,flagVirtCol,flagVirtDep
      !WRITE(*,*) 'iPcl,iRow,iCol,iDep',iPcl,iRow,iCol,iDep
      !WRITE(*,*) 'iPcl,icg,nCells,GlobReg',iPcl,icg,pGrid%nCells,pRegion%iRegionGlobal
      !WRITE(*,*) 'radLoc,theLoc,zLoc=',radLoc,theLoc,zLoc

      icgVirt = icg - pGrid%nCells
      iBorder = pGrid%vc2border(1,icgVirt)

      ! Subbu - Print out info if iBorder exceeds limits
      IF (iBorder .GT. pGrid%nBorders .OR. iBorder .LE. 0) THEN
        WRITE(*,*) ' ************* ERROR : iBorder is outisde bounds  ***************'
        WRITE(*,*) 'Error in file :',__FILE__
        WRITE(*,*) 'Line # :',__LINE__
        STOP
        !WRITE(*,*) 'iBorder=',iBorder
        !WRITE(*,*) 'r,t,z=',radLoc,theLoc,zLoc
        !WRITE(*,*) 'flags-Dep,Row,Col=',flagVirtDep,flagVirtRow,flagVirtCol
        !WRITE(*,*) 'Lowlim,UppLim,theMinExt,dthe=',pGrid%theLowLim,pGrid%theUppLim,pGrid%theMinVirtExt,pGrid%dthe
        !WRITE(*,*) 'R,C,D=',iRow,iCol,iDep
        !WRITE(*,*) 'icgVirt,Reg,icg=',icgVirt,pRegion%iRegionGlobal,icg
        !WRITE(*,*) 'Sgn,SgntheVirtMinExt=',Sgn,pGrid%SgntheMinVirtExt
        !icg = icg + 1
        !icgVirt = icg - pGrid%nCells
        !iBorder = pGrid%vc2border(1,icgVirt)
        !IF (iBorder .GT. pGrid%nBorders) THEN
        !  icg = icg - 2
        !  icgVirt = icg - pGrid%nCells
        !  iBorder = pGrid%vc2border(1,icgVirt)
        !END IF
      END IF
      ! Subbu - Print out info if iBorder exceeds limits

      pBorder => pGrid%borders(iBorder)
      pBorder%nPclsSend = pBorder%nPclsSend + 1
      IF ( pBorder%nPclsSend > pBorder%nPclsSendMax ) &
         CALL RFLU_MPI_RecreateBufferIPclSend(pRegion,pBorder)

      pBorder%iPclSend(1,pBorder%nPclsSend) = iPcl
      pBorder%iPclSend(2,pBorder%nPclsSend) = pGrid%vc2border(2,icgVirt)
      pPlag%aiv(AIV_PLAG_STATUS,iPcl) = PLAG_STATUS_COMM
    END IF 

100 pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = pRegion%grid%cellGlob2Loc(2,icg)

200 pPlag%cv(CV_PLAG_XPOS,iPcl) = xLoc
    pPlag%cv(CV_PLAG_YPOS,iPcl) = yLoc
    pPlag%cv(CV_PLAG_ZPOS,iPcl) = zLoc
     
    pPlag%cvOld(CV_PLAG_XPOS,iPcl) = xLocOld
    pPlag%cvOld(CV_PLAG_YPOS,iPcl) = yLocOld
    pPlag%cvOld(CV_PLAG_ZPOS,iPcl) = zLocOld 

 END DO ! iPcl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

 END SUBROUTINE PLAG_RFLU_FindCellsHardCode

! Subbu - End subroutine hardcode pcl tracking




! ******************************************************************************
!
! Purpose: Determine cells which contain particles using Lohners approach.
!
! Description: Follow particle path by passing particle to face which failed
!   in-cell test by largest amount.
!
! Input: 
!  pRegion      Pointer to region
!
! Output: None.
!
! Notes: 
!   1. See R. Lohner and J. Ambrosiano, A Vectorized Particle Tracer for 
!      Unstructured Grids, J. Comp. Phys., Vol. 91, pp. 22-31, 1990.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsLohner(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  LOGICAL :: testInCell
  CHARACTER(CHRLEN) :: errorString
  INTEGER :: c1,c1k,c2,c2k,errorFlag,icg,icgOut,ifg,iloc,iPcl,loopCounter
  REAL(RFREAL) :: xLoc,yLoc,zLoc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsLohner',__FILE__)

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

! ******************************************************************************
! Loop over particles
! ******************************************************************************

  DO iPcl = 1,pPlag%nPcls
    loopCounter = 0 ! Reset loop counter 

    xLoc = pPlag%cv(CV_PLAG_XPOS,iPcl)
    yLoc = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLoc = pPlag%cv(CV_PLAG_ZPOS,iPcl)  
                   
    icg = pPlag%aivOld(AIV_PLAG_ICELLS,iPcl)

! ==============================================================================
!   Loop until distance travelled along trajectory consumed 
! ==============================================================================        
        
    infLoop: DO     
      loopCounter = loopCounter + 1  
        
! ------------------------------------------------------------------------------
!     Find appropriate intersection and associated face  
! ------------------------------------------------------------------------------
                    
      CALL RFLU_ICT_TestInCellLohner(pRegion,xLoc,yLoc,zLoc,icg,testInCell, &
                                     iloc,ifg)                                                        
                                   
! ------------------------------------------------------------------------------
!     Check whether have remaining total distance 
! ------------------------------------------------------------------------------

! --- In-cell test successful --------------------------------------------------

      IF ( testInCell .EQV. .TRUE. ) THEN 
        pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg

        EXIT infLoop
        
! --- In-cell test failed ------------------------------------------------------
        
      ELSE 

! ----- Intersect interior face

        IF ( iloc == 0 ) THEN           
          c1 = pGrid%f2c(1,ifg)
          c2 = pGrid%f2c(2,ifg)
          
          c1k = RFLU_GetGlobalCellKind(global,pGrid,c1)
          c2k = RFLU_GetGlobalCellKind(global,pGrid,c2)

          SELECT CASE ( RFLU_GetFaceKind(global,c1k,c2k) ) 
            CASE ( FACE_KIND_AA ) ! Actual-actual face
              IF ( c1 == icg ) THEN 
                icg = c2
              ELSE 
                icg = c1
              END IF ! c1 
! TO DO
!            CASE ( FACE_KIND_AV ) ! Actual-virtual face
! Count particle intersections for each buffer           
! END TO DO
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! RFLU_GetFaceKind                    
                  
! ----- Intersect boundary face          
          
        ELSE ! Boundary face
          pPatch => pRegion%patches(iloc)

          SELECT CASE ( pPatch%bcType )             
            CASE DEFAULT 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)           
          END SELECT ! pPatch%bcType
        END IF ! iLoc 
      END IF ! testInCell
      
! ------------------------------------------------------------------------------
!     Guard against infinite loop 
! ------------------------------------------------------------------------------  
  
      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter
                    
    END DO infLoop     
  END DO ! iPcl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsLohner







! ******************************************************************************
!
! Purpose: Determine cells which contain particles using Octree approach.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!
! Output: None.
!
! Notes: 
!   1. It is not enough to simply get closest cell from Octree, because 
!      proximity to cell centroid does not guarantee that the specified 
!      location is actually inside the cell.
!   2. This routine is suitable to determine the cells which contain an 
!      initial distribution of particles, but should not be used to track 
!      moving particles, because it will be inefficient on account of it
!      not using knowledge of past position. It will also overwrite entry 
!      containing initial region in particle data structure. 
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsOct(pRegion)

  USE RFLU_ModGeometryTools, ONLY: RFLU_TestInBoundBox

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  LOGICAL :: foundFlag
  CHARACTER(CHRLEN) :: errorString
  INTEGER :: ccSize,errorFlag,icg,iPcl
  INTEGER, DIMENSION(:), ALLOCATABLE :: cc
  REAL(RFREAL) :: delFrac,xDel,xLoc,xMax,xMin,yDel,yLoc,yMax,yMin,zDel,zLoc, & 
                  zMax,zMin
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsOct',__FILE__)

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

  delFrac = 0.01_RFREAL
  
  ccSize = MIN(100,pGrid%nCells) ! Must be larger than unity 

  ALLOCATE(cc(ccSize),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cc')
  END IF ! global%error

! ******************************************************************************  
! Build bounding box
! ******************************************************************************

  xMin = MINVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
  xMax = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
  yMin = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
  yMax = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
  zMin = MINVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))
  zMax = MAXVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))

  xDel = xMax - xMin 
  yDel = yMax - yMin
  zDel = zMax - zMin

  xMin = xMin - delFrac*xDel
  xMax = xMax + delFrac*xDel 
  yMin = yMin - delFrac*yDel
  yMax = yMax + delFrac*yDel 
  zMin = zMin - delFrac*zDel
  zMax = zMax + delFrac*zDel  

! ******************************************************************************  
! Build Octree with cell centroids: NOTE only for actual cells
! ******************************************************************************

  CALL RFLU_CreateOctree(global,pGrid%nCells)

  CALL RFLU_BuildOctree(pGrid%cofg(XCOORD,1:pGrid%nCells), & 
                        pGrid%cofg(YCOORD,1:pGrid%nCells), & 
                        pGrid%cofg(ZCOORD,1:pGrid%nCells), & 
                        xMin,xMax,yMin,yMax,zMin,zMax)

! ******************************************************************************
! Loop over particles
! ******************************************************************************

  DO iPcl = 1,pPlag%nPcls

! ==============================================================================  
!   Test particle location against bounding box of partition
! ==============================================================================  
  
    xLoc = pPlag%cv(CV_PLAG_XPOS,iPcl)
    yloc = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLoc = pPlag%cv(CV_PLAG_ZPOS,iPcl)
  
    IF ( RFLU_TestInBoundBox(global,xLoc,yLoc,zLoc,xMin,xMax,yMin,yMax, & 
                             zMin,zMax) .EQV. .TRUE. ) THEN  
                               
! ------------------------------------------------------------------------------  
!     Query octree to get closest cells  
! ------------------------------------------------------------------------------ 
    
      CALL RFLU_QueryOctree(xLoc,yLoc,zLoc,ccSize,cc)
        
! ------------------------------------------------------------------------------
!     Test cells obtained from Octree if they contain specified location
! ------------------------------------------------------------------------------        
           
      foundFlag = .FALSE.     
           
      cellLoop: DO icg = 1,ccSize           
        IF ( RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,cc(icg)) .EQV. .TRUE. ) THEN
          pPlag%aiv(AIV_PLAG_REGINI,iPcl) = pRegion%iRegionGlobal
          pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = cc(icg)

          foundFlag = .TRUE.

          EXIT cellLoop
        END IF ! RFLU_ICT_TestInCell        
      END DO cellLoop
      
! ------------------------------------------------------------------------------
!     Check whether particle was found
! ------------------------------------------------------------------------------        
            
      IF ( foundFlag .EQV. .FALSE. ) THEN
        WRITE(errorString,'(I6)') iPcl 
        CALL ErrorStop(global,ERR_PLAG_PCL_NOT_FOUND,__LINE__, &
                       TRIM(errorString))
      END IF ! foundFlag

! ==============================================================================  
!   Particle located outside bounding box of partition
! ==============================================================================  

    ELSE   
      WRITE(errorString,'(I6)') iPcl 
      CALL ErrorStop(global,ERR_PLAG_PCL_NOT_FOUND,__LINE__,TRIM(errorString))      
    END IF ! RFLU_TestInBoundNox 
  END DO ! iPcl

! ******************************************************************************
! Destroy Octree and deallocate temporary memory
! ******************************************************************************

  CALL RFLU_DestroyOctree(global)

  DEALLOCATE(cc,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cc')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsOct








! ******************************************************************************
!
! Purpose: Kernel for Octree search
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to region
!   xLoc,yLoc,zLoc      Location
!
! Output:
!   icgOut              Cell containing location
!
! Notes: 
!   1. It is not enough to simply get closest cell from Octree, because 
!      proximity to cell centroid does not guarantee that the specified 
!      location is actually inside the cell.
!   2. This routine is suitable to determine the cell containing a particle 
!      after the latters position has been updated. 
!   3. This routine cannot be used in parallel runs and for particles which 
!      leave the domain or bounce from solid walls.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsOctKernel(pRegion,xLoc,yLoc,zLoc,xMin,xMax,&
                                        yMin,yMax,zMin,zMax,icgOut)

  USE RFLU_ModGeometryTools, ONLY: RFLU_TestInBoundBox
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(INOUT)   :: icgOut
  REAL(RFREAL), INTENT(IN) :: xLoc,xMax,xMin,yLoc,yMax,yMin,zLoc,zMax,zMin
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  LOGICAL :: foundFlag
  INTEGER :: ccSize,errorFlag,icgLoop,iPcl
  INTEGER, DIMENSION(:), ALLOCATABLE :: cc
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsOctKernel',__FILE__)
  
! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid

!  ccSize = MIN(100,pGrid%nCells) ! Must be larger than unity 

! TEMPORARY
  ccSize = MIN(50,pGrid%nCells) ! Must be larger than unity 
! END TEMPORARY

  ALLOCATE(cc(ccSize),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cc')
  END IF ! global%error

! ******************************************************************************
! Test particle location against bounding box of partition
! ******************************************************************************

  IF ( RFLU_TestInBoundBox(global,xLoc,yLoc,zLoc,xMin,xMax,yMin,yMax, & 
                             zMin,zMax) .EQV. .TRUE. ) THEN  

! ============================================================================== 
!  Query octree to get closest cells  
! ============================================================================== 
    
   CALL RFLU_QueryOctree(xLoc,yLoc,zLoc,ccSize,cc)

! ============================================================================== 
!  Test cells obtained from Octree if they contain specified location
! ============================================================================== 
! ------------------------------------------------------------------------------ 
           
   foundFlag = .FALSE.     

   cellLoop: DO icgLoop = 1,ccSize           
      IF ( RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,cc(icgLoop)) .EQV. .TRUE. ) THEN
        icgOut = cc(icgLoop)

        foundFlag = .TRUE.

        EXIT cellLoop
      END IF ! RFLU_ICT_TestInCell        
    END DO cellLoop

! ============================================================================== 
!   Check whether particle was found 
! ============================================================================== 
            
    IF ( foundFlag .EQV. .FALSE. ) THEN
      icgOut = CRAZY_VALUE_INT
    END IF ! foundFlag

! ******************************************************************************
! Particle located outside bounding box of partition
! ******************************************************************************

  ELSE   
    icgOut = CRAZY_VALUE_INT    
  END IF ! RFLU_TestInBoundNox

! ******************************************************************************
! Deallocate temporary memory
! ******************************************************************************

  DEALLOCATE(cc,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cc')
  END IF ! global%error
    
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsOctKernel








! ******************************************************************************
!
! Purpose: Determine cells which contain particles using modified Octree 
!   approach.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!
! Output: None.
!
! Notes: 
!   1. It is not enough to simply get closest cell from Octree, because 
!      proximity to cell centroid does not guarantee that the specified 
!      location is actually inside the cell.
!   2. This routine is suitable to determine the cell containing a particle 
!      after the latters position has been updated. 
!   3. This routine cannot be used in parallel runs and for particles which 
!      leave the domain or bounce from solid walls.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsOctMod(pRegion)

  USE RFLU_ModGeometryTools, ONLY: RFLU_TestInBoundBox

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  LOGICAL :: foundFlag
  CHARACTER(CHRLEN) :: errorString
  INTEGER :: ccSize,errorFlag,icg,iPcl
  INTEGER, DIMENSION(:), ALLOCATABLE :: cc
  REAL(RFREAL) :: delFrac,xDel,xLoc,xMax,xMin,yDel,yLoc,yMax,yMin,zDel,zLoc, & 
                  zMax,zMin
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsOctMod',__FILE__)

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

  delFrac = 0.01_RFREAL
  
!  ccSize = MIN(100,pGrid%nCells) ! Must be larger than unity 

! TEMPORARY
  ccSize = MIN(50,pGrid%nCells) ! Must be larger than unity 
! END TEMPORARY

  ALLOCATE(cc(ccSize),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cc')
  END IF ! global%error

! ******************************************************************************  
! Build bounding box
! ******************************************************************************

  xMin = MINVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
  xMax = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
  yMin = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
  yMax = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
  zMin = MINVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))
  zMax = MAXVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))

  xDel = xMax - xMin 
  yDel = yMax - yMin
  zDel = zMax - zMin

  xMin = xMin - delFrac*xDel
  xMax = xMax + delFrac*xDel 
  yMin = yMin - delFrac*yDel
  yMax = yMax + delFrac*yDel 
  zMin = zMin - delFrac*zDel
  zMax = zMax + delFrac*zDel  

! ******************************************************************************  
! Build Octree with cell centroids: NOTE only for actual cells
! ******************************************************************************

  CALL RFLU_CreateOctree(global,pGrid%nCells)

  CALL RFLU_BuildOctree(pGrid%cofg(XCOORD,1:pGrid%nCells), & 
                        pGrid%cofg(YCOORD,1:pGrid%nCells), & 
                        pGrid%cofg(ZCOORD,1:pGrid%nCells), & 
                        xMin,xMax,yMin,yMax,zMin,zMax)

! ******************************************************************************
! Loop over particles
! ******************************************************************************

  DO iPcl = 1,pPlag%nPcls

! ==============================================================================  
!   Get particle cell and position
! ==============================================================================  

    icg  = pPlag%aivOld(AIV_PLAG_ICELLS,iPcl) ! NOTE OLD cell index
  
    xLoc = pPlag%cv(CV_PLAG_XPOS,iPcl) ! NOTE already NEW particle position
    yloc = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLoc = pPlag%cv(CV_PLAG_ZPOS,iPcl)

! ==============================================================================  
!   Check if particle still in same cell 
! ==============================================================================  
    
    IF ( RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,icg) .EQV. .FALSE. ) THEN    
 
! ------------------------------------------------------------------------------  
!     Test particle location against bounding box of partition
! ------------------------------------------------------------------------------  
  
      IF ( RFLU_TestInBoundBox(global,xLoc,yLoc,zLoc,xMin,xMax,yMin,yMax, & 
                               zMin,zMax) .EQV. .TRUE. ) THEN  
                               
! ----- Query octree to get closest cells -------------------------------------- 
    
        CALL RFLU_QueryOctree(xLoc,yLoc,zLoc,ccSize,cc)
        
! ----- Test cells obtained from Octree if they contain specified location -----
           
        foundFlag = .FALSE.     

        cellLoop: DO icg = 1,ccSize           
          IF ( RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,cc(icg)) .EQV. .TRUE. ) THEN
            pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = cc(icg)

            foundFlag = .TRUE.

            EXIT cellLoop
          END IF ! RFLU_ICT_TestInCell        
        END DO cellLoop
      
! ----- Check whether particle was found ---------------------------------------
            
        IF ( foundFlag .EQV. .FALSE. ) THEN
          WRITE(errorString,'(I6)') iPcl 
          CALL ErrorStop(global,ERR_PLAG_PCL_NOT_FOUND,__LINE__, &
                         TRIM(errorString))
        END IF ! foundFlag

! ------------------------------------------------------------------------------  
!     Particle located outside bounding box of partition
! ------------------------------------------------------------------------------  

      ELSE   

! TEMPORARY
        WRITE(*,*) 'timeCurrent,iPcl,xLocOld,yLocOld,zLocOld,xLoc,yLoc,zLoc = ',&
        global%currentTime,iPcl,pPlag%cvOld(CV_PLAG_XPOS:CV_PLAG_ZPOS,iPcl),&
        xLoc,yLoc,zLoc
! END TEMPORARY

        WRITE(errorString,'(I6)') iPcl 
        CALL ErrorStop(global,ERR_PLAG_PCL_NOT_FOUND,__LINE__,TRIM(errorString))      
      END IF ! RFLU_TestInBoundNox

! ==============================================================================  
!   Particle still in same cell, copy cell index 
! ==============================================================================  
    
    ELSE   
      pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg
    END IF ! RFLU_ICT_TestInCell 
  END DO ! iPcl

! ******************************************************************************
! Destroy Octree and deallocate temporary memory
! ******************************************************************************

  CALL RFLU_DestroyOctree(global)

  DEALLOCATE(cc,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cc')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsOctMod







! ******************************************************************************
!
! Purpose: Determine cells which contain particles by following trajectory 
!   using fast algorithm.
!
! Description: Follow particle path by intersecting particle trajectory with
!   faces and taking appropriate action determined by type of intersected face.
!
! Input: 
!  pRegion      Pointer to region
!  iPclBeg      Beginning particle index
!  iPclEnd      Ending particle index
!
! Output: None.
!
! Notes:
!   1. Need to properly handle particle deletion from outflow or motion
!      between domain during an intermediate RK-stage as its final 
!      position at (n+1) step might not be outside the region.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsTrajFast(pRegion,iPclBeg,iPclEnd)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER :: iPclBeg,iPclEnd
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  LOGICAL :: inCellCheckFlag
  CHARACTER(CHRLEN) :: errorString
  INTEGER :: c1,c1k,c2,c2k,errorFlag,iBorder,icg,ifg,ifgOut,ifl,iloc,ilocOut, & 
             iPatch,iPcl,loopCounter
  REAL(RFREAL) :: dist,distTot,distTotCutoff,eps,fnx,fny,fnz,iMagTraj,theta, &
                  xLoc,xLocNew,xLocOld,xTraj,yLoc,yLocNew,yLocOld,yTraj,zLoc, &
                  zLocNew,zLocOld,zTraj
  TYPE(t_border), POINTER :: pBorder
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsTrajFast',__FILE__)

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

  eps = EPSILON(1.0_RFREAL)
  distTotCutoff = 10*EPSILON(1.0_RFREAL) ! NOTE must be less than tolerICT

! ******************************************************************************
! Loop over particles
! ******************************************************************************

  DO iPcl = iPclBeg,iPclEnd
    loopCounter = 0 ! Reset loop counter 

! ==============================================================================  
!   Set distance travelled and compute trajectory. NOTE that cannot compute 
!   trajectory (unit) vector from distTot because distTot diminishes as the 
!   particle travels along the trajectory but xLocNew and xLocOld (as well as
!   y and z components) do not change.
! ==============================================================================  

    distTot = pPlag%arv(ARV_PLAG_DISTOT,iPcl)

    IF ( distTot < distTotCutOff ) THEN  
      CYCLE 
    END IF ! distTot

    xLocNew = pPlag%cv(CV_PLAG_XPOS,iPcl)
    yLocNew = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLocNew = pPlag%cv(CV_PLAG_ZPOS,iPcl)  

    xLocOld = pPlag%cvOld(CV_PLAG_XPOS,iPcl)
    ylocOld = pPlag%cvOld(CV_PLAG_YPOS,iPcl)
    zLocOld = pPlag%cvOld(CV_PLAG_ZPOS,iPcl)
        
    xTraj = xLocNew - xLocOld
    yTraj = yLocNew - yLocOld     
    zTraj = zLocNew - zLocOld
             
    iMagTraj = 1.0_RFREAL/(SQRT(xTraj*xTraj + yTraj*yTraj + zTraj*zTraj) + eps)
    
    xTraj = iMagTraj*xTraj
    yTraj = iMagTraj*yTraj
    zTraj = iMagTraj*zTraj    
    
! ------------------------------------------------------------------------------
!   Set variables for trajectory loop
! ------------------------------------------------------------------------------    
    
    xLoc = xLocNew - distTot*xTraj
    yLoc = yLocNew - distTot*yTraj
    zLoc = zLocNew - distTot*zTraj

    icg = pPlag%aivOld(AIV_PLAG_ICELLS,iPcl)

! ==============================================================================
!   Loop until distance travelled along trajectory consumed 
! ==============================================================================        
        
    trajLoop: DO     
      loopCounter = loopCounter + 1  
        
! ------------------------------------------------------------------------------
!     Find appropriate intersection and associated face  
! ------------------------------------------------------------------------------
                    
      CALL RFLU_ComputeLineCellXSectFast(pRegion,xLoc,yLoc,zLoc,xTraj,yTraj, &
                                         zTraj,icg,dist,iloc,ifg)

! ------------------------------------------------------------------------------
!     Update total distance travelled                             
! ------------------------------------------------------------------------------
                                   
      distTot = distTot - dist

! ------------------------------------------------------------------------------
!     Check whether have remaining total distance 
! ------------------------------------------------------------------------------

! --- No distance remaining: Found cell containing new location ----------------

      IF ( distTot <= 0.0_RFREAL ) THEN 
        pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg

        EXIT trajLoop
        
! --- Distance remaining: Keep searching ---------------------------------------        
        
      ELSE 

! ----- Trajectory intersects interior face

        IF ( iloc == 0 ) THEN           
          c1 = pGrid%f2c(1,ifg)
          c2 = pGrid%f2c(2,ifg)
          
          c1k = RFLU_GetGlobalCellKind(global,pGrid,c1)
          c2k = RFLU_GetGlobalCellKind(global,pGrid,c2)

          SELECT CASE ( RFLU_GetFaceKind(global,c1k,c2k) ) 
            CASE ( FACE_KIND_AA ) ! Actual-actual face
              IF ( c1 == icg ) THEN 
                icg = c2
              ELSE 
                icg = c1
              END IF ! c1                                    
            CASE ( FACE_KIND_AV ) ! Actual-virtual face
              ifl = ifg - pGrid%nFaces + pGrid%nFacesAV

              iBorder = pGrid%avf2b(1,ifl)

              pBorder => pGrid%borders(iBorder)
              
              pBorder%nPclsSend = pBorder%nPclsSend + 1

              IF ( pBorder%nPclsSend > pBorder%nPclsSendMax ) &
                CALL RFLU_MPI_RecreateBufferIPclSend(pRegion,pBorder)
          
              pBorder%iPclSend(1,pBorder%nPclsSend) = iPcl
              pBorder%iPclSend(2,pBorder%nPclsSend) = pGrid%avf2b(2,ifl)            

              pPlag%aiv(AIV_PLAG_STATUS,iPcl) = PLAG_STATUS_COMM
              pPlag%arv(ARV_PLAG_DISTOT,iPcl) = distTot

              IF ( c1 == icg ) THEN 
                icg = c2
              ELSE 
                icg = c1
              END IF ! c1
              
              pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg
              
	      iPatch = pGrid%avf2p(ifl)
	      
	      IF ( iPatch /= CRAZY_VALUE_INT ) THEN ! Transform coordinates
	        pPatch => pRegion%patches(iPatch)
		
		IF ( pPatch%bcType == BC_PERIODIC ) THEN 
		  CALL RFLU_RELP_TransformVector(pRegion,pPatch,xLocOld, &
		                                 yLocOld,zLocOld)
		  CALL RFLU_RELP_TransformVector(pRegion,pPatch,xLocNew, &
		                                 yLocNew,zLocNew)						 
		ELSE
		  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
		END IF ! pPatch%bcType
	      END IF ! iPatch	      
	      
              EXIT trajLoop                  
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! RFLU_GetFaceKind                    
                  
! ----- Trajectory intersects boundary face          
          
        ELSE ! Boundary face
          pPatch => pRegion%patches(iloc)
          
          fnx = pPatch%fn(XCOORD,ifg)
          fny = pPatch%fn(YCOORD,ifg)
          fnz = pPatch%fn(ZCOORD,ifg) 
                      
          theta = xTraj*fnx + yTraj*fny + zTraj*fnz 
          
          SELECT CASE ( pPatch%bcType )             
            CASE ( BC_SLIPWALL )
              IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
                CALL PLAG_GatherSurfStats(pRegion,pPlag,pPatch%statsPlag,&
                                          ifg,iPcl,theta)
              END IF ! pPatch%plotStatsFlag
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE ( BC_NOSLIPWALL_HFLUX, &
		   BC_NOSLIPWALL_TEMP   ) 
              IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
                CALL PLAG_GatherSurfStats(pRegion,pPlag,pPatch%statsPlag,&
                                          ifg,iPcl,theta)
              END IF ! pPatch%plotStatsFlag
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE ( BC_INJECTION ) 
              IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
                CALL PLAG_GatherSurfStats(pRegion,pPlag,pPatch%statsPlag,&
                                          ifg,iPcl,theta)
              END IF ! pPatch%plotStatsFlag
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE ( BC_OUTFLOW )
              IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
                CALL PLAG_GatherSurfStats(pRegion,pPlag,pPatch%statsPlag,&
                                          ifg,iPcl,theta)
              END IF ! pPatch%plotStatsFlag
              pPlag%aiv(AIV_PLAG_STATUS,iPcl) = PLAG_STATUS_DELETE
              EXIT trajLoop            
            CASE ( BC_FARFIELD )
              pPlag%aiv(AIV_PLAG_STATUS,iPcl) = PLAG_STATUS_DELETE
              EXIT trajLoop 
            CASE ( BC_SYMMETRY )
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE ( BC_VIRTUAL )
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE DEFAULT 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)           
          END SELECT ! pPatch%bcType

        END IF ! iLoc 
      END IF ! distTot
    
! ------------------------------------------------------------------------------
!     Guard against infinite loop 
! ------------------------------------------------------------------------------  
  
      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
             'Infinite loop encountered in particle cell search algorithm'
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Module: Lagrangian Particle (PLAG).'

        WRITE(STDOUT,'(A,3X,A,1X,1PE12.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime

        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                         pRegion%iRegionGlobal
        WRITE(STDOUT,'(A,6X,A,11(1X,A))') SOLVER_NAME,'#', &
                                         '    iPcl     ', &
                                         '    idIni    ', &
                                         '    RegIni   ', &
                                         '    icg      ', &
                                         '  x-location ', &
                                         '  y-location ', &
                                         '  z-Location ', &
                                         '   Energy    ', &
                                         '   Diameter  '

        WRITE(STDOUT,'(A,4X,4(1X,I16),6(1X,E13.6))') SOLVER_NAME,iPcl, &
                        pPlag%aivOld(AIV_PLAG_PIDINI,iPcl),           &
                        pPlag%aivOld(AIV_PLAG_REGINI,iPcl),           &
                        icg,xLoc,yLoc,zLoc,                           &
                        pPlag%cv(CV_PLAG_ENER,iPcl),                  &
                        pPlag%dv(DV_PLAG_DIAM,iPcl)

        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter        
    END DO trajLoop 

! ==============================================================================
!   Check that kept particles are indeed located in new cells
! ==============================================================================        

    IF ( (global%checkLevel == CHECK_HIGH) .AND. &
	 (pPlag%aiv(AIV_PLAG_STATUS,iPcl) == PLAG_STATUS_KEEP) ) THEN
      CALL RFLU_ICT_TestInCellFancy(pRegion,xLocNew,yLocNew,zLocNew,icg, &
				    inCellCheckFlag,ilocOut,ifgOut)

      IF ( inCellCheckFlag .EQV. .FALSE. ) THEN 
	WRITE(STDERR,'(A,1X,A,1X,I6)') SOLVER_NAME,'Particle index:',iPcl
	WRITE(STDERR,'(A,1X,A,1X,I6)') SOLVER_NAME,'Cell which failed test:', & 
						   icg
	WRITE(STDERR,'(A,1X,A,2(1X,I6))') SOLVER_NAME, &
					  'Face which failed test:', & 
					  ilocOut,ifgOut								 
	WRITE(STDERR,'(A,1X,A,3(1X,E23.16))') SOLVER_NAME, &
					      'Particle old location:', &
					      xLocOld,yLocOld,zLocOld
	WRITE(STDERR,'(A,1X,A,3(1X,E23.16))') SOLVER_NAME, &
					      'Particle new location:', &
					      xLocNew,yLocNew,zLocNew

	WRITE(errorString,'(I6)') iPcl
	CALL ErrorStop(global,ERR_PLAG_PCL_NOT_FOUND,__LINE__, & 
		       TRIM(errorString))	    
      END IF ! inCellCheckFlag
    END IF ! global%checkLevel    
    
! ==============================================================================
!   Store new position. IMPORTANT: Need to store new position because it may 
!   have been reflected and old position because it may have been transformed.
! ==============================================================================

    pPlag%cv(CV_PLAG_XPOS,iPcl) = xLocNew
    pPlag%cv(CV_PLAG_YPOS,iPcl) = yLocNew
    pPlag%cv(CV_PLAG_ZPOS,iPcl) = zLocNew
     
    pPlag%cvOld(CV_PLAG_XPOS,iPcl) = xLocOld
    pPlag%cvOld(CV_PLAG_YPOS,iPcl) = yLocOld
    pPlag%cvOld(CV_PLAG_ZPOS,iPcl) = zLocOld       
  END DO ! iPcl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsTrajFast







! ******************************************************************************
!
! Purpose: Determine cells which contain particles by following trajectory 
!   using safe algorithm.
!
! Description: Follow particle path by intersecting particle trajectory with
!   faces and taking appropriate action determined by type of intersected face.
!
! Input: 
!  pRegion      Pointer to region
!  iPclBeg      Beginning particle index
!  iPclEnd      Ending particle index
!
! Output: None.
!
! Notes:
!   1. Need to properly handle particle deletion from outflow or motion
!      between domain during an intermediate RK-stage as its final 
!      position at (n+1) step might not be outside the region.
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_FindCellsTrajSafe(pRegion,iPclBeg,iPclEnd)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER :: iPclBeg,iPclEnd
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  LOGICAL :: inCellCheckFlag
  CHARACTER(CHRLEN) :: errorString
  INTEGER :: c1,c1k,c2,c2k,errorFlag,iBorder,icg,ifg,ifgOut,ifl,iloc, &
             ilocOut,iPatch,iPcl,loopCounter
  REAL(RFREAL) :: dist,distTot,distTotCutoff,eps,fnx,fny,fnz,iMagTraj,theta, &
                  xLoc,xLocNew,xLocOld,xTraj,yLoc,yLocNew,yLocOld,yTraj,zLoc, &
                  zLocNew,zLocOld,zTraj
  TYPE(t_border), POINTER :: pBorder
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************
  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_FindCellsTrajSafe',__FILE__)

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

  eps = EPSILON(1.0_RFREAL)
  distTotCutoff = 10*eps ! NOTE must be less than tolerICT

! ******************************************************************************
! Loop over particles
! ******************************************************************************

  DO iPcl = iPclBeg,iPclEnd
    loopCounter = 0 ! Reset loop counter 

! ==============================================================================  
!   Set distance travelled and compute trajectory. NOTE that cannot compute 
!   trajectory (unit) vector from distTot because distTot diminishes as the 
!   particle travels along the trajectory but xLocNew and xLocOld (as well as
!   y and z components) do not change.
! ==============================================================================  

    distTot = pPlag%arv(ARV_PLAG_DISTOT,iPcl)
    
    IF ( distTot < distTotCutoff ) THEN  
      CYCLE 
    END IF ! distTot

    xLocNew = pPlag%cv(CV_PLAG_XPOS,iPcl)
    yLocNew = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zLocNew = pPlag%cv(CV_PLAG_ZPOS,iPcl)  

    xLocOld = pPlag%cvOld(CV_PLAG_XPOS,iPcl)
    ylocOld = pPlag%cvOld(CV_PLAG_YPOS,iPcl)
    zLocOld = pPlag%cvOld(CV_PLAG_ZPOS,iPcl)
    
    xTraj = xLocNew - xLocOld
    yTraj = yLocNew - yLocOld      
    zTraj = zLocNew - zLocOld

    iMagTraj = 1.0_RFREAL/(SQRT(xTraj*xTraj + yTraj*yTraj + zTraj*zTraj) + eps)
    
    xTraj = iMagTraj*xTraj
    yTraj = iMagTraj*yTraj
    zTraj = iMagTraj*zTraj    
    
! ------------------------------------------------------------------------------
!   Set variables for trajectory loop
! ------------------------------------------------------------------------------    
    
    xLoc = xLocNew - distTot*xTraj
    yLoc = yLocNew - distTot*yTraj
    zLoc = zLocNew - distTot*zTraj
    
    icg = pPlag%aivOld(AIV_PLAG_ICELLS,iPcl)

! ==============================================================================
!   Loop until distance travelled along trajectory consumed 
! ==============================================================================        
        
    trajLoop: DO     
      loopCounter = loopCounter + 1  
        
! ------------------------------------------------------------------------------
!     Find appropriate intersection and associated face  
! ------------------------------------------------------------------------------
      CALL RFLU_ComputeLineCellXSectSafe(pRegion,xLoc,yLoc,zLoc,xTraj,yTraj, &
                                         zTraj,icg,dist,iloc,ifg)  
                                                                                              
! ------------------------------------------------------------------------------
!     Update total distance travelled                             
! ------------------------------------------------------------------------------
                                   
      distTot = distTot - dist
      
! ------------------------------------------------------------------------------
!     Check whether have remaining total distance 
! ------------------------------------------------------------------------------

! --- No distance remaining: Found cell containing new location ----------------

      IF ( distTot <= 0.0_RFREAL ) THEN       
        pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg
	
        EXIT trajLoop
        
! --- Distance remaining: Keep searching ---------------------------------------        
        
      ELSE 

! ----- Trajectory intersects interior face

        IF ( iloc == 0 ) THEN          
          c1 = pGrid%f2c(1,ifg)
          c2 = pGrid%f2c(2,ifg)
          
          c1k = RFLU_GetGlobalCellKind(global,pGrid,c1)
          c2k = RFLU_GetGlobalCellKind(global,pGrid,c2)

          SELECT CASE ( RFLU_GetFaceKind(global,c1k,c2k) ) 
            CASE ( FACE_KIND_AA ) ! Actual-actual face
              IF ( c1 == icg ) THEN 
                icg = c2
              ELSE 
                icg = c1
              END IF ! c1         
                          
            CASE ( FACE_KIND_AV ) ! Actual-virtual face
              ifl = ifg - pGrid%nFaces + pGrid%nFacesAV

              iBorder = pGrid%avf2b(1,ifl)

              pBorder => pGrid%borders(iBorder)
              
              pBorder%nPclsSend = pBorder%nPclsSend + 1
                
              IF ( pBorder%nPclsSend > pBorder%nPclsSendMax ) &
                CALL RFLU_MPI_RecreateBufferIPclSend(pRegion,pBorder)
              
              pBorder%iPclSend(1,pBorder%nPclsSend) = iPcl
              pBorder%iPclSend(2,pBorder%nPclsSend) = pGrid%avf2b(2,ifl)            

              pPlag%aiv(AIV_PLAG_STATUS,iPcl) = PLAG_STATUS_COMM
              pPlag%arv(ARV_PLAG_DISTOT,iPcl) = distTot
              
              IF ( c1 == icg ) THEN 
                icg = c2
              ELSE 
                icg = c1
              END IF ! c1
                          
              pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg  
              
	      iPatch = pGrid%avf2p(ifl)
	      
	      IF ( iPatch /= CRAZY_VALUE_INT ) THEN ! Transform coordinates
	        pPatch => pRegion%patches(iPatch)
		
		IF ( pPatch%bcType == BC_PERIODIC ) THEN 
		  CALL RFLU_RELP_TransformVector(pRegion,pPatch,xLocOld, &
		                                 yLocOld,zLocOld)
		  CALL RFLU_RELP_TransformVector(pRegion,pPatch,xLocNew, &
		                                 yLocNew,zLocNew)						 
		ELSE
		  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
		END IF ! pPatch%bcType
	      END IF ! iPatch

              EXIT trajLoop                             
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! RFLU_GetFaceKind                    
                  
! ----- Trajectory intersects boundary face          
          
        ELSE ! Boundary face
          pPatch => pRegion%patches(iloc)
          
          fnx = pPatch%fn(XCOORD,ifg)
          fny = pPatch%fn(YCOORD,ifg)
          fnz = pPatch%fn(ZCOORD,ifg) 
                      
          theta = xTraj*fnx + yTraj*fny + zTraj*fnz 
          
          SELECT CASE ( pPatch%bcType )             
            CASE ( BC_SLIPWALL )
              IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
                CALL PLAG_GatherSurfStats(pRegion,pPlag,pPatch%statsPlag,&
                                          ifg,iPcl,theta)
              END IF ! pPatch%plotStatsFlag
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE ( BC_NOSLIPWALL_HFLUX, &
		   BC_NOSLIPWALL_TEMP   ) 
              IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
                CALL PLAG_GatherSurfStats(pRegion,pPlag,pPatch%statsPlag,&
                                          ifg,iPcl,theta)
              END IF ! pPatch%plotStatsFlag
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE ( BC_INJECTION ) 
              IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
                CALL PLAG_GatherSurfStats(pRegion,pPlag,pPatch%statsPlag,&
                                          ifg,iPcl,theta)
              END IF ! pPatch%plotStatsFlag
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE ( BC_OUTFLOW )
              IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
                CALL PLAG_GatherSurfStats(pRegion,pPlag,pPatch%statsPlag,&
                                          ifg,iPcl,theta)
              END IF ! pPatch%plotStatsFlag
              pPlag%aiv(AIV_PLAG_STATUS,iPcl) = PLAG_STATUS_DELETE
              EXIT trajLoop            
            CASE ( BC_FARFIELD )
              pPlag%aiv(AIV_PLAG_STATUS,iPcl) = PLAG_STATUS_DELETE
              EXIT trajLoop 
            CASE ( BC_SYMMETRY )
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE ( BC_VIRTUAL )
              CALL PLAG_ReflectParticleData(pPatch,pPlag,ifg,iPcl,xLocOld, &
                                            yLocOld,zLocOld,xLocNew,yLocNew, &
                                            zLocNew,xTraj,yTraj,zTraj)
            CASE DEFAULT 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)           
          END SELECT ! pPatch%bcType

        END IF ! iLoc 
      END IF ! distTot
    
! ------------------------------------------------------------------------------
!     Guard against infinite loop 
! ------------------------------------------------------------------------------  

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
             'Infinite loop encountered in particle cell search algorithm'
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Module: Lagrangian Particle (PLAG).'        
          
        WRITE(STDOUT,'(A,3X,A,1X,1PE12.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime              
                                            
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal  
        WRITE(STDOUT,'(A,6X,A,11(1X,A))') SOLVER_NAME,'#', &
                                         '    iPcl     ', &
                                         '    idIni    ', &
                                         '    RegIni   ', &
                                         '    icg      ', &
                                         '  x-location ', &
                                         '  y-location ', &
                                         '  z-Location ', &
                                         '   Energy    ', &
                                         '   Diameter  '       

        WRITE(STDOUT,'(A,4X,4(1X,I16),6(1X,E13.6))') SOLVER_NAME,iPcl, & 
                        pPlag%aivOld(AIV_PLAG_PIDINI,iPcl),           &
                        pPlag%aivOld(AIV_PLAG_REGINI,iPcl),           &
                        icg,xLoc,yLoc,zLoc,                           &
                        pPlag%cv(CV_PLAG_ENER,iPcl),                  &
                        pPlag%dv(DV_PLAG_DIAM,iPcl)

        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter        
    END DO trajLoop 

! ==============================================================================
!   Check that kept particles are indeed located in new cells
! ==============================================================================        

    IF ( (global%checkLevel == CHECK_HIGH) .AND. &
	 (pPlag%aiv(AIV_PLAG_STATUS,iPcl) == PLAG_STATUS_KEEP) ) THEN
      CALL RFLU_ICT_TestInCellFancy(pRegion,xLocNew,yLocNew,zLocNew,icg, &
				    inCellCheckFlag,ilocOut,ifgOut)

      IF ( inCellCheckFlag .EQV. .FALSE. ) THEN 
	WRITE(STDERR,'(A,1X,A,1X,I6)') SOLVER_NAME,'Particle index:',iPcl
	WRITE(STDERR,'(A,1X,A,1X,I6)') SOLVER_NAME,'Cell which failed test:', & 
						   icg
	WRITE(STDERR,'(A,1X,A,2(1X,I6))') SOLVER_NAME, &
					  'Face which failed test:', & 
					  ilocOut,ifgOut								 
	WRITE(STDERR,'(A,1X,A,3(1X,E23.16))') SOLVER_NAME, &
					      'Particle old location:', &
					      xLocOld,yLocOld,zLocOld
	WRITE(STDERR,'(A,1X,A,3(1X,E23.16))') SOLVER_NAME, &
					      'Particle new location:', &
					      xLocNew,yLocNew,zLocNew

	WRITE(errorString,'(I6)') iPcl
	CALL ErrorStop(global,ERR_PLAG_PCL_NOT_FOUND,__LINE__, & 
		       TRIM(errorString))	    
      END IF ! inCellCheckFlag
    END IF ! global%checkLevel    
    
! ==============================================================================
!   Store new position. IMPORTANT: Need to store new position because it may 
!   have been reflected and old position because it may have been transformed.
! ==============================================================================

    pPlag%cv(CV_PLAG_XPOS,iPcl) = xLocNew
    pPlag%cv(CV_PLAG_YPOS,iPcl) = yLocNew
    pPlag%cv(CV_PLAG_ZPOS,iPcl) = zLocNew
     
    pPlag%cvOld(CV_PLAG_XPOS,iPcl) = xLocOld
    pPlag%cvOld(CV_PLAG_YPOS,iPcl) = yLocOld
    pPlag%cvOld(CV_PLAG_ZPOS,iPcl) = zLocOld                  	      
  END DO ! iPcl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_FindCellsTrajSafe





! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE PLAG_RFLU_ModFindCells

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_ModFindCells.F90,v $
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
! Revision 1.3  2007/09/04 13:20:46  haselbac
! Removed BC_NOSLIPWALL from CASE statements, caused compile errors on UF cluster
!
! Revision 1.2  2007/08/07 21:54:01  fnajjar
! Removed statements with BC_RANGE since obsolete for RocfluMP
!
! Revision 1.1  2007/04/09 18:50:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:34  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.14  2006/08/18 21:11:35  fnajjar
! Enabled periodicity, cosmetics
!
! Revision 1.13  2006/05/05 02:15:49  haselbac
! Bug fixes: Incorrect computation of traj, missing EXIT, wrong IF
!
! Revision 1.12  2006/05/02 17:46:23  fnajjar
! Allowed surface statistics to be gathered on active patches
!
! Revision 1.11  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.10  2005/12/24 21:39:47  haselbac
! Adapted to changes in ICT
!
! Revision 1.9  2005/12/19 16:47:45  fnajjar
! Added verbosity around error trap for infinite loop counter
!
! Revision 1.8  2005/12/14 21:21:21  fnajjar
! Added call for dynamic allocation of iPclSend
!
! Revision 1.7  2005/12/02 20:07:50  fnajjar
! Added particle reflection for BC_VIRTUAL
!
! Revision 1.6  2005/11/12 00:34:08  fnajjar
! Bug fix for iMagTraj in Fast and Safe
!
! Revision 1.5  2005/09/20 15:46:35  fnajjar
! Fixed bug in TrajFast and added error trap for iPclSend overflow
!
! Revision 1.4  2005/07/18 20:47:41  fnajjar
! Aligned fast routine with safe adding proper construct for trajectory
!
! Revision 1.3  2005/05/18 22:17:50  fnajjar
! Added PLAG_RFLU_ComputeDistTot, changed comp of distTot and xyzLoc, 
! infrast for comm
!
! Revision 1.2  2005/05/02 21:59:38  haselbac
! Preparation for parallelization of FindCells routines, cosmetics
!
! Revision 1.1  2005/04/27 14:55:55  fnajjar
! Initial import
!
! ******************************************************************************

