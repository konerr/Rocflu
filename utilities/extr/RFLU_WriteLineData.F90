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
! Purpose: Extract data along line in quadrilateral grid. 
!
! Description: None.
!
! Input: 
!   regions		Region data
!   iRegStart		Starting region index
!   icgStart  		Starting cell index
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_WriteLineData.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_WriteLineData(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModMPI
  USE ModParameters

  USE ModBuildFileNames, ONLY: BuildFileNamePlain

  USE RFLU_ModPlottingVars
  
  USE ModInterfaces, ONLY: RFLU_AllocateReadComputeVars, & 
                           RFLU_DeallocateReadComputeVars
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName,RCSIdentString
  CHARACTER(2*CHRLEN) :: header
  INTEGER :: errorFlag,icg,icl,iPv,iPv2,iRegNew,iRegOld,iVar,iXSect,iXv, & 
             nVars,nXv
  REAL(RFREAL) :: dist,distLoc,distTot,nx,ny,rvn,rvt,tx,ty,wght,xBeg,xEnd, & 
                  xLoc,yBeg,yEnd,yLoc,zBeg,zEnd,zLoc
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: var
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pXv
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGridSerial
  TYPE(t_region), POINTER :: pRegion,pRegionSerial

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_WriteLineData.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_WriteLineData',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing line data...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pRegionSerial => regions(0)
  pGridSerial   => pRegionSerial%grid

  xBeg = global%extrXCoordBeg
  yBeg = global%extrYCoordBeg
  zBeg = global%extrZCoordBeg

  xEnd = global%extrXCoordEnd
  yEnd = global%extrYCoordEnd
  zEnd = global%extrZCoordEnd

  distTot = SQRT((xEnd-xBeg)**2 + (yEnd-yBeg)**2 + (zEnd-zBeg)**2)

  tx = (xEnd-xBeg)/distTot
  ty = (yEnd-yBeg)/distTot

  nx = -ty
  ny =  tx 

! ******************************************************************************
! Allocate memory and open file 
! ******************************************************************************

  CALL RFLU_CountPlottingVars(pRegionSerial)
  CALL RFLU_CreatePlottingVarMaps(pRegionSerial)
  CALL RFLU_BuildPlottingVarMaps(pRegionSerial)

  nVars = pRegionSerial%mixtInput%nCv & 
        + pRegionSerial%mixtInput%nDv & 
        + pRegionSerial%plot%nPv
        
  ALLOCATE(var(nVars),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'var')
  END IF ! global%error

  CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.lin',iFileName)

  OPEN(IF_EXTR_DATA1,FILE=iFileName,FORM='FORMATTED',STATUS='UNKNOWN', &
       IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! =============================================================================
! Write header
! =============================================================================

  WRITE(header,'(12(1X,A))') 'x','y','z','s', &
                             'r','run','rut','rw','rE', & 
                             'p','T','a'

  DO iPv = 1,pRegionSerial%plot%nPv
    iPv2 = pRegionSerial%plot%pvi2pv(iPv)

    WRITE(header,'(A,1X,A)') TRIM(header), &
                             TRIM(pRegionSerial%plot%pvNameShort(iPv2))
  END DO ! iPv 

  WRITE(IF_EXTR_DATA1,'(A,1X,A)') 'VARIABLES=',TRIM(header)
  WRITE(IF_EXTR_DATA1,'(A,1X,I6)') 'ZONE I=',pGridSerial%nXSect

! ******************************************************************************
! Loop over intersection locations
! ******************************************************************************

  iRegOld = CRAZY_VALUE_INT ! Initial value important

  DO iXSect = 1,pGridSerial%nXSect
    iRegNew = pGridSerial%xSectList(1,iXSect)

    pRegion => regions(iRegNew)

! ==============================================================================
!   Allocate, read, and compute read data if necessary
! ==============================================================================

    IF ( (iXSect == 1) .OR. (iRegNew /= iRegOld) ) THEN       
      CALL RFLU_AllocateReadComputeVars(pRegion)
    END IF ! iXSect

! ==============================================================================
!   Get intersection location
! ==============================================================================

    xLoc = pGridSerial%xSectGeom(1,iXSect)
    yLoc = pGridSerial%xSectGeom(2,iXSect)
    zLoc = pGridSerial%xSectGeom(3,iXSect)

    distLoc = SQRT((xLoc-xBeg)**2 + (yLoc-yBeg)**2 + (zLoc-zBeg)**2)
    dist    = distLoc/distTot

    IF ( (global%myProcid == MASTERPROC) .AND. & 
         (global%verbLevel > VERBOSE_LOW ) ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I4)') SOLVER_NAME,'Intersection:',iXSect
      WRITE(STDOUT,'(A,5X,A,1X,I5.5)') SOLVER_NAME,'Region:',iRegNew
      WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME,'Distance:',dist
    END IF ! global%verbLevel

! ==============================================================================
!   Interpolate to intersection location 
! ==============================================================================

    DO iVar = 1,nVars
      var(iVar) = 0.0_RFREAL
    END DO ! iVar
       
    DO icl = 1,4
      icg  = pGridSerial%xSectList(1+icl,iXSect)
      wght = pGridSerial%xSectGeom(3+icl,iXSect)

      iVar = 0 

! ------------------------------------------------------------------------------
!     Conserved variables 
! ------------------------------------------------------------------------------

      pXv => pRegion%mixt%cv
      nXv =  UBOUND(pXv,1)

      DO iXv = 1,nXv
        iVar = iVar + 1

        var(iVar) = var(iVar) + wght*pXv(iXv,icg)
      END DO ! iXv

! ------------------------------------------------------------------------------
!     Dependent variables 
! ------------------------------------------------------------------------------

      pXv => pRegion%mixt%dv
      nXv =  UBOUND(pXv,1)

      DO iXv = 1,nXv
        iVar = iVar + 1

        var(iVar) = var(iVar) + wght*pXv(iXv,icg)
      END DO ! iXv

! ------------------------------------------------------------------------------
!     Plotting variables 
! ------------------------------------------------------------------------------

      pXv => pRegion%plot%pv
      nXv =  UBOUND(pXv,1)

      DO iXv = 1,nXv
        iVar = iVar + 1

        var(iVar) = var(iVar) + wght*pXv(iXv,icg)
      END DO ! iXv  
    END DO ! icl

! ==============================================================================
!  Transform x- and y-momenta into normal and tangential components (wrt 
!  extraction line). NOTE must be done here (i.e. outside icl loop because 
!  otherwise transform repeatedly.
! ==============================================================================

   rvn = var(CV_MIXT_XMOM)*nx + var(CV_MIXT_YMOM)*ny
   rvt = var(CV_MIXT_XMOM)*tx + var(CV_MIXT_YMOM)*ty

   var(CV_MIXT_XMOM) = rvn
   var(CV_MIXT_YMOM) = rvt
 
! ==============================================================================
!   Write to file. NOTE fixed format.
! ==============================================================================

    WRITE(IF_EXTR_DATA1,'(30(1X,E13.6))') xLoc,yLoc,zLoc,dist, & 
                                          (var(iVar),iVar=1,nVars) 

! ==============================================================================
!   Deallocate memory if necessary
! ==============================================================================

    IF ( iXSect == pGridSerial%nXSect ) THEN 
      CALL RFLU_DeallocateReadComputeVars(pRegion)
    ELSE 
      iRegOld = iRegNew
      iRegNew = pGridSerial%xSectList(1,iXSect+1)

      IF ( iRegNew /= iRegOld ) THEN 
        CALL RFLU_DeallocateReadComputeVars(pRegion)
      END IF ! iRegNew
    END IF ! iXSect
  END DO ! iXSect

! ******************************************************************************
! Deallocate memory and close file 
! ******************************************************************************

  CALL RFLU_DestroyPlottingVarMaps(pRegionSerial)

  DEALLOCATE(var,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'var')
  END IF ! global%error

  CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( (global%myProcid == MASTERPROC) .AND. & 
       (global%verbLevel > VERBOSE_NONE) ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing line data done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_WriteLineData

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_WriteLineData.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.5  2008/12/06 08:43:55  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:07  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2008/01/08 19:15:50  haselbac
! Bug fix for normal and tangential velocities
!
! Revision 1.2  2007/12/05 13:26:31  haselbac
! Added writing of dv and pv, normal and tangential vel, file header
!
! Revision 1.1  2007/11/27 13:17:26  haselbac
! Initial revision
!
!******************************************************************************

