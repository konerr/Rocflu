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
! Purpose: Initialize particle solution in serial region for 1d cases.
!
! Description: None.
!
! Input: 
!   pRegion	Pointer to serial region
!
! Output: None.
!
! Notes: 
!   1. Assume grid spacing is constant!
!   2. Assume cells are numbered in ascending order in direction of increasing
!      x-coordinate.
!
! ******************************************************************************
!
! $Id: PLAG_RFLU_InitSolSerial_1D.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_InitSolSerial_1D(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY: t_plag
  USE ModParameters

  USE RFLU_ModFaceList, ONLY: RFLU_CreateCell2FaceList, & 
                              RFLU_BuildCell2FaceList, & 
                              RFLU_DestroyCell2FaceList    
  USE RFLU_ModGeometryTools, ONLY: RFLU_GetGridSpacingCartesian
  USE RFLU_ModInCellTest, ONLY: RFLU_ICT_TestInCell

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

  LOGICAL :: foundFlag
  CHARACTER(CHRLEN) :: errorString,RCSIdentString
  INTEGER :: dIcg,dIcgDel,dIcgMax,dIcgMin,errorFlag,icg,icg2,iPcl
  REAL(RFREAL) :: dx,idx,xMin,xPcl,yPcl,zPcl
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: PLAG_RFLU_InitSolSerial_1D.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_InitSolSerial_1D', &
                        __FILE__)
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing particle solution '// & 
                             'for serial region in 1d...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                     pRegion%iRegionGlobal
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and values
! ******************************************************************************
  
  pGrid => pRegion%grid
  pPlag => pRegion%plag  

! ******************************************************************************
! Get geometric data
! ******************************************************************************

  xMin = MINVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))

  CALL RFLU_GetGridSpacingCartesian(pRegion,XCOORD,dx)
  idx = 1.0_RFREAL/dx
       
! ******************************************************************************
! Build cell-to-face list (needed for in-cell test)
! ******************************************************************************

  CALL RFLU_CreateCell2FaceList(pRegion)
  CALL RFLU_BuildCell2FaceList(pRegion)

! ******************************************************************************
! Loop over particles in serial region
! ******************************************************************************

  DO iPcl = 1,pPlag%nPcls    

! ==============================================================================
!   Get particle location and determine cell which should contain this particle
! ==============================================================================

    xPcl = pPlag%cv(CV_PLAG_XPOS,iPcl)
    yPcl = pPlag%cv(CV_PLAG_YPOS,iPcl)
    zPcl = pPlag%cv(CV_PLAG_ZPOS,iPcl)

    icg = INT(1.0_RFREAL + idx*(xPcl-xMin))

! ==============================================================================
!   Carry out in-cell test
! ==============================================================================

    IF ( RFLU_ICT_TestInCell(pRegion,xPcl,yPcl,zPcl,icg) .EQV. .TRUE. ) THEN
      pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg
      pPlag%aiv(AIV_PLAG_REGINI,iPcl) = pRegion%iRegionGlobal

! ==============================================================================
!   If the in-cell test failed, then the only reason should be that particle is
!   located right on a face and may be regarded as being in adjacent cell due 
!   to machine precision issues, so check those cells, too.
! ==============================================================================

    ELSE 
      IF ( icg == 1 ) THEN 
        dIcgMin = 1
        dIcgMax = 1
        dIcgDel = 1 
      ELSE IF ( icg == pGrid%nCells ) THEN 
        dIcgMin = -1
        dIcgMax = -1
        dIcgDel =  1 
      ELSE 
        dIcgMin = -1
        dIcgMax =  1
        dIcgDel =  2 
      END IF ! icg 

      foundFlag = .FALSE.

      DO dIcg = dIcgMin,dIcgMax,dIcgDel
        icg2 = icg + dIcg

        IF ( RFLU_ICT_TestInCell(pRegion,xPcl,yPcl,zPcl,icg2) .EQV. .TRUE. ) THEN
          pPlag%aiv(AIV_PLAG_ICELLS,iPcl) = icg2
          pPlag%aiv(AIV_PLAG_REGINI,iPcl) = pRegion%iRegionGlobal

          foundFlag = .TRUE.
        END IF ! RFLU_ICT_TestInCell
      END DO ! dIcg

! ------------------------------------------------------------------------------
!     If that check fails, too, then there must be an error
! ------------------------------------------------------------------------------

      IF ( foundFlag .EQV. .FALSE. ) THEN 
        WRITE(errorString,'(I6)') iPcl
        CALL ErrorStop(global,ERR_PLAG_PCL_NOT_FOUND,__LINE__,TRIM(errorString))
      END IF ! foundFlag 
    END IF ! RFLU_ICT_TestInCell
  END DO ! iPcl
    
! ******************************************************************************
! Destroy cell-to-face list
! ******************************************************************************

  CALL RFLU_DestroyCell2FaceList(pRegion)        

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Initializing particle solution for serial region in 1d done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_InitSolSerial_1D

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_InitSolSerial_1D.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:51  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2008/04/05 18:26:30  haselbac
! Made determination of cell containing particle more general
!
! Revision 1.1  2007/04/09 18:50:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2007/03/27 00:22:01  haselbac
! Fixed indentation
!
! Revision 1.1  2007/03/15 21:58:29  haselbac
! Initial revision
!
! ******************************************************************************

