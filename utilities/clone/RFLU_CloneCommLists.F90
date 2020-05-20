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
! Purpose: Clone communication lists.
!
! Description: None.
!
! Input:
!   pRegionSerial    	Pointer to serial region
!   pRegion     	Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_CloneCommLists.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_CloneCommLists(pRegionSerial,pRegion)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModBorder, ONLY: t_border
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  USE RFLU_ModCommLists

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================
   
  TYPE(t_region), POINTER :: pRegion,pRegionSerial 
   
! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iBorder,icl,ivl
  TYPE(t_border), POINTER :: pBorder,pBorderSerial
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid,pGridSerial
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_CloneCommLists.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global
  
  CALL RegisterFunction(global,'RFLU_CloneCommLists',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Cloning communication lists...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                     pRegion%iRegionGlobal       
  END IF ! global%verbLevel

! ******************************************************************************
! 
! ******************************************************************************

  pGrid       => pRegion%grid
  pGridSerial => pRegionSerial%grid

  pGrid%nBorders = pGridSerial%nBorders
  
  CALL RFLU_COMM_CreateBorders(pRegion,CREATE_BORDERS_MODE_DIM_UNKNOWN)
  
  DO iBorder = 1,pGrid%nBorders
    pBorder       => pGrid%borders(iBorder)
    pBorderSerial => pGridSerial%borders(iBorder)
  
    pBorder%nCellsSend  = pBorderSerial%nCellsSend    
    pBorder%nCellsRecv  = pBorderSerial%nCellsRecv
    pBorder%nVertSend   = pBorderSerial%nVertSend   
    pBorder%nVertRecv   = pBorderSerial%nVertRecv
    pBorder%nVertShared = pBorderSerial%nVertShared           
  END DO ! iBorder

  IF ( pRegion%iRegionGlobal == 1 ) THEN 
    pGrid%borders(1)%iRegionGlobal = global%nRegionsLocal
    pGrid%borders(2)%iRegionGlobal = 2    
  ELSE IF ( pRegion%iRegionGlobal == global%nRegionsLocal ) THEN 
    pGrid%borders(1)%iRegionGlobal = global%nRegionsLocal-1
    pGrid%borders(2)%iRegionGlobal = 1
  ELSE 
    pGrid%borders(1)%iRegionGlobal = pRegion%iRegionGlobal-1
    pGrid%borders(2)%iRegionGlobal = pRegion%iRegionGlobal+1  
  END IF ! pRegion%iRegionGlobal

  pGrid%borders(1)%iBorder = 2
  pGrid%borders(2)%iBorder = 1 

  CALL RFLU_COMM_CreateCommLists(pRegion)
  
  DO iBorder = 1,pGrid%nBorders
    pBorder       => pGrid%borders(iBorder)
    pBorderSerial => pGridSerial%borders(iBorder)
  
    DO icl = 1,pBorder%nCellsSend
      pBorder%icgSend(icl) = pBorderSerial%icgSend(icl)
    END DO ! icl
    
    DO icl = 1,pBorder%nCellsRecv
      pBorder%icgRecv(icl) = pBorderSerial%icgRecv(icl)
    END DO ! icl
    
    DO ivl = 1,pBorder%nVertSend
      pBorder%ivgSend(ivl) = pBorderSerial%ivgSend(ivl)
    END DO ! ivl  
    
    DO ivl = 1,pBorder%nVertRecv
      pBorder%ivgRecv(ivl) = pBorderSerial%ivgRecv(ivl)
    END DO ! ivl  
    
    DO ivl = 1,pBorder%nVertShared
      pBorder%ivgShared(ivl) = pBorderSerial%ivgShared(ivl)
    END DO ! ivl                       
  END DO ! iBorder  

! ******************************************************************************
! End
! ******************************************************************************
 
  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Cloning communication lists done.'
  END IF ! global%verbLevel  
 
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CloneCommLists


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CloneCommLists.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:54  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:06  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/08/07 17:13:59  haselbac
! Initial revision
!
! ******************************************************************************

