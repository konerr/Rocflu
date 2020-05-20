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
! Purpose: Write file with special cells.
!
! Description: None.
!
! Input: 
!   pRegion              Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_WritePostInfo.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_WritePostInfo(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModMPI
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: ics,iFile,ifs
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_WritePostInfo.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_WritePostInfo',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing post-processor info...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                     pRegion%iRegionGlobal
  END IF ! global%myProcid

! ******************************************************************************
! Set grid pointer and variables
! ******************************************************************************

  pGrid => pRegion%grid

  iFile = IF_POSTINFO

! ******************************************************************************
! Write to file
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,5X,A,1X,I3)') SOLVER_NAME,'Number of special cells:', &
                                   pGrid%nCellsSpecial
    WRITE(STDOUT,'(A,5X,A,1X,I3)') SOLVER_NAME,'Number of special faces:', &
                                   pGrid%nFacesSpecial                                   
  END IF ! global%myProcid

  WRITE(iFile,'(A,1X,I5.5)') '# Region',pRegion%iRegionGlobal
  WRITE(iFile,'(L1)') pRegion%postActiveFlag
  
  IF ( pRegion%postActiveFlag .EQV. .TRUE. ) THEN   
    WRITE(iFile,*) pGrid%nCellsSpecial

    DO ics = 1,pGrid%nCellsSpecial
      WRITE(iFile,*) pGrid%cellsSpecial(ics)
    END DO ! ics
    
    WRITE(iFile,*) pGrid%nFacesSpecial
    
    DO ifs = 1,pGrid%nFacesSpecial
      WRITE(iFile,*) pGrid%facesSpecial(1:2,ifs)
    END DO ! ifs    
  END IF ! pRegion%postActiveFlag

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing post-processor info done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_WritePostInfo

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_WritePostInfo.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:57  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:11  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:57:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.5  2004/09/27 01:40:14  haselbac
! Added writing of special faces
!
! Revision 1.4  2004/03/23 03:18:46  haselbac
! Changed format statements
!
! Revision 1.3  2004/03/15 21:10:09  haselbac
! Fixed alignment
!
! Revision 1.2  2003/08/07 15:35:50  haselbac
! Changed var name, fixed small bug
!
! Revision 1.1.1.1  2003/06/04 22:31:20  haselbac
! Initial revision
!
! ******************************************************************************

