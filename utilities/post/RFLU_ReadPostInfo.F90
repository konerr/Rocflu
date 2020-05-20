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
! Purpose: Read file with with post-processor information.
!
! Description: None.
!
! Input: 
!   pRegion              Pointer to region data
!   readMode             Reading mode
!
! Output: None.
!
! Notes: 
!   1. Need to have two reading modes because for the activation flag to be
!      useful, need to do this before any quantities are read or allocated, 
!      but this means that nCellsSpecial will be overwritten when grid is 
!      created. So read twice, first for activation flag only, and then after 
!      grid is created, read nCellsSpecial and the actual indices of the 
!      special cells.
!
! ******************************************************************************
!
! $Id: RFLU_ReadPostInfo.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ReadPostInfo(pRegion,readMode)

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

  INTEGER, INTENT(IN) :: readMode
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: dummyLogical
  CHARACTER(CHRLEN) :: dummyString,dummyString2,iRegionString,RCSIdentString
  INTEGER :: dummyInteger,dummyInteger2,ics,iFile,ifs,indx,iReg,iRegionGlobal, &
             nCellsSpecial
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ReadPostInfo.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ReadPostInfo',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading post-processor info...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                     pRegion%iRegionGlobal 
  END IF ! global%verbLevel

! ******************************************************************************
! Set grid pointer and variables
! ******************************************************************************

  pGrid => pRegion%grid

  iFile = IF_POSTINFO

! ******************************************************************************
! Rewind file - must be done because stays open while read several times
! ******************************************************************************

  REWIND(iFile)

! ******************************************************************************
! Read from file
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    READ(iFile,'(A)') dummyString

    WRITE(iRegionString,'(I5.5)') pRegion%iRegionGlobal
    indx = INDEX(dummyString,TRIM(iRegionString))
    
    IF ( indx /= 0 ) THEN
      dummyString2 = dummyString(indx:indx+LEN_TRIM(iRegionString)-1)
      READ(dummyString2,*) iRegionGlobal
    ELSE 
      iRegionGlobal = CRAZY_VALUE_INT ! Anything but pRegion%iRegionGlobal
    END IF ! indx
        
! ==============================================================================
!   Found entry for my region
! ==============================================================================

    IF ( iRegionGlobal == pRegion%iRegionGlobal ) THEN
    
! ------------------------------------------------------------------------------
!     Read activation flag only
! ------------------------------------------------------------------------------
    
      IF ( readMode == INFOFILE_READMODE_FLAG ) THEN
        READ(iFile,'(L1)') pRegion%postActiveFlag
        
        IF ( pRegion%postActiveFlag .EQV. .TRUE. ) THEN
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Region is active.'
          END IF ! global%verbLevel          
                         
          READ(iFile,*) dummyInteger ! Special cells

          DO ics = 1,dummyInteger
            READ(iFile,*) dummyInteger2
          END DO ! ics
          
          READ(iFile,*) dummyInteger ! Special faces

          DO ifs = 1,dummyInteger
            READ(iFile,*) dummyInteger2
          END DO ! ifs          
        ELSE 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Region is not active.'
          END IF ! global%verbLevel         
        END IF ! pRegion%postActiveFlag   

! ------------------------------------------------------------------------------
!     Read special cells and faces only
! ------------------------------------------------------------------------------

      ELSE IF ( readMode == INFOFILE_READMODE_DATA ) THEN
        READ(iFile,'(L1)') dummyLogical
            
        IF ( dummyLogical .EQV. .TRUE. ) THEN 
          READ(iFile,*) pGrid%nCellsSpecial

          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME, &
                                           'Number of special cells:', &
                                           pGrid%nCellsSpecial
          END IF ! global%verbLevel

          DO ics = 1,pGrid%nCellsSpecial
            READ(iFile,*) pGrid%cellsSpecial(ics)
          END DO ! ics 
          
          READ(iFile,*) pGrid%nFacesSpecial

          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME, &
                                           'Number of special faces:', &
                                           pGrid%nFacesSpecial
          END IF ! global%verbLevel

          DO ifs = 1,pGrid%nFacesSpecial
            READ(iFile,*) pGrid%facesSpecial(1,ifs), & 
                          pGrid%facesSpecial(2,ifs)
          END DO ! ifs                                        
        ELSE 
          pGrid%nCellsSpecial = 0
          pGrid%nFacesSpecial = 0          
        END IF ! dummyLogical
      ELSE 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! readMode
      
! ==============================================================================
!   Found entry for other region(s) - skip by reading dummy variables
! ==============================================================================

    ELSE 
      READ(iFile,'(L1)') dummyLogical
    
      IF ( dummyLogical .EQV. .TRUE. ) THEN 
        READ(iFile,*) dummyInteger ! Special cells

        DO ics = 1,dummyInteger
          READ(iFile,*) dummyInteger2
        END DO ! ics 
        
        READ(iFile,*) dummyInteger ! Special faces

        DO ifs = 1,dummyInteger
          READ(iFile,*) dummyInteger2
        END DO ! ifs          
      END IF ! dummyLogical    
    END IF ! iRegionGlobal
  END DO ! iReg

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading post-processor info done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ReadPostInfo

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ReadPostInfo.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:58  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:13  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:58:09  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.5  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.4  2004/09/27 01:43:50  haselbac
! Modified to read info about special faces
!
! Revision 1.3  2004/03/23 03:17:34  haselbac
! Changed format statements
!
! Revision 1.2  2003/08/07 15:37:16  haselbac
! Changed var names
!
! Revision 1.1  2003/06/04 22:44:14  haselbac
! Initial revision
!
! ******************************************************************************

