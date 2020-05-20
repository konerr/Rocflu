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
! Purpose: Pick regions.
!
! Description: None.
!
! Input: 
!   regions		Region data
!
! Output: None.
!
! Notes: 
!   1. If selecting some regions, initialize postActiveFlag to FALSE only if 
!      coordinate-activation not selected.
!
!******************************************************************************
!
! $Id: RFLU_PickRegionsManual.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_PickRegionsManual(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModMPI
  USE ModParameters
  
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

  CHARACTER :: infoType
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iReg
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PickRegionsManual.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_PickRegionsManual',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Picking regions manually...'
  END IF ! global%verbLevel

! ******************************************************************************
! Enter information
! ******************************************************************************

  WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Enter information on regions:'
  WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'a - Pick all regions'   
  WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'s - Pick some regions'
  WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'n - Pick no regions'
  WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Enter information type:'
  READ(STDIN,'(A)') infoType
   
  SELECT CASE ( infoType ) 
  
! ==============================================================================
!   All regions
! ==============================================================================
    
    CASE ( 'a' ) 
      DO iReg = 1,global%nRegionsLocal
        regions(iReg)%postActiveFlag = .TRUE.
      END DO ! iReg    

! ==============================================================================
!   Some regions
! ==============================================================================

    CASE ( 's' ) 
      IF ( global%pickCoordFlag .EQV. .FALSE. ) THEN 
        DO iReg = 1,global%nRegionsLocal
          regions(iReg)%postActiveFlag = .FALSE.
        END DO ! iReg         
      END IF ! global%pickCoordFlag 
            
      DO     
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
                                 'Enter global region index (< 1 to exit):'       
        READ(STDIN,*) iReg
        
        IF ( iReg < 1 ) THEN 
          EXIT
        ELSE IF ( iReg > global%nRegionsLocal ) THEN
          global%warnCounter = global%warnCounter + 1         
         
          WRITE(STDOUT,'(A,5X,A,1X,A,I7,1X,A)') SOLVER_NAME, & 
                       '*** WARNING *** Invalid input.', & 
                       'There are only ',global%nRegionsLocal,'regions.'
        ELSE IF ( regions(iReg)%postActiveFlag .EQV. .TRUE. ) THEN 
          global%warnCounter = global%warnCounter + 1         
         
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, & 
                       '*** WARNING *** Region already active.'
        END IF ! iReg
        
        regions(iReg)%postActiveFlag = .TRUE.
      END DO ! <empty>

! ==============================================================================
!   No regions
! ==============================================================================
      
    CASE ( 'n' ) 

! ==============================================================================
!   Default
! ==============================================================================
      
    CASE DEFAULT 
      global%warnCounter = global%warnCounter + 1     
    
      WRITE(STDOUT,'(A,5X,A,1X,A)') SOLVER_NAME, &
                                    '*** WARNING *** Invalid input.', & 
                                    'No regions selected.' 
  END SELECT ! infoType
  
! ******************************************************************************
! End
! ******************************************************************************

  IF ( (global%myProcid == MASTERPROC) .AND. & 
       (global%verbLevel > VERBOSE_NONE) ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Picking regions manually done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PickRegionsManual

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PickRegionsManual.F90,v $
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
! Revision 1.2  2005/12/10 23:30:01  haselbac
! Renamed postCoordFlag
!
! Revision 1.1  2003/08/07 15:15:04  haselbac
! Used to be called RFLU_PickRegions
!
! Revision 1.2  2003/07/22 02:08:33  haselbac
! Added global%warnCounter
!
! Revision 1.1.1.1  2003/06/04 22:31:20  haselbac
! Initial revision
!
!******************************************************************************

