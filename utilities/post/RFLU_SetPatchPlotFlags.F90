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
! Purpose: Set patch plotting flags.
!
! Description: None.
!
! Input: 
!   pRegion      Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_SetPatchPlotFlags.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_SetPatchPlotFlags(pRegion)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global 
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModParameters
  USE ModMPI

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion
  
! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: bcVirtCntr,iPatch
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_SetPatchPlotFlags.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SetPatchPlotFlags',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Setting patch plot flags...'
  END IF ! global%verbLevel

! ******************************************************************************  
! Loop over patches
! ******************************************************************************

  bcVirtCntr = 0 

  SELECT CASE ( pRegion%mixtInput%dimens ) 
    CASE ( 1,2 ) 
      IF ( global%postPlotPatchFlag .EQV. .TRUE. ) THEN 
        DO iPatch = 1,pRegion%grid%nPatches
          pPatch => pRegion%patches(iPatch)

          IF ( pPatch%bcType == BC_VIRTUAL ) THEN 
            bcVirtCntr = bcVirtCntr + 1
          
            IF ( bcVirtCntr == 1 ) THEN 
              pPatch%plotFlag = .TRUE.
            ELSE 
              pPatch%plotFlag = .FALSE.
            END IF ! bcVirtCntr
          ELSE 
            pPatch%plotFlag = .FALSE.
          END IF ! pPatch%bcType
        END DO ! iPatch
      ELSE 
        DO iPatch = 1,pRegion%grid%nPatches
          pPatch => pRegion%patches(iPatch)

          pPatch%plotFlag = .TRUE. 
        END DO ! iPatch
      END IF ! global%postPlotPatchFlag
    CASE ( 3 ) 
      DO iPatch = 1,pRegion%grid%nPatches
        pPatch => pRegion%patches(iPatch)

        pPatch%plotFlag = .TRUE.
        ! Subbu - Debug
        pPatch%plotFlag = .FALSE.
        ! Subbu - End Debug
      END DO ! iPatch      
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixtInput%dimens
 
! ******************************************************************************
! Write info
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Patch plot flag information:'
    WRITE(STDOUT,'(A,5X,A,1X,A,1X,A)') SOLVER_NAME,'Local','Global','Flag' 

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

      WRITE(STDOUT,'(A,4X,I4,2X,I4,5X,L1)') SOLVER_NAME,iPatch, & 
                                            pPatch%iPatchGlobal, & 
                                            pPatch%plotFlag
    END DO ! iPatch
  END IF ! global%verbLevel

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Setting patch plot flags done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SetPatchPlotFlags

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_SetPatchPlotFlags.F90,v $
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
! Revision 1.5  2007/02/27 13:25:04  haselbac
! Enabled 1d computations
!
! Revision 1.4  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.3  2005/10/28 19:18:51  haselbac
! Added use of postPlotPatchFlag in setting pPatch%plotFlag
!
! Revision 1.2  2005/09/21 19:42:18  haselbac
! Added writing of info
!
! Revision 1.1  2005/08/09 01:14:30  haselbac
! Initial revision
!
! ******************************************************************************

