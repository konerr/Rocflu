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
! Purpose: Impose patch coordinates.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: 
!   1. This routine is of interest if the grid is known not to satisfy exactly
!      certain coordinate values, which may affect mass conservation tests
!      because the volume is not conserved. 
!   2. This routine will have to be hard-coded for each case.
!   3. Routine coded for flexibility, not efficiency...
!
!******************************************************************************
!
! $Id: RFLU_USER_EnforcePatchCoords.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_USER_EnforcePatchCoords(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters

  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region), POINTER :: pRegion

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: ibv,iPatch,ivg
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_USER_EnforcePatchCoords.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_USER_EnforcePatchCoords', &
                        __FILE__)

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Enforcing patch coordinates...'
  END IF ! global%myProcid

! *****************************************************************************
! Set pointers
! *****************************************************************************

  pGrid => pRegion%grid

! *****************************************************************************
! Select case-dependent boundary patch deformation
! *****************************************************************************

  SELECT CASE ( TRIM(global%casename) ) 
    
! -----------------------------------------------------------------------------
!   Endburner problem (old, edge length 0.1)
! -----------------------------------------------------------------------------

    CASE ( "endburner3pt","endburner5pt","endburner9pt" )      
      DO iPatch=1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

        IF ( global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME,'Patch:',iPatch       
        END IF ! global%myProcid 

        DO ibv = 1,pPatch%nBVert
          ivg = pPatch%bv(ibv)

          IF ( iPatch == 1 ) THEN 
            pGrid%xyz(ZCOORD,ivg) = 0.0_RFREAL
          ELSE IF ( iPatch == 2 ) THEN 
            pGrid%xyz(XCOORD,ivg) = 0.0_RFREAL
          ELSE IF ( iPatch == 3 ) THEN 
            pGrid%xyz(YCOORD,ivg) = 0.0_RFREAL
          ELSE IF ( iPatch == 4 ) THEN 
            pGrid%xyz(YCOORD,ivg) = 0.1_RFREAL
          ELSE IF ( iPatch == 5 ) THEN 
            pGrid%xyz(ZCOORD,ivg) = 0.1_RFREAL
          ELSE IF ( iPatch == 6 ) THEN 
            pGrid%xyz(XCOORD,ivg) = 0.1_RFREAL
          ELSE 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END IF ! iPatch          
        END DO ! ibv
      END DO ! iPatch

! -----------------------------------------------------------------------------
!   Endburner problem (new, edge length 1.0)
! -----------------------------------------------------------------------------

    CASE ( "endburner3ptnew" )    
      DO iPatch=1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

        IF ( global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME,'Patch:',iPatch       
        END IF ! global%myProcid 

        DO ibv = 1,pPatch%nBVert
          ivg = pPatch%bv(ibv)

          IF ( iPatch == 1 ) THEN 
            pGrid%xyz(ZCOORD,ivg) = 1.0_RFREAL
          ELSE IF ( iPatch == 2 ) THEN 
            pGrid%xyz(YCOORD,ivg) = 0.0_RFREAL
          ELSE IF ( iPatch == 3 ) THEN 
            pGrid%xyz(XCOORD,ivg) = 1.0_RFREAL
          ELSE IF ( iPatch == 4 ) THEN 
            pGrid%xyz(XCOORD,ivg) = 0.0_RFREAL
          ELSE IF ( iPatch == 5 ) THEN 
            pGrid%xyz(YCOORD,ivg) = 1.0_RFREAL
          ELSE IF ( iPatch == 6 ) THEN 
            pGrid%xyz(ZCOORD,ivg) = 0.0_RFREAL
          ELSE 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END IF ! iPatch          
        END DO ! ibv
      END DO ! iPatch

! -----------------------------------------------------------------------------
!   Endburner problem (angled, edge length 0.1)
! -----------------------------------------------------------------------------

    CASE ( "endburner3pt_angled" )      
      DO iPatch=1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

        IF ( global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME,'Patch:',iPatch       
        END IF ! global%myProcid 

        DO ibv = 1,pPatch%nBVert
          ivg = pPatch%bv(ibv)

          IF ( iPatch == 1 ) THEN 
            pGrid%xyz(ZCOORD,ivg) = 0.0_RFREAL
          ELSE IF ( iPatch == 2 ) THEN 
            pGrid%xyz(XCOORD,ivg) = 0.0_RFREAL
          ELSE IF ( iPatch == 3 ) THEN 
            pGrid%xyz(YCOORD,ivg) = 0.0_RFREAL
          ELSE IF ( iPatch == 4 ) THEN 
            pGrid%xyz(YCOORD,ivg) = 0.1_RFREAL
          ELSE IF ( iPatch == 5 ) THEN 
            pGrid%xyz(ZCOORD,ivg) = 0.1_RFREAL
          ELSE IF ( iPatch == 6 ) THEN 
            pGrid%xyz(XCOORD,ivg) = 0.1_RFREAL
          ELSE 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END IF ! iPatch          
        END DO ! ibv
      END DO ! iPatch

! ----------------------------------------------------------------------------
!   Default - reaching this here is no error
! ----------------------------------------------------------------------------

    CASE DEFAULT
      IF ( global%verbLevel > VERBOSE_LOW ) THEN
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Nothing to be done.'
      END IF ! global%myProcid       
  END SELECT ! global%casename

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A,1X,A)') SOLVER_NAME,'Enforcing patch coordinates', &
                                  'done.'
  END IF ! global%myProcid

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_USER_EnforcePatchCoords

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_USER_EnforcePatchCoords.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:56  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:10  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:56:52  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2005/04/15 15:09:19  haselbac
! Initial revision
!
! Revision 1.3  2003/03/31 16:22:20  haselbac
! Added CASE for endburner3pt_angled
!
! Revision 1.2  2003/03/20 20:02:32  haselbac
! Modified RegFun call to avoid probs with
! long __FILE__ names
!
! Revision 1.1  2003/02/01 01:14:42  haselbac
! Initial revision
! 
!******************************************************************************

