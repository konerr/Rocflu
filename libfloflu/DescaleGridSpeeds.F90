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
! Purpose: Descale grid speeds.
!
! Description: None. 
!
! Input: 
!   region     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Grid speeds need to be descaled because they are scaled during the 
!      Runge-Kutta stages, but for consistent restarts, they need to be 
!      written out consistent with their values when they were used for the 
!      time-step computation.
!
!******************************************************************************
!
! $Id: DescaleGridSpeeds.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE DescaleGridSpeeds( region )

  USE ModDataTypes
  USE ModGrid, ONLY      : t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY  : t_patch
  USE ModGlobal, ONLY    : t_global
  USE ModError
  USE ModParameters

  USE RFLU_ModGrid

  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region) :: region

! =============================================================================
! Locals
! =============================================================================

  INTEGER :: ifc, iPatch, irk, irkStep
  REAL(RFREAL) :: scaleFactor, term
  REAL(RFREAL) :: ark(5), grk(5)
  REAL(RFREAL), DIMENSION(:), POINTER :: gs
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  global => region%global

  CALL RegisterFunction(global,'DescaleGridSpeeds',__FILE__)
  
! *****************************************************************************
! Set pointers and variables
! *****************************************************************************

  irkStep = region%irkStep

  ark(:) = region%mixtInput%ark(:)
  grk(:) = region%mixtInput%grk(:)

! *****************************************************************************
! Descale grid speeds
! *****************************************************************************

! =============================================================================
! Determine scaling factor (compare with ScaleGridSpeeds.F90)
! ============================================================================= 

  term = 0.0_RFREAL

  DO irk = 1,global%nrkSteps-1
    term = term + grk(irk)/ark(irk)
  END DO ! irk

  scaleFactor = 1.0_RFREAL/(1.0_RFREAL/ark(global%nrkSteps) - term)

! =============================================================================
! Interior faces
! ============================================================================= 

  gs => region%grid%gs

  DO ifc = 1,region%grid%nFaces
    gs(ifc) = scaleFactor*gs(ifc)
  END DO ! ifc

! =============================================================================
! Patch faces
! ============================================================================= 

  DO iPatch = 1,region%grid%nPatches
    pPatch => region%patches(iPatch)
    gs     => pPatch%gs

    DO ifc = 1,pPatch%nBFaces
      gs(ifc) = scaleFactor*gs(ifc)      
    END DO ! ifc
  END DO ! iPatch

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE DescaleGridSpeeds

!******************************************************************************
!
! RCS Revision history:
!
! $Log: DescaleGridSpeeds.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:24  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 16:48:29  haselbac
! Initial revision after changing case
!
! Revision 1.7  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.4  2003/10/01 23:52:09  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.3  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.2  2003/03/15 16:14:30  haselbac
! Changed loop limit
!
! Revision 1.1  2003/02/25 21:42:24  haselbac
! Initial revision
!
!******************************************************************************

