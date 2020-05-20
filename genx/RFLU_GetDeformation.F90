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
! Purpose: Obtain deformation to boundary nodes from GENX.
!
! Description: None.
!
! Input: 
!   region	grid dimensions and topology
!
! Output: None.
!
! Notes: 
!   1. Get deformation for ALL boundaries, not just coupled ones!
!
! ******************************************************************************
!
! $Id: RFLU_GetDeformation.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_GetDeformation(region)

  USE ModDataTypes
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  USE ModMPI
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region) :: region

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iPatch,ibv
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_GetDeformation.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction(global,'RFLU_GetDeformation',__FILE__)
  
! ******************************************************************************
! Loop over ALL boundaries 
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A,1X,I3)') SOLVER_NAME, &
                                   'Getting displacements from GENX...'
  END IF ! global%myProcid
  
  DO iPatch=1,region%grid%nPatches
    pPatch => region%patches(iPatch)

    DO ibv = 1,pPatch%nBVert
      pPatch%dXyz(XCOORD,ibv) = pPatch%duAlp(XCOORD,ibv)
      pPatch%dXyz(YCOORD,ibv) = pPatch%duAlp(YCOORD,ibv)
      pPatch%dXyz(ZCOORD,ibv) = pPatch%duAlp(ZCOORD,ibv)
    END DO ! ibv
    
    DO ibv = pPatch%nBVert+1,pPatch%nBVertTot
      pPatch%dXyz(XCOORD,ibv) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      pPatch%dXyz(YCOORD,ibv) = REAL(CRAZY_VALUE_INT,KIND=RFREAL) 
      pPatch%dXyz(ZCOORD,ibv) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    END DO ! ibv

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      IF ( pPatch%nBVert > 0 ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME,'Patch:',iPatch       
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Minimum/maximum values:'       
        WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'dXyz.x:', & 
              MINVAL(pPatch%dXyz(XCOORD,1:pPatch%nBVert)), & 
              MAXVAL(pPatch%dXyz(XCOORD,1:pPatch%nBVert))
        WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'dXyz.y:', & 
              MINVAL(pPatch%dXyz(YCOORD,1:pPatch%nBVert)), & 
              MAXVAL(pPatch%dXyz(YCOORD,1:pPatch%nBVert))
        WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'dXyz.z:', & 
              MINVAL(pPatch%dXyz(ZCOORD,1:pPatch%nBVert)), & 
              MAXVAL(pPatch%dXyz(ZCOORD,1:pPatch%nBVert))
      END IF ! pPatch%nBVert
    END IF ! global%myProcid
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A,1X,I3)') SOLVER_NAME, &
                                   'Getting displacements from GENX done.'
  END IF ! global%myProcid

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_GetDeformation

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_GetDeformation.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:30  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:47:51  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:58:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.8  2005/06/09 20:17:06  haselbac
! Cosmetics only
!
! Revision 1.7  2004/10/19 19:35:06  haselbac
! Cosmetics only
!
! Revision 1.6  2003/08/20 02:09:58  haselbac
! Changed verbosity conditions to reduce solver output in GENx runs
!
! Revision 1.5  2003/05/13 23:42:02  haselbac
! Added init of dummy vertices
!
! Revision 1.4  2003/04/12 21:33:22  haselbac
! Made parallel: Added MASTERPROC checks
!
! Revision 1.3  2002/10/27 18:43:21  haselbac
! Mainly cosmetic changes to output
!
! Revision 1.2  2002/10/07 14:09:32  haselbac
! Changed ModBndpatch to ModBndPatch
!
! Revision 1.1  2002/10/05 18:28:18  haselbac
! Initial revision
!
! ******************************************************************************

