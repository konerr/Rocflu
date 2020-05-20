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
! Purpose: Wrapper routine for moving grid.
!
! Description: None. 
!
! Input: 
!   regions     Region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_MoveGridWrapper.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_MoveGridWrapper(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModDataStruct, ONLY: t_region  
  USE ModMPI
  USE ModParameters

#ifdef GENX
  USE RFLU_ModGENXTools, ONLY: RFLU_GENX_MoveGrid, &
                               RFLU_GENX_ConstrainDisp
#endif

  USE ModInterfaces, ONLY: RFLU_MoveGridDisp, & 
                           RFLU_MoveGridXyz

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global 
  TYPE(t_region), POINTER :: pRegion
  INTEGER :: iReg

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_MoveGridWrapper.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_MoveGridWrapper',__FILE__)

! ******************************************************************************
! Call appropriate grid motion routine
! ******************************************************************************

#ifdef GENX
  DO iReg = 1,global%nRegionsLocal
     pRegion => regions(iReg)
     IF(global%cnstrCaseRad > 0.0) THEN
        CALL RFLU_GENX_ConstrainDisp(pRegion)
     ENDIF
  END DO ! iReg
#endif
  IF ( regions(1)%mixtInput%moveGridType == MOVEGRID_TYPE_DISP ) THEN 
    CALL RFLU_MoveGridDisp(regions)  
  ELSE IF ( regions(1)%mixtInput%moveGridType == MOVEGRID_TYPE_XYZ ) THEN 
    CALL RFLU_MoveGridXyz(regions,MOVEGRID_CONTEXT_MOVESMOOTH)
#ifdef GENX    
  ELSE IF ( regions(1)%mixtInput%moveGridType == MOVEGRID_TYPE_GENX ) THEN 
    CALL RFLU_GENX_MoveGrid(regions)
#endif    
  ELSE 
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END IF ! regions(1)%mixtInput%moveGridType

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_MoveGridWrapper

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_MoveGridWrapper.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2009/08/28 18:29:48  mtcampbe
! RocfluMP integration with Rocstar and some makefile tweaks.  To build
! Rocstar with new Rocflu:
! make ROCFLU=RocfluMP
! To build Rocstar with the new RocfluND:
! make ROCFLU=RocfluMP HYPRE=/the/hypre/install/path
!
! Revision 1.3  2008/12/06 08:43:48  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:57  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:01  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2004/10/19 19:29:22  haselbac
! Added GENX grid motion option, cosmetics
!
! Revision 1.1  2003/03/31 16:05:39  haselbac
! Initial revision
!
! ******************************************************************************

