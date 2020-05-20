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
! Purpose: Allocate memory for Lagrangian particles sources 
!          related to time stepping.
!
! Description: none.
! Input: 
!   pRegion     Region pointer
!   pPlag       Pointer to particle data structure
!
! Output: None.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_INRT_AllocMemTStep.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_INRT_AllocMemTStep( pRegion,pPlag )

  USE ModDataTypes
  USE ModError
  USE ModParameters  
  USE ModGlobal, ONLY: t_global 
  USE ModDataStruct, ONLY: t_region
  USE ModPartLag, ONLY: t_plag  
  USE ModMPI
   
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion
  TYPE(t_plag), POINTER :: pPlag

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iPcl,iVar,maxDisEdges,nPclsMax,errorFlag

  TYPE(t_global),        POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_INRT_AllocMemTStep.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_INRT_AllocMemTStep',__FILE__)

! ******************************************************************************
! Set pointers and variables 
! ******************************************************************************

  nPclsMax = pRegion%plag%nPclsMax
  maxDisEdges = pRegion%inrtInput%maxDisEdges

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  IF (maxDisEdges > 0) THEN
    ALLOCATE( pPlag%inrtSources(maxDisEdges,nPclsMax),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN 
      CALL ErrorStop( global,ERR_ALLOCATE,__LINE__,'pPlag%inrtSources')
    END IF ! global%error 

  ELSE
    NULLIFY( pPlag%inrtSources )

  END IF ! maxDisEdges

! ******************************************************************************
! Initialize memory
! ******************************************************************************

  DO iPcl = 1,nPclsMax
    DO iVar = 1, maxDisEdges
      pPlag%inrtSources(iVar,iPcl) = 0.0_RFREAL
    END DO ! iVar
  END DO ! iPcl
    
! ******************************************************************************
! End
! ******************************************************************************

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_INRT_AllocMemTStep

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_INRT_AllocMemTStep.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:50  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:50:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.6  2007/03/06 23:13:13  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.5  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.4  2004/08/18 15:19:20  fnajjar
! Bug fix for updated definition of nPclsTot
!
! Revision 1.3  2004/07/26 18:47:05  fnajjar
! Cleaned up header comments to properly reflect calling sequence
!
! Revision 1.2  2004/07/26 18:43:16  fnajjar
! Removed temporary debugging io
!
! Revision 1.1  2004/07/26 17:06:54  fnajjar
! initial import
!
!******************************************************************************

