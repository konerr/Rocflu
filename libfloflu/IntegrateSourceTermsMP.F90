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
! Purpose: Integrate source terms using operator-splitting method.
!
! Description: None.
!
! Input: 
!   regions             Data for all grid regions
!
! Output: None.
!
! Notes: 
!   1. This routine is work in progress... 
!
! ******************************************************************************
!
! $Id: IntegrateSourceTermsMP.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE IntegrateSourceTermsMP(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters

#ifdef RFLU
#ifdef SPEC
  USE SPEC_RFLU_ModChemistry, ONLY: SPEC_RFLU_IntegrateChemSrcTerm
  USE ModInterfaces, ONLY: RFLU_SetVarsWrapper
#endif  
#endif

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

  TYPE(t_region), POINTER :: regions(:)
#ifdef RFLU
  TYPE(t_region), POINTER :: pRegion  
#endif

! ==============================================================================
! Arguments
! ==============================================================================

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
#ifdef RFLU
  INTEGER :: iReg
#endif  
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: IntegrateSourceTermsMP.F90,v $ $Revision: 1.1.1.1 $'

! ******************************************************************************
! Set pointers and variables
! ****************************************************************************** 

  global => regions(1)%global

  CALL RegisterFunction(global,'IntegrateSourceTermsMP',__FILE__)

#ifdef RFLU
#ifdef SPEC  
! ******************************************************************************
! Integrate chemistry source term and update dependent variables
! ******************************************************************************  

! TEMPORARY - At present, do not use operator-split integration of chemistry
!             source terms based on modifications by Luca. This means that 
!             call source term residual function directly in SourceTermsMP.F90
!  DO iReg = 1,global%nRegionsLocal
!    pRegion => regions(iReg)
!      
!    IF ( pRegion%specInput%sourceFlag .EQV. .TRUE. ) THEN 
!      CALL SPEC_RFLU_IntegrateChemSrcTerm(pRegion)
!      CALL RFLU_SetVarsWrapper(pRegion,1,pRegion%grid%nCellsTot)
!    END IF ! pRegion%specInput%sourceType
!  END DO ! iReg
! END TEMPORARY
#endif
#endif

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE IntegrateSourceTermsMP

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: IntegrateSourceTermsMP.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.3  2005/06/06 14:22:09  haselbac
! Adapted to Lucas changes
!
! Revision 1.2  2005/04/15 15:06:02  haselbac
! Adapted call to RFLU_SetVarsWrapper
!
! Revision 1.1  2004/12/01 16:48:43  haselbac
! Initial revision after changing case
!
! Revision 1.3  2004/11/14 20:30:23  haselbac
! Bug fix: Moved USE inside ifdef RFLU
!
! Revision 1.2  2004/11/14 19:35:06  haselbac
! Replaced call to UpdateDependentVarsMP by RFLU_SetVarsWrapper, cosmetics
!
! Revision 1.1  2004/04/01 21:22:17  haselbac
! Initial revision
!
! ******************************************************************************

