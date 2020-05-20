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
! Purpose: Read in user input related to mixture.
!
! Description: None.
!
! Input: 
!   regions     Data for regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ReadMixtureSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReadMixtureSection(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  USE ModDataStruct, ONLY: t_region  
  
  USE ModInterfaces, ONLY: ReadSection  
  
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

  INTEGER :: iReg,nVals

#ifdef RFLU
  INTEGER, PARAMETER :: NVALS_MAX = 2
#endif

  CHARACTER(10) :: keys(NVALS_MAX)
  LOGICAL :: defined(NVALS_MAX)
  REAL(RFREAL) :: vals(NVALS_MAX)
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction(global,'ReadMixtureSection',__FILE__)

! specify keywords and search for them

  nVals = NVALS_MAX

#ifdef RFLU
  keys(1) = 'FROZENFLAG'
  keys(2) = 'GASMODEL'     
#endif

#ifdef RFLU
  CALL ReadSection(global,IF_INPUT,nVals,keys,vals,defined)

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)  
    IF ( defined(1) .EQV. .TRUE. ) THEN 
      IF ( NINT(vals(1)) == 1 ) THEN 
        regions(iReg)%mixtInput%frozenFlag = .TRUE.
      ELSE 
        regions(iReg)%mixtInput%frozenFlag = .FALSE.
      END IF ! NINT(vals(1))
    END IF ! defined               
  END DO ! iReg

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)  
    IF ( defined(2) .EQV. .TRUE. ) THEN 
      regions(iReg)%mixtInput%gasModel = NINT(vals(2))
    END IF ! defined               
  END DO ! iReg
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE ReadMixtureSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadMixtureSection.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.2  2005/11/10 22:20:06  fnajjar
! ACH: Added frozenFlag
!
! Revision 1.1  2005/10/31 19:23:53  haselbac
! Initial revision
!
! ******************************************************************************

