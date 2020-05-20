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
! Purpose: set dummycells to zero.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: regions%levels%mixt%rhs = residual in dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ZeroDummyCellsMP.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ZeroDummyCellsMP( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE ModInterfaces, ONLY: RFLU_ZeroVirtualCellVars

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), TARGET :: region

! ... local variables
  LOGICAL :: specUsed

  REAL(RFREAL), POINTER :: rhs(:,:), rhsPeul(:,:)

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'ZeroDummyCellsMP',__FILE__ )

! get dimensions and pointers =================================================

  specUsed = global%specUsed

  pRegion => region

! zero out residuals in dummy cells -------------------------------------------

  CALL RFLU_ZeroVirtualCellVars(pRegion,pRegion%mixt%rhs)

#ifdef SPEC
  IF ( specUsed .EQV. .TRUE. ) THEN
    CALL RFLU_ZeroVirtualCellVars(pRegion,pRegion%spec%rhs)
  END IF ! specUsed
#endif

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE ZeroDummyCellsMP

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ZeroDummyCellsMP.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2006/01/12 02:38:30  wasistho
! added modelClass condition in turb
!
! Revision 1.2  2005/05/16 20:40:32  haselbac
! Renamed RFLU_ZeroDummyCells, now use pRegion
!
! Revision 1.1  2004/12/01 16:52:32  haselbac
! Initial revision after changing case
!
! Revision 1.14  2004/07/28 15:29:18  jferry
! created global variable for spec use
!
! Revision 1.13  2004/03/19 02:40:15  wasistho
! renamed TURB_RFLO_RansZeroDummyCells to TURB_RansZeroDummyCells
!
! Revision 1.12  2004/03/11 03:33:06  wasistho
! changed rocturb nomenclature
!
! Revision 1.11  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.10  2004/03/02 21:50:44  jferry
! Corrected typo
!
! Revision 1.9  2003/11/25 21:01:46  haselbac
! Added support for rocspecies, changed TURB call
!
! Revision 1.8  2003/11/20 16:40:36  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/10/16 20:14:10  wasistho
! turbulence zero dummy cells computed through one routine
!
! Revision 1.4  2003/10/03 20:13:21  wasistho
! initial installation of turbModel SA and DES
!
! Revision 1.3  2003/10/01 23:52:10  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.2  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.1  2003/03/28 19:48:44  fnajjar
! Initial import for RocfluidMP
!
!******************************************************************************

