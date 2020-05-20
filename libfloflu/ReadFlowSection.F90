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
! Purpose: read in user input related to flow model.
!
! Description: none.
!
! Input: user input file.
!
! Output: regions = flow model, moving grid.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ReadFlowSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadFlowSection( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  INTEGER :: nVals

  INTEGER :: iReg
  INTEGER, PARAMETER :: NVALS_MAX = 3

  CHARACTER(10) :: keys(NVALS_MAX)

  LOGICAL :: defined(NVALS_MAX)

  REAL(RFREAL) :: vals(NVALS_MAX)

!******************************************************************************

  CALL RegisterFunction( regions(1)%global,'ReadFlowSection',__FILE__ )

! specify keywords and search for them

  nVals = NVALS_MAX

  keys(1) = 'MODEL'
  keys(2) = 'MOVEGRID'
  keys(3) = 'PBAFLAG'

  CALL ReadSection( regions(1)%global,IF_INPUT,nVals,keys(1:nVals),vals(1:nVals), &
                    defined(1:nVals) )
  
  IF ( defined(1) ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%flowModel = NINT(vals(1))
    END DO ! iReg
  END IF ! defined
  
  IF ( defined(2) ) THEN
    IF ( NINT(vals(2)) == 0 ) THEN
      DO iReg = LBOUND(regions,1),UBOUND(regions,1)      
        regions(iReg)%mixtInput%moveGrid = .FALSE.
      END DO ! iReg
    ELSEIF ( NINT(vals(2)) == 1 ) THEN      
      DO iReg = LBOUND(regions,1),UBOUND(regions,1)
        regions(iReg)%mixtInput%moveGrid = .TRUE.
      END DO ! iReg
    END IF ! NINT
  END IF ! defined 

  IF ( defined(3) ) THEN
    IF ( NINT(vals(3)) == 0 ) THEN
      DO iReg = LBOUND(regions,1),UBOUND(regions,1)      
        regions(iReg)%global%pbaFlag = .FALSE.
      END DO ! iReg
    ELSEIF ( NINT(vals(3)) == 1 ) THEN      
      DO iReg = LBOUND(regions,1),UBOUND(regions,1)
        regions(iReg)%global%pbaFlag     = .TRUE.
        regions(iReg)%global%pbaBurnFlag = .TRUE.
      END DO ! iReg
    END IF ! NINT
  END IF ! defined 

! finalize

  CALL DeregisterFunction( regions(1)%global )

END SUBROUTINE ReadFlowSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadFlowSection.F90,v $
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
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 16:50:12  haselbac
! Initial revision after changing case
!
! Revision 1.16  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.13  2003/05/16 02:27:43  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.12  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.11  2003/03/15 16:22:59  haselbac
! Added KIND qualifyers
!
! Revision 1.10  2002/10/16 21:09:43  haselbac
! Fixed bug in RFLU code segment
!
! Revision 1.9  2002/09/09 14:01:16  haselbac
! mixtInput now under regions
!
! Revision 1.8  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.7  2002/03/26 19:01:46  haselbac
! Added ROCFLU functionality
!
! Revision 1.6  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.5  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.4  2002/01/23 03:51:24  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.3  2002/01/11 17:18:31  jblazek
! Updated description of I/O variables.
!
! Revision 1.2  2001/12/22 00:09:38  jblazek
! Added routines to store grid and solution.
!
! Revision 1.1  2001/12/07 16:54:31  jblazek
! Added files to read user input.
!
!******************************************************************************

