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
! Purpose: Read in user input related to the grid motion scheme.
!
! Description: None.
!
! Input: from file.
!
! Output:
!   regions = region data 
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ReadGridMotionSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReadGridMotionSection(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModInterfaces, ONLY: ReadSection
  USE ModError
  USE ModParameters

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: iReg
  INTEGER, PARAMETER :: NVALS_MAX = 3 
  INTEGER :: nVals

  CHARACTER(10) :: keys(NVALS_MAX)

  LOGICAL :: defined(NVALS_MAX)

  REAL(RFREAL) :: vals(NVALS_MAX)

! ******************************************************************************
! Start, specify keywords and search for them
! ******************************************************************************

  keys(1) = 'TYPE'
  keys(2) = 'NITER'
  keys(3) = 'SFACT'

  nVals = 3

  CALL RegisterFunction(regions(1)%global,'ReadGridMotionSection',__FILE__)

  CALL ReadSection(regions(1)%global,IF_INPUT,nVals,keys,vals,defined) 
  
  IF ( defined(1) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%moveGridType = NINT(vals(1))
    END DO ! iReg    
  END IF ! defined  
  
  IF ( defined(2) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%moveGridNIter = NINT(vals(2))
    END DO ! iReg
  END IF ! defined
  
  IF ( defined(3) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)      
      regions(iReg)%mixtInput%moveGridSFact = vals(3)
    END DO ! iReg
  END IF ! defined 

  CALL DeregisterFunction( regions(1)%global )

! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE ReadGridMotionSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadGridMotionSection.F90,v $
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
! Revision 1.18  2006/03/18 13:28:03  wasistho
! added orthDir and orthWghtX,Y,Z
!
! Revision 1.17  2006/03/08 06:38:41  wasistho
! added moveGridSiter and Viter
!
! Revision 1.16  2006/03/02 03:52:51  wasistho
! bug fixed, defined(11) was defined(10)
!
! Revision 1.15  2006/03/02 01:26:12  wasistho
! split movegrid_epde to elglobal and elframe
!
! Revision 1.14  2006/02/08 04:02:35  wasistho
! added movegrid_epde
!
! Revision 1.13  2005/10/28 22:46:34  wasistho
! added orthocell
!
! Revision 1.12  2005/10/28 05:41:48  wasistho
! read FOMS scheme
!
! Revision 1.11  2005/08/28 23:47:55  wasistho
! added orthoWght for block orthogonality of RFLO global-gridmotion
!
! Revision 1.10  2005/08/18 19:47:46  wasistho
! added moveGridNsmatch
!
! Revision 1.9  2005/06/23 05:50:21  wasistho
! changed NEIGHBOUR to NEIGHBOR
!
! Revision 1.8  2005/06/23 03:32:57  wasistho
! fixed bug in reading NEIGHBOUR
!
! Revision 1.7  2005/06/23 01:35:49  wasistho
! added input parameter NEIGHBOUR for rocflo
!
! Revision 1.6  2005/06/04 01:02:40  wasistho
! distinguished to AMPLIFX,Y,Z
!
! Revision 1.5  2005/06/02 22:57:53  wasistho
! added moveGridAmplif and moveGridPower
!
! Revision 1.4  2005/06/02 03:21:54  wasistho
! shuffle MoveGridVms with MoveGridFrame
!
! Revision 1.3  2005/05/28 21:24:16  wasistho
! added moveGridFrame
!
! Revision 1.2  2005/05/21 04:16:21  wasistho
! added MOVEGRID_VMS in Rocflo part
!
! Revision 1.1  2004/12/01 16:50:23  haselbac
! Initial revision after changing case
!
! Revision 1.13  2004/10/19 19:25:41  haselbac
! Cosmetics only
!
! Revision 1.12  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.9  2003/08/27 21:33:58  haselbac
! Changed logic so vars not overwritten if not present
!
! Revision 1.8  2003/08/25 21:51:24  jblazek
! Full version of global grid motion scheme.
!
! Revision 1.7  2003/08/11 21:51:17  jblazek
! Added basic global grid smoothing scheme.
!
! Revision 1.6  2003/05/16 02:27:43  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.5  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.4  2003/03/31 16:31:56  haselbac
! Added reading of new argument
!
! Revision 1.3  2003/03/15 16:26:10  haselbac
! Added KIND qualifyer
!
! Revision 1.2  2003/02/07 23:10:43  haselbac
! Fixed stupid mistake in setting of sFact
!
! Revision 1.1  2003/01/28 16:13:38  haselbac
! Initial revision
!
! ******************************************************************************

