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
! Purpose: read in user input related to reference values.
!
! Description: none.
!
! Input: user input file.
!
! Output: global = reference variables (viscosity set in DerivedInputValues).
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ReadReferenceSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadReferenceSection( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  INTEGER :: nVals
  INTEGER, PARAMETER :: NVALS_MAX = 11

  CHARACTER(10) :: keys(NVALS_MAX)

  LOGICAL :: defined(NVALS_MAX)

  REAL(RFREAL) :: vals(NVALS_MAX)

!******************************************************************************

  CALL RegisterFunction( global,'ReadReferenceSection',__FILE__ )

! specify keywords and search for them

  nVals = NVALS_MAX

  keys( 1) = 'ABSVEL'
  keys( 2) = 'PRESS'
  keys( 3) = 'DENS'
  keys( 4) = 'CP'
  keys( 5) = 'GAMMA'
  keys( 6) = 'LENGTH'
  keys( 7) = 'RENUM'
  keys( 8) = 'PRLAM'
  keys( 9) = 'PRTURB'
  keys(10) = 'SCNLAM'
  keys(11) = 'SCNTURB'
  
  CALL ReadSection( global,IF_INPUT,nVals,keys,vals,defined )

  IF (defined( 1)) global%refVelocity = ABS(vals( 1))
  IF (defined( 2)) global%refPressure = ABS(vals( 2))
  IF (defined( 3)) global%refDensity  = ABS(vals( 3))
  IF (defined( 4)) global%refCp       = ABS(vals( 4))
  IF (defined( 5)) global%refGamma    = ABS(vals( 5))
  IF (defined( 6)) global%refLength   = ABS(vals( 6))
  IF (defined( 7)) global%refReNum    = ABS(vals( 7))
  IF (defined( 8)) global%prLam       = ABS(vals( 8))
  IF (defined( 9)) global%prTurb      = ABS(vals( 9))
  IF (defined(10)) global%scnLam      = ABS(vals(10))
  IF (defined(11)) global%scnTurb     = ABS(vals(11))

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadReferenceSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadReferenceSection.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:32  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:16:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/11/28 23:17:44  mparmar
! Changed variable name refREnum to refReNum
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 16:50:48  haselbac
! Initial revision after changing case
!
! Revision 1.9  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.5  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.4  2002/04/11 18:46:16  haselbac
! Use parameter to specify size of keys,vals,defined
!
! Revision 1.3  2002/01/11 17:18:31  jblazek
! Updated description of I/O variables.
!
! Revision 1.2  2001/12/22 00:09:38  jblazek
! Added routines to store grid and solution.
!
! Revision 1.1  2001/12/07 16:54:32  jblazek
! Added files to read user input.
!
!******************************************************************************

