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
! Purpose: read in user input related to rockets
!
! Description: none.
!
! Input: user input file.
!
! Output: global = constraint params
!
! Notes: none.
!
!******************************************************************************


SUBROUTINE ReadRocketSection( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
  USE ModMPI

  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  INTEGER :: nVals
  INTEGER, PARAMETER :: NVALS_MAX = 15

  CHARACTER(10) :: keys(NVALS_MAX)
  LOGICAL :: hasACase
  LOGICAL :: defined(NVALS_MAX)
  REAL(RFREAL) :: vals(NVALS_MAX)

!******************************************************************************

  CALL RegisterFunction( global,'ReadRocketSection',__FILE__ )
! specify keywords and search for them

  nVals = NVALS_MAX

  keys(1)  =  'CASERAD'
  keys(2)  =  'HEADEND'
  keys(3)  =  'AFTEND'
  keys(4)  =  'COORDL'
  keys(5)  =  'TOL1'
  keys(6)  =  'TOL2'
  keys(7)  =  'ELLIPSL'
  keys(8)  =  'ELLIPST'
  keys(9)  =  'NOZY'
  keys(10) =  'LMIN'
  keys(11) =  'LMAX'
  keys(12) =  'T1MIN'
  keys(13) =  'T1MAX'
  keys(14) =  'T2MIN'
  keys(15) =  'T2MAX'
  
  CALL ReadSection( global,IF_INPUT,nVals,keys,vals,defined )

  IF (defined( 1)) THEN
     hasACase = .TRUE.
     global%cnstrCaseRad        = vals( 1)
  ELSE
     global%cnstrCaseRad = -999999_RFREAL
  ENDIF

  IF (defined( 2)) THEN
     global%cnstrHeadEnd = vals( 2)
  ELSE
     global%cnstrHeadEnd = -9999999_RFREAL
  ENDIF

  IF (defined( 3)) THEN
     global%cnstrAftEnd = vals( 3)
  ELSE
     global%cnstrAftEnd = 9999999_RFREAL
  ENDIF

  IF (defined( 4)) THEN
     IF(vals( 4) == 1.0) THEN
        global%cnstrCoordL  = XCOORD
        global%cnstrCoordT1 = YCOORD
        global%cnstrCoordT2 = ZCOORD
     ELSE IF (vals(4) == 2.0) THEN
        global%cnstrCoordL  = YCOORD
        global%cnstrCoordT1 = ZCOORD
        global%cnstrCoordT2 = XCOORD
     ELSE
        global%cnstrCoordL  = ZCOORD
        global%cnstrCoordT1 = XCOORD
        global%cnstrCoordT2 = YCOORD
     ENDIF
  ELSE
     global%cnstrCoordL  = XCOORD
     global%cnstrCoordT1 = YCOORD
     global%cnstrCoordT2 = ZCOORD
  ENDIF

  IF (defined( 5)) THEN
     global%cnstrTol1 = vals( 5)
  ELSE
     global%cnstrTol1 = 0.000001_RFREAL
  ENDIF
  IF (defined( 6)) THEN
     global%cnstrTol2 = vals(6)
  ELSE
     global%cnstrTol2 = 0.00000001_RFREAL
  ENDIF

  IF (defined( 7)) THEN
     global%cnstrEllipsL   = vals(7)
  ELSE
     global%cnstrEllipsL = 1.0_RFREAL
  ENDIF

  IF (defined( 8)) THEN
     global%cnstrEllipsT   = vals(8)
  ELSE
     global%cnstrEllipsT = global%cnstrCaseRad
  ENDIF

  IF (defined( 9)) THEN
     global%cnstrNozY   = vals(9)
  ELSE
     global%cnstrNozY = 1.0_RFREAL
  ENDIF
  IF (defined( 10)) THEN
     global%cnstrLMinPlane   = vals(10)
  ELSE
     global%cnstrLMinPlane = -9999999_RFREAL
  ENDIF
  IF (defined( 11)) THEN
     global%cnstrLMaxPlane   = vals(11)
  ELSE
     global%cnstrLMaxPlane = 9999999_RFREAL
  ENDIF
  IF (defined( 12)) THEN
     global%cnstrT1MinPlane   = vals(12)
  ELSE
     global%cnstrT1MinPlane = -9999999_RFREAL
  ENDIF
  IF (defined( 13)) THEN
     global%cnstrT1MaxPlane   = vals(13)
  ELSE
     global%cnstrT1MaxPlane = 9999999_RFREAL
  ENDIF
  IF (defined( 14)) THEN
     global%cnstrT2MinPlane   = vals(14)
  ELSE
     global%cnstrT2MinPlane = -9999999_RFREAL
  ENDIF
  IF (defined( 15)) THEN
     global%cnstrT2MaxPlane   = vals(15)
  ELSE
     global%cnstrT2MaxPlane = 9999999_RFREAL
  ENDIF

     
! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadRocketSection


