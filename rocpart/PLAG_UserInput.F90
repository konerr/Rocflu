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
! Purpose: read user input and store it in the data structure.
!          Check user input.
!
! Description: none.
!
! Input: regions = dimensions and topology (finest grid).
!
! Output: regions = dimensions, topology and user input on all grid levels.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_UserInput.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_UserInput( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE PLAG_ModInterfaces, ONLY : PLAG_CheckUserInput, &
                                 PLAG_DerivedInputValues, &
                                 PLAG_InitInputValues, &
                                 PLAG_ReadInputFile
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE (t_region), POINTER :: regions(:)

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString, msg

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_UserInput.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction( global, 'PLAG_UserInput',__FILE__ )

! Initialize parameters -------------------------------------------------------

  CALL PLAG_InitInputValues( regions )

! Read user input (global%plagUsed is set here) -------------------------------

  CALL PLAG_ReadInputFile( regions )

  IF (global%plagUsed) THEN

! Set derived user input ------------------------------------------------------

    CALL PLAG_DerivedInputValues( regions )

! Check user input ------------------------------------------------------------

    CALL PLAG_CheckUserInput( regions )

  END IF ! plagUsed

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_UserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_UserInput.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:52  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:50:27  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:34  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 20:58:19  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.4  2004/02/26 21:02:25  haselbac
! Added call to set derived input values
!
! Revision 1.3  2003/11/21 22:42:16  fnajjar
! Added plagActive
!
! Revision 1.2  2003/04/14 18:58:30  fnajjar
! Added PLAG_InitInputValues call
!
! Revision 1.1  2002/10/25 14:20:32  f-najjar
! Initial Import of Rocpart
!
!******************************************************************************

