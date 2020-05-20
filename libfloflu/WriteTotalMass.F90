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
! Purpose: Write total mass and related info to file.
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
!******************************************************************************
!
! $Id: WriteTotalMass.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE WriteTotalMass(regions)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModMPI
  USE ModParameters
  USE ModDataStruct, ONLY: t_region
  IMPLICIT NONE

! ... parameters

  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction(global,'WriteTotalMass',__FILE__)

! steady flow -----------------------------------------------------------------

  IF ( global%flowType == FLOW_STEADY .AND. & 
       global%myProcid == MASTERPROC ) THEN
    WRITE(IF_MASS,'(I6,4(1X,E23.16))') global%currentIter,global%totalMass, &
                                       global%massIn,global%massOut, & 
                                       global%totalVol

! unsteady flow ---------------------------------------------------------------

  ELSE IF ( global%flowType == FLOW_UNSTEADY .AND. & 
            global%myProcid == MASTERPROC ) THEN
    WRITE(IF_MASS,'(5(1X,E23.16))') global%currentTime,global%totalMass, &
                                    global%massIn,global%massOut, & 
                                    global%totalVol
  END IF ! global%flowType

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction(global)

END SUBROUTINE WriteTotalMass

!******************************************************************************
!
! RCS Revision history:
!
! $Log: WriteTotalMass.F90,v $
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
! Revision 1.2  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 16:52:30  haselbac
! Initial revision after changing case
!
! Revision 1.4  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.3  2002/11/15 21:29:33  haselbac
! Deleted RFLU stuff (moved elsewhere), now only write
!
! Revision 1.2  2002/11/15 14:09:25  haselbac
! Changed output format
!
! Revision 1.1  2002/11/08 21:55:48  haselbac
! Initial revision
!
!******************************************************************************

