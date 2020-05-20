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
! Purpose: read in user input related to multigrid and successive
!          grid refinement.
!
! Description: none.
!
! Input: user input file.
!
! Output: global = start level, cycle type, iterations before refinement.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ReadMultigridSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadMultigridSection( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  CHARACTER(10) :: keys(3)

  LOGICAL :: defined(3)

  REAL(RFREAL) :: vals(3)

!******************************************************************************

  CALL RegisterFunction( global,'ReadMultigridSection',__FILE__ )

! specify keywords and search for them

  keys(1) = 'START'
  keys(2) = 'CYCLE'
  keys(3) = 'REFINE'
  
  CALL ReadSection( global,IF_INPUT,3,keys,vals,defined )

  IF (defined(1)) global%startLevel = MAX(1,INT(ABS(vals(1))))
  IF (defined(2)) THEN
                                       global%cycleType = MGCYCLE_NO
    IF (vals(2)>0.9 .AND. vals(2)<1.1) global%cycleType = MGCYCLE_V
    IF (vals(2)>1.9 .AND. vals(2)<2.1) global%cycleType = MGCYCLE_W
  ENDIF
  IF (defined(3)) global%refineIter = INT(ABS(vals(3))+0.5_RFREAL)

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadMultigridSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadMultigridSection.F90,v $
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
! Revision 1.1  2004/12/01 16:50:33  haselbac
! Initial revision after changing case
!
! Revision 1.8  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/08/11 21:51:17  jblazek
! Added basic global grid smoothing scheme.
!
! Revision 1.4  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.3  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/01/11 17:18:31  jblazek
! Updated description of I/O variables.
!
! Revision 1.1  2002/01/02 16:00:03  jblazek
! Added input for multigrid parameters.
!
!******************************************************************************

