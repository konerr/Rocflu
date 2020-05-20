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
! Purpose: read in user input related to formats of grid and solution file.
!
! Description: none.
!
! Input: user input file.
!
! Output: global = format of grid and solution file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ReadFormatsSection.F90,v 1.2 2015/07/23 23:11:18 brollin Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadFormatsSection( global )

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
  INTEGER, PARAMETER :: NVALS_MAX = 3

  CHARACTER(10) :: keys(NVALS_MAX)

  LOGICAL :: defined(NVALS_MAX)

  REAL(RFREAL) :: vals(NVALS_MAX)

!******************************************************************************

  CALL RegisterFunction( global,'ReadFormatsSection',__FILE__ )

! specify keywords and search for them

  nVals = NVALS_MAX
  
  keys(1) = 'GRID'
  keys(2) = 'SOLUTION'

#ifdef RFLU
  keys(3) = 'GRIDSRC'
#endif   
  
  CALL ReadSection( global,IF_INPUT,nVals,keys(1:nVals),vals(1:nVals), & 
                    defined(1:nVals) )

  IF (defined(1)) THEN
!                                       global%gridFormat = FORMAT_ASCII
!    IF (vals(1)>0.9 .AND. vals(1)<1.1) global%gridFormat = FORMAT_BINARY
!    IF (vals(1) > 1.9)                 global%gridFormat = FORMAT_HDF
  ! BBR - begin 
  ! allow for choice of endianness
                                       global%gridFormat = FORMAT_ASCII
    IF (vals(1)>0.9 .AND. vals(1)<1.1) global%gridFormat = FORMAT_BINARY
    IF (vals(1)>1.9 .AND. vals(1)<2.1) global%gridFormat = FORMAT_HDF
    IF (vals(1)>2.9 .AND. vals(1)<3.1) global%gridFormat = FORMAT_BINARY_L
    IF (vals(1)>3.9)                   global%gridFormat = FORMAT_BINARY_B
  ! BBR - end
  ENDIF
  IF (defined(2)) THEN
!                                       global%solutFormat = FORMAT_ASCII
!    IF (vals(2)>0.9 .AND. vals(2)<1.1) global%solutFormat = FORMAT_BINARY
!    IF (vals(2) > 1.9)                 global%solutFormat = FORMAT_HDF
  ! BBR - begin 
  ! allow for choice of endianness
                                       global%solutFormat = FORMAT_ASCII
    IF (vals(2)>0.9 .AND. vals(2)<1.1) global%solutFormat = FORMAT_BINARY
    IF (vals(2)>1.9 .AND. vals(2)<2.1) global%solutFormat = FORMAT_HDF
    IF (vals(2)>2.9 .AND. vals(2)<3.1) global%solutFormat = FORMAT_BINARY_L
    IF (vals(2)>3.9)                   global%solutFormat = FORMAT_BINARY_B
  ! BBR - end

  ENDIF

#ifdef RFLU
  IF (defined(3)) THEN
    global%gridSource = NINT(vals(3))
  ENDIF ! defined  
#endif 

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadFormatsSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadFormatsSection.F90,v $
! Revision 1.2  2015/07/23 23:11:18  brollin
! 1) The pressure coefficient of the  collision model has been changed back to its original form
! 2) New options in the format of the grid and solutions have been added. Now the user can choose the endianness, and convert from one to the over in rfluconv.
! 3) The solutions are now stored in folders named by timestamp or iteration number.
! 4) The address enty in the hashtable has been changed to an integer(8) for cases when the grid becomes very large.
! 5) RFLU_WritePM can now compute PM2 on the fly for the Macroscale problem
!
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
! Revision 1.1  2004/12/01 16:50:18  haselbac
! Initial revision after changing case
!
! Revision 1.11  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.8  2003/05/16 02:27:43  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.7  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.6  2003/03/15 16:23:51  haselbac
! Added KIND qualifyer
!
! Revision 1.5  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.4  2002/03/26 19:02:34  haselbac
! Added ROCFLU functionality
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

