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
! Purpose: Read in user input related to preprocessor.
!
! Description: None.
!
! Input: 
!   global      Pointer to global type
!
! Output: None.
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: ReadPrepSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReadPrepSection(global)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
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

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER, PARAMETER :: NVALS_MAX = 1  

  LOGICAL :: defined(NVALS_MAX)
  CHARACTER(10) :: keys(NVALS_MAX)
  CHARACTER(CHRLEN) :: RCSIdentString    
  INTEGER :: nVals  
  REAL(RFREAL) :: vals(NVALS_MAX)

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: ReadPrepSection.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction(global,'ReadPrepSection',__FILE__)

! ******************************************************************************
! Specify keywords and search for them
! ******************************************************************************

#ifdef RFLU
  nVals = NVALS_MAX

  keys(1) = 'PARTMODE'

  CALL ReadSection(global,IF_INPUT,nVals,keys,vals,defined)

  IF ( defined(1) .EQV. .TRUE. ) THEN
    IF ( NINT(vals(1)) == PARTITION_MODE_PROPER ) THEN 
      global%prepPartMode = PARTITION_MODE_PROPER
    ELSE 
      global%prepPartMode = PARTITION_MODE_IMPOSED
    END IF ! NINT(vals)
  END IF ! defined
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE ReadPrepSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadPrepSection.F90,v $
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
! Revision 1.2  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 16:50:41  haselbac
! Initial revision after changing case
!
! Revision 1.7  2004/10/19 19:25:47  haselbac
! Removed SURFFLAG, cosmetics
!
! Revision 1.6  2003/08/07 15:29:01  haselbac
! Changed var names
!
! Revision 1.5  2003/05/16 02:27:43  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.4  2003/05/07 00:21:40  haselbac
! Changed logic for check
!
! Revision 1.3  2003/05/01 14:06:40  haselbac
! Fixed bug: Did not update properly
!
! Revision 1.2  2003/04/29 21:47:37  haselbac
! Added SURFFLAG
!
! Revision 1.1  2003/04/29 14:55:48  haselbac
! Initial revision
!
! ******************************************************************************

