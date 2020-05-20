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
! Purpose: Read in user input related to calculation of forces.
!
! Description: None.
!
! Input: 
!    global     Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ReadForcesSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReadForcesSection( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================   
! Locals
! ==============================================================================    

  INTEGER, PARAMETER :: NVALS_MAX = 7  
  LOGICAL :: defined(NVALS_MAX)
  CHARACTER(10) :: keys(NVALS_MAX)
  INTEGER :: nVals   
  REAL(RFREAL) :: vals(NVALS_MAX)

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction( global,'ReadForcesSection',__FILE__ )

! ******************************************************************************
! Specify keywords and search for them
! ******************************************************************************

  nVals = NVALS_MAX

  keys(1) = 'FLAG'
  keys(2) = 'REFLENGTH'
  keys(3) = 'REFAREA'
  keys(4) = 'REFXCOORD'
  keys(5) = 'REFYCOORD'
  keys(6) = 'REFZCOORD'
  keys(7) = 'PATCHFLAG'    
  
  CALL ReadSection(global,IF_INPUT,nVals,keys,vals,defined)

! ******************************************************************************
! Set variables
! ******************************************************************************

  IF ( defined(1) .EQV. .TRUE. ) THEN 
    IF ( NINT(vals(1)) == 1 ) THEN 
      global%forceFlag = .TRUE.
    ELSE 
      global%forceFlag = .FALSE.
    END IF ! NINT
  END IF ! defined

  IF ( defined(2) .EQV. .TRUE. ) THEN 
    global%forceRefLength = vals(2)
  END IF ! defined
  
  IF ( defined(3) .EQV. .TRUE. ) THEN 
    global%forceRefArea = vals(3)
  END IF ! defined
  
  IF ( defined(4) .EQV. .TRUE. ) THEN 
    global%forceRefXCoord = vals(4)
  END IF ! defined 
  
  IF ( defined(5) .EQV. .TRUE. ) THEN 
    global%forceRefYCoord = vals(5)
  END IF ! defined 
  
  IF ( defined(6) .EQV. .TRUE. ) THEN 
    global%forceRefZCoord = vals(6)
  END IF ! defined  
  
  IF ( defined(7) .EQV. .TRUE. ) THEN
    IF ( NINT(vals(1)) == 1 ) THEN 
      global%patchCoeffFlag = .TRUE.
    ELSE 
      global%patchCoeffFlag = .FALSE.
    END IF ! NINT
  END IF ! defined        

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE ReadForcesSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadForcesSection.F90,v $
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
! Revision 1.5  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.4  2006/03/10 01:46:33  wasistho
! read acBndBox coordinates for Rocflo
!
! Revision 1.3  2006/03/09 20:47:48  wasistho
! prepared for aerodyn.coeffs calc in RFLO
!
! Revision 1.2  2005/08/09 00:52:39  haselbac
! Added reading of PATCHFLAG
!
! Revision 1.1  2004/12/01 16:50:15  haselbac
! Initial revision after changing case
!
! Revision 1.9  2004/06/16 20:00:13  haselbac
! Added RFLU code
!
! Revision 1.8  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.4  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
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
! ******************************************************************************

