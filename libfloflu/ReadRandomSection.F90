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
! Purpose: read in user input related to random number generation
!
! Description: none.
!
! Input: user input file.
!
! Output: global%randSeedOffset
!
! Notes:
!
! * The value of global%randSeedOffset should be 0 typically (the default).
!   It can be set to 1, 2, 3, etc. to perform runs in the same configuration,
!   but with the random number generators (in each region) seeded differently.
!
! * Example input section:
!
!-----
! # RANDOM
! SEED_OFFSET  0  ! Offset for seed of RNG (default = 0, otherwise: 1,2,3,etc)
! #
!-----
!   
!******************************************************************************
!
! $Id: ReadRandomSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadRandomSection( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  INTEGER, PARAMETER :: NKEYS_MAX = 10

  CHARACTER(CHRLEN)     :: keys(NKEYS_MAX)

  INTEGER :: iKeySeedOffset,iKeySeedType,nKeys

  LOGICAL :: defined(NKEYS_MAX)

  REAL(RFREAL) :: vals(NKEYS_MAX)

!******************************************************************************

  CALL RegisterFunction( global,'ReadRandomSection',__FILE__ )

! begin -----------------------------------------------------------------------

! define keys

  iKeySeedOffset = 1
  iKeySeedType   = 2
  nKeys = 2

  keys(iKeySeedOffset) = 'SEEDOFFSET'
  keys(iKeySeedType)   = 'SEEDTYPE'
  
  IF (nKeys > NKEYS_MAX) CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

  CALL ReadSection( global,IF_INPUT,nKeys,keys,vals,defined )

  IF (defined(iKeySeedOffset)) global%randSeedOffset = &
    ABS(NINT(vals(iKeySeedOffset)))

  IF (defined(iKeySeedType)) global%randSeedType = &
    ABS(NINT(vals(iKeySeedType)))

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE ReadRandomSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadRandomSection.F90,v $
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
! Revision 1.2  2005/12/01 17:08:07  fnajjar
! Added reading of randSeedType and code cleanup to remove underscore from OFFSET
!
! Revision 1.1  2004/12/01 16:50:46  haselbac
! Initial revision after changing case
!
! Revision 1.1  2003/11/21 22:33:10  fnajjar
! Updated Random Number Generator
!
!******************************************************************************

