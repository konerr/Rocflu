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
! Purpose: read in user input related to viscosity model.
!
! Description: none.
!
! Input: user input file.
!
! Output: regions = flow model, moving grid.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ReadViscositySection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadViscositySection( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  INTEGER, PARAMETER :: NKEYS = 4
  CHARACTER(10) :: keys(NKEYS)
  INTEGER :: iReg
 
  LOGICAL :: defined(NKEYS)
 
  REAL(RFREAL) :: vals(NKEYS)

!******************************************************************************

  CALL RegisterFunction( regions(1)%global,'ReadViscositySection',__FILE__ )

! specify keywords and search for them

  keys(1) = 'MODEL'
  keys(2) = 'VISCOSITY'
  keys(3) = 'REFTEMP'
  keys(4) = 'SUTHCOEF'

  CALL ReadSection( regions(1)%global,IF_INPUT,NKEYS,keys,vals,defined ) 
  
  IF ( defined(1) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%viscModel = NINT(vals(1))
    END DO ! iReg
  END IF ! defined
  
  IF ( defined(2) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)      
      regions(iReg)%mixtInput%refVisc = ABS(vals(2))
    END DO ! iReg
  ELSE 
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)      
      regions(iReg)%mixtInput%refVisc = REAL(CRAZY_VALUE_INT,RFREAL)
    END DO ! iReg    
  END IF ! defined 
  
  IF ( defined(3) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)      
      regions(iReg)%mixtInput%refViscTemp = ABS(vals(3))
    END DO ! iReg
  END IF ! defined 
  
  IF ( defined(4) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)      
      regions(iReg)%mixtInput%suthCoef = ABS(vals(4))
    END DO ! iReg
  END IF ! defined 

! finalize

  CALL DeregisterFunction( regions(1)%global )

END SUBROUTINE ReadViscositySection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadViscositySection.F90,v $
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
! Revision 1.2  2007/11/28 23:17:47  mparmar
! Changed variable name from refTemp to refViscTemp
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 16:50:59  haselbac
! Initial revision after changing case
!
! Revision 1.8  2003/12/04 03:23:04  haselbac
! Added parameter, setting of refVisc if undefined, fixed bug
!
! Revision 1.7  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.4  2003/05/20 22:16:16  jblazek
! Corrected bug in viscosity model input.
!
! Revision 1.3  2003/05/16 02:27:43  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.2  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.1  2003/04/10 23:28:27  fnajjar
! Initial import for viscosity models
!
!******************************************************************************

