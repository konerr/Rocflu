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
! Purpose: Check input values for species.
!
! Description: None.
!
! Input: 
!   regions             Data associated with regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: SPEC_CheckUserInput.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_CheckUserInput(regions)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters
          
  USE ModSpecies, ONLY: t_spec_type          
            
  USE ModInterfaces, ONLY: MixtPerf_R_M            
              
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iReg,iSpec
  REAL(RFREAL) :: cp,g,gc
  TYPE(t_global), POINTER :: global
  TYPE(t_spec_type), POINTER :: pSpecType  

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_CheckUserInput.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'SPEC_CheckUserInput',__FILE__)

! ******************************************************************************
! Check input values
! ******************************************************************************

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)

! ==============================================================================
!   Gas properties
! ==============================================================================

    DO iSpec = 1,regions(iReg)%specInput%nSpecies
      pSpecType => regions(iReg)%specInput%specType(iSpec)
      
      gc = MixtPerf_R_M(pSpecType%pMaterial%molw)
      cp = pSpecType%pMaterial%spht
      g  = cp/(cp-gc)
      
      IF ( (g < 1.0_RFREAL) .OR. (g > 4.0_RFREAL) ) THEN 
        CALL ErrorStop(global,ERR_SPEC_PROPS_INVALID,__LINE__)
      END IF ! pSpecType%pMaterial%spht
    END DO ! iSpec

! ==============================================================================
!   EEM velocity
! ==============================================================================

    SELECT CASE ( regions(iReg)%mixtInput%indSd )
      CASE ( 0 )
      CASE ( 1 )
        IF ( regions(iReg)%mixtInput%computeTv .EQV. .FALSE. ) THEN
          CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__, &
                         'Attempting to use EEM without viscosity')
        END IF ! regions(iReg)%mixtInput%computeTv

        IF ( regions(iReg)%mixtInput%moveGrid .EQV. .TRUE. ) THEN
          CALL ErrorStop(global,ERR_UNKNOWN_OPTION,__LINE__, &
                         'EEM not yet implemented with moving grid')
        END IF ! regions(iReg)%mixtInput%moveGrid
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! regions(iReg)%mixtInput%indSd
  END DO ! iReg

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_CheckUserInput

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_CheckUserInput.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:53  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:51:22  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.3  2005/11/10 02:33:38  haselbac
! Clean-up, added checks for gas properties
!
! Revision 1.2  2004/07/30 22:47:37  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
! ******************************************************************************

