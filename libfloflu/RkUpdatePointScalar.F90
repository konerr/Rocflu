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
! Purpose: Updates particle velocity with classical 4-stage Runge-Kutta method.
!
! Description: None.
!
! Input: 
!   region	Region data
!   iStage	Runge-Kutta stage
!   ivBeg	Beginning index for variable update
!   ivEnd	Ending index for variable update
!   var		Conserved variables
!   varOld	Old conserved variables
!   rhs		Residual
!   rhsSum	Residual sum
!
!
! Output: 
!   var		Variables
!   rhsSum	Residual sum
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RkUpdatePointScalar.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RkUpdatePointScalar(region,iStage,ivBeg,ivEnd,var,varOld,rhs,rhsSum)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  
  IMPLICIT NONE

! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region) :: region
  INTEGER, INTENT(IN) :: iStage,ivBeg,ivEnd
  REAL(RFREAL), DIMENSION(:), POINTER :: var,varOld,rhs,rhsSum

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iv
  REAL(RFREAL) :: fac
  REAL(RFREAL) :: ark(5),grk(5)
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RkUpdatePointScalar.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction(global,'RkUpdatePointScalar',__FILE__)

! *****************************************************************************
! Set pointers and variables
! *****************************************************************************

  ark(:) = region%mixtInput%ark(:)
  grk(:) = region%mixtInput%grk(:)

! *****************************************************************************
! Update
! *****************************************************************************

  fac = ark(iStage)*global%dtMin

  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
    fac = ark(iStage)*global%dtMin*global%refVelocity/global%refLength
  END IF ! global%solverType

  SELECT CASE ( global%rkScheme ) 
  CASE ( RK_SCHEME_4_CLASSICAL ) 
    IF ( iStage == 1 ) THEN
      DO iv = ivBeg,ivEnd                               
        var(iv)     = varOld(iv) + fac*rhs(iv)
        rhsSum(iv)  = rhs(iv)
      END DO ! iv
    ELSE IF ( iStage == global%nrkSteps ) THEN
      DO iv = ivBeg,ivEnd                               
        var(iv)     = varOld(iv) + fac*(rhs(iv) + rhsSum(iv))
      END DO ! iv
    ELSE
      DO iv = ivBeg,ivEnd                               
        var(iv)     = varOld(iv) + fac*rhs(iv)
        rhsSum(iv)  = rhsSum(iv) + grk(iStage)*rhs(iv)
      END DO ! iv
    END IF ! iStage
  CASE ( RK_SCHEME_3_WRAY ) 
    IF ( iStage == 1 ) THEN
      DO iv = ivBeg,ivEnd                               
        var(iv)     = varOld(iv) + fac*rhs(iv)
        rhsSum(iv)  = rhs(iv)
      END DO ! iv
    ELSE IF ( iStage == 2 ) THEN
      DO iv = ivBeg,ivEnd                               
        var(iv)     = varOld(iv) + fac*(rhs(iv) - grk(iStage)*rhsSum(iv))
        rhsSum(iv)  = rhs(iv)
      END DO ! iv
    ELSE
      DO iv = ivBeg,ivEnd                               
        var(iv)     = varOld(iv) + fac*(rhs(iv) - grk(iStage)*rhsSum(iv))
      END DO ! iv
    END IF ! iStage      
  CASE ( RK_SCHEME_1_EULER ) 
    DO iv = ivBeg,ivEnd                               
      var(iv)     = varOld(iv) + fac*rhs(iv)
      rhsSum(iv)  = rhs(iv)
    END DO ! iv
  CASE DEFAULT
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%rkScheme

! *****************************************************************************
! End
! *****************************************************************************

 CALL DeregisterFunction(global)

END SUBROUTINE RkUpdatePointScalar

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RkUpdatePointScalar.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2009/07/08 19:11:26  mparmar
! Non-dimensionalized time for SOLV_IMPLICIT_HM
!
! Revision 1.3  2008/12/06 08:43:32  mtcampbe
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
! Revision 1.1  2006/08/19 15:37:27  mparmar
! Initial revision
!
! ******************************************************************************

