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
! Purpose: Initialize Runge-Kutta scheme.
!
! Description: None.
!
! Input: 
!   region      Region data
!   iStage      Runge-Kutta stage
!   icBeg       Beginning index for cell update
!   icEnd       Ending index for cell update
!   ivBeg       Beginning index for variable update
!   ivEnd       Ending index for variable update
!   cv          Conserved variables
!   cvOld       Old conserved variables
!   diss        Dissipation residual
!
! Output: 
!   cvOld       Old conserved variables
!   diss        Dissipation residual
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RkInitGeneric.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RkInitGeneric(region,iStage,icBeg,icEnd,ivBeg,ivEnd,cv,cvOld,diss)

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

  INTEGER, INTENT(IN) :: icBeg,icEnd,iStage,ivBeg,ivEnd
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,cvOld,diss
  TYPE(t_region) :: region

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: ic,iv
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RkInitGeneric.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction(global,'RkInitGeneric',__FILE__)

! *****************************************************************************
! Initialize Runge-Kutta scheme
! *****************************************************************************

  SELECT CASE ( global%rkScheme ) 
    CASE ( RK_SCHEME_4_CLASSICAL ) 
      IF ( iStage == 1 ) THEN
        DO ic = icBeg,icEnd
          DO iv = ivBeg,ivEnd
            cvOld(iv,ic) = cv(iv,ic)
            diss(iv,ic)  = 0.0_RFREAL
          END DO ! iv
        END DO ! ic
      ELSE
        DO ic = icBeg,icEnd
          DO iv = ivBeg,ivEnd
            diss(iv,ic) = 0.0_RFREAL
          END DO ! iv      
        END DO ! ic
      END IF ! iStage
    CASE ( RK_SCHEME_3_WRAY ) 
      DO ic = icBeg,icEnd
        DO iv = ivBeg,ivEnd
          cvOld(iv,ic) = cv(iv,ic)
          diss(iv,ic)  = 0.0_RFREAL
        END DO ! iv
      END DO ! ic   
    CASE ( RK_SCHEME_1_EULER ) 
      DO ic = icBeg,icEnd
        DO iv = ivBeg,ivEnd
          cvOld(iv,ic) = cv(iv,ic)
          diss(iv,ic)  = 0.0_RFREAL
        END DO ! iv
      END DO ! ic   
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%rkScheme

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RkInitGeneric

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RkInitGeneric.F90,v $
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
! Revision 1.1  2004/12/01 16:51:02  haselbac
! Initial revision after changing case
!
! Revision 1.2  2004/11/17 16:23:40  haselbac
! Added initialization for RK3
!
! Revision 1.1  2003/11/25 21:01:50  haselbac
! Initial revision
!
!******************************************************************************

