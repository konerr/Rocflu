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
! Purpose: Updates solution with classical 4-stage Runge-Kutta method.
!
! Description: None.
!
! Input: 
!   region      Region data
!   varType     Variable type to be updated
!   iStage      Runge-Kutta stage
!   icBeg       Beginning index for cell update
!   icEnd       Ending index for cell update
!   ivBeg       Beginning index for variable update
!   ivEnd       Ending index for variable update
!   cv          Conserved variables
!   cvOld       Old conserved variables
!   rhs         Residual
!   rhsSum      Residual sum
!
! Output: 
!   cv          Conserved variables
!   rhsSum      Residual sum
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RkUpdateGeneric.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RkUpdateGeneric(region,varType,iStage,icBeg,icEnd,ivBeg,ivEnd, &
                           cv,cvOld,rhs,rhsSum)

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
  INTEGER, INTENT(IN) :: icBeg,icEnd,iStage,ivBeg,ivEnd,varType
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,cvOld,rhs,rhsSum

! =============================================================================
! Locals
! =============================================================================

  LOGICAL :: moveGrid
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: ic,iv
  REAL(RFREAL) :: adtv,fac,volRat  
  REAL(RFREAL) :: ark(5),grk(5)
  REAL(RFREAL), DIMENSION(:), POINTER :: vol,volOld
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RkUpdateGeneric.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction(global,'RkUpdateGeneric',__FILE__)

! *****************************************************************************
! Set pointers and variables
! *****************************************************************************

  ark(:) = region%mixtInput%ark(:)
  grk(:) = region%mixtInput%grk(:)

! =============================================================================
! Set volume(s) when updating cell-based variables
! =============================================================================

  IF ( varType == VAR_TYPE_CELL ) THEN 
    moveGrid = region%mixtInput%moveGrid

    vol => region%grid%vol

    IF ( moveGrid .EQV. .TRUE. ) THEN 
      volOld => region%gridOld%vol
    END IF ! moveGrid
  END IF ! varType

! *****************************************************************************
! Update
! *****************************************************************************

  fac = ark(iStage)*global%dtMin

  SELECT CASE ( varType ) 
  
! =============================================================================
!   Update cell-based variables, for which we need the volume (and volume 
!   ratio for moving grid computations). 
! =============================================================================
  
    CASE ( VAR_TYPE_CELL )   
        
! -----------------------------------------------------------------------------
!     Update for moving grid
! -----------------------------------------------------------------------------

      IF ( moveGrid .EQV. .TRUE. ) THEN
        SELECT CASE ( global%rkScheme ) 
          CASE ( RK_SCHEME_4_CLASSICAL ) 
            IF ( iStage == 1 ) THEN
              DO ic = icBeg,icEnd
                adtv   = fac/vol(ic)
                volRat = volOld(ic)/vol(ic)

                DO iv = ivBeg,ivEnd
                  cv(iv,ic)     = volRat*cvOld(iv,ic) - adtv*rhs(iv,ic)
                  rhsSum(iv,ic) = rhs(iv,ic)
                END DO ! iv
              END DO ! ic      
            ELSE IF ( iStage == global%nrkSteps ) THEN
              DO ic = icBeg,icEnd
                adtv   = fac/vol(ic)
                volRat = volOld(ic)/vol(ic)

                DO iv = ivBeg,ivEnd                                       
                  cv(iv,ic) = volRat*cvOld(iv,ic) & 
                            - adtv*(rhs(iv,ic)+rhsSum(iv,ic))
                END DO ! iv
              END DO ! ic      
            ELSE
              DO ic = icBeg,icEnd
                adtv   = fac/vol(ic)
                volRat = volOld(ic)/vol(ic)

                DO iv = ivBeg,ivEnd
                  cv(iv,ic)     = volRat*cvOld(iv,ic) - adtv*rhs(iv,ic)
                  rhsSum(iv,ic) = rhsSum(iv,ic) + grk(iStage)*rhs(iv,ic)
                END DO ! iv
              END DO ! ic
            END IF ! iStage
          CASE ( RK_SCHEME_3_WRAY ) 
            IF ( iStage == 1 ) THEN
              DO ic = icBeg,icEnd
                adtv   = fac/vol(ic)
                volRat = volOld(ic)/vol(ic)

                DO iv = ivBeg,ivEnd
                  cv(iv,ic)     = volRat*cvOld(iv,ic) - adtv*rhs(iv,ic)
                  rhsSum(iv,ic) = rhs(iv,ic)
                END DO ! iv
              END DO ! ic      
            ELSE IF ( iStage == 2 ) THEN
              DO ic = icBeg,icEnd
                adtv   = fac/vol(ic)
                volRat = volOld(ic)/vol(ic)

                DO iv = ivBeg,ivEnd                                       
                  cv(iv,ic)     = volRat*cvOld(iv,ic) & 
                                - adtv*(rhs(iv,ic) - grk(iStage)*rhsSum(iv,ic))
                  rhsSum(iv,ic) = rhs(iv,ic)              
                END DO ! iv
              END DO ! ic      
            ELSE
              DO ic = icBeg,icEnd
                adtv   = fac/vol(ic)
                volRat = volOld(ic)/vol(ic)

                DO iv = ivBeg,ivEnd
                  cv(iv,ic) = volRat*cvOld(iv,ic) & 
                            - adtv*(rhs(iv,ic) - grk(iStage)*rhsSum(iv,ic))
                END DO ! iv
              END DO ! ic
            END IF ! iStage      
          CASE ( RK_SCHEME_1_EULER ) 
            DO ic = icBeg,icEnd
              adtv   = fac/vol(ic)
              volRat = volOld(ic)/vol(ic)

              DO iv = ivBeg,ivEnd
                cv(iv,ic)     = volRat*cvOld(iv,ic) - adtv*rhs(iv,ic)
              END DO ! iv
            END DO ! ic      
          CASE DEFAULT 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! global%rkScheme

! =============================================================================
!     Update for non-moving grid
! =============================================================================

      ELSE 
        SELECT CASE ( global%rkScheme ) 
          CASE ( RK_SCHEME_4_CLASSICAL ) 
            IF ( iStage == 1 ) THEN
              DO ic = icBeg,icEnd
                adtv = fac/vol(ic) 

                DO iv = ivBeg,ivEnd                               
                  cv(iv,ic)     = cvOld(iv,ic) - adtv*rhs(iv,ic)
                  rhsSum(iv,ic) = rhs(iv,ic)
                END DO ! iv
              END DO ! ic
            ELSE IF ( iStage == global%nrkSteps ) THEN
              DO ic = icBeg,icEnd
                adtv = fac/vol(ic)

                DO iv = ivBeg,ivEnd
                  cv(iv,ic) = cvOld(iv,ic) - adtv*(rhs(iv,ic) + rhsSum(iv,ic))
                END DO ! iv
              END DO !ic      
            ELSE
              DO ic = icBeg,icEnd
                adtv = fac/vol(ic) 

                DO iv = ivBeg,ivEnd          
                  cv(iv,ic)     = cvOld(iv,ic) - adtv*rhs(iv,ic)
                  rhsSum(iv,ic) = rhsSum(iv,ic) + grk(iStage)*rhs(iv,ic)
                END DO ! iv
              END DO !ic      
            END IF ! iStage
          CASE ( RK_SCHEME_3_WRAY ) 
            IF ( iStage == 1 ) THEN
              DO ic = icBeg,icEnd
                adtv = fac/vol(ic) 

                DO iv = ivBeg,ivEnd                               
                  cv(iv,ic)     = cvOld(iv,ic) - adtv*rhs(iv,ic)
                  rhsSum(iv,ic) = rhs(iv,ic)
                END DO ! iv
              END DO ! ic
            ELSE IF ( iStage == 2 ) THEN
              DO ic = icBeg,icEnd
                adtv = fac/vol(ic)

                DO iv = ivBeg,ivEnd
                  cv(iv,ic)     = cvOld(iv,ic) & 
                                - adtv*(rhs(iv,ic) - grk(iStage)*rhsSum(iv,ic))
                  rhsSum(iv,ic) = rhs(iv,ic)              
                END DO ! iv
              END DO !ic      
            ELSE
              DO ic = icBeg,icEnd
                adtv = fac/vol(ic) 

                DO iv = ivBeg,ivEnd          
                  cv(iv,ic) = cvOld(iv,ic) & 
                            - adtv*(rhs(iv,ic) - grk(iStage)*rhsSum(iv,ic))
                END DO ! iv
              END DO !ic      
            END IF ! iStage      
          CASE ( RK_SCHEME_1_EULER ) 
            DO ic = icBeg,icEnd
              adtv = fac/vol(ic) 

              DO iv = ivBeg,ivEnd                               
                cv(iv,ic)     = cvOld(iv,ic) - adtv*rhs(iv,ic)
              END DO ! iv
            END DO ! ic
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! global%rkScheme
      END IF ! moveGrid

! =============================================================================
!   Update point-based variables, for which we DO NOT need the volume (and 
!   volume ratio for moving grid computations). 
! =============================================================================
  
    CASE ( VAR_TYPE_POINT )
      adtv = fac

      SELECT CASE ( global%rkScheme ) 
        CASE ( RK_SCHEME_4_CLASSICAL ) 
          IF ( iStage == 1 ) THEN
            DO ic = icBeg,icEnd
              DO iv = ivBeg,ivEnd                               
                cv(iv,ic)     = cvOld(iv,ic) - adtv*rhs(iv,ic)
                rhsSum(iv,ic) = rhs(iv,ic)
              END DO ! iv
            END DO ! ic
          ELSE IF ( iStage == global%nrkSteps ) THEN
            DO ic = icBeg,icEnd
              DO iv = ivBeg,ivEnd
                cv(iv,ic) = cvOld(iv,ic) - adtv*(rhs(iv,ic) + rhsSum(iv,ic))
              END DO ! iv
            END DO !ic      
          ELSE
            DO ic = icBeg,icEnd
              DO iv = ivBeg,ivEnd          
                cv(iv,ic)     = cvOld(iv,ic) - adtv*rhs(iv,ic)
                rhsSum(iv,ic) = rhsSum(iv,ic) + grk(iStage)*rhs(iv,ic)
              END DO ! iv
            END DO !ic      
          END IF ! iStage
        CASE ( RK_SCHEME_3_WRAY ) 
          IF ( iStage == 1 ) THEN
            DO ic = icBeg,icEnd
              DO iv = ivBeg,ivEnd                               
                cv(iv,ic)     = cvOld(iv,ic) - adtv*rhs(iv,ic)
                rhsSum(iv,ic) = rhs(iv,ic)
              END DO ! iv
            END DO ! ic
          ELSE IF ( iStage == 2 ) THEN
            DO ic = icBeg,icEnd
              DO iv = ivBeg,ivEnd
                cv(iv,ic)     = cvOld(iv,ic) & 
                              - adtv*(rhs(iv,ic) - grk(iStage)*rhsSum(iv,ic))
                rhsSum(iv,ic) = rhs(iv,ic)              
              END DO ! iv
            END DO !ic      
          ELSE
            DO ic = icBeg,icEnd
              DO iv = ivBeg,ivEnd          
                cv(iv,ic) = cvOld(iv,ic) & 
                          - adtv*(rhs(iv,ic) - grk(iStage)*rhsSum(iv,ic))
              END DO ! iv
            END DO !ic      
          END IF ! iStage      
        CASE ( RK_SCHEME_1_EULER ) 
          DO ic = icBeg,icEnd
            DO iv = ivBeg,ivEnd                               
              cv(iv,ic)     = cvOld(iv,ic) - adtv*rhs(iv,ic)
            END DO ! iv
          END DO ! ic
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! global%rkScheme
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)       
  END SELECT ! varType

! *****************************************************************************
! End
! *****************************************************************************

 CALL DeregisterFunction(global)

END SUBROUTINE RkUpdateGeneric

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RkUpdateGeneric.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
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
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 16:51:08  haselbac
! Initial revision after changing case
!
! Revision 1.2  2004/11/17 16:24:30  haselbac
! Added varType and RK3
!
! Revision 1.1  2003/11/25 21:01:50  haselbac
! Initial revision
!
!******************************************************************************

