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
! Purpose: find new values for TBCs and reset BCs accordingly.
!
! Description: none.
!
! Input: region index, substep time and timestep, and whether at final substep
!
! Output: modifies TBC and BC data in regions 
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: UpdateTbc.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE UpdateTbc(region,t,dt,final)

  USE ModDataTypes
  USE ModBndPatch
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  
  USE ModInterfaces, ONLY : UpdateTbcPiecewise, &
                            UpdateTbcSinusoidal, &
                            UpdateTbcStochastic, &
                            UpdateTbcWhitenoise
                            
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  LOGICAL, INTENT(IN) :: final
  REAL(RFREAL), INTENT(IN) :: t, dt
  TYPE(t_region), INTENT(INOUT) :: region

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: finiteDt
  INTEGER :: iPatch, n, ftype
  REAL(RFREAL), POINTER :: vals(:,:)
  TYPE(t_patch),  POINTER :: patch
  TYPE(t_tbcvalues), POINTER :: tbc, tbcs(:)
  TYPE(t_bcvalues), POINTER :: bc
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => region%global

  CALL RegisterFunction(global,'UpdateTbc',__FILE__)

  finiteDt = (dt > 1.E-10_RFREAL*global%dtMin)

! ******************************************************************************
! Loop over patches
! ******************************************************************************

  DO iPatch=1,region%grid%nPatches
    patch => region%patches(iPatch)

! ==============================================================================
!   Loop over modules
! ==============================================================================

    DO ftype = 1,FTYPE_MAX
      SELECT CASE ( ftype )
        CASE ( FTYPE_MIXT )
          bc => patch%mixt
#ifdef SPEC
        CASE ( FTYPE_SPEC )
          IF ( global%specUsed .EQV. .FALSE. ) THEN
            CYCLE
          END IF ! global%specUsed
          bc => patch%spec
#endif
        CASE DEFAULT
          CYCLE
      END SELECT ! ftype

      IF ( ASSOCIATED(bc%vals) .NEQV. .TRUE. ) THEN 
        CYCLE
      ELSE
        vals => bc%vals
      END IF ! ASSOCIATED
      
! ------------------------------------------------------------------------------
!     Loop over variables
! ------------------------------------------------------------------------------

      DO n = 1,bc%nData
        tbc => bc%tbcs(n)

        IF ( tbc%tbcType /= TBC_NONE ) THEN

! ------- Check if TBC is on ---------------------------------------------------

          IF ( t <= tbc%params(TBCDAT_ONTIME) ) THEN
            vals(n,:) = tbc%mean
          ELSE IF ( t >= tbc%params(TBCDAT_OFFTIME) ) THEN
            tbc%tbcType = TBC_NONE
            vals(n,:) = tbc%mean
          ELSE

! --------- Update TBC data

            SELECT CASE ( tbc%tbcType )
              CASE ( TBC_SINUSOIDAL )
                IF (finiteDt) THEN
                  CALL UpdateTbcSinusoidal(global,tbc,t)
                  vals(n,:) = tbc%mean*(1.0_RFREAL + tbc%svals(TBCSTO_VAL))
                END IF
              CASE ( TBC_STOCHASTIC )
                IF ( finiteDt ) THEN
                  CALL UpdateTbcStochastic(region,tbc,dt)
                  vals(n,:) = tbc%mean*tbc%bvals(TBCSTO_FACTOR,:)
                END IF
              CASE ( TBC_WHITENOISE )
                IF ( (tbc%switches(TBCSWI_SUBSTEP) == TBCOPT_STEP .AND. final) .OR. &
                     (tbc%switches(TBCSWI_SUBSTEP) == TBCOPT_SUBSTEP .AND. finiteDt) ) THEN
                  CALL UpdateTbcWhitenoise(region,tbc)
                  vals(n,:) = tbc%mean*(1.0_RFREAL + tbc%bvals(TBCSTO_VAL,:))
                END IF
              CASE ( TBC_PIECEWISE )
                CALL UpdateTbcPiecewise(global,tbc,t)
                vals(n,:) = tbc%mean*tbc%svals(TBCSTO_VAL)
              CASE DEFAULT
                CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END SELECT ! tbc%tbctype 
          END IF ! t 
        END IF ! tbc%tbctype 
      END DO ! n
    END DO ! ftype
  END DO ! iPatch

! ******************************************************************************
! Finalize
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE UpdateTbc

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: UpdateTbc.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:33  mtcampbe
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
! Revision 1.6  2006/08/19 15:38:36  mparmar
! Renamed patch variables
!
! Revision 1.5  2006/05/20 19:10:25  fnajjar
! Fixed bug to include a CYCLE statement when specUsed is not active
!
! Revision 1.4  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.3  2006/01/07 04:49:36  wasistho
! cycled turb and radi
!
! Revision 1.2  2005/04/27 02:07:23  haselbac
! Cosmetics only
!
! Revision 1.1  2004/12/01 16:52:03  haselbac
! Initial revision after changing case
!
! Revision 1.14  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.11  2003/06/21 20:39:36  haselbac
! Added ifdefs as workaround for IBM problems, removd finiteDt if for pw linear
!
! Revision 1.10  2003/06/10 22:54:42  jferry
! Added Piecewise TBC
!
! Revision 1.9  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.8  2003/02/17 19:31:12  jferry
! Implemented portable random number generator ModRandom
!
! Revision 1.7  2002/10/12 19:10:30  haselbac
! Added check for association status of vals (does not work otherwise)
!
! Revision 1.6  2002/09/27 00:57:09  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.5  2002/09/25 18:29:57  jferry
! simplified TBC parameter lists
!
! Revision 1.4  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.3  2002/09/18 21:50:49  jferry
! Streamlined inelegant coding
!
! Revision 1.2  2002/09/18 15:25:30  jferry
! Fixed RFLU compilation bug
!
! Revision 1.1  2002/09/17 13:42:59  jferry
! Added Time-dependent boundary conditions
!
! ******************************************************************************

