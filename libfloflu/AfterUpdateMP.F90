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
! Purpose: check and/or modify conserved variable fields for RocfluidMP
!
! Description: none.
!
! Input: pRegion = data of current region,
!
! Output: pRegion%levels%*%cv = modified solution
!
! Notes: none.
!
!******************************************************************************
!
! $Id: AfterUpdateMP.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
!******************************************************************************

SUBROUTINE AfterUpdateMP( pRegion,istage )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters
#ifdef INRT
  USE INRT_ModParameters
  ! Subbu - Add ModInteract module
  USE ModInteract
  ! Subbu - End Add ModInteract module
#endif

  USE ModInterfaces, ONLY : RFLU_CheckPositivityWrapper, &
                            RFLU_CheckValidityWrapper,   &
                            RFLU_EnforceBoundsWrapper

#ifdef INRT
  USE ModInterfacesInteract, ONLY : INRT_SetParticleTemp,       &
                                    INRT_VaporEnergyConversion, &
                                    INRT_BurnStatusUpdate
#endif

#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY: PLAG_RFLU_ShiftUnsteadyData
  USE PLAG_ModCheckVars, ONLY: PLAG_CheckValidity
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER    :: pRegion
  INTEGER,        INTENT(IN) :: istage

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  LOGICAL :: finalStage

  TYPE(t_global), POINTER :: global

  ! Subbu - Type t_inrt_interact
#ifdef INRT
  TYPE(t_inrt_interact), POINTER :: inrt
#endif
  REAL(RFREAL) :: tBeg,tEnd
  ! Subbu - Type t_inrt_interact
!******************************************************************************

  RCSIdentString = &
    '$RCSfile: AfterUpdateMP.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction( global,'AfterUpdateMP',__FILE__ )

  finalStage = (istage == global%nrkSteps)

! check positivity ------------------------------------------------------------

  CALL RFLU_CheckValidityWrapper(pRegion)
  CALL RFLU_EnforceBoundsWrapper(pRegion)
  CALL RFLU_CheckPositivityWrapper(pRegion)

#ifdef PLAG
  IF (global%plagUsed .AND. global%inrtUsed .AND. finalStage) THEN
    ! Subbu - Turn of PLAG_RFLU_ShiftUnsteady data for viscous force computation
    !IF (1==2) THEN 
    inrt => pRegion%inrtInput%inrts(INRT_TYPE_DRAG)
    !IF(inrt%switches(INRT_SWI_DRAG_UNSTEADY) == INRT_DRAG_UNSTEADY_USE) THEN
    IF(inrt%switches(INRT_SWI_DRAG_VU) == INRT_DRAG_VISCUNST_USE) THEN
      CALL PLAG_RFLU_ShiftUnsteadyData(pRegion)
    END IF
    !END IF
    ! Subbu - End Turn of PLAG_RFLU_ShiftUnsteady data for viscous force computation
  END IF ! plagUsed
#endif  

#ifdef INRT
  IF (global%inrtUsed .AND. finalStage) THEN

    IF (pRegion%inrtInput%inrts(INRT_TYPE_BURNING)%used) THEN
      CALL INRT_SetParticleTemp(pRegion)
      CALL INRT_VaporEnergyConversion(pRegion)
    END IF ! burning Used

    IF (pRegion%inrtInput%inrts(INRT_TYPE_BURNING)%used) THEN
      CALL INRT_BurnStatusUpdate(pRegion)
    END IF ! burning used

  END IF ! inrtUsed
#endif

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE AfterUpdateMP

!******************************************************************************
!
! RCS Revision history:
!
! $Log: AfterUpdateMP.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:24  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/03/27 15:17:26  haselbac
! Fixed screwed-up last check-in
!
! Revision 1.3  2006/03/26 20:21:07  haselbac
! Changed to wrappers bcos of GL model
!
! Revision 1.2  2005/12/01 21:50:11  fnajjar
! Added call to PLAG_CheckValidity
!
! Revision 1.1  2004/12/01 16:47:41  haselbac
! Initial revision after changing case
!
! Revision 1.5  2004/08/04 00:26:49  wasistho
! call rflo_checkvalidity only at finalstage
!
! Revision 1.4  2004/07/26 19:02:36  wasistho
! add RFLO_CheckValidity
!
! Revision 1.3  2004/03/25 21:14:20  jferry
! changed AfterUpdate to call most subroutines only after final RK stage
!
! Revision 1.2  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.1  2004/03/02 21:47:28  jferry
! Added After Update interactions
!
!******************************************************************************

