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
! Purpose: Write Integral quantities history.
!
! Description: None.
!
! Input: 
!   global      Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_WriteInteg.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_WriteInteg(regions)

  USE ModParameters
  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModMPI
  USE RFLU_ModPlottingVars

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: regions(:)
  TYPE(t_region), POINTER :: pRegion
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,i,iReg
  REAL(RFREAL) :: KinGas,KinPcls,KinGas_l,KinPcls_l
  REAL(RFREAL) :: cntg_l,cntg,cntp_l,cntp 

#ifdef PLAG
  LOGICAL :: plagFlag
  INTEGER :: iLocUp,iLocVp,iLocWp,iLocVFp,iLocYp
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pPv
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_WriteInteg',__FILE__)

! ******************************************************************************
! Compute integral quantities
! ******************************************************************************

    KinGas_l = 0.0_RFREAL;KinGas = 0.0_RFREAL;
    cntg_l = 0.0_RFREAL;cntg = 0.0_RFREAL;
    KinPcls_l = 0.0_RFREAL;KinPcls = 0.0_RFREAL;
    cntp_l = 0.0_RFREAL;cntp = 0.0_RFREAL;

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)

      DO i = 1,pRegion%grid%nCells
       KinGas_l = KinGas_l + 0.5_RFREAL*pRegion%mixt%cv(CV_MIXT_DENS,i)*SQRT( &
         (pRegion%mixt%cv(CV_MIXT_XMOM,i)/pRegion%mixt%cv(CV_MIXT_DENS,i))**2 &
       + (pRegion%mixt%cv(CV_MIXT_YMOM,i)/pRegion%mixt%cv(CV_MIXT_DENS,i))**2 &
       + (pRegion%mixt%cv(CV_MIXT_ZMOM,i)/pRegion%mixt%cv(CV_MIXT_DENS,i))**2)
         cntg_l = cntg_l + 1.0_RFREAL
      END DO
     END DO !iReg

     CALL MPI_Allreduce(KinGas_l,KinGas,1,MPI_RFREAL,MPI_SUM, &
                         global%mpiComm,global%mpierr )
     CALL MPI_Allreduce(cntg_l,cntg,1,MPI_RFREAL,MPI_SUM, &
                         global%mpiComm,global%mpierr )

#ifdef PLAG

    CALL RFLU_CountPlottingVars(pRegion)
    CALL RFLU_CreatePlottingVarMaps(pRegion)
    CALL RFLU_BuildPlottingVarMaps(pRegion)
    CALL RFLU_CreatePlottingVars(pRegion)
    CALL RFLU_ComputePlottingVarsWrapper(pRegion)

    IF ( ASSOCIATED(pRegion%plot%pv) .EQV. .TRUE. ) THEN
      pPv => pRegion%plot%pv

    END IF ! ASSOCIATED

  plagFlag = .FALSE.

  IF ( (global%plagUsed .EQV. .TRUE.) .AND. &
       (global%postLag2EulFlag .EQV. .TRUE.) ) THEN
    iLocUp = pRegion%plot%pv2pvi(PV_PLAG_XVEL)
    iLocVp = pRegion%plot%pv2pvi(PV_PLAG_YVEL)
    iLocWp = pRegion%plot%pv2pvi(PV_PLAG_ZVEL)
    iLocYp = pRegion%plot%pv2pvi(PV_PLAG_MFRC)
    iLocVFp = pRegion%plot%pv2pvi(PV_PLAG_VFRC)
  END IF


    IF ( (iLocUp /= CRAZY_VALUE_INT) .AND. &
         (iLocVp /= CRAZY_VALUE_INT) .AND. &
         (iLocWp /= CRAZY_VALUE_INT) .AND. &
         (iLocYp /= CRAZY_VALUE_INT) .AND. &
         (iLocVFp /= CRAZY_VALUE_INT) ) THEN
      plagFlag = .TRUE.
    END IF ! iLocUp

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      DO i = 1,pRegion%grid%nCells
         KinPcls_l = KinPcls_l + 0.5_RFREAL*SQRT(pPv(iLocUp,i)**2 &
           + pPv(iLocVp,i)**2 +pPv(iLocWp,i)**2)
         cntp_l = cntp + 1.0_RFREAL
      END DO
    END DO ! iReg

    CALL MPI_Allreduce(KinPcls_l,KinPcls,1,MPI_RFREAL,MPI_SUM, &
                       global%mpiComm,global%mpierr )
    CALL MPI_Allreduce(cntp_l,cntp,1,MPI_RFREAL,MPI_SUM, &
                       global%mpiComm,global%mpierr )
#endif


! ******************************************************************************
! Write integral quantities data
! ******************************************************************************

  IF ( (global%flowType == FLOW_UNSTEADY) .AND. & 
            (global%myProcid==MASTERPROC) ) THEN

    WRITE(IF_INTEG,'(1PE12.5,2E13.4)') global%currentTime &
                   ,KinGas/MAX(1.0_RFREAL,cntg) &
                   ,KinPcls/MAX(1.0_RFREAL,cntp)
  END IF ! global%flowType

! ******************************************************************************
! Finish
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_WriteInteg

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_WriteInteg.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:48  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:57  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:01  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/04/07 15:19:22  haselbac
! Removed tabs
!
! Revision 1.1  2005/04/15 15:07:10  haselbac
! Initial revision
!
! ******************************************************************************

