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
! Purpose: update conserved variable fields for RocfluidMP framework.
!
! Description: none.
!
! Input: region = data of current region,
!        iReg   = index of currect region,
!        istage = RK current stage.
!
! Output: region%levels%*%cv,rhs = new solution after one time step.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RkUpdateMP.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RKUpdateMP( region,iReg,istage )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModMixture,    ONLY : t_mixt
  USE ModSpecies,    ONLY : t_spec
  USE ModError
  USE ModParameters
  USE ModBndPatch, ONLY: t_patch

  USE RFLU_ModNSCBC, ONLY: RFLU_NSCBC_DecideHaveNSCBC

  USE ModInterfaces, ONLY: RkUpdateGeneric, &
                           RkUpdatePointScalar

#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY : PLAG_RkUpdateWrapper
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region
  INTEGER, INTENT(IN) :: iReg, istage

! ... local variables
  INTEGER :: ibc, iec

  LOGICAL :: plagUsed, specUsed

  TYPE(t_mixt),   POINTER :: mixt
  TYPE(t_spec),   POINTER :: spec
  TYPE(t_global), POINTER :: global

  INTEGER :: iPatch
  TYPE(t_patch), POINTER :: pPatch

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RKUpdateMP',__FILE__ )

! set flags ===================================================================

  plagUsed = global%plagUsed
  specUsed = global%specUsed

  region%irkStep = istage

! get dimensions and pointers -------------------------------------------------

  ibc = 1
  iec = region%grid%nCellsTot
  mixt => region%mixt
  spec => region%spec

! update solution for flow itself and all multiphysics modules ----------------

  IF ( region%mixtInput%frozenFlag .EQV. .FALSE. ) THEN 
! DEBUG: Manoj: 2012-04-24: why explosive region rhs is not zero          
IF (1==2) THEN
    WRITE(*,*) "=========================================================="
    WRITE(*,*) "time= ",global%currentTimeRK
    WRITE(*,'(A,3(2X,E24.8))') "311: ",spec%rhs(1,311), spec%rhs(2,311), spec%rhs(3,311)
    WRITE(*,'(A,3(2X,E24.8))') "411: ",spec%rhs(1,411), spec%rhs(2,411), spec%rhs(3,411)
    WRITE(*,'(A,3(2X,E24.8))') "511: ",spec%rhs(1,511), spec%rhs(2,511), spec%rhs(3,511)
    WRITE(*,'(A,3(2X,E24.8))') "611: ",spec%rhs(1,611), spec%rhs(2,611), spec%rhs(3,611)
    WRITE(*,*) "------------------------"
    WRITE(*,'(A,3(2X,E24.8))') "311: ",mixt%rhs(1,311), mixt%rhs(2,311), mixt%rhs(5,311)
    WRITE(*,'(A,3(2X,E24.8))') "411: ",mixt%rhs(1,411), mixt%rhs(2,411), mixt%rhs(5,411)
    WRITE(*,'(A,3(2X,E24.8))') "511: ",mixt%rhs(1,511), mixt%rhs(2,511), mixt%rhs(5,511)
    WRITE(*,'(A,3(2X,E24.8))') "611: ",mixt%rhs(1,611), mixt%rhs(2,611), mixt%rhs(5,611)
    WRITE(*,*) "------------------------"
    WRITE(*,'(A,3(2X,E24.8))') "312: ",spec%rhs(1,312), spec%rhs(2,312), spec%rhs(3,312)
    WRITE(*,'(A,3(2X,E24.8))') "412: ",spec%rhs(1,412), spec%rhs(2,412), spec%rhs(3,412)
    WRITE(*,'(A,3(2X,E24.8))') "512: ",spec%rhs(1,512), spec%rhs(2,512), spec%rhs(3,512)
    WRITE(*,'(A,3(2X,E24.8))') "612: ",spec%rhs(1,612), spec%rhs(2,612), spec%rhs(3,612)
    WRITE(*,*) "------------------------"
    WRITE(*,'(A,3(2X,E24.8))') "312: ",mixt%rhs(1,312), mixt%rhs(2,312), mixt%rhs(5,312)
    WRITE(*,'(A,3(2X,E24.8))') "412: ",mixt%rhs(1,412), mixt%rhs(2,412), mixt%rhs(5,412)
    WRITE(*,'(A,3(2X,E24.8))') "512: ",mixt%rhs(1,512), mixt%rhs(2,512), mixt%rhs(5,512)
    WRITE(*,'(A,3(2X,E24.8))') "612: ",mixt%rhs(1,612), mixt%rhs(2,612), mixt%rhs(5,612)
    WRITE(*,*) "=========================================================="
END IF ! 1==2
!    WRITE(*,*) "Stopping here ..."
!    STOP
! END DEBUG
    CALL RkUpdateGeneric(region,VAR_TYPE_CELL,istage,ibc,iec,1,CV_MIXT_NEQS, &
                         mixt%cv,mixt%cvOld,mixt%rhs,mixt%rhsSum)
  END IF ! region%mixtInput%frozenFlag

! update particle velocity for moving reference frame formulation  ------------
  IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
    IF ( global%mvfAccFlag .EQV. .FALSE. ) THEN
      CALL RkUpdatePointScalar(region,istage,1,3,region%mvfVel, &
                               region%mvfVelOld,region%mvfAcc,region%mvfAccSum)
    END IF ! global%mvfAccFlag

    CALL RkUpdatePointScalar(region,istage,1,3,region%mvfLoc,region%mvfLocOld, &
                             region%mvfVel,region%mvfVelSum)
  END IF ! global%mvFrameFlag

! loop over patches in a region and update boundary variables -------------
  DO iPatch = 1,region%grid%nPatches
    pPatch => region%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      CALL RkUpdateGeneric(region,VAR_TYPE_POINT,istage,1,pPatch%nBFaces, &
                           1,CV_MIXT_NEQS,pPatch%mixt%cv,pPatch%mixt%cvOld, &
                           pPatch%mixt%rhs,pPatch%mixt%rhsSum)

    END IF ! pPatch%bcKind
  END DO ! iPatch

#ifdef SPEC
  IF ( specUsed ) THEN
    CALL RkUpdateGeneric(region,VAR_TYPE_CELL,istage,ibc,iec,1, & 
                         region%specInput%nSpecies,spec%cv,spec%cvOld, &
                         spec%rhs,spec%rhsSum)
  END IF ! specUsed
#endif

#ifdef PLAG
  IF ( plagUsed ) THEN
    CALL PLAG_RkUpdateWrapper(region,iReg,istage)

  END IF ! plagUsed
#endif

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE RKUpdateMP

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RkUpdateMP.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.5  2008/12/06 08:43:32  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:16:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2008/05/29 01:35:09  mparmar
! Removed Setting of reference frame velocity and update of BC from here
!
! Revision 1.2  2007/06/18 17:41:30  mparmar
! Added calls to update moving reference frame variables
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.7  2006/08/19 15:38:33  mparmar
! Added update of Runge-Kutta scheme for boundary arrays
!
! Revision 1.6  2006/02/13 21:00:58  wasistho
! added ifdef PEUL
!
! Revision 1.5  2005/11/10 22:20:38  fnajjar
! ACH: Added IF on frozenFlag
!
! Revision 1.4  2005/04/06 02:16:49  wasistho
! mv call to PERI_CoMeanCorrection to UpdateBoundaryConditionsMP
!
! Revision 1.3  2005/03/11 04:22:57  wasistho
! commented PERI_coMeanCorrection temporarily while in testing
!
! Revision 1.2  2005/03/07 05:05:04  wasistho
! install hybrid DESSA turbulence model
!
! Revision 1.1  2004/12/01 16:51:11  haselbac
! Initial revision after changing case
!
! Revision 1.21  2004/11/17 23:44:27  wasistho
! used generic RK-update for rocturb
!
! Revision 1.20  2004/11/17 16:25:11  haselbac
! Added varType to interface of RkUpdateGeneric
!
! Revision 1.19  2004/07/28 15:29:18  jferry
! created global variable for spec use
!
! Revision 1.18  2004/07/23 22:43:15  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.17  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.16  2004/03/03 23:55:08  jferry
! Made module calls more uniform
!
! Revision 1.15  2004/03/02 21:50:30  jferry
! Changed rkInit and rkUpdate routines to call generic procedures
!
! Revision 1.14  2004/02/26 21:01:45  haselbac
! Removed ifdef RFLO around PLAG_rkUpdateWrapper
!
! Revision 1.13  2004/02/02 22:48:32  haselbac
! Added ifdef RFLO - temporary measure
!
! Revision 1.12  2003/11/25 21:01:44  haselbac
! Added rocspecies support with rkUpdateGeneric routine
!
! Revision 1.11  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.8  2003/09/26 21:45:34  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.7  2003/08/28 20:32:10  wasistho
! excluced ModInterfacesTurbulence,Radiation,Periodic from ModInterfaces
!
! Revision 1.6  2003/08/14 01:46:30  wasistho
! fixed ifdef around TURB_solutionUpdate
!
! Revision 1.5  2003/08/06 15:53:09  wasistho
! added vorticities computation
!
! Revision 1.4  2003/07/08 21:21:36  jblazek
! Modified start up procedure for dual-time stepping.
!
! Revision 1.3  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.2  2003/04/05 01:58:45  wasistho
! install ROCPERI
!
! Revision 1.1  2003/03/28 19:47:19  fnajjar
! Initial import for RocfluidMP
!
!******************************************************************************

