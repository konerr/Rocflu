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
! Purpose: set initial solution field for first RK stage.
!
! Description: none.
!
! Input: region = data of current region,
!        istage  = current RK stage.
!
! Output: region%levels%*%cvOld,diss = initializations for new time step.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RkInitMP.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RKInitMP( region, istage )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModMixture,    ONLY : t_mixt
  USE ModSpecies,    ONLY : t_spec
  USE ModError
  USE ModParameters
  USE ModBndPatch, ONLY: t_patch 

  USE ModInterfaces, ONLY : RkInitGeneric, &
                            RkInitPointScalar, &
                            RkInitSD

  USE RFLU_ModNSCBC, ONLY: RFLU_NSCBC_DecideHaveNSCBC

#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY : PLAG_RkInit
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region
  INTEGER,        INTENT(IN)            :: istage

! ... local variables
  INTEGER :: ibc, iec

  LOGICAL :: moveGrid, specUsed, plagUsed

  TYPE(t_mixt),   POINTER :: mixt
  TYPE(t_spec),   POINTER :: spec
  TYPE(t_global), POINTER :: global

  INTEGER :: iPatch
  TYPE(t_patch), POINTER :: pPatch

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RKInitMP',__FILE__ )

! set flags ===================================================================

  plagUsed = global%plagUsed
  specUsed = global%specUsed

  moveGrid = region%mixtInput%moveGrid

! get dimensions and pointers -------------------------------------------------

  ibc = 1
  iec = region%grid%nCellsTot
  mixt => region%mixt
  spec => region%spec

! initialize for flow itself and all multiphysics modules ---------------------

  CALL RkInitGeneric(region,istage,ibc,iec,1,CV_MIXT_NEQS, &
                     mixt%cv,mixt%cvOld,mixt%diss)

! initialize particle velocity for moving reference frame formulation  --------
  IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
    CALL RkInitPointScalar(region,istage,1,3,region%mvfVel,region%mvfVelOld)

    CALL RkInitPointScalar(region,istage,1,3,region%mvfLoc,region%mvfLocOld)
  END IF ! global%mvFrameFlag

! loop over patches in a region and initialize boundary variables -------------
  DO iPatch = 1,region%grid%nPatches
    pPatch => region%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      CALL RkInitGeneric(region,istage,1,pPatch%nBFaces,1,CV_MIXT_NEQS, &
                         pPatch%mixt%cv,pPatch%mixt%cvOld,pPatch%mixt%rhs)
    END IF ! pPatch%bcKind
  END DO ! iPatch

#ifdef SPEC
  IF ( specUsed ) THEN
    CALL RkInitGeneric(region,istage,ibc,iec,1,region%specInput%nSpecies, &
                       spec%cv,spec%cvOld,spec%diss)
  END IF ! specUsed
#endif
  SELECT CASE ( region%mixtInput%indSd )

  CASE (0)
    CALL RkInitSD(region,  0,  1,1,SD_ZMOM-SD_XMOM+1,mixt%sd)

  CASE (1)
    CALL RkInitSD(region,ibc,iec,1,SD_ZMOM-SD_XMOM+1,mixt%sd)

  CASE DEFAULT
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)

  END SELECT ! indSd

#ifdef PLAG
  IF ( plagUsed ) THEN
    CALL PLAG_RkInit(region,istage)
  END IF ! plagUsed
#endif


! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE rkInitMP

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RkInitMP.F90,v $
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
! Revision 1.2  2007/06/18 17:40:56  mparmar
! Added initialization of moving reference frame variables
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/08/19 15:38:31  mparmar
! Added initialization of Runge-Kutta scheme for boundary arrays
!
! Revision 1.3  2006/02/13 21:00:51  wasistho
! added ifdef PEUL
!
! Revision 1.2  2005/03/31 16:30:04  haselbac
! Replaced nSd by parameters
!
! Revision 1.1  2004/12/01 16:51:04  haselbac
! Initial revision after changing case
!
! Revision 1.19  2004/07/30 22:47:33  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.18  2004/07/28 15:29:18  jferry
! created global variable for spec use
!
! Revision 1.17  2004/07/23 22:43:15  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.16  2004/04/14 02:06:42  haselbac
! No longer call ScaleGridSpeeds for RFLU
!
! Revision 1.15  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.14  2004/03/03 23:55:08  jferry
! Made module calls more uniform
!
! Revision 1.13  2004/03/02 21:50:30  jferry
! Changed rkInit and rkUpdate routines to call generic procedures
!
! Revision 1.12  2004/02/26 21:13:06  wasistho
! changed TURB_ransRkInit to TURB_rkInit
!
! Revision 1.11  2004/02/26 21:01:43  haselbac
! Removed ifdef RFLO around PLAG_rkInitMP
!
! Revision 1.10  2004/02/02 22:48:21  haselbac
! Added ifdef RFLO - temporary measure
!
! Revision 1.9  2003/11/25 21:01:43  haselbac
! Added rocspecies support with rkInitGeneric routine
!
! Revision 1.8  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/10/03 20:12:53  wasistho
! initial installation of turbModel SA and DES
!
! Revision 1.4  2003/10/01 23:52:09  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.3  2003/09/26 21:45:34  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.2  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.1  2003/03/28 19:46:59  fnajjar
! Initial import for RocfluidMP
!
!******************************************************************************

