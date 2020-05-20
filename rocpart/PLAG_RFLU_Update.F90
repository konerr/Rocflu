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
! Purpose: Update step for PLAG module.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
!   iStage      Current RK stage
!
! Output: None.
!
! Notes: 
!   1. This corresponds to Part IV Steps 10d and Part V in RocfluidMP 
!      framework description.
!
! ******************************************************************************
!
! $Id: PLAG_RFLU_Update.F90,v 1.2 2015/12/18 23:15:41 rahul Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_Update(pRegion,iStage)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModPartLag, ONLY: t_plag
  USE ModError
  USE ModMPI

  USE ModParameters
  USE PLAG_ModParameters         

  USE ModInterfaces, ONLY: RFLU_DecidePrint, & 
                           RkUpdateGeneric

  USE PLAG_ModInterfaces, ONLY: PLAG_CalcBreakup, &
                                PLAG_CalcRhsPosition, &
                                PLAG_ContinuousRandomWalk, &
                                PLAG_RFLU_InjectionDriver, &
                                PLAG_SetDependentVarsOld, &
                                PLAG_UpdateDataStruct

  USE PLAG_ModReallocateMemory
  USE PLAG_RFLU_ModFindCells
  USE PLAG_RFLU_ModComm
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(IN) :: iStage
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  LOGICAL :: doPrint
  INTEGER :: iReg,ic
  TYPE(t_global), POINTER :: global
  TYPE(t_plag),   POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_Update.F90,v $ $Revision: 1.2 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_Update',__FILE__)  

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pPlag => pRegion%plag

  iReg  = pRegion%iRegionGlobal

! ******************************************************************************
! Initialize send counters even for null nPcls
! ******************************************************************************
  CALL PLAG_RFLU_InitSendCounters(pRegion)     
 

! ******************************************************************************
! Main algorithm for discrete particle evolution
! ******************************************************************************
! Updates Eulerian volume fraction at each RK step. Rahul
 IF (pRegion%mixtInput%spaceDiscr .EQ. DISCR_UPW_AUSMPLUSUP) THEN
  DO ic = 1,pRegion%grid%nCellsTot
    pRegion%plag%VFracEOld(1,ic) = pRegion%plag%vFracE(1,ic)
  END DO 
 END IF
! End - Rahul
  IF ( pPlag%nPcls > 0 ) THEN

! ============================================================================== 
!   Calculate RHS for position vector
! ============================================================================== 

    CALL PLAG_CalcRhsPosition(pRegion)

! ============================================================================== 
!   Invoke RK Update
! ============================================================================== 
IF(1==2) THEN
    CALL RkUpdateGeneric(pRegion,VAR_TYPE_POINT,iStage,1,pPlag%nPcls,1, & 
                         pPlag%nCv,pPlag%cv,pPlag%cvOld,pPlag%rhs,pPlag%rhsSum)
END IF
! ============================================================================== 
!   Continuous Random Walk implementation for particle velocity fluctuation                         
! ============================================================================== 

    IF ( (pRegion%plagInput%flagCRW .EQV. .TRUE.) .AND. &
         (iStage == 1) ) THEN
      CALL PLAG_SetDependentVarsOld(pRegion)
    END IF ! pRegion%plagInput%flagCRW

    IF ( (pRegion%plagInput%flagCRW .EQV. .TRUE.) .AND. &
         (iStage == global%nrkSteps) ) THEN
! DEBUG: Manoj: 2012-06-25         
!IF (1==2) THEN         
      CALL PLAG_ContinuousRandomWalk(pRegion)
!END IF ! 1==2
! END DEBUG
    END IF ! pRegion%plagInput%flagCRW

! TEMPORARY
! ============================================================================== 
!   Create New location for particle search testing  
! ============================================================================== 
!
!    CALL PLAG_CreateLocation(pRegion)
! END TEMPORARY

! ============================================================================== 
!   Determine total distance travelled
! ============================================================================== 

    CALL PLAG_RFLU_ComputeDistTot(pRegion)      

! ============================================================================== 
!   Invoke particle location search. NOTE brute-force and Octree cannot deal 
!   with particles hitting the wall and bouncing. They should only be used for
!   testing purposes.
! ============================================================================== 

    SELECT CASE ( pRegion%plagInput%findPclMethod )
      CASE ( FIND_PCL_METHOD_TRAJ_FAST ) 
        CALL PLAG_RFLU_FindCellsTrajFast(pRegion,1,pPlag%nPcls)      
      CASE ( FIND_PCL_METHOD_TRAJ_SAFE ) 
        CALL PLAG_RFLU_FindCellsTrajSafe(pRegion,1,pPlag%nPcls)  
      CASE ( FIND_PCL_METHOD_BRUTE )
        CALL PLAG_RFLU_FindCellsBruteMod(pRegion)
      ! Subbu - Call hard coded pcl tracking routine
      CASE ( FIND_PCL_METHOD_HARDCODE)
        CALL PLAG_RFLU_FindCellsHardCode(pRegion,1,pPlag%nPcls)
      ! Subbu - End Call hard coded pcl tracking routine
      CASE ( FIND_PCL_METHOD_OCT ) 
        CALL PLAG_RFLU_FindCellsOctMod(pRegion)
      CASE ( FIND_PCL_METHOD_LOHNER ) 
        CALL PLAG_RFLU_FindCellsLohner(pRegion)
      CASE DEFAULT 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%plagInput%findPclMethod

! ============================================================================== 
!   Update particle data structure 
!    Note: done for serial computations, else for parallel runs
!          update performed in PLAG_RFLU_CommDriver
! ============================================================================== 
    
    IF ( global%nRegions == 1 ) &
      CALL PLAG_UpdateDataStruct(pRegion)   
    END IF ! nPcls

! ******************************************************************************
! Invoke injection algorithm at last RK stage 
! ******************************************************************************

  IF ( iStage == global%nrkSteps ) THEN 
    CALL PLAG_RFLU_InjectionDriver(pRegion)
  END IF ! iStage

! ******************************************************************************
! Print number of particles in the domain 
! ******************************************************************************

! TEMPORARY
!  IF ( iStage == global%nrkSteps ) THEN 
!    doPrint = RFLU_DecidePrint(global)
!    IF ( (global%verbLevel > VERBOSE_NONE) .AND. (doPrint .EQV. .TRUE.)  ) & 
!      WRITE(STDOUT,'(A,I4,I8)') 'iReg: nPcls = ',iReg,pPlag%nPcls
!  END IF ! iStage
! END TEMPORARY

! ******************************************************************************
! Invoke breakup algorithm at last RK stage
! ******************************************************************************

  IF ( pRegion%plagInput%breakupModel > PLAG_BREAKUP_NOMODEL .AND. &
       iStage == global%nrkSteps ) THEN 
    CALL PLAG_CalcBreakup(pRegion,iReg)
  END IF ! pRegion%plagInput%breakupModel

! ******************************************************************************
! Invoke memory reallocation at last RK stage and for non null size
! ******************************************************************************

  IF ( pPlag%nPcls > 0 .AND. iStage == global%nrkSteps ) THEN 
    CALL PLAG_ReallocMemWrapper(pRegion)   
  END IF ! pPlag%nPcls

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLU_Update

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_Update.F90,v $
! Revision 1.2  2015/12/18 23:15:41  rahul
! Added Eulerian particle volume fraction update at each RK step. This is
! required for multiphase AUSM+up scheme.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:52  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:50:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:34  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.20  2007/03/20 22:04:07  fnajjar
! Turned off printing of number of particles for each region separately
!
! Revision 1.19  2006/05/20 19:08:24  fnajjar
! Fixed bug by moving initialization of nPclsSend out IF nPcls Statement
!
! Revision 1.18  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.17  2005/12/19 16:47:12  fnajjar
! Added if statement to call PLAG_UpdateDataStruc for single region
!
! Revision 1.16  2005/05/18 22:21:29  fnajjar
! Added init of send counters, comp of distTot
!
! Revision 1.15  2005/05/02 21:59:02  haselbac
! Adapted to changes in interface of PLAG_RFLU_FindCellsTraj routines
!
! Revision 1.14  2005/04/27 14:57:58  fnajjar
! Included module call to PLAG_RFLU_ModFindCells
!
! Revision 1.13  2005/03/11 02:28:16  haselbac
! Adapted to new PLAG_RFLU_FindCellsTrajXYZ routines
!
! Revision 1.12  2005/01/01 21:32:21  haselbac
! Added routine to find particle cells using Aptes method
!
! Revision 1.11  2004/11/17 16:42:41  haselbac
! Replaced RKUpdate call with call to RkUpdateGeneric
!
! Revision 1.10  2004/11/06 21:08:41  fnajjar
! Bug fix: added interface entries for particle-cell search routines
!
! Revision 1.9  2004/11/06 18:48:44  fnajjar
! Added doPrint statement to decide printout for nPcls
!
! Revision 1.8  2004/11/05 21:50:09  fnajjar
! Updated calls to particle-cell search routines
!
! Revision 1.7  2004/10/11 22:10:01  haselbac
! Renamed procedures, changed IF to SELECT
!
! Revision 1.6  2004/10/08 22:11:20  haselbac
! Added brute-force search option for finding particle positions
!
! Revision 1.5  2004/07/29 16:32:22  fnajjar
! Removed temporary io flag for cleaner version
!
! Revision 1.4  2004/07/29 16:05:36  fnajjar
! Removed temporary io file after memory reallocation
!
! Revision 1.3  2004/07/28 18:58:54  fnajjar
! Included dynamic memory reallocation call
!
! Revision 1.2  2004/05/05 20:59:55  fnajjar
! Clean up; removed output statements
!
! Revision 1.1  2004/03/26 21:33:24  fnajjar
! Initial import for RFLU-specific update routine
!
! ******************************************************************************

