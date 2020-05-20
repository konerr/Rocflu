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
! Purpose: Write Prediction Metrics' history.
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
! $Id: RFLU_WritePM.F90,v 1.4 2015/12/19 01:16:49 rahul Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_WritePM(regions)

  USE ModParameters
  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModMPI
  USE RFLU_ModPlottingVars
  USE ModMixture, ONLY: t_mixt_input

#ifdef PLAG
  USE ModPartLag, ONLY: t_plag
  USE PLAG_ModParameters
#endif

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: regions(:)
  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion
  TYPE(t_mixt_input), POINTER :: pMixtInput
#ifdef PLAG
  TYPE(t_plag),   POINTER :: pPlag
#endif
! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,i,iReg
  REAL(RFREAL) :: pm1,pm1_l
  REAL(RFREAL) :: cntpm1_l,cntpm1
  REAL(RFREAL) :: eps,pout,minp,maxp,pressure
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: globalValsReal,localValsReal
#ifdef PLAG
  INTEGER :: iPcl
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: globalVF,localVF
  REAL(RFREAL) :: maxVF,radius,roc
  REAL(RFREAL) :: pm2,pm2_l
  REAL(RFREAL) :: cntpm2_l,cntpm2
#endif
! ******************************************************************************
! Start
! ******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_WritePM',__FILE__)

! ******************************************************************************
! Compute Predictive Metric
! ******************************************************************************

  pMixtInput => regions(1)%mixtInput

  pm1 = 0.0_RFREAL; pm1_l = 0.0_RFREAL
  cntpm1_l = 0.0_RFREAL; cntpm1 = 0.0_RFREAL
#ifdef PLAG
  pm2 = 0.0_RFREAL; pm2_l = 0.0_RFREAL
  cntpm2_l = 0.0_RFREAL; cntpm2 = 0.0_RFREAL
#endif

  ALLOCATE(localValsReal(0:global%nRegions),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'localValsReal')
  END IF ! global%error

  ALLOCATE(globalValsReal(0:global%nRegions),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'globalValsReal')
  END IF ! global%error

#ifdef PLAG
  ALLOCATE(localVF(0:global%nRegions),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'localVF')
  END IF ! global%error

  ALLOCATE(globalVF(0:global%nRegions),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'globalVF')
  END IF ! global%error
#endif

  DO iReg = 0,global%nRegions
    localValsReal(iReg) = HUGE(1.0_RFREAL)
#ifdef PLAG 
    localVF(iReg) = -HUGE(1.0_RFREAL)
#endif
  END DO ! iReg

  DO iReg = 1,global%nRegionsLocal
    pRegion => regions(iReg)

    localValsReal(pRegion%iRegionGlobal) =  &
                    MINVAL(pRegion%mixt%dv(DV_MIXT_PRES,:))
#ifdef PLAG
    pPlag => pRegion%plag
    localVF(pRegion%iRegionGlobal) = &
                    MAXVAL(pPlag%vFracL(1,:))
#endif
  END DO ! iReg

  CALL MPI_Allreduce(localValsReal(0:global%nRegions), &
                    globalValsReal(0:global%nRegions),global%nRegions+1, &
                    MPI_RFREAL,MPI_MIN,global%mpiComm,errorFlag)

  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
  END IF ! global%errorFlag

  pout = MINVAL(globalValsReal)

#ifdef PLAG
  CALL MPI_Allreduce(localVF(0:global%nRegions), &
                    globalVF(0:global%nRegions),global%nRegions+1, &
                    MPI_RFREAL,MPI_MAX,global%mpiComm,errorFlag)

  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
  END IF ! global%errorFlag

  maxVF = MAXVAL(globalVF)
#endif

  DO iReg = 1,global%nRegionsLocal
    pRegion => regions(iReg)

    minp = MINVAL(pRegion%mixt%dv(DV_MIXT_PRES,:))
    maxp = MAXVAL(pRegion%mixt%dv(DV_MIXT_PRES,:)) 

    IF(ABS(minp-pout)/pout.LT.0.1_RFREAL .AND. maxp.GT.3.0_RFREAL*minp)THEN
       DO i = 1,pRegion%grid%nCellsTot
          pressure = pRegion%mixt%dv(DV_MIXT_PRES,i)
          IF(pressure .GE. 0.99_RFREAL*maxp)THEN
            pm1_l = pm1_l + SQRT(pRegion%grid%xyz(XCOORD,i)**2 &
                          + pRegion%grid%xyz(YCOORD,i)**2)
            cntpm1_l = cntpm1_l + 1.0_RFREAL
          END IF
       END DO !i
    END IF !pout

#ifdef PLAG
    pPlag => pRegion%plag
    DO iPcl = 1, pRegion%plag%nPcls
       radius = SQRT(pPlag%cv(CV_PLAG_XPOS,iPcl)**2 &
                   +pPlag%cv(CV_PLAG_YPOS,iPcl)**2)
       roc = ABS(pPlag%vFracL(1,iPcl)-maxVF)/maxVF
       IF ( ((0.09 .LT. roc) .AND. (roc .LT. 0.1)) &
          .AND. (radius .GT. 0.05) ) THEN
          pm2_l = pm2_l + radius
          cntpm2_l = cntpm2_l + 1.0_RFREAL
       END IF
    END DO
#endif

  END DO !iReg

  CALL MPI_Allreduce(pm1_l,pm1,1,MPI_RFREAL,MPI_SUM, &
                         global%mpiComm,global%mpierr )
  CALL MPI_Allreduce(cntpm1_l,cntpm1,1,MPI_RFREAL,MPI_SUM, &
                         global%mpiComm,global%mpierr )

#ifdef PLAG
  CALL MPI_Allreduce(pm2_l,pm2,1,MPI_RFREAL,MPI_SUM, &
                         global%mpiComm,global%mpierr )
  CALL MPI_Allreduce(cntpm2_l,cntpm2,1,MPI_RFREAL,MPI_SUM, &
                         global%mpiComm,global%mpierr )
#endif
! ******************************************************************************
! Write PM data
! ******************************************************************************
  IF ( (global%flowType == FLOW_UNSTEADY) .AND. & 
            (global%myProcid==MASTERPROC) ) THEN
     IF ( global%plagUsed .EQV. .TRUE. ) THEN
#ifdef PLAG
       WRITE(IF_PM,'(1PE12.5,2X,E13.4,2X,E13.4)') &
       global%currentTime,pm1/MAX(1.0_RFREAL,cntpm1), &
       MAX(0.05_RFREAL,pm2/MAX(1.0_RFREAL,cntpm2))
#endif
     ELSE
       WRITE(IF_PM,'(1PE12.5,2X,E13.4)') & 
       global%currentTime,pm1/MAX(1.0_RFREAL,cntpm1)
     ENDIF
  END IF ! global%flowType

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_WritePM

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_WritePM.F90,v $
! Revision 1.4  2015/12/19 01:16:49  rahul
! Bug fix: Added PLAG definitions to pm2, pm2_l, cntpm2_l, cntpm2. Without
! this code will not compile for gas only cases.
!
! Revision 1.3  2015/08/12 19:41:42  brollin
! Updating module declaration in rfluinit.F90
!
! Revision 1.2  2015/07/23 23:11:19  brollin
! 1) The pressure coefficient of the  collision model has been changed back to its original form
! 2) New options in the format of the grid and solutions have been added. Now the user can choose the endianness, and convert from one to the over in rfluconv.
! 3) The solutions are now stored in folders named by timestamp or iteration number.
! 4) The address enty in the hashtable has been changed to an integer(8) for cases when the grid becomes very large.
! 5) RFLU_WritePM can now compute PM2 on the fly for the Macroscale problem
!
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

