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
! Purpose: Obtain boundary values from Genx.
!
! Description: None.
!
! Input: 
!   region	Dimensions and topology
!
! Output: 
!   region%levels%patch 	Input values for boundary conditions
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_GetBoundaryValues.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_GetBoundaryValues(region)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters
  USE ModMPI
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region) :: region

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: errorString,RCSIdentString
  INTEGER :: ifl,iPatch,boCount,boTot
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pVals
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_GetBoundaryValues.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction(global,'RFLU_GetBoundaryValues',__FILE__)

! ******************************************************************************
! Loop over ALL boundaries 
! ******************************************************************************


  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A,1X,I3)') SOLVER_NAME, &
                                   'Getting values from GENX buffers...'
  END IF ! global%myProcid

  boCount = 0
  DO iPatch=1,region%grid%nPatches
    pPatch  => region%patches(iPatch)
 
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%checkLevel == CHECK_HIGH .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN 
      IF ( pPatch%bcCoupled == BC_NOT_COUPLED ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A)') SOLVER_NAME,'Patch:',iPatch, & 
	                                    '(not interacting)'
      ELSE IF ( pPatch%bcCoupled == BC_NOT_BURNING ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A)') SOLVER_NAME,'Patch:',iPatch, & 
	                                    '(interacting)'      
      ELSE IF ( pPatch%bcCoupled == BC_BURNING ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A)') SOLVER_NAME,'Patch:',iPatch, & 
	                                    '(interacting and burning)'          
      END IF ! pPatch%bcCoupled        
    END IF ! global%myProcid 
           
! ==============================================================================
!   Coupled burning boundary
! ==============================================================================
           
    IF ( pPatch%bcCoupled == BC_BURNING )THEN    
      pVals => pPatch%mixt%vals
      DO ifl = 1,pPatch%nBFaces      
! TZ DOES NOT WORK PROPERLY WITH THE ROCMAN EVALUATED MDOT (leaving this in for now)
!        pPatch%mdotAlp(ifl) = global%tzRhos*global%tzA*pPatch%pf(ifl)**global%tzN
! Rocman3 has been modified to always pass in the nominal mdots, so we trust the mdot
! coming from Rocman3 now.
!
! We've the BFLAG off for faces that have burned out.  We need to zero these mdots
! in the fluid RK substeps.  The following check accomplishes mdot zeroing for 
! burned out faces.
        IF(pPatch%bflag(ifl) .ne. 0) THEN
           pVals(BCDAT_INJECT_MFRATE,ifl) = pPatch%mdotAlp(ifl)
        ELSE
	   boCount = boCount + 1
           pVals(BCDAT_INJECT_MFRATE,ifl) = 0.0
           pPatch%mdotAlp(ifl) = 0.0
        ENDIF
        pVals(BCDAT_INJECT_TEMP  ,ifl) = pPatch%tflmAlp(ifl)
      END DO ! ifl

      DO ifl = pPatch%nBFaces+1,pPatch%nBFacesTot
        pVals(BCDAT_INJECT_MFRATE,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pVals(BCDAT_INJECT_TEMP  ,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      END DO ! ifl

      IF ( global%myProcid == MASTERPROC .AND. & 
           global%checkLevel == CHECK_HIGH .AND. & 
           global%verbLevel > VERBOSE_LOW ) THEN 
         IF ( pPatch%nBFaces > 0 ) THEN
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Minimum/maximum values:'      
            WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'mdotAlp:    ', & 
                 MINVAL(pPatch%mdotAlp(1:pPatch%nBFaces)), & 
                 MAXVAL(pPatch%mdotAlp(1:pPatch%nBFaces))
            WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'TflmAlp:    ', & 
                 MINVAL(pPatch%tflmAlp(1:pPatch%nBFaces)), & 
                 MAXVAL(pPatch%tflmAlp(1:pPatch%nBFaces))
         END IF ! pPatch%nBFaces
      END IF ! global%myProcid       
      
      IF ( MINVAL(pPatch%mdotAlp(1:pPatch%nBFaces)) < 0.0_RFREAL ) THEN 
        WRITE(errorString,*) 'Patch:',pPatch%iPatchGlobal
        CALL ErrorStop(global,ERR_MDOT_NEGATIVE,__LINE__,TRIM(errorString))
      END IF ! MINVAL
    
      IF ( MINVAL(pPatch%tflmAlp(1:pPatch%nBFaces)) < 0.0_RFREAL ) THEN 
        WRITE(errorString,*) 'Patch:',pPatch%iPatchGlobal
        CALL ErrorStop(global,ERR_TFLM_NEGATIVE,__LINE__,TRIM(errorString))
      END IF ! MINVAL

! ==============================================================================
!   Coupled non-burning boundary. NOTE only copy boundary temperature for 
!   isothermal no-slip walls.
! ==============================================================================    
    
    ELSE IF ( pPatch%bcCoupled == BC_NOT_BURNING ) THEN
      IF ( pPatch%bcType == BC_NOSLIPWALL_TEMP ) THEN  
	pVals => pPatch%mixt%vals

	DO ifl = 1,pPatch%nBFaces      
          pVals(BCDAT_NOSLIP_T,ifl) = pPatch%tbAlp(ifl)
	END DO ! ifl

	DO ifl = pPatch%nBFaces+1,pPatch%nBFacesTot
          pVals(BCDAT_NOSLIP_T,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
	END DO ! ifl

	IF ( global%myProcid == MASTERPROC .AND. & 
             global%checkLevel == CHECK_HIGH .AND. & 
             global%verbLevel > VERBOSE_LOW ) THEN 
          IF ( pPatch%nBFaces > 0 ) THEN
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Minimum/maximum values:'      
            WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'TbAlp:    ', & 
                  MINVAL(pPatch%tbAlp(1:pPatch%nBFaces)), & 
                  MAXVAL(pPatch%tbAlp(1:pPatch%nBFaces))
          END IF ! pPatch%nBFaces
	END IF ! global%myProcid       

	IF ( MINVAL(pPatch%tbAlp(1:pPatch%nBFaces)) < 0.0_RFREAL ) THEN 
          WRITE(errorString,*) 'Patch:',pPatch%iPatchGlobal
          CALL ErrorStop(global,ERR_TB_NEGATIVE,__LINE__,TRIM(errorString))
	END IF ! MINVAL 
      END IF ! pPatch%bcType   
    
! ==============================================================================
!   Noncoupled boundary
! ==============================================================================    

    ELSE IF ( pPatch%bcCoupled == BC_NOT_COUPLED ) THEN   

! ==============================================================================
!   Impossible combination - defensive programming
! ==============================================================================    
      
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
    END IF ! pPatch%bcCoupled  
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%checkLevel == CHECK_HIGH .AND. boCount > 0) THEN
     WRITE(STDOUT,'(A,3X,A,1X,I9,1X,I9)') SOLVER_NAME, &
          'Local burned out faces:',boCount,global%myProcid
  ENDIF

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN 
     WRITE(STDOUT,'(A,1X,A,1X,I3)') SOLVER_NAME, &
           'Getting values from GENX buffers done.'
  END IF ! global%myProcid

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_GetBoundaryValues

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_GetBoundaryValues.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2009/08/28 18:29:38  mtcampbe
! RocfluMP integration with Rocstar and some makefile tweaks.  To build
! Rocstar with new Rocflu:
! make ROCFLU=RocfluMP
! To build Rocstar with the new RocfluND:
! make ROCFLU=RocfluMP HYPRE=/the/hypre/install/path
!
! Revision 1.3  2008/12/06 08:43:30  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:47:51  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:58:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.13  2006/08/19 15:37:53  mparmar
! Renamed patch variables
!
! Revision 1.12  2005/10/14 13:59:46  haselbac
! Added setting of boundary temperature for non-burning patches
!
! Revision 1.11  2004/10/19 19:23:04  haselbac
! Removed setting of relative velocity on injecting boundaries, cosmetics
!
! Revision 1.10  2003/08/20 02:09:58  haselbac
! Changed verbosity conditions to reduce solver output in GENx runs
!
! Revision 1.9  2003/05/13 23:41:35  haselbac
! Added init of dummy faces and check for temp
!
! Revision 1.8  2003/05/01 20:44:16  haselbac
! Added check for negative mdot
!
! Revision 1.7  2003/05/01 14:06:11  haselbac
! Removed warning for non-burning boundaries
!
! Revision 1.6  2003/04/12 21:33:01  haselbac
! Made parallel: Added MASTERPROC checks
!
! Revision 1.5  2003/01/28 18:22:28  haselbac
! Cosmetics only
!
! Revision 1.4  2002/10/12 19:08:45  haselbac
! Cosmetic changes only
!
! Revision 1.3  2002/10/12 14:48:15  haselbac
! Corrected bcFlag bug, some cosmetic changes
!
! Revision 1.2  2002/10/07 14:09:32  haselbac
! Changed ModBndpatch to ModBndPatch
!
! Revision 1.1  2002/10/05 18:28:18  haselbac
! Initial revision
!
! ******************************************************************************

