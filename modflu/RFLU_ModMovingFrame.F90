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
! Purpose: Collection of routines to compute particle acceleration, velocity.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModMovingFrame.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModMovingFrame

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModError
  USE ModMPI

  USE RFLU_ModForcesMoments
  USE ModInterfaces, ONLY: MixtPerf_G_CpR, &
                           MixtPerf_R_M

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModMovingFrame.F90,v $ $Revision: 1.1.1.1 $'

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: RFLU_MVF_ComputeAcceleration, &  
            RFLU_MVF_CreatePatchVelAccel, & 
            RFLU_MVF_DestroyPatchVelAccel, &           
            RFLU_MVF_InitPatchVelAccel, & 
            RFLU_MVF_ModifyFlowField, &
            RFLU_MVF_ReadPatchVelAccel, &
            RFLU_MVF_SetVelocity, &
            RFLU_MVF_UpdateBC, &
            RFLU_MVF_WritePatchVelAccel

! ==============================================================================
! Private functions
! ==============================================================================

  PRIVATE :: RFLU_MVF_ClosePatchVelAccel, &
             RFLU_MVF_OpenPatchVelAccel

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS






! ******************************************************************************
!
! Purpose: Close mvf restart file.
!
! Description: None.
!
! Input: 
!   global              Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MVF_ClosePatchVelAccel(global)
  
    IMPLICIT NONE

! *****************************************************************************
!   Declarations and definitions
! *****************************************************************************

! =============================================================================
!   Arguments
! =============================================================================

    TYPE(t_global), POINTER :: global

! =============================================================================
!   Locals
! =============================================================================

    CHARACTER(CHRLEN) :: iFileName
    INTEGER :: errorFlag
  
! *****************************************************************************
!   Start
! *****************************************************************************

    CALL RegisterFunction(global,'RFLU_MVF_ClosePatchVelAccel',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing mvf restart file...'
    END IF ! global%verbLevel

! =============================================================================
!   Close file
! =============================================================================

    WRITE(iFileName,'(A,I4.4)') TRIM(global%outDir)//TRIM(global%casename)// &
                                '.vel_',global%iPatchGlobalMvFrame

    CLOSE(IF_VELACCEL,IOSTAT=errorFlag)                        
    global%error = errorFlag
    IF ( global%error /= 0 ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error

! *****************************************************************************
!   End
! *****************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing mvf restart file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MVF_ClosePatchVelAccel







! ******************************************************************************
!
! Purpose: Calculate force and acceleration on a particle represented by a patch.
!
! Description: None.
!
! Input:
!   regions             Data for all grid regions
!
! Output:
!   regions             Data for all grid regions
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MVF_ComputeAcceleration(regions)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: iPatch,iReg,ix
    REAL(RFREAL) :: accX,accY,accZ,deltaT,mvfMass,refForce,timeAccEnd, &
                    timeAccStart
    TYPE(t_global), POINTER :: global
    TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => regions(1)%global

    CALL RegisterFunction(global,'RFLU_MVF_ComputeAcceleration',__FILE__)

! ******************************************************************************
!   Compute force coefficient on given patch. NOTE now we compute forces on all 
!   patches, although they would only be needed for specified patch. May need to
!   be modified in the future for efficiency reasons.
! ******************************************************************************

    IF ( global%mvfAccFlag .EQV. .FALSE. ) THEN
      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_ComputeLocalForcesMoments(pRegion)
      END DO ! iReg

      CALL RFLU_ComputeGlobalForcesMoments(regions)
    END IF ! global%mvfAccFlag

! ******************************************************************************
!   Compute forces and acceleration
! ******************************************************************************

    iPatch   = global%iPatchGlobalMvFrame
    mvfMass  = global%mvfMass
    refForce = 0.5_RFREAL*global%refDensity*global%refVelocity & 
                         *global%refVelocity*global%forceRefArea

    timeAccStart = global%mvfAccTs
    timeAccEnd   = global%mvfAccTe
    accX         = global%mvfAccX
    accY         = global%mvfAccY
    accZ         = global%mvfAccZ

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)

      IF ( global%mvfAccFlag .EQV. .TRUE. ) THEN
        IF ( global%currentTime >= timeAccStart .AND. &
             global%currentTime < timeAccEnd ) THEN
          deltaT = global%currentTime-timeAccStart

          SELECT CASE (global%mvfAccType)
            CASE (ACC_CONSTANT)
              pRegion%mvfAcc(XCOORD) = accX
              pRegion%mvfAcc(YCOORD) = accY
              pRegion%mvfAcc(ZCOORD) = accZ
            CASE (ACC_SINUSOIDAL)
              pRegion%mvfAcc(XCOORD) = accX*SIN(global%mvfOmega*deltaT &
                                               +global%mvfPhase)
              pRegion%mvfAcc(YCOORD) = accY*SIN(global%mvfOmega*deltaT &
                                               +global%mvfPhase)
              pRegion%mvfAcc(ZCOORD) = accZ*SIN(global%mvfOmega*deltaT &
                                               +global%mvfPhase)
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pPatch%bcType
        ELSE
          pRegion%mvfAcc(XCOORD) = 0.0_RFREAL
          pRegion%mvfAcc(YCOORD) = 0.0_RFREAL
          pRegion%mvfAcc(ZCOORD) = 0.0_RFREAL
        END IF ! currentTime
      ELSE
        pRegion%mvfAcc(XCOORD) = pRegion%forceCoeffsGlobal(XCOORD,COMP_PRES, &
                                                           iPatch)
        pRegion%mvfAcc(YCOORD) = pRegion%forceCoeffsGlobal(YCOORD,COMP_PRES, &
                                                           iPatch)
        pRegion%mvfAcc(ZCOORD) = pRegion%forceCoeffsGlobal(ZCOORD,COMP_PRES, &
                                                           iPatch)

        pRegion%mvfAcc(XCOORD) = pRegion%mvfAcc(XCOORD)*refForce/mvfMass
        pRegion%mvfAcc(YCOORD) = pRegion%mvfAcc(YCOORD)*refForce/mvfMass
        pRegion%mvfAcc(ZCOORD) = pRegion%mvfAcc(ZCOORD)*refForce/mvfMass
      END IF ! global%mvfAccFlag

      IF (pRegion%mixtInput%axiFlag .EQV. .TRUE. ) THEN
        pRegion%mvfAcc(YCOORD) = 0.0_RFREAL
        pRegion%mvfAcc(ZCOORD) = 0.0_RFREAL
      END IF ! pRegion%mixtInput%axiFlag

      IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN      
        DO ix = XCOORD,ZCOORD
          pRegion%mvfAcc(ix)    = pRegion%mvfAcc(ix)*global%refLength &
                                       /(global%refVelocity*global%refVelocity)
        END DO ! ix 
      END IF ! global%solverType          
    END DO ! iReg

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MVF_ComputeAcceleration   


  





! *******************************************************************************
!
! Purpose: Create particle acceleration and velocity variables.
!
! Description: None.
!
! Input:
!   pRegion		Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MVF_CreatePatchVelAccel(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
  
    CALL RegisterFunction(global,'RFLU_MVF_CreatePatchVelAccel',__FILE__)
    
! ******************************************************************************
!   Allocate memory
! ******************************************************************************
        
    ALLOCATE(pRegion%mvfLoc(XCOORD:ZCOORD),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%mvfLoc')
    END IF ! global%error  
          
    ALLOCATE(pRegion%mvfLocOld(XCOORD:ZCOORD),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%mvfLocOld')
    END IF ! global%error  
         
    ALLOCATE(pRegion%mvfVel(XCOORD:ZCOORD),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%mvfVel')
    END IF ! global%error  
          
    ALLOCATE(pRegion%mvfVelOld(XCOORD:ZCOORD),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%mvfVelOld')
    END IF ! global%error  
          
    ALLOCATE(pRegion%mvfVelSum(XCOORD:ZCOORD),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%mvfVelSum')
    END IF ! global%error  

    ALLOCATE(pRegion%mvfAcc(XCOORD:ZCOORD),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%mvfAcc')
    END IF ! global%error  
    
    IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN      
      ALLOCATE(pRegion%mvfAccOld(XCOORD:ZCOORD),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%mvfAccOld')
      END IF ! global%error  
    END IF ! global%solverType          

    ALLOCATE(pRegion%mvfAccSum(XCOORD:ZCOORD),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%mvfAccSum')
    END IF ! global%error  

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MVF_CreatePatchVelAccel
 







! *******************************************************************************
!
! Purpose: Destroy particle acceleration and velocity variables.
!
! Description: None.
!
! Input:
!   pRegion		Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MVF_DestroyPatchVelAccel(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag
    TYPE(t_global), POINTER :: global    
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_MVF_DestroyPatchVelAccel',__FILE__)
    
! ******************************************************************************
!   Deallocate memory
! ******************************************************************************
        
    DEALLOCATE(pRegion%mvfLoc,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mvfLoc')
    END IF ! global%error  
          
    DEALLOCATE(pRegion%mvfLocOld,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mvfLocOld')
    END IF ! global%error  
          
    DEALLOCATE(pRegion%mvfVel,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mvfVel')
    END IF ! global%error  
          
    DEALLOCATE(pRegion%mvfVelOld,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mvfVelOld')
    END IF ! global%error  
          
    DEALLOCATE(pRegion%mvfVelSum,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mvfVelSum')
    END IF ! global%error  
          
    DEALLOCATE(pRegion%mvfAcc,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mvfAcc')
    END IF ! global%error  
          
    IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN      
      DEALLOCATE(pRegion%mvfAccOld,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mvfAccOld')
      END IF ! global%error  
    END IF ! global%solverType          

    DEALLOCATE(pRegion%mvfAccSum,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mvfAccSum')
    END IF ! global%error  
          
! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    NULLIFY(pRegion%mvfLoc)
    NULLIFY(pRegion%mvfLocOld)
    NULLIFY(pRegion%mvfVel)
    NULLIFY(pRegion%mvfVelOld)
    NULLIFY(pRegion%mvfVelSum)
    NULLIFY(pRegion%mvfAcc)
    IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN      
      NULLIFY(pRegion%mvfAccOld)
    END IF ! global%solverType          
    NULLIFY(pRegion%mvfAccSum)

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MVF_DestroyPatchVelAccel






! *******************************************************************************
!
! Purpose: Initialize particle acceleration and velocity variables.
!
! Description: None.
!
! Input:
!   pRegion		Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MVF_InitPatchVelAccel(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: ix
    TYPE(t_global), POINTER :: global
   
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_MVF_InitPatchVelAccel',__FILE__)
    
! ******************************************************************************
!   Initialize memory
! ******************************************************************************

    DO ix = XCOORD,ZCOORD
      pRegion%mvfLoc(ix)    = 0.0_RFREAL
      pRegion%mvfLocOld(ix) = 0.0_RFREAL
      pRegion%mvfVel(ix)    = 0.0_RFREAL
      pRegion%mvfVelOld(ix) = 0.0_RFREAL
      pRegion%mvfVelSum(ix) = 0.0_RFREAL
      pRegion%mvfAcc(ix)    = 0.0_RFREAL
      pRegion%mvfAccSum(ix) = 0.0_RFREAL 
    END DO ! ix 
    
    pRegion%mvfLoc(XCOORD) = global%mvfLocInitX
    pRegion%mvfLoc(YCOORD) = global%mvfLocInitY
    pRegion%mvfLoc(ZCOORD) = global%mvfLocInitZ

    pRegion%mvfVel(XCOORD) = global%mvfVelInitX
    pRegion%mvfVel(YCOORD) = global%mvfVelInitY
    pRegion%mvfVel(ZCOORD) = global%mvfVelInitZ

    IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN      
      DO ix = XCOORD,ZCOORD
        pRegion%mvfLoc(ix)    = pRegion%mvfLoc(ix)/global%refLength
        pRegion%mvfLocOld(ix) = pRegion%mvfLocOld(ix)/global%refLength
        pRegion%mvfVel(ix)    = pRegion%mvfVel(ix)/global%refVelocity
        pRegion%mvfVelOld(ix) = pRegion%mvfVelOld(ix)/global%refVelocity
        pRegion%mvfVelSum(ix) = pRegion%mvfVelSum(ix)/global%refVelocity
        pRegion%mvfAcc(ix)    = pRegion%mvfAcc(ix)*global%refLength &
                                       /(global%refVelocity*global%refVelocity)
        pRegion%mvfAccSum(ix) = pRegion%mvfAccSum(ix)*global%refLength & 
                                       /(global%refVelocity*global%refVelocity)
      END DO ! ix 
    END IF ! global%solverType          
 
! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MVF_InitPatchVelAccel
 






! ******************************************************************************
!
! Purpose: Modify flow field data to moving reference frame attached to particle
!
! Description: None.
!
! Input:
!   pRegion             Data for a region
!
! Output:
!   pRegion             Data for a region
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MVF_ModifyFlowField(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: icg
    REAL(RFREAL) :: Eo,ir,r,u,uOld,v,vOld,w,wOld,p 
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MVF_ModifyFlowField',__FILE__)

! ******************************************************************************
!   Modify velocity field and keep the pressure constant 
! ******************************************************************************

    DO icg=1,pRegion%grid%nCellsTot
      r  = pRegion%mixt%cv(CV_MIXT_DENS,icg)
      ir = 1.0_RFREAL/r

      uOld = pRegion%mixt%cv(CV_MIXT_XMOM,icg)*ir
      vOld = pRegion%mixt%cv(CV_MIXT_YMOM,icg)*ir
      wOld = pRegion%mixt%cv(CV_MIXT_ZMOM,icg)*ir

      u = uOld - pRegion%mvfVel(XCOORD)
      v = vOld - pRegion%mvfVel(YCOORD)
      w = wOld - pRegion%mvfVel(ZCOORD)

      pRegion%mixt%cv(CV_MIXT_XMOM,icg) = r*u
      pRegion%mixt%cv(CV_MIXT_YMOM,icg) = r*v
      pRegion%mixt%cv(CV_MIXT_ZMOM,icg) = r*w

      Eo = pRegion%mixt%cv(CV_MIXT_ENER,icg)*ir

      Eo = Eo + 0.5_RFREAL*(u*u + v*v + w*w - uOld*uOld - vOld*vOld - wOld*wOld) 

      pRegion%mixt%cv(CV_MIXT_ENER,icg) = r*Eo 
    END DO ! icg

! ******************************************************************************
!   End
! * *****************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MVF_ModifyFlowField









! ******************************************************************************
!
! Purpose: Open mvf restart file.
!
! Description: None.
!
! Input: 
!   global              Pointer to global data
!   filePosition        Position at which file to be opened (start or end)
!
! Output: 
!   fileExists          Logical indicating whether file exists
!
! Notes: 
!   1. The filePosition parameter is needed because the mvf restart file
!      is opened with two goals. The first is to open it with the goal of
!      getting the last output time and corresponding patch velocities. 
!      The second is to be able to write to it by appending additional lines.
!   2. The fileExists parameter is needed because on if the file exists
!      on restarting the code, it will need to be read to determine the 
!      last output iteration or time.
!
! ******************************************************************************

  SUBROUTINE RFLU_MVF_OpenPatchVelAccel(global,filePosition,fileExists)

    USE ModBuildFileNames, ONLY: BuildFileNamePlain

! *****************************************************************************
!   Declarations and definitions
! *****************************************************************************

! =============================================================================
!   Arguments
! =============================================================================

    LOGICAL, INTENT(OUT) :: fileExists
    INTEGER, INTENT(IN) :: filePosition
    TYPE(t_global), POINTER :: global

! =============================================================================
!   Locals
! =============================================================================

    CHARACTER(CHRLEN) :: iFileName
    INTEGER :: errorFlag

! *****************************************************************************
!   Start
! *****************************************************************************

    CALL RegisterFunction(global,'RFLU_MVF_OpenPatchVelAccel',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Opening mvf restart file...'
    END IF ! global%verbLevel

! =============================================================================
!   Open file
! =============================================================================

    WRITE(iFileName,'(A,I4.4)') TRIM(global%outDir)//TRIM(global%casename)// &
                                '.vel_',global%iPatchGlobalMvFrame

    INQUIRE(FILE=iFileName,EXIST=fileExists)

    IF ( fileExists .EQV. .TRUE. ) THEN
      IF ( filePosition == FILE_POSITION_START ) THEN
        OPEN(IF_VELACCEL,FILE=iFileName,FORM='FORMATTED',STATUS='OLD', &
             IOSTAT=errorFlag)
      ELSE IF ( filePosition == FILE_POSITION_END ) THEN
        OPEN(IF_VELACCEL,FILE=iFileName,FORM='FORMATTED',STATUS='OLD', &
             POSITION='APPEND',IOSTAT=errorFlag)
      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! filePosition
    ELSE
      OPEN(IF_VELACCEL,FILE=iFileName,FORM='FORMATTED',STATUS='NEW', &
           IOSTAT=errorFlag)
    END IF ! file

    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error

! *****************************************************************************
!   End
! *****************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Opening mvf restart file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MVF_OpenPatchVelAccel







! ******************************************************************************
!
! Purpose: Read mvf restart info file to get last particle velocity.
!
! Description: None.
!
! Input: 
!   global              Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MVF_ReadPatchVelAccel(pRegion)

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    LOGICAL :: fileExists
    INTEGER :: dummyInteger,errorFlag
    REAL(RFREAL) :: accx,accy,accz,dummyRFReal,locx,locy,locz,mvfTime,velx, &
                    vely,velz
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MVF_ReadPatchVelAccel',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Reading mvf restart file...'
    END IF ! global%verbLevel

! ==============================================================================
!   Open file
! ==============================================================================

    CALL RFLU_MVF_OpenPatchVelAccel(global,FILE_POSITION_START,fileExists)
  
! ==============================================================================
!   Set or read restart iteration or time
! ==============================================================================

    IF ( fileExists .EQV. .TRUE. ) THEN
      DO 
        READ(IF_VELACCEL,*,IOSTAT=errorFlag) mvfTime,locx,locy,locz,velx,vely, &
                                             velz,accx,accy,accz
        
        IF ( errorFlag /= ERR_NONE ) THEN 
          EXIT
        ELSE
! TEMPORARY : Need to have proper logic to read vel file
!          IF ( mvfTime > global%currentTime ) THEN
!             EXIT
!          ELSE 
          pRegion%mvfLoc(XCOORD) = locx
          pRegion%mvfLoc(YCOORD) = locy
          pRegion%mvfLoc(ZCOORD) = locz
          pRegion%mvfVel(XCOORD) = velx
          pRegion%mvfVel(YCOORD) = vely
          pRegion%mvfVel(ZCOORD) = velz
          pRegion%mvfAcc(XCOORD) = accx
          pRegion%mvfAcc(YCOORD) = accy
          pRegion%mvfAcc(ZCOORD) = accz
!        END IF ! mvfTime
! END TEMPORARY
        END IF ! errorFlag
      END DO ! <empty>
    END IF ! fileExists

! ==============================================================================
!   Write info
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A,1X,E24.16)') SOLVER_NAME,'Restart time:',mvfTime
      WRITE(STDOUT,'(A,3X,A,1X,E24.16)') SOLVER_NAME,'VelocityX:', &
                                          pRegion%mvfVel(XCOORD)
      WRITE(STDOUT,'(A,3X,A,1X,E24.16)') SOLVER_NAME,'VelocityY:', &
                                          pRegion%mvfVel(YCOORD)
      WRITE(STDOUT,'(A,3X,A,1X,E24.16)') SOLVER_NAME,'VelocityZ:', &
                                          pRegion%mvfVel(ZCOORD)
      WRITE(STDOUT,'(A,3X,A,1X,E24.16)') SOLVER_NAME,'AccelerationX:', &
                                          pRegion%mvfAcc(XCOORD)
      WRITE(STDOUT,'(A,3X,A,1X,E24.16)') SOLVER_NAME,'AccelerationY:', &
                                          pRegion%mvfAcc(YCOORD)
      WRITE(STDOUT,'(A,3X,A,1X,E24.16)') SOLVER_NAME,'AccelerationZ:', &
                                          pRegion%mvfAcc(ZCOORD)
    END IF ! global%myProcid

! ==============================================================================
!   Close file
! ==============================================================================

    CALL RFLU_MVF_ClosePatchVelAccel(global)        

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Reading mvf restart file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MVF_ReadPatchVelAccel







! ******************************************************************************
!
! Purpose: Force a predecided velocity profile on MVF velocity between given
!          time interval 
!
! Description: None.
!
! Input:
!   pRegion             Data for a region
!
! Output:
!   pRegion             Data for a region
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MVF_SetVelocity(region)


    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region) :: region 

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: ix
    REAL(RFREAL) :: accX,accY,accZ,deltaT,omega,phase,timeAccStart,timeAccEnd
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => region%global

    CALL RegisterFunction(global,'RFLU_MVF_SetVelocity',__FILE__)

    timeAccStart = global%mvfAccTs
    timeAccEnd   = global%mvfAccTe
    accX         = global%mvfaccX
    accY         = global%mvfaccY
    accZ         = global%mvfaccZ
    omega        = global%mvfOmega
    phase        = global%mvfPhase

! ******************************************************************************
!   Enforce particle velocity to velocity function value
! ******************************************************************************

    IF ( global%currentTime < timeAccStart ) THEN
      region%mvfVel(XCOORD) = 0.0_RFREAL 
      region%mvfVel(YCOORD) = 0.0_RFREAL 
      region%mvfVel(ZCOORD) = 0.0_RFREAL 
    ELSE IF ( global%currentTime >= timeAccStart .AND. &
              global%currentTime < timeAccEnd ) THEN
      deltaT = global%currentTime-timeAccStart

      SELECT CASE (global%mvfAccType)
        CASE (ACC_CONSTANT)
          region%mvfVel(XCOORD) = global%mvfVelInitX + accX*deltaT
          region%mvfVel(YCOORD) = global%mvfVelInitY + accY*deltaT
          region%mvfVel(ZCOORD) = global%mvfVelInitZ + accZ*deltaT
        CASE (ACC_SINUSOIDAL)
          region%mvfVel(XCOORD) = global%mvfVelInitX &
                                + (accX/omega)*COS(phase) &
                                - (accX/omega)*COS(omega*deltaT+phase) 
          region%mvfVel(YCOORD) = global%mvfVelInitY &
                                + (accY/omega)*COS(phase) &
                                - (accY/omega)*COS(omega*deltaT+phase) 
          region%mvfVel(ZCOORD) = global%mvfVelInitZ &
                                + (accZ/omega)*COS(phase) &
                                - (accZ/omega)*COS(omega*deltaT+phase) 
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%bcType
    ELSE
      deltaT = timeAccEnd-timeAccStart

      SELECT CASE (global%mvfAccType)
        CASE (ACC_CONSTANT)
          region%mvfVel(XCOORD) = global%mvfVelInitX + accX*deltaT
          region%mvfVel(YCOORD) = global%mvfVelInitY + accY*deltaT
          region%mvfVel(ZCOORD) = global%mvfVelInitZ + accZ*deltaT
        CASE (ACC_SINUSOIDAL)
          region%mvfVel(XCOORD) = global%mvfVelInitX &
                                + (accX/omega)*COS(phase) &
                                - (accX/omega)*COS(omega*deltaT+phase) 
          region%mvfVel(YCOORD) = global%mvfVelInitY &
                                + (accY/omega)*COS(phase) &
                                - (accY/omega)*COS(omega*deltaT+phase) 
          region%mvfVel(ZCOORD) = global%mvfVelInitZ &
                                + (accZ/omega)*COS(phase) &
                                - (accZ/omega)*COS(omega*deltaT+phase) 
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%bcType
    END IF ! currentTime

    IF (region%mixtInput%axiFlag .EQV. .TRUE. ) THEN
      region%mvfVel(YCOORD) = 0.0_RFREAL
      region%mvfVel(ZCOORD) = 0.0_RFREAL
    END IF ! region%mixtInput%axiFlag

    IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN      
      DO ix = XCOORD,ZCOORD
        region%mvfVel(ix) = region%mvfVel(ix)/global%refVelocity
      END DO ! ix 
    END IF ! global%solverType          

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MVF_SetVelocity








! ******************************************************************************
!
! Purpose: Update boundary data as moving reference velocity changes.
!
! Description: None.
!
! Input:
!   pRegion             Data for a region
!
! Output:
!   pRegion             Data for a region
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MVF_UpdateBC(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region) :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: c1,distrib,ifl,indCp,indMol,iPatch,ix
    REAL(RFREAL) :: aoa,aos,cp,ggas,mf,mm,Mmvf,rgas,SpeedOfSound,tf,VelX,VelY, &
                    VelZ,VelM
    REAL(RFREAL), DIMENSION(:,:), POINTER :: gv
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MVF_UpdateBC',__FILE__)

! ******************************************************************************
!   Enforce particle velocity to velocity function value
! ******************************************************************************

    indCp   = pRegion%mixtInput%indCp
    indMol  = pRegion%mixtInput%indMol

    gv  => pRegion%mixt%gv

    IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN      
      DO ix = XCOORD,ZCOORD
        pRegion%mvfVel(ix) = pRegion%mvfVel(ix)*global%refVelocity
      END DO ! ix 
    END IF ! global%solverType          

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)
   
      distrib = pPatch%mixt%distrib

      SELECT CASE( pPatch%bcType )
        CASE (BC_SLIPWALL)
        CASE (BC_NOSLIPWALL_HFLUX)
        CASE (BC_NOSLIPWALL_TEMP)
        CASE (BC_INFLOW_TOTANG)
        CASE (BC_INFLOW_VELTEMP)
          DO ifl = 1,pPatch%nBFaces
            c1   = pPatch%bf2c(ifl)
  
            VelX = pPatch%mixt%ivtU - pRegion%mvfVel(XCOORD)
            VelY = pPatch%mixt%ivtV - pRegion%mvfVel(YCOORD)
            VelZ = pPatch%mixt%ivtW - pRegion%mvfVel(ZCOORD)

            pPatch%mixt%vals(BCDAT_INFLOW_U,distrib*ifl) = VelX
            pPatch%mixt%vals(BCDAT_INFLOW_V,distrib*ifl) = VelY
            pPatch%mixt%vals(BCDAT_INFLOW_W,distrib*ifl) = VelZ
          END DO ! ifl
        CASE (BC_OUTFLOW)
        CASE (BC_FARFIELD)
          DO ifl = 1,pPatch%nBFaces
            c1   = pPatch%bf2c(ifl)

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              mm = gv(GV_MIXT_MOL,indMol*c1)
              cp = gv(GV_MIXT_CP ,indCp *c1)

              rgas = MixtPerf_R_M(mm)
              ggas = global%refGamma
            ELSE
              mm = gv(GV_MIXT_MOL,indMol*c1)
              cp = gv(GV_MIXT_CP ,indCp *c1)

              rgas = MixtPerf_R_M(mm)
              ggas = MixtPerf_G_CpR(cp,rgas)
            END IF ! solverType

            mf   = pPatch%mixt%ffMachno
            aoa  = pPatch%mixt%ffAttack
            aos  = pPatch%mixt%ffSlip

            tf   = pPatch%mixt%vals(BCDAT_FARF_TEMP  ,distrib*ifl)

            SpeedOfSound = SQRT(ggas*rgas*tf) 

            VelX = mf*SpeedOfSound*COS(aoa)*COS(aos) - pRegion%mvfVel(XCOORD)
            VelY = mf*SpeedOfSound*SIN(aoa)*COS(aos) - pRegion%mvfVel(YCOORD)
            VelZ = mf*SpeedOfSound*SIN(aos)          - pRegion%mvfVel(ZCOORD)

            VelM = SQRT(VelX*VelX + VelY*VelY + VelZ*VelZ)

            Mmvf = VelM/SpeedOfSound 
        
            pPatch%mixt%vals(BCDAT_FARF_MACH,distrib*ifl) = Mmvf
            pPatch%mixt%vals(BCDAT_FARF_ATTACK,distrib*ifl) = ATAN2(VelY,VelX)
            pPatch%mixt%vals(BCDAT_FARF_SLIP,distrib*ifl) = ATAN2(VelZ, &
                                                   SQRT(VelX*VelX+VelY*VelY))
          END DO ! ifl
        CASE (BC_INJECTION)
        CASE (BC_PERIODIC, &
              BC_SYMMETRY, &
              BC_VIRTUAL)
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%bcType
    END DO ! iPatch

    IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN      
      DO ix = XCOORD,ZCOORD
        pRegion%mvfVel(ix) = pRegion%mvfVel(ix)/global%refVelocity
      END DO ! ix 
    END IF ! global%solverType          

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MVF_UpdateBC






! *******************************************************************************
!
! Purpose: Write Patch (particle) velocity and acceleration to file.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine should only be called for the master process. It will not 
!      do anything if called by any other process (safeguard).
!
! ******************************************************************************

  SUBROUTINE RFLU_MVF_WritePatchVelAccel(pRegion)

    IMPLICIT NONE
      
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
  
! ==============================================================================
!   Arguments
! ==============================================================================
  
    TYPE(t_region), POINTER :: pRegion
  
! ==============================================================================
!   Locals
! ==============================================================================
    
    LOGICAL :: fileExists
    CHARACTER(CHRLEN) :: iFileName
    INTEGER :: errorFlag,iFile
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
  
    CALL RegisterFunction(global,'RFLU_MVF_WritePatchVelAccel',__FILE__)

! ******************************************************************************
!   Write velocity and Accelerations
! ******************************************************************************

    iFile = IF_VELACCEL

    IF ( global%myProcid == MASTERPROC ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN   
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing patch velocity and '// &
                                 'Acceleration...' 
      END IF ! global%verbLevel

      IF ( global%verbLevel > VERBOSE_LOW ) THEN   
        WRITE(STDOUT,'(A,3X,A,1X,I2)') SOLVER_NAME,'Patch:', &
                                       global%iPatchGlobalMvFrame
      END IF ! global%verbLevel

! ------------------------------------------------------------------------------
!     Open file
! ------------------------------------------------------------------------------

      CALL RFLU_MVF_OpenPatchVelAccel(global,FILE_POSITION_END,fileExists)

! ------------------------------------------------------------------------------
!     Write to file
! ------------------------------------------------------------------------------

      WRITE(iFile,1007,IOSTAT=errorFlag) global%currentTime, &
                     pRegion%mvfLoc(XCOORD),pRegion%mvfLoc(YCOORD), &
                     pRegion%mvfLoc(ZCOORD),pRegion%mvfVel(XCOORD), &
                     pRegion%mvfVel(YCOORD),pRegion%mvfVel(ZCOORD), &
                     pRegion%mvfAcc(XCOORD),pRegion%mvfAcc(YCOORD), &
                     pRegion%mvfAcc(ZCOORD)
 
1007 FORMAT(1PE24.16,9E24.16) 

! ------------------------------------------------------------------------------
!     Close file
! ------------------------------------------------------------------------------

      CALL RFLU_MVF_ClosePatchVelAccel(global)

      IF ( global%verbLevel > VERBOSE_NONE ) THEN       
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing patch velocity and '// & 
                                             'Acceleration done.'             
      END IF ! global%verbLevel
    END IF ! global%myProcid

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_MVF_WritePatchVelAccel






! ******************************************************************************
! End  
! ******************************************************************************

END MODULE RFLU_ModMovingFrame

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModMovingFrame.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:36  brollin
! New Stable version
!
! Revision 1.6  2009/09/28 14:21:50  mparmar
! Minor bug fixes
!
! Revision 1.5  2009/07/08 19:11:57  mparmar
! Adapted for SOLV_IMPLICIT_HM
!
! Revision 1.4  2008/12/06 08:43:41  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:16:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2008/05/29 01:35:29  mparmar
! Added capability for different types of motion and generic BC treatment
!
! Revision 1.1  2007/06/18 17:32:48  mparmar
! Initial revision
!
! ******************************************************************************

