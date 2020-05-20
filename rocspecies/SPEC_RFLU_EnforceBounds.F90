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
! Purpose: Enforce bounds on species.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
! 
! Output: None.
!
! Notes: 
!   1. At present, only enforce minimum bound. In future, could also enforce
!      maximum bound on each individual component, and that the sum is equal
!      to the mixture density.
!
!******************************************************************************
!
! $Id: SPEC_RFLU_EnforceBounds.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE SPEC_RFLU_EnforceBounds(pRegion)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region 
  USE ModMPI
  
  USE ModInterfaces, ONLY: RFLU_DecidePrint
    
  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Parameters
! ==============================================================================   

  TYPE(t_region), POINTER :: pRegion 

! ==============================================================================
! Local variables
! ==============================================================================
 
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: icg,iSpec,negValCntr
  REAL(RFREAL) :: negValMin
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvMixt,pCvSpec
  TYPE(t_grid), POINTER :: pGrid  
  TYPE(t_global), POINTER :: global
  
  RCSIdentString = '$RCSfile: SPEC_RFLU_EnforceBounds.F90,v $ $Revision: 1.1.1.1 $'
  
! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_EnforceBounds',__FILE__)

! ==============================================================================
! Set pointers and variables
! ==============================================================================

  pCvMixt => pRegion%mixt%cv
  pCvSpec => pRegion%spec%cv
  pGrid   => pRegion%grid
  
  negValCntr = 0
  negValMin  = HUGE(1.0_RFREAL)
  
! ******************************************************************************
! Enforce bounds
! ******************************************************************************

  DO icg = 1,pGrid%nCells
    DO iSpec = 1,pRegion%specInput%nSpecies
      IF ( pCvSpec(iSpec,icg) < 0.0_RFREAL ) THEN
        negValCntr = negValCntr + 1 
        negValMin = MIN(negValMin,pCvSpec(iSpec,icg))        
        pCvSpec(iSpec,icg) = 0.0_RFREAL                   
! DEBUG: Manoj: 2012-05-01
!WRITE(*,*) "-ve spec mass frac:   icg=",icg,"  iSpec=",iSpec
! END DEBUG
      END IF ! pCvSpec
    END DO ! iSpec
  END DO ! icg
           
! ******************************************************************************
! Write information
! ******************************************************************************
           
  IF ( (global%verbLevel > VERBOSE_LOW) .AND. &
       (global%myProcid == MASTERPROC) ) THEN 
    IF ( (RFLU_DecidePrint(global) .EQV. .TRUE.) .AND. (negValCntr > 0) ) THEN
      global%warnCounter = global%warnCounter + 1
     
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
            '*** WARNING *** Detected negative species mass fractions!'
            
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal 
      IF ( global%flowType == FLOW_UNSTEADY ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                            global%currentTime
      END IF ! global%flowType
      
      WRITE(STDOUT,'(A,3X,A,1X,I1)') SOLVER_NAME,'Runge-Kutta stage:', &
                                     pRegion%irkStep      

      WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME, & 
            'Number of locations with negative species mass fractions:', &
            negValCntr
      WRITE(STDOUT,'(A,3X,A,1X,E23.16)') SOLVER_NAME, &
            'Largest negative species mass fraction:',negValMin
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
            'Negative species mass fractions set to zero.'            
    END IF ! RFLU_DecidePrint
  END IF ! global%verbLevel           
           
! ******************************************************************************
! End
! ******************************************************************************  
      
  CALL DeregisterFunction(global)  
    
END SUBROUTINE SPEC_RFLU_EnforceBounds


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: SPEC_RFLU_EnforceBounds.F90,v $
!   Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
!   merged rocflu micro and macro
!
!   Revision 1.1.1.1  2014/07/15 14:31:38  brollin
!   New Stable version
!
!   Revision 1.3  2008/12/06 08:43:53  mtcampbe
!   Updated license.
!
!   Revision 1.2  2008/11/19 22:17:05  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.1  2007/04/09 18:51:23  haselbac
!   Initial revision after split from RocfloMP
!
!   Revision 1.1  2007/04/09 18:01:50  haselbac
!   Initial revision after split from RocfloMP
!
!   Revision 1.1  2004/01/29 22:59:29  haselbac
!   Initial revision
!
! ******************************************************************************

