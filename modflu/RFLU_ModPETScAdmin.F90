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
! Purpose: Suite of routines to handle interaction with PETSc.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModPETScAdmin.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModPETScAdmin

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  IMPLICIT NONE
   
  PRIVATE
  PUBLIC ::  RFLU_PETSC_Finalize, & 
             RFLU_PETSC_Init
      
#include "include/finclude/petsc.h"    
         
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
     
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModPETScAdmin.F90,v $ $Revision: 1.1.1.1 $'               
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  




! ******************************************************************************
!
! Purpose: Finalize PETSc.
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

  SUBROUTINE RFLU_PETSC_Finalize(global)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_global), POINTER :: global
    
! ==============================================================================
!   Locals      
! ==============================================================================

    INTEGER :: errorFlag

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_PETSC_Finalize',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Finalizing PETSc...' 
    END IF ! global%verbLevel

! ******************************************************************************
!   Initialize PETSc
! ******************************************************************************

    CALL PetscFinalize(errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error    

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Finalizing PETSc done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PETSC_Finalize






! ******************************************************************************
!
! Purpose: Initialize PETSc.
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

  SUBROUTINE RFLU_PETSC_Init(global)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_global), POINTER :: global
    
! ==============================================================================
!   Locals      
! ==============================================================================

    INTEGER :: errorFlag

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_PETSC_Init',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing PETSc...' 
    END IF ! global%verbLevel

! ******************************************************************************
!   Initialize PETSc
! ******************************************************************************

    CALL PetscInitialize(PETSC_NULL_CHARACTER,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_PETSC_OUTPUT,__LINE__)
    END IF ! global%error    

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing PETSc done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_PETSC_Init
  
  
  

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModPETScAdmin


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModPETScAdmin.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:43  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:41  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/19 15:40:45  haselbac
! Initial revision
!
! ******************************************************************************

