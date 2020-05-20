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
!*******************************************************************************
!
! Purpose: Suite of routines to build boundary-face gradient access list.
!
! Description: None.
!
! Notes: 
!   1. The routine which builds the access lists MUST be called after the face 
!      lists have been built. 
!
!*******************************************************************************
!
! $Id: RFLU_ModBFaceGradAccessList.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
!*******************************************************************************

MODULE RFLU_ModBFaceGradAccessList

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_NullifyBFaceGradAccessList, & 
            RFLU_CreateBFaceGradAccessList, & 
            RFLU_BuildBFaceGradAccessList, & 
            RFLU_DestroyBFaceGradAccessList
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModBFaceGradAccessList.F90,v $ $Revision: 1.1.1.1 $' 
              
! *****************************************************************************
! Routines
! *****************************************************************************

  CONTAINS
  


! ******************************************************************************
!   Nullify boundary-face gradient access list
! ******************************************************************************
  
    SUBROUTINE RFLU_NullifyBFaceGradAccessList(pRegion)

      IMPLICIT NONE

! ******************************************************************************
!     Declarations and definitions
! ******************************************************************************
      
! ==============================================================================
!     Arguments
! ==============================================================================
  
      TYPE(t_region), POINTER :: pRegion  
  
! ==============================================================================
!     Locals
! ==============================================================================  

      INTEGER :: errorFlag,iPatch
      TYPE(t_grid), POINTER :: pGrid
      TYPE(t_patch), POINTER :: pPatch
      TYPE(t_global), POINTER :: global

! ******************************************************************************
!     Start
! ******************************************************************************      
        
      global => pRegion%global  
        
      CALL RegisterFunction(global,'RFLU_NullifyBFaceGradAccessList',__FILE__) 
        
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN                 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
              'Nullifying boundary-face gradient access list...'
      END IF ! global%verbLevel

! ******************************************************************************
!     Set grid pointer
! ******************************************************************************

      pGrid => pRegion%grid     
                
! ******************************************************************************
!     Loop over patches
! ******************************************************************************
          
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

! TEMPORARY : removing usage of bf2bg from everywhere
!        NULLIFY(pPatch%bf2bg)
!        NULLIFY(pPatch%bf2bgTot)
      END DO ! iPatch          
        
! ******************************************************************************
!     End
! ******************************************************************************

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
              'Nullifying boundary-face gradient access list done.'
      END IF ! global%verbLevel  
  
      CALL DeregisterFunction(global)  
  
    END SUBROUTINE RFLU_NullifyBFaceGradAccessList




! ******************************************************************************
!   Create boundary-face gradient access list
! ******************************************************************************
  
    SUBROUTINE RFLU_CreateBFaceGradAccessList(pRegion)

      IMPLICIT NONE

! ******************************************************************************
!     Declarations and definitions
! ******************************************************************************
      
! ==============================================================================
!     Arguments
! ==============================================================================
  
      TYPE(t_region), POINTER :: pRegion  
  
! ==============================================================================
!     Locals
! ==============================================================================  

      INTEGER :: errorFlag,iPatch
      TYPE(t_grid), POINTER :: pGrid
      TYPE(t_patch), POINTER :: pPatch
      TYPE(t_global), POINTER :: global

! ******************************************************************************
!     Start
! ******************************************************************************      
        
      global => pRegion%global  
        
      CALL RegisterFunction(global,'RFLU_CreateBFaceGradAccessList',__FILE__) 
        
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN                 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
              'Creating boundary-face gradient access list...'
      END IF ! global%verbLevel

! ******************************************************************************
!     Set grid pointer
! ******************************************************************************

      pGrid => pRegion%grid     
                
! ******************************************************************************
!     Loop over patches
! ******************************************************************************
          
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

! TEMPORARY : removing usage of bf2bg from everywhere
!        ALLOCATE(pPatch%bf2bg(BF2BG_BEG:BF2BG_END),STAT=errorFlag)
!        global%error = errorFlag   
!        IF ( global%error /= ERR_NONE ) THEN 
!          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2bg')          
!        END IF ! global%error
        
! TEMPORARY : removing usage of bf2bg from everywhere
!        ALLOCATE(pPatch%bf2bgTot(BF2BG_BEG:BF2BG_END),STAT=errorFlag)
!        global%error = errorFlag   
!        IF ( global%error /= ERR_NONE ) THEN 
!          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2bgTot')          
!        END IF ! global%error        
      END DO ! iPatch        
        
! ******************************************************************************
!     End
! ******************************************************************************

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
              'Creating boundary-face gradient access list done.'
      END IF ! global%verbLevel  
  
      CALL DeregisterFunction(global)  
  
    END SUBROUTINE RFLU_CreateBFaceGradAccessList
    



  
! ******************************************************************************
!   Build boundary-face gradient access list
! ****************************************************************************** 
    
    SUBROUTINE RFLU_BuildBFaceGradAccessList(pRegion)

      IMPLICIT NONE

! ******************************************************************************
!     Declarations and definitions
! ******************************************************************************
      
! ==============================================================================
!     Arguments
! ==============================================================================

      TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!     Locals
! ==============================================================================

      INTEGER :: bf2bgLast,bf2bgTotLast,iPatch
      TYPE(t_global), POINTER :: global
      TYPE(t_grid), POINTER :: pGrid
      TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!     Start
! ******************************************************************************

      global => pRegion%global

      CALL RegisterFunction(global,'RFLU_BuildBFaceGradAccessList',__FILE__)

      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel > VERBOSE_NONE ) THEN    
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
              'Building boundary-face gradient access lists...'
      END IF ! global%verbLevel

! ******************************************************************************
!     Set grid pointer
! ******************************************************************************

      pGrid => pRegion%grid  

! ==============================================================================
!     Loop over patches
! ==============================================================================

      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

        IF ( iPatch /= 1 ) THEN    
! TEMPORARY : removing usage of bf2bg from everywhere
!          pPatch%bf2bg(BF2BG_BEG) = bf2bgLast + 1        
!          pPatch%bf2bg(BF2BG_END) = pPatch%bf2bg(BF2BG_BEG) & 
!                                  + pPatch%nBFaces - 1

! TEMPORARY : removing usage of bf2bg from everywhere
!          pPatch%bf2bgTot(BF2BG_BEG) = bf2bgTotLast + 1        
!          pPatch%bf2bgTot(BF2BG_END) = pPatch%bf2bgTot(BF2BG_BEG) & 
!                                     + pPatch%nBFacesTot - 1
        ELSE             
! TEMPORARY : removing usage of bf2bg from everywhere
!          pPatch%bf2bg(BF2BG_BEG) = 1 
!          pPatch%bf2bg(BF2BG_END) = pPatch%nBFaces
 
! TEMPORARY : removing usage of bf2bg from everywhere         
!          pPatch%bf2bgTot(BF2BG_BEG) = 1 
!          pPatch%bf2bgTot(BF2BG_END) = pPatch%nBFacesTot              
        END IF ! iPatch

! TEMPORARY : removing usage of bf2bg from everywhere
!        bf2bgLast    = pPatch%bf2bg(BF2BG_END) 
!        bf2bgTotLast = pPatch%bf2bgTot(BF2BG_END)            
      END DO ! iPatch

! ******************************************************************************
!     End
! ******************************************************************************

      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel > VERBOSE_NONE ) THEN  
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
              'Building boundary-face gradient access lists done.'
      END IF ! global%verbLevel

      CALL DeregisterFunction(global)

    END SUBROUTINE RFLU_BuildBFaceGradAccessList
    
      
      
  
 


! ******************************************************************************
!   Destroy boundary-face gradient access list
! ******************************************************************************
  
    SUBROUTINE RFLU_DestroyBFaceGradAccessList(pRegion)

      IMPLICIT NONE

! ******************************************************************************
!     Declarations and definitions
! ******************************************************************************
      
! ==============================================================================
!     Arguments
! ==============================================================================
  
      TYPE(t_region), POINTER :: pRegion  
  
! ==============================================================================
!     Locals
! ==============================================================================  

      INTEGER :: errorFlag,iPatch
      TYPE(t_grid), POINTER :: pGrid
      TYPE(t_patch), POINTER :: pPatch
      TYPE(t_global), POINTER :: global

! ******************************************************************************
!     Start
! ******************************************************************************      
        
      global => pRegion%global  
        
      CALL RegisterFunction(global,'RFLU_DestroyBFaceGradAccessList',__FILE__) 
        
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN                 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
              'Destroying boundary-face gradient access list...'
      END IF ! global%verbLevel

! ******************************************************************************
!     Set grid pointer
! ******************************************************************************

      pGrid => pRegion%grid     
                
! ******************************************************************************
!     Loop over patches
! ******************************************************************************
        
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)
      
! TEMPORARY : removing usage of bf2bg from everywhere
!        DEALLOCATE(pPatch%bf2bg,STAT=errorFlag)
!        global%error = errorFlag   
!        IF ( global%error /= ERR_NONE ) THEN 
!          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%bf2bg')          
!        END IF ! global%error
 
! TEMPORARY : removing usage of bf2bg from everywhere       
!        DEALLOCATE(pPatch%bf2bgTot,STAT=errorFlag)
!        global%error = errorFlag   
!        IF ( global%error /= ERR_NONE ) THEN 
!          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%bf2bgTot')          
!        END IF ! global%error        
      END DO ! iPatch          
        
! ******************************************************************************
!     Nullify memory
! ******************************************************************************

      CALL RFLU_NullifyBFaceGradAccessList(pRegion)

! ******************************************************************************
!     End
! ******************************************************************************

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
              'Destroying boundary-face gradient access list done.'
      END IF ! global%verbLevel  
  
      CALL DeregisterFunction(global)  
  
    END SUBROUTINE RFLU_DestroyBFaceGradAccessList
 




! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModBFaceGradAccessList


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLU_ModBFaceGradAccessList.F90,v $
!   Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
!   merged rocflu micro and macro
!
!   Revision 1.1.1.1  2014/07/15 14:31:37  brollin
!   New Stable version
!
!   Revision 1.3  2008/12/06 08:43:39  mtcampbe
!   Updated license.
!
!   Revision 1.2  2008/11/19 22:16:53  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.1  2007/04/09 18:49:23  haselbac
!   Initial revision after split from RocfloMP
!
!   Revision 1.1  2007/04/09 18:00:39  haselbac
!   Initial revision after split from RocfloMP
!
!   Revision 1.5  2006/08/19 15:39:21  mparmar
!   Removed bf2bg,bf2bgTot
!
!   Revision 1.4  2004/05/25 01:34:16  haselbac
!   Added code for bf2bgTot array; needed for Rocturb Sij access
!
!   Revision 1.3  2004/01/22 16:03:58  haselbac
!   Made contents of modules PRIVATE, only procs PUBLIC, to avoid errors on ALC and titan
!
!   Revision 1.2  2003/12/09 03:58:01  haselbac
!   Bug fix for building of bf2bg list (no virtual faces)
!
!   Revision 1.1  2003/11/03 03:51:54  haselbac
!   Initial revision
!
! ******************************************************************************

