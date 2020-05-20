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
! Purpose: Collection of routines for particle statistics on Eulerian grid.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: PLAG_ModStats.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE PLAG_ModStats

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY: t_plag 
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModMPI
  USE PLAG_ModParameters

  IMPLICIT NONE
  
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: PLAG_ModStats.F90,v $ $Revision: 1.1.1.1 $'        

! ==============================================================================
! Public functions
! ==============================================================================
  
  PUBLIC :: PLAG_CreateStat,  & 
            PLAG_DestroyStat, & 
            PLAG_InitStat

! ==============================================================================
! Private functions
! ==============================================================================

!  PRIVATE :: 

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  
  
  
  
  
  
  
  

  





! *******************************************************************************
!
! Purpose: Create statistics.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPlag               Pointer to plag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_CreateStat(pRegion,pPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion

    TYPE(t_plag), POINTER :: pPlag

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,ibc,iec,nTav
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_CreateStat',__FILE__)

    nTav =global%plagNStat     

    pGrid => pRegion%grid
    ibc = 1
    iec = pGrid%nCells
  
! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    ALLOCATE(pPlag%tav(nTav,ibc:iec),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%tav')
    END IF ! global%error 

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_CreateStat 






! *******************************************************************************
!
! Purpose: Destroy statistics.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPlag               Pointer to plag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_DestroyStat(pRegion,pPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion

    TYPE(t_plag), POINTER :: pPlag
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag
    TYPE(t_global), POINTER :: global    

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_DestroyStat',__FILE__)

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    DEALLOCATE(pPlag%tav,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPlag%tav')
    END IF ! global%error
  
! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    NULLIFY(pPlag%tav)
 
! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_DestroyStat  





! *******************************************************************************
!
! Purpose: Initialize statistics.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPlag               Pointer to plag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_InitStat(pRegion,pPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
     
    TYPE(t_region), POINTER :: pRegion

    TYPE(t_plag), POINTER :: pPlag
  
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,ibc,iec,iCell,iVar,nVars
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'PLAG_InitStat',__FILE__)

    nVars = global%plagNStat      

    pGrid => pRegion%grid
    ibc = 1
    iec = pGrid%nCells

! ******************************************************************************
!   Initialize memory
! ******************************************************************************
        
    DO iCell = ibc,iec
      DO iVar = 1, nVars
        pPlag%tav(iVar,iCell) = 0.0_RFREAL
      END DO ! iVar                  
    END DO ! iCell                      

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_InitStat








END MODULE PLAG_ModStats

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModStats.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:51  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/04/16 23:21:41  fnajjar
! Deleted all RFLO-pertinent calls and statements
!
! Revision 1.1  2007/04/09 18:50:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.1  2005/01/08 20:44:32  fnajjar
! Initial import for PLAG statistics
!
! ******************************************************************************

