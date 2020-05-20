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
! Purpose: Collection of routines for plotting.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: PLAG_ModPlotting.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
! ******************************************************************************

MODULE PLAG_ModPlotting

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
    RCSIdentString = '$RCSfile: PLAG_ModPlotting.F90,v $ $Revision: 1.1.1.1 $'        

! ==============================================================================
! Public functions
! ==============================================================================
  
  PUBLIC :: PLAG_CreatePlotFlags, &
            PLAG_DestroyPlotFlags, &
            PLAG_SetPlotFlags

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
! Purpose: Create plot flag.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_CreatePlotFlags(pRegion)

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
    TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
!   Start, set pointers
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'PLAG_CreatePlotFlags',__FILE__)

    pPlag => pRegion%plag

    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating particle plot flags...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    IF ( pPlag%nPcls > 0 ) THEN 
      ALLOCATE(pPlag%plotFlag(pPlag%nPcls),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%plotFlag')
      END IF ! global%error
    ELSE 
      NULLIFY(pPlag%plotFlag)
    END IF ! pPlag%nPcls 

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating particle plot flags done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_CreatePlotFlags

  
  
  
! *******************************************************************************
!
! Purpose: Destroy plot flag.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_DestroyPlotFlags(pRegion)

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
    TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
!   Start, set pointers
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'PLAG_DestroyPlotFlags',__FILE__)

    pPlag => pRegion%plag

    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying particle plot flags...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    IF ( pPlag%nPcls > 0 ) THEN 
      DEALLOCATE(pPlag%plotFlag,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPlag%plotFlag')
      END IF ! global%error
    END IF ! pPlag%nPcls 

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying particle plot '// & 
                               'flags done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_DestroyPlotFlags








! *******************************************************************************
!
! Purpose: Set plot flag.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_SetPlotFlags(pRegion)

    USE ModRandom, ONLY: Rand1Uniform

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

    INTEGER :: emptyLoopCounter,errorFlag,iPcl,iPcl2,nPclsPlot
    REAL(RFREAL) :: rn
    TYPE(t_global), POINTER :: global
    TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
!   Start, set pointers
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'PLAG_SetPlotFlags',__FILE__)

    pPlag => pRegion%plag

    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Setting particle plot flags...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel

! ******************************************************************************
!   Set nPclsPlot and set plot flags 
! ******************************************************************************

    pPlag%nPclsPlot = NINT(global%postPlagFrac*pPlag%nPcls)

    DO iPcl = 1,pPlag%nPcls
      pPlag%plotFlag(iPcl) = .FALSE.
    END DO ! iPcl    

    emptyLoopCounter = 0
    nPclsPlot = 0

    emptyLoop: DO 
      emptyLoopCounter = emptyLoopCounter + 1 

      rn = Rand1Uniform(pRegion%randData)

      iPcl2 = 1 + NINT(rn*(pPlag%nPcls-1))

      IF ( pPlag%plotFlag(iPcl2) .EQV. .FALSE. ) THEN 
        pPlag%plotFlag(iPcl2) = .TRUE.

        nPclsPlot = nPclsPlot + 1

        IF ( nPclsPlot == pPlag%nPclsPlot ) THEN 
          EXIT emptyLoop
        END IF ! nPclsPlot
      END IF ! pPlag%plotFlag

      IF ( emptyLoopCounter == 100*pPlag%nPcls ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! emptyLoopCounter
    END DO emptyLoop

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Setting particle plot flags done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_SetPlotFlags






END MODULE PLAG_ModPlotting

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModPlotting.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:51  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/12 18:06:59  haselbac
! Initial revision
!
! ******************************************************************************

