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
! Purpose: Display minimum and maximum values of species state vector for 
!   given local domain.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
! 
! Output: N/A.
!
! Notes: None. 
!
!******************************************************************************
!
! $Id: SPEC_RFLU_PrintFlowInfo.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE SPEC_RFLU_PrintFlowInfo(pRegion)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region 
    
  USE ModInterfaces, ONLY: RFLU_PrintLocInfo   

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
  INTEGER :: errorFlag,iSpec,nSpecies,uppLim
  INTEGER :: dummy(1) 
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: loc
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
  TYPE(t_grid), POINTER :: pGrid  
  TYPE(t_global), POINTER :: global
  
  RCSIdentString = '$RCSfile: SPEC_RFLU_PrintFlowInfo.F90,v $ $Revision: 1.1.1.1 $'
  
! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_PrintFlowInfo',__FILE__)

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing species information...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                     pRegion%iRegionGlobal 
    IF ( global%flowType == FLOW_UNSTEADY ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                          global%currentTime
    END IF ! global%flowType
  END IF ! global%verbLevel

! ==============================================================================
! Set pointers and variables
! ==============================================================================

  pCv   => pRegion%spec%cv
  pGrid => pRegion%grid

! ==============================================================================
! Allocate memory for location array
! ==============================================================================

  nSpecies = pRegion%specInput%nSpecies
  
  ALLOCATE(loc(nSpecies,MIN_VAL:MAX_VAL),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'loc')
  END IF ! global%error

! ==============================================================================
! Set upper limit: useful for checking of cell/dummy cell values
! ==============================================================================

  uppLim = pGrid%nCells ! Only internal cells
!  uppLim = pGrid%nCellsTot ! Also dummy cells

! ==============================================================================
! Find locations of extrema: NOTE Asinine coding needed because of poor 
! FORTRAN interface for MINLOC and MAXLOC functions... 
! ==============================================================================

  DO iSpec = 1,nSpecies
    dummy              = MINLOC(pCv(iSpec,1:uppLim))
    loc(iSpec,MIN_VAL) = dummy(1)

    dummy              = MAXLOC(pCv(iSpec,1:uppLim))
    loc(iSpec,MAX_VAL) = dummy(1)
  END DO ! iSpec
  
! ==============================================================================
! Print locations of extrema
! ==============================================================================  
  
  DO iSpec = 1,nSpecies
    WRITE(STDOUT,'(A,3X,A,2(1X,E23.16),2(1X,I9))') & 
                  SOLVER_NAME,'Density (kg/m^3):', & 
                  MINVAL(pCv(iSpec,1:uppLim)),MAXVAL(pCv(iSpec,1:uppLim)), & 
                  loc(iSpec,MIN_VAL),loc(iSpec,MAX_VAL)
  END DO ! iSpec  
    
! ==============================================================================
! Print out locations of cells at which extrema occur
! ==============================================================================  
     
  IF ( global%verbLevel /= VERBOSE_LOW ) THEN   
    CALL RFLU_PrintLocInfo(pRegion,loc,nSpecies,LOCINFO_MODE_SILENT, &
                           OUTPUT_MODE_MASTER_ONLY)
  END IF ! global%verbLevel
      
! ==============================================================================
! Deallocate memory for location array
! ==============================================================================
  
  DEALLOCATE(loc,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'loc')
  END IF ! global%error      
      
! ******************************************************************************
! End
! ******************************************************************************  
  
  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing species information done.'
  END IF ! global%verbLevel  
    
  CALL DeregisterFunction(global)  
    
END SUBROUTINE SPEC_RFLU_PrintFlowInfo


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: SPEC_RFLU_PrintFlowInfo.F90,v $
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
!   Revision 1.1  2003/11/25 21:08:37  haselbac
!   Initial revision
!
! ******************************************************************************

