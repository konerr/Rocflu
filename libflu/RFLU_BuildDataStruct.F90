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
! Purpose: Prepare data structure for ROCFLU.
!
! Description: None.
!
! Input: 
!   global              Pointer to global type
!
! Output: None.
!
! Notes: 
!   1. Allocation for regions on each level starts from 0, which is used 
!      to store undecomposed grid on finest level
!   2. The allocation of the regions was removed from this routine and put 
!      into separate routine because compilation of this routine with 
!      optimization took a very long time and it was determined that the 
!      compilation of the allocation statement was responsible. This was
!      likely due to compiler trying to optimize the allocation as it was 
!      inside a DO-loop. Putting allocation into separate routine (without 
!      DO-loop) will break optimization. 
!
! ******************************************************************************
!
! $Id: RFLU_BuildDataStruct.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_BuildDataStruct(global,levels)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModError
  USE ModMPI

  USE ModInterfaces, ONLY: RFLU_CreateRegions

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global
  TYPE(t_level), DIMENSION(:), POINTER :: levels
  TYPE(t_region), DIMENSION(:), POINTER :: pRegions

! ==============================================================================
! Local variables
! ==============================================================================

  INTEGER :: errorFlag,iLev,iReg,iRegLow,iRegUpp
  CHARACTER(CHRLEN) :: RCSIdentString
  
! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_BuildDataStruct.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction(global,'RFLU_BuildDataStruct',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building data structure...'
  END IF ! global%verbLevel
  
  ALLOCATE(levels(global%nLevels),STAT=errorFlag)
  global%error = errorFlag    
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'levels')
  END IF ! global%error
  
! ******************************************************************************
! Loop over levels and allocate memory for regions
! ******************************************************************************
  
  DO iLev = 1,global%nLevels       
    CALL RFLU_CreateRegions(global,iLev,levels)
  END DO ! iLev    

! ******************************************************************************
! Set global pointer and postActiveFlag
! ******************************************************************************

  DO iLev = 1,global%nLevels
    pRegions => levels(iLev)%regions

    DO iReg = LBOUND(pRegions,1),UBOUND(pRegions,1)
      pRegions(iReg)%pRegion => pRegions(iReg)    
      pRegions(iReg)%global  => global
      
      pRegions(iReg)%postActiveFlag = .TRUE.
    END DO ! iReg
  END DO ! iLev
    
! ******************************************************************************
! End
! ******************************************************************************
      
  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building data structure done.'    
  END IF ! global%verbLevel
  
  CALL DeregisterFunction(global)  
    
END SUBROUTINE RFLU_BuildDataStruct

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_BuildDataStruct.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:34  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:52  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:48  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.6  2004/10/19 19:24:06  haselbac
! Cosmetics only
!
! Revision 1.5  2003/08/07 15:30:09  haselbac
! Changed var name
!
! Revision 1.4  2003/07/22 01:54:26  haselbac
! Added setting of region%pRegion
!
! Revision 1.3  2003/06/09 14:02:40  haselbac
! Rewrite to give faster compilation with optimization
!
! Revision 1.2  2003/06/04 21:58:30  haselbac
! Added setting of activeFlag (for rflupost)
!
! Revision 1.1  2003/01/28 15:53:31  haselbac
! Initial revision
!
! Revision 1.11  2002/10/27 18:52:28  haselbac
! Removed tabs
!
! Revision 1.10  2002/10/08 15:48:56  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.9  2002/10/05 18:47:44  haselbac
! Integration into GENX: Now also pass levels, some cosmetics
!
! Revision 1.8  2002/09/09 14:15:01  haselbac
! global now under regions
!
! Revision 1.7  2002/07/25 14:42:23  haselbac
! Now only write out messages for MASTERPROC
!
! Revision 1.6  2002/06/17 13:31:22  haselbac
! Prefixed SOLVER_NAME to all screen output
!
! Revision 1.5  2002/06/14 20:08:36  haselbac
! ModLocal deleted, so local% becomes global%
!
! Revision 1.4  2002/05/04 16:12:04  haselbac
! Clean up and added type local
!
! Revision 1.3  2002/03/26 18:54:50  haselbac
! Added IF statement for verbosity level check
!
! Revision 1.2  2002/03/01 16:06:02  haselbac
! Deleted previous region and renamed domain to region
!
! Revision 1.1  2002/02/08 14:54:30  haselbac
! Moved to libflu from rfluprep
!
! Revision 1.2  2002/01/14 20:26:44  haselbac
! Added RCS revision history section...
!
! ******************************************************************************


