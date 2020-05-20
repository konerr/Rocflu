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
! Purpose: Suite of routines to manage time.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModTime.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModTime

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  IMPLICIT NONE
   
  PRIVATE
  PUBLIC :: RFLU_GetTimeRK, & 
            RFLU_SetTimeRK
     
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
     
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModTime.F90,v $ $Revision: 1.1.1.1 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
 
 
 



! ******************************************************************************
!
! Purpose: Get time within RK stages.
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
    
  FUNCTION RFLU_GetTimeRK(global)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Function   
! ==============================================================================

    REAL(RFREAL) :: RFLU_GetTimeRK

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    RFLU_GetTimeRK = global%currentTimeRK
     
! ******************************************************************************
!   End
! ******************************************************************************

  END FUNCTION RFLU_GetTimeRK





! ******************************************************************************
!
! Purpose: Set time within RK stages.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   iStage      Index of RK stage
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_SetTimeRK(pRegion,iStage)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    INTEGER, INTENT(IN) :: iStage
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_SetTimeRK',__FILE__)  

! ******************************************************************************
!   
! ******************************************************************************

    global%currentTimeRK = & 
      global%currentTime + global%dtMin*(pRegion%mixtInput%trk(iStage) & 
                                       - pRegion%mixtInput%trk(     1))

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_SetTimeRK

 
 
 
 
! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModTime


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModTime.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:46  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:58  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:27  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:42  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.1  2005/03/31 17:22:54  haselbac
! Initial revision
!
! ******************************************************************************

