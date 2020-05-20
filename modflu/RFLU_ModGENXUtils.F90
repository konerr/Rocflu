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
! Purpose: Utility routines for Rocflu-GENX interface. 
!
! Description: None
!
! Notes: 
!   1. This module must only be compiled when GENX=1. 
!
! ******************************************************************************
!
! $Id: RFLU_ModGENXUtils.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModGENXUtils

  USE ModDataTypes
  USE ModParameters

  IMPLICIT NONE

  INCLUDE 'roccomf90.h'

  PRIVATE
  PUBLIC :: RFLU_GENX_BuildPaneId, &  
            RFLU_GENX_BuildTimeString, &
            RFLU_GENX_GetFileStemSurf, &
            RFLU_GENX_GetFileStemVol, &            
            RFLU_GENX_GetFileStemGridSurf, &            
            RFLU_GENX_GetFileStemGridVol, &
            RFLU_GENX_GetFileStemGSpSurf, &            
            RFLU_GENX_GetFileStemGSpVol, &
            RFLU_GENX_GetFileStemMixtSurf, &            
            RFLU_GENX_GetFileStemMixtVol, & 
            RFLU_GENX_SetCnstrType     

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================  
! Private
! ==============================================================================  

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModGENXUtils.F90,v $ $Revision: 1.1.1.1 $'       
  
! ==============================================================================  
! Public
! ==============================================================================  
  
  INTEGER, PARAMETER, PUBLIC :: GENX_TIME_STRING_LEN = 9

! ******************************************************************************
! Contained routines
! ******************************************************************************

  CONTAINS 
  
  


! ******************************************************************************
!
! Purpose: Build unique identifier for each pane.
!
! Description: None.
!
! Input:
!   iRegion     Region index
!   iPatch      Patch index
!
! Output: 
!   paneId      Pane identifier
!
! Notes: 
!   1. Pane id must not be zero, so add unit offset.
!   2. If the mapping is changed, the inverse mapping in the routine
!      RFLU_GENX_GetPatchId must also be changed.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_BuildPaneId(iRegion,iPatch,paneId)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    INTEGER, INTENT(IN) :: iPatch,iRegion
    INTEGER, INTENT(OUT) :: paneId
    
! ******************************************************************************
!   Start
! ******************************************************************************
    
    paneId = 1 + iRegion*REGION_INDEX_OFFSET + iPatch

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_BuildPaneId 






! ******************************************************************************
!
! Purpose: Build unique identifier for each pane.
!
! Description: None.
!
! Input:
!   iRegion     Region index
!   iPatch      Patch index
!
! Output: 
!   paneId      Pane identifier
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_BuildTimeString(time,timeString)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    REAL(RFREAL), INTENT(IN) :: time
    CHARACTER(GENX_TIME_STRING_LEN), INTENT(OUT) :: timeString

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(15) :: tempString
    
! ******************************************************************************
!   Start
! ******************************************************************************
    
    WRITE(tempString,'(E12.6)') 1.0E9_RFREAL*time
    
    timeString = tempString(11:12)//'.'//tempString(3:8)    

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_BuildTimeString 





! ******************************************************************************
!
! Purpose: Get file stem for surface files.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   fileStem    File stem 
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_GetFileStemSurf(fileStem)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    CHARACTER(*), INTENT(OUT) :: fileStem

! ******************************************************************************
!   Start
! ******************************************************************************
    
    fileStem = 'ifluid'

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_GetFileStemSurf







! ******************************************************************************
!
! Purpose: Get file stem for volume files.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   fileStem    File stem 
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_GetFileStemVol(fileStem)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    CHARACTER(*), INTENT(OUT) :: fileStem

! ******************************************************************************
!   Start
! ******************************************************************************
    
    fileStem = 'fluid'

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_GetFileStemVol









! ******************************************************************************
!
! Purpose: Get file stem for surface grid file.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   fileStem    File stem 
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_GetFileStemGridSurf(fileStem)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    CHARACTER(*), INTENT(OUT) :: fileStem

! ******************************************************************************
!   Start
! ******************************************************************************
    
    CALL RFLU_GENX_GetFileStemSurf(fileStem)
    
    fileStem = TRIM(fileStem)//'-grid'

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_GetFileStemGridSurf








! ******************************************************************************
!
! Purpose: Get file stem for volume grid file.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   fileStem    File stem 
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_GetFileStemGridVol(fileStem)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    CHARACTER(*), INTENT(OUT) :: fileStem

! ******************************************************************************
!   Start
! ******************************************************************************
    
    CALL RFLU_GENX_GetFileStemVol(fileStem)    
    
    fileStem = TRIM(fileStem)//'-grid'

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_GetFileStemGridVol









! ******************************************************************************
!
! Purpose: Get file stem for surface grid speed file.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   fileStem    File stem 
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_GetFileStemGSpSurf(fileStem)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    CHARACTER(*), INTENT(OUT) :: fileStem

! ******************************************************************************
!   Start
! ******************************************************************************
    
    CALL RFLU_GENX_GetFileStemSurf(fileStem)
    
    fileStem = TRIM(fileStem)//'-gsp'

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_GetFileStemGSpSurf








! ******************************************************************************
!
! Purpose: Get file stem for volume grid speed file.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   fileStem    File stem 
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_GetFileStemGSpVol(fileStem)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    CHARACTER(*), INTENT(OUT) :: fileStem

! ******************************************************************************
!   Start
! ******************************************************************************
    
    CALL RFLU_GENX_GetFileStemVol(fileStem)    
    
    fileStem = TRIM(fileStem)//'-gsp'

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_GetFileStemGSpVol







! ******************************************************************************
!
! Purpose: Get file stem for mixture solution file.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   fileStem    File stem 
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_GetFileStemMixtSurf(fileStem)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    CHARACTER(*), INTENT(OUT) :: fileStem

! ******************************************************************************
!   Start
! ******************************************************************************
    
    CALL RFLU_GENX_GetFileStemSurf(fileStem)     
    
    fileStem = TRIM(fileStem)//'-mixt'

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_GetFileStemMixtSurf









! ******************************************************************************
!
! Purpose: Get file stem for mixture solution file.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   fileStem    File stem 
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_GetFileStemMixtVol(fileStem)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    CHARACTER(*), INTENT(OUT) :: fileStem

! ******************************************************************************
!   Start
! ******************************************************************************
    
    CALL RFLU_GENX_GetFileStemVol(fileStem)     
    
    fileStem = TRIM(fileStem)//'-mixt'

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_GetFileStemMixtVol






! ******************************************************************************
!
! Purpose: Set constraint type.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   movePatchDir        Direction in which patch is allowed to move.
!
! Notes: None.
!
! ******************************************************************************

  INTEGER FUNCTION RFLU_GENX_SetCnstrType(movePatchDir)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    INTEGER, INTENT(IN) :: movePatchDir

! ******************************************************************************
!   Start
! ******************************************************************************
    
    SELECT CASE ( movePatchDir ) 
      CASE ( MOVEPATCH_DIR_NONE ) 
        RFLU_GENX_SetCnstrType =    2
      CASE ( MOVEPATCH_DIR_X ) 
        RFLU_GENX_SetCnstrType =  120
      CASE ( MOVEPATCH_DIR_Y ) 
        RFLU_GENX_SetCnstrType =  121
      CASE ( MOVEPATCH_DIR_Z ) 
        RFLU_GENX_SetCnstrType =  122   
      CASE ( MOVEPATCH_DIR_XY ) 
        RFLU_GENX_SetCnstrType = -122
      CASE ( MOVEPATCH_DIR_XZ ) 
        RFLU_GENX_SetCnstrType = -121
      CASE ( MOVEPATCH_DIR_YZ ) 
        RFLU_GENX_SetCnstrType = -120
      CASE ( MOVEPATCH_DIR_XYZ ) 
        RFLU_GENX_SetCnstrType =    0                   
      CASE DEFAULT 
        RFLU_GENX_SetCnstrType =    2
    END SELECT ! movePatchDir

! ******************************************************************************
!   End
! ******************************************************************************
  
  END FUNCTION RFLU_GENX_SetCnstrType






END MODULE RFLU_ModGENXUtils

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModGENXUtils.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:40  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:40  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.2  2005/06/09 20:21:17  haselbac
! Added function to set constraint type
!
! Revision 1.1  2004/10/19 19:27:04  haselbac
! Initial revision
!
! ******************************************************************************

