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
! Purpose: Pick special faces.
!
! Description: None.
!
! Input:
!   pRegion		Pointer to region data
! 
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_PickSpecialFaces.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PickSpecialFaces(pRegion)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError  
  USE ModBndPatch, ONLY: t_patch
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModSortSearch

  USE RFLU_ModGrid

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion
 
! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER :: infoType,stencilType
  CHARACTER(CHRLEN) :: RCSIdentString  
  INTEGER :: errorFlag,faceIndx,iFacesSpecial,patchIndx
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PickSpecialFaces.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_PickSpecialFaces', &
                        __FILE__)

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Picking special faces...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                     pRegion%iRegionGlobal  
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and initialize
! ******************************************************************************

  pGrid => pRegion%grid

  iFacesSpecial = 0 
  pGrid%facesSpecial(1:2,1:NFACES_SPECIAL_MAX) = 0

! ******************************************************************************
! Get information from user
! ******************************************************************************

  WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter information on special faces:' 
  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'b - boundary face'    
  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'i - interior face'      
  WRITE(STDOUT,'(A,7X,A)') SOLVER_NAME,'q - quit'

! ******************************************************************************
! Set up infinite loop
! ******************************************************************************

  DO
  
! ==============================================================================
!   Enter information type
! ==============================================================================
   
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Enter information type:'
    READ(STDIN,'(A)') infoType
  
    SELECT CASE ( infoType )

! ------------------------------------------------------------------------------    
!     Boundary face
! ------------------------------------------------------------------------------

      CASE ( 'b' ) 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter patch index:'
        READ(STDIN,*,IOSTAT=errorFlag) patchIndx

        IF ( errorFlag /= ERR_NONE ) THEN 
          global%warnCounter = global%warnCounter + 1         
        
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.'
          CYCLE
        END IF ! errorFlag        
          
        IF ( patchIndx > 0 .AND. patchIndx <= pGrid%nPatches ) THEN
          pPatch => pRegion%patches(patchIndx)
         
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter face index:'
          READ(STDIN,*,IOSTAT=errorFlag) faceIndx           

          IF ( errorFlag /= ERR_NONE ) THEN
            global%warnCounter = global%warnCounter + 1           
           
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, & 
                                     '*** WARNING *** Invalid input.'
            CYCLE
          END IF ! errorFlag 

          IF ( faceIndx > 0 .AND. faceIndx <= pPatch%nBFacesTot ) THEN
            IF ( iFacesSpecial == NFACES_SPECIAL_MAX ) THEN 
              CALL ErrorStop(global,ERR_NFACES_SPECIAL_MAX,__LINE__)
            END IF ! iFacesSpecial         
           
            iFacesSpecial = iFacesSpecial + 1          
            pGrid%facesSpecial(1,iFacesSpecial) = patchIndx
            pGrid%facesSpecial(2,iFacesSpecial) = faceIndx            
  
            WRITE(STDOUT,'(A,5X,A,1X,I8)') SOLVER_NAME,'Added face:',faceIndx        
          ELSE 
            global%warnCounter = global%warnCounter + 1           
          
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, & 
                                     '*** WARNING *** Invalid input.' 
            CYCLE 
          END IF ! faceIndx        
        ELSE 
          global%warnCounter = global%warnCounter + 1         
        
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.' 
          CYCLE 
        END IF ! patchIndx
   
! ------------------------------------------------------------------------------
!     Interior face     
! ------------------------------------------------------------------------------

      CASE ( 'i' ) 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Enter interior face index:'
        READ(STDIN,*,IOSTAT=errorFlag) faceIndx
 
        IF ( errorFlag /= ERR_NONE ) THEN
          global%warnCounter = global%warnCounter + 1         
         
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.'
          CYCLE
        END IF ! errorFlag         
        
        IF ( faceIndx > 0 .AND. faceIndx <= pGrid%nFacesTot ) THEN
          IF ( iFacesSpecial == NFACES_SPECIAL_MAX ) THEN ! NOTE
            CALL ErrorStop(global,ERR_NFACES_SPECIAL_MAX,__LINE__)
          END IF ! iFacesSpecial       
         
          iFacesSpecial = iFacesSpecial + 1          
          pGrid%facesSpecial(1,iFacesSpecial) = 0
          pGrid%facesSpecial(2,iFacesSpecial) = faceIndx
                                        
          WRITE(STDOUT,'(A,5X,A,1X,I8)') SOLVER_NAME,'Added face:',faceIndx   
        ELSE 
          global%warnCounter = global%warnCounter + 1         
        
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.' 
          CYCLE 
        END IF ! faceIndx        
            
! ------------------------------------------------------------------------------
!     Quit       
! ------------------------------------------------------------------------------
      
      CASE ( 'q' )
        EXIT
        
! ------------------------------------------------------------------------------
!     Default        
! ------------------------------------------------------------------------------
        
      CASE DEFAULT 
        global%warnCounter = global%warnCounter + 1 
      
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'*** WARNING *** Invalid input.' 
        CYCLE 
    END SELECT  
  END DO ! <empty> 

! ******************************************************************************
! Set number of special faces
! ******************************************************************************

  pGrid%nFacesSpecial = iFacesSpecial

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Picking special faces done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PickSpecialFaces

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PickSpecialFaces.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:57  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:11  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:57:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/09/27 02:02:23  haselbac
! Initial revision
!
! ******************************************************************************

