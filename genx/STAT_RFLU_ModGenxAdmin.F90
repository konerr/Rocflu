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
! Purpose: Collection of routines related to GENX interaction.
!
! Description: None
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: STAT_RFLU_ModGenxAdmin.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE STAT_RFLU_ModGenxAdmin

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid 
  USE ModMPI 

  USE RFLU_ModGENXUtils

  IMPLICIT NONE

  INCLUDE 'roccomf90.h'

  PRIVATE

  PUBLIC :: STAT_RFLU_GenxCreateAttr, &
            STAT_RFLU_GenxRegisterData, &
            STAT_RFLU_GenxGetData

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================  
! Private
! ==============================================================================  

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: STAT_RFLU_ModGenxAdmin.F90,v $ $Revision: 1.1.1.1 $'       
  
! ==============================================================================  
! Public
! ==============================================================================  
  

! ******************************************************************************
! Contained routines
! ******************************************************************************

  CONTAINS 
  
  
  
! ******************************************************************************
!
! Purpose: Create new attributes for time-averaged statistics.
!
! Description: None.
!
! Input:
!   pRegion		Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE STAT_RFLU_GenxCreateAttr(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Loop
! ==============================================================================  
    INTEGER :: iStat

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: winv
    CHARACTER(CHRLEN), POINTER :: statNm(:,:,:)
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global
  
! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE STAT_RFLU_GenxCreateAttr
  
  
  
  
  


! ******************************************************************************
!
! Purpose: Register statistics data.
!
! Description: None.
!
! Input:
!   pRegion		Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE STAT_RFLU_GenxRegisterData(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Loop
! ==============================================================================  
    INTEGER :: iStat

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: winv
    CHARACTER(CHRLEN), POINTER :: statNm(:,:,:)
    INTEGER :: paneId, ilb
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering statistics data...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Create new attributes
! ******************************************************************************

! ==============================================================================  
!   Volume
! ==============================================================================  

    winv = global%volWinName
    IF ((global%flowType == FLOW_UNSTEADY).AND.(global%doStat == ACTIVE)) THEN

! ------------------------------------------------------------------------------
!     Global
! ------------------------------------------------------------------------------
      CALL COM_new_attribute( TRIM(winv)//'.tStat','w',COM_DOUBLE,1,'s' )
      CALL COM_set_array( TRIM(winv)//'.tStat',0, global%integrTime )

! ------------------------------------------------------------------------------
!     Mixture
! ------------------------------------------------------------------------------
      IF (global%mixtNStat > 0) THEN
        statNm => global%mixtStatNm
        DO iStat=1,global%mixtNStat
          CALL COM_new_attribute( TRIM(winv)//'.'//TRIM(statNm(1,1,iStat)),'e',&
                                  COM_DOUBLE,1,TRIM(statNm(1,2,iStat)) )
        ENDDO
      ENDIF

#ifdef TURB
      IF ((global%turbActive .EQV. .true.) .AND. (global%turbNStat > 0)) THEN
        statNm => global%turbStatNm
        DO iStat=1,global%turbNStat
          CALL COM_new_attribute( TRIM(winv)//'.'//TRIM(statNm(1,1,iStat)),'e',&
                                  COM_DOUBLE,1,TRIM(statNm(1,2,iStat)) )
        ENDDO
      ENDIF  ! turbNStat
#endif
    ENDIF   ! unsteady and dostat

    
! ******************************************************************************
!   Register data
! ******************************************************************************

! ==============================================================================  
!   Volume
! ==============================================================================  

    winv = global%volWinName       

    CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,0,paneId) 

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'Window name:',TRIM(winv)
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Pane id:', paneId           
    END IF ! global%verbLevel
   
    IF ((global%flowType==FLOW_UNSTEADY) .AND. (global%doStat==ACTIVE)) THEN
      ilb = 1
! ------------------------------------------------------------------------------
!     Mixture statistics
! ------------------------------------------------------------------------------   
      IF (global%mixtNStat > 0) THEN
        statNm => global%mixtStatNm
        DO iStat=1,global%mixtNStat
          CALL COM_set_array( TRIM(winv)//'.'//TRIM(statNm(1,1,iStat)),paneId,&
             pRegion%mixt%tav(iStat,ilb), global%mixtNStat)
        ENDDO
      ENDIF  ! mixtNStat

! ------------------------------------------------------------------------------
!     Turbulence statistics
! ------------------------------------------------------------------------------   
#ifdef TURB
      IF ((global%turbActive .EQV. .true.) .AND. (global%turbNStat > 0)) THEN
        statNm => global%turbStatNm
        DO iStat=1,global%turbNStat
          CALL COM_set_array( TRIM(winv)//'.'//TRIM(statNm(1,1,iStat)),paneId,&
             pRegion%turb%tav(iStat,ilb), global%turbNStat)
        ENDDO
      ENDIF  ! turbNStat
#endif   
    ENDIF    ! flowType and doStat

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering statistics data done.'
    END IF ! global%verbLevel 
  
  END SUBROUTINE STAT_RFLU_GenxRegisterData








! ******************************************************************************
!
! Purpose: Get statistics data through Roccom.
!
! Description: None.
!
! Input:
!   pRegion		Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE STAT_RFLU_GenxGetData(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Loop
! ==============================================================================  
    INTEGER :: iStat

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: winName,winNameIn
    CHARACTER(CHRLEN), POINTER :: statNm(:,:,:)
    INTEGER :: handleIn,handleObtain,handleOut
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global

    handleObtain = global%handleObtain

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting statistics data...'
    END IF ! global%verbLevel 
        
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
    END IF ! global%myProcid

! ******************************************************************************
!   Get data
! ******************************************************************************

    winNameIn = global%volWinNameInput
    winName   = global%volWinName  
    
! ==============================================================================
!   Global
! ==============================================================================   

    handleIn = COM_get_attribute_handle_const(TRIM(winNameIn)//'.tStat')
    handleOut = COM_get_attribute_handle(TRIM(winName)//'.tStat')
    CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    
! ==============================================================================
!   Statistics variables
! ==============================================================================   
    
! ------------------------------------------------------------------------------
!   Mixture statistics
! ------------------------------------------------------------------------------

    IF (global%mixtNStat > 0) THEN
      statNm => global%mixtStatNm
      DO iStat=1,global%mixtNStat

        handleIn = COM_get_attribute_handle_const(TRIM(winNameIn)//'.'// &
                                                  TRIM(statNm(1,1,iStat)))
        handleOut = COM_get_attribute_handle(TRIM(winName)//'.'// &
                                             TRIM(statNm(1,1,iStat)))
        CALL COM_call_function(handleObtain,2,handleIn,handleOut)

      ENDDO
    ENDIF  ! mixtNStat

    
! ------------------------------------------------------------------------------
!   Turbulence statistics
! ------------------------------------------------------------------------------

#ifdef TURB
    IF (global%turbNStat > 0) THEN
      statNm => global%turbStatNm
      DO iStat=1,global%turbNStat

        handleIn = COM_get_attribute_handle_const(TRIM(winNameIn)//'.'// &
                                                  TRIM(statNm(1,1,iStat)))
        handleOut = COM_get_attribute_handle(TRIM(winName)//'.'// &
                                             TRIM(statNm(1,1,iStat)))
        CALL COM_call_function(handleObtain,2,handleIn,handleOut)

      ENDDO
    ENDIF  ! turbNStat
#endif   
    
! ==============================================================================
!   Set variables back to accumulated mode
! ==============================================================================   
    
! ------------------------------------------------------------------------------
!   Mixture statistics
! ------------------------------------------------------------------------------

    IF (global%mixtNStat > 0) THEN
      DO iStat=1,global%mixtNStat
        pRegion%mixt%tav(iStat,:) = &
        pRegion%mixt%tav(iStat,:)*global%integrTime
      ENDDO
    ENDIF ! mixtNstat
    
! ------------------------------------------------------------------------------
!   Turbulence statistics
! ------------------------------------------------------------------------------

#ifdef TURB
    IF ((global%turbActive .EQV. .true.) .AND. (global%turbNStat > 0)) THEN
      DO iStat=1,global%turbNStat
        pRegion%turb%tav(iStat,:) = &
        pRegion%turb%tav(iStat,:)*global%integrTime
      ENDDO
    ENDIF ! turbNstat
#endif

! ******************************************************************************
!   End
! ******************************************************************************
 
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting statistics data done.'
    END IF ! global%verbLevel   
  
  END SUBROUTINE STAT_RFLU_GenxGetData






! ******************************************************************************
! End Module
! ******************************************************************************




END MODULE STAT_RFLU_ModGenxAdmin

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: STAT_RFLU_ModGenxAdmin.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:30  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:47:51  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:58:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.6  2006/02/05 03:46:00  wasistho
! added global%turbActive
!
! Revision 1.5  2006/01/10 06:32:39  wasistho
! mixt%tav to turb%tav within turb stats
!
! Revision 1.4  2006/01/10 05:00:15  wasistho
! added GenxGetData
!
! Revision 1.3  2006/01/04 20:06:13  wasistho
! modified data registration
!
! Revision 1.2  2006/01/03 09:51:43  wasistho
! get rid of ifdef rflu
!
! Revision 1.1  2006/01/03 06:34:41  wasistho
! initial import
!
!
! ******************************************************************************

