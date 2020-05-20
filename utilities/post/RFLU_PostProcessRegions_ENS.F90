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
! Purpose: Wrapper for postprocessing results without merging the regions 
!   specific to EnSight.
!
! Description: None.
!
! Input: 
!   levels      Level data
!
! Output: None.
!
! Notes: 
!   1. This routine is specific to EnSight because the requirement that data be
!      written out for a certain number of servers means sharing code with 
!      routine for Tecplot becomes cumbersome.
!
! ******************************************************************************
!
! $Id: RFLU_PostProcessRegions_ENS.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PostProcessRegions_ENS(levels)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global 
  USE ModError
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModParameters

  USE RFLU_ModDimensions
  USE RFLU_ModENSIGHT
  USE RFLU_ModPlottingVars

  USE ModInterfaces, ONLY: RFLU_ClosePostInfo, &
                           RFLU_CreateGrid, & 
                           RFLU_DestroyGrid, & 
                           RFLU_OpenPostInfo, &
                           RFLU_PostProcessRegionsCommon1, & 
                           RFLU_PostProcessRegionsCommon2, &
                           RFLU_ReadPostInfo
    
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_level), POINTER :: levels(:)
  
! ==============================================================================
! Local variables
! ==============================================================================

  LOGICAL :: postInfoFileExists
  CHARACTER(CHRLEN) :: casename,RCSIdentString
  INTEGER :: errorFlag,iReg,iServer,iServerTarget
  TYPE(t_region), POINTER :: pRegion
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PostProcessRegions_ENS.F90,v $ $Revision: 1.1.1.1 $'

  global => levels(1)%regions(1)%global

  CALL RegisterFunction(global,'RFLU_PostProcessRegions_ENS',__FILE__)

! ******************************************************************************  
! Open post info file
! ******************************************************************************

  CALL RFLU_OpenPostInfo(global,FILE_STATUS_OLD,postInfoFileExists) 

! ******************************************************************************  
! Loop over servers
! ******************************************************************************

  DO iServer = 1,global%postNServers
    pRegion => levels(1)%regions(1)

    CALL RFLU_CountPlottingVars(pRegion)
    CALL RFLU_CreatePlottingVarMaps(pRegion)
    CALL RFLU_BuildPlottingVarMaps(pRegion)
    CALL RFLU_PrintPlottingVarsInfo(pRegion)

    CALL RFLU_ENS_BuildDataInfo(pRegion,iServer) 

    CALL RFLU_DestroyPlottingVarMaps(pRegion)

    CALL RFLU_ENS_WriteFileCase(global,iServer)
    CALL RFLU_ENS_InitPartNumber(global)

    CALL RFLU_ENS_OpenFileGeometry(global,iServer)
    CALL RFLU_ENS_OpenFileFlowWrapper(global,iServer)

! ==============================================================================  
!   Open post info
! ==============================================================================  

    CALL RFLU_OpenPostInfo(global,FILE_STATUS_OLD,postInfoFileExists) 

! ==============================================================================  
!   Loop over regions
! ==============================================================================  

    DO iReg = 1,global%nRegionsLocal
      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)
      ELSE 
        pRegion => levels(1)%regions(iReg)
      END IF ! global%nRegions

      IF ( postInfoFileExists .EQV. .TRUE. ) THEN 
        CALL RFLU_ReadPostInfo(pRegion,INFOFILE_READMODE_FLAG)    
      END IF ! postInfoFileExists  

! ------------------------------------------------------------------------------
!     Check for activation status
! ------------------------------------------------------------------------------

      IF ( pRegion%postActiveFlag .EQV. .TRUE. ) THEN 
        iServerTarget = RFLU_ENS_MapRegion2Server(iReg,global%nRegions, & 
                                                  global%postNServers)
            
! ----- If target server not equal to this server, write empty parts -----------            
            
        IF ( iServerTarget /= iServer ) THEN 
          CALL RFLU_ReadDimensionsWrapper(pRegion)
          CALL RFLU_CreateGrid(pRegion)
        
          CALL RFLU_ENS_StorePartNumber(pRegion%global)   
          CALL RFLU_ENS_WriteGridWrapper(pRegion,.TRUE.)    
          CALL RFLU_ENS_WriteFlowWrapper(pRegion,.TRUE.)
                   
          CALL RFLU_DestroyGrid(pRegion)                       

! ----- If target server equal to this server, write parts ---------------------

        ELSE
          CALL RFLU_PostProcessRegionsCommon1(pRegion,postInfoFileExists)
           
          CALL RFLU_ENS_StorePartNumber(pRegion%global)   
          CALL RFLU_ENS_WriteGridWrapper(pRegion)    
          CALL RFLU_ENS_WriteFlowWrapper(pRegion)                       

          CALL RFLU_PostProcessRegionsCommon2(pRegion)
        END IF ! ! iServerTarget
      END IF ! pRegion%postActiveFlag                       
    END DO ! iReg

! ==============================================================================  
!   Close files
! ==============================================================================  
   
    IF ( postInfoFileExists .EQV. .TRUE. ) THEN 
      CALL RFLU_ClosePostInfo(global) 
    END IF ! postInfoFileExists    

    CALL RFLU_ENS_CloseFileGeometry(global)
    CALL RFLU_ENS_CloseFileFlowWrapper(global)     
    CALL RFLU_ENS_DestroyDataInfo(global)   
  END DO ! iServer
       
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PostProcessRegions_ENS

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PostProcessRegions_ENS.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:58  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:12  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:58:09  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2007/03/19 21:46:05  haselbac
! Bug fix: Enabled use of plotting variables, done as part of changes to pv
!
! Revision 1.2  2006/01/24 21:21:15  mparmar
! Now pass pRegion instead of global into RFLU_ENS_BuildDataInfo
!
! Revision 1.1  2005/10/05 20:23:34  haselbac
! Initial revision
!
! ******************************************************************************

