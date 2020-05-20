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
! Purpose: Read grid file and convert if necessary.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ReadConvGridWrapper.F90,v 1.2 2015/07/23 23:11:19 brollin Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ReadConvGridWrapper(pRegion)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region

  USE RFLU_ModCENTAUR
  USE RFLU_ModCOBALT
  USE RFLU_ModGAMBIT
  USE RFLU_ModMESH3D
  USE RFLU_ModOMG
  USE RFLU_ModSTARCD
  USE RFLU_ModTETMESH
  USE RFLU_ModVGRIDns

#ifdef GENX
  USE RFLU_ModGENXIO, ONLY: RFLU_GENX_DecideReadFile, & 
                            RFLU_GENX_GetGrid
#endif 

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================
   
  TYPE(t_region), POINTER :: pRegion 
   
! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ReadConvGridWrapper.F90,v $ $Revision: 1.2 $'

  global => pRegion%global
  
  CALL RegisterFunction(global,'RFLU_ReadConvGridWrapper',__FILE__)

! ******************************************************************************
! Read grid file and convert if necessary 
! ******************************************************************************

  IF ( global%gridSource == GRID_SRC_CENTAUR_ASCII ) THEN    
    CALL RFLU_ReadGridCENTAURASCII(pRegion)
    CALL RFLU_ConvCENTAUR2ROCFLU(pRegion)
!  ELSE IF ( global%gridSource == GRID_SRC_CENTAUR_BINARY ) THEN    
! BBR - begin - choice of endianness
   ELSE IF ( global%gridSource == GRID_SRC_CENTAUR_BINARY .OR. & 
             global%gridSource == GRID_SRC_CENTAUR_BINARY_L .OR. &
             global%gridSource == GRID_SRC_CENTAUR_BINARY_B ) THEN
! BBR - end 
    CALL RFLU_ReadGridCENTAURBinary(pRegion)
    CALL RFLU_ConvCENTAUR2ROCFLU(pRegion)
  ELSE IF ( global%gridSource == GRID_SRC_VGRIDNS ) THEN 
    CALL RFLU_ReadGridVGRIDns(pRegion)
    CALL RFLU_ConvVGRIDns2ROCFLU(pRegion)
  ELSE IF ( global%gridSource == GRID_SRC_MESH3D ) THEN
    CALL RFLU_ReadGridMESH3D(pRegion)
    CALL RFLU_ConvMESH3D2ROCFLU(pRegion) 
  ELSE IF ( global%gridSource == GRID_SRC_TETMESH ) THEN 
    CALL RFLU_ReadGridTETMESH(pRegion)
    CALL RFLU_ConvTETMESH2ROCFLU(pRegion)
  ELSE IF ( global%gridSource == GRID_SRC_COBALT ) THEN 
    CALL RFLU_ReadGridCOBALT(pRegion)
    CALL RFLU_ConvCOBALT2ROCFLU(pRegion)
  ELSE IF ( global%gridSource == GRID_SRC_GAMBIT_NEUTRAL ) THEN 
    CALL RFLU_ReadGridGAMBITNeutral(pRegion)
    CALL RFLU_ConvGAMBIT2ROCFLU(pRegion)    
  ELSE IF ( global%gridSource == GRID_SRC_STARCD ) THEN
    CALL RFLU_ReadGridSTARCD(pRegion)
    CALL RFLU_ConvSTARCD2ROCFLU(pRegion) 
  ELSE IF ( global%gridSource == GRID_SRC_OMG ) THEN 
    CALL RFLU_ReadGridOMGASCII(pRegion)
    CALL RFLU_ConvOMG2ROCFLU(pRegion) 
  ELSE 
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
  END IF ! global%gridSource

! ******************************************************************************
! End
! ******************************************************************************
 
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ReadConvGridWrapper


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ReadConvGridWrapper.F90,v $
! Revision 1.2  2015/07/23 23:11:19  brollin
! 1) The pressure coefficient of the  collision model has been changed back to its original form
! 2) New options in the format of the grid and solutions have been added. Now the user can choose the endianness, and convert from one to the over in rfluconv.
! 3) The solutions are now stored in folders named by timestamp or iteration number.
! 4) The address enty in the hashtable has been changed to an integer(8) for cases when the grid becomes very large.
! 5) RFLU_WritePM can now compute PM2 on the fly for the Macroscale problem
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.5  2010/05/24 16:08:41  haselbac
! Added calls for OMG format
!
! Revision 1.4  2008/12/06 08:43:56  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:10  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2008/02/09 23:03:19  haselbac
! Added calls to StarCD routines
!
! Revision 1.1  2007/04/09 18:56:52  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2005/04/15 15:09:18  haselbac
! Initial revision
!
! Revision 1.2  2004/11/03 15:06:30  haselbac
! Added GAMBIT grid conversion option
!
! Revision 1.1  2004/10/19 19:30:36  haselbac
! Initial revision
!
! ******************************************************************************

