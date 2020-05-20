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
! Purpose: Read in user input related to boundary conditions for species 
!   (done on all processors).
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
! $Id: SPEC_RFLU_ReadBcInputFile.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_RFLU_ReadBcInputFile(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModError
  USE ModParameters
  USE ModMPI
  
  USE ModBuildFileNames, ONLY: BuildFileNamePlain
  
  USE ModInterfacesSpecies, ONLY: SPEC_RFLU_ReadBcFarfSection, & 
                                  SPEC_RFLU_ReadBcInflowSection, & 
                                  SPEC_RFLU_ReadBcInjectSection, & 
                                  SPEC_RFLU_ReadBcSectionDummy
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName,RCSIdentString
  CHARACTER(256) :: line
  INTEGER :: errorFlag,iPatch,loopCounter
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_grid) :: grid
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_RFLU_ReadBcInputFile.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_ReadBcInputFile',__FILE__)

! ******************************************************************************
! Open file
! ******************************************************************************

  CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.bc',iFileName)
  
  OPEN(IF_INPUT,FILE=iFileName,FORM='FORMATTED',STATUS='OLD',IOSTAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! Read file looking for keywords
! ******************************************************************************

  loopCounter = 0

  KeyWordLoop: DO
    READ(IF_INPUT,'(A256)',IOSTAT=errorFlag) line
    
    IF ( errorFlag > 0 ) THEN ! Error occurred
      CALL ErrorStop(global,ERR_FILE_READ,__LINE__,'File: '//TRIM(iFileName))
    ELSE IF ( errorFlag < 0 ) THEN ! Encountered end of file
      EXIT KeyWordLoop
    END IF ! errorFlag    
    
    SELECT CASE( TRIM(line) )
      CASE ('# BC_SLIPW')
        CALL SPEC_RFLU_ReadBcSectionDummy(pRegion)
      CASE ('# BC_NOSLIP')
        CALL SPEC_RFLU_ReadBcSectionDummy(pRegion)      
! TEMPORARY - Keep this for backward compatibility
      CASE ('# BC_INFLOW')
        CALL SPEC_RFLU_ReadBcInflowSection(pRegion)
! END TEMPORARY    
      CASE ('# BC_INFLOW_TOTANG')
        CALL SPEC_RFLU_ReadBcInflowSection(pRegion)
      CASE ('# BC_INFLOW_VELTEMP')
        CALL SPEC_RFLU_ReadBcInflowSection(pRegion)             
      CASE ('# BC_OUTFLOW')
        CALL SPEC_RFLU_ReadBcSectionDummy(pRegion)      
      CASE ('# BC_FARF')
        CALL SPEC_RFLU_ReadBcFarfSection(pRegion)
      CASE ('# BC_INJECT')
        CALL SPEC_RFLU_ReadBcInjectSection(pRegion)
! Subbu - Correction - Added BC_VIRTUAL 
      CASE ('# BC_VIRTUAL')
        CALL SPEC_RFLU_ReadBcSectionDummy(pRegion)
! Subbu - End Correction
      CASE ('# END')
        EXIT
    END SELECT ! TRIM(line)

    loopCounter = loopCounter + 1 ! Prevent infinite loop
    IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
      CALL ErrorStop(global,ERR_INFINITE_LOOP ,__LINE__)
    END IF ! loopCounter
  END DO KeyWordLoop

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_INPUT,IOSTAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_RFLU_ReadBcInputFile

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_RFLU_ReadBcInputFile.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:53  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:05  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:51:23  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:50  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.5  2006/04/07 15:19:25  haselbac
! Removed tabs
!
! Revision 1.4  2005/04/27 02:13:37  haselbac
! Added routine to read INFLOW_VELTEMP
!
! Revision 1.3  2005/03/31 17:20:03  haselbac
! Fixed problem with crashes by reading with dummy routine
!
! Revision 1.2  2004/06/16 20:01:23  haselbac
! Added use of ModBuildFileNames, cosmetics
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
! ******************************************************************************

