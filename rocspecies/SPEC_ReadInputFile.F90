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
! Purpose: Read user input, store it in the data structure and check.
!
! Description: None.
!
! Input:
!   regions                Data associated with regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: SPEC_ReadInputFile.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_ReadInputFile(regions)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters

  USE ModBuildFileNames, ONLY: BuildFileNamePlain

  USE SPEC_ModInterfaces, ONLY: SPEC_InitInputValuesSpecType, &
                                SPEC_ReadSpecSection, &
                                SPEC_ReadSpecTypeSection
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName,RCSIdentString
  CHARACTER(256) :: line
  INTEGER :: errorFlag,iPass,iReg,iSpec,iSpec2,iSpecType,nSpecies
  LOGICAL :: usedSomewhere,unusedSomewhere
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_ReadInputFile.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'SPEC_ReadInputFile',__FILE__)

! ******************************************************************************
! Read user input for species
! ******************************************************************************

! ==============================================================================
! Build file name
! ==============================================================================

  CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.inp',iFileName)

! ==============================================================================
! Look for keywords and read appropriate sections, read in two passes
! ==============================================================================

  Outer: DO iPass = 1,2
    iSpecType = 0

! ------------------------------------------------------------------------------
!   Open file
! ------------------------------------------------------------------------------

    OPEN(IF_INPUT,FILE=iFileName,FORM='FORMATTED',STATUS='UNKNOWN', &
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= 0 ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error

! ------------------------------------------------------------------------------
!   In each pass, look for keywords
! ------------------------------------------------------------------------------

    Inner: DO
      READ(IF_INPUT,'(A256)',IOSTAT=errorFlag) line

      IF ( errorFlag > 0 ) THEN ! Error occurred
        CALL ErrorStop(global,ERR_FILE_READ,__LINE__,'File: '//TRIM(iFileName))
      ELSE IF ( errorFlag < 0 ) THEN ! Encountered end of file
        EXIT Inner
      END IF ! errorFlag

! --- First pass: Read only SPECIES section ------------------------------------

      IF ( iPass == 1 ) THEN
        SELECT CASE( TRIM(line) )
          CASE ( '# SPECIES' )
            CALL SPEC_ReadSpecSection(regions)
        END SELECT

! --- Second pass: Read only SPECIES_TYPE sections -----------------------------

      ELSE IF ( iPass == 2 ) THEN
        SELECT CASE( TRIM(line) )
          CASE ( '# SPECIES_TYPE' )
            iSpecType = iSpecType + 1
            CALL SPEC_ReadSpecTypeSection(regions,iSpecType)
        END SELECT

! --- Invalid pass number (defensive coding) -----------------------------------

      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! iPass
    END DO Inner

! ------------------------------------------------------------------------------
!   End of pass 1: Allocate memory for species-specific input and initialize
! ------------------------------------------------------------------------------

    IF ( iPass == 1 ) THEN
      DO iReg = LBOUND(regions,1),UBOUND(regions,1)
        nSpecies = regions(iReg)%specInput%nSpecies

        IF ( nSpecies > 0 ) THEN
          ALLOCATE(regions(iReg)%specInput%specType(nSpecies),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, &
                           'regions%specInput%specType')
          END IF ! global%error

          DO iSpec = 1,nSpecies
            CALL SPEC_InitInputValuesSpecType(regions(iReg),iSpec)
          END DO ! iSpec
        ELSE
          NULLIFY(regions(iReg)%specInput%specType)
        END IF ! nSpecies
      END DO ! iReg
    END IF ! iPass

! ------------------------------------------------------------------------------
!   Check that have enough SPECIES_TYPE sections if species USED
! ------------------------------------------------------------------------------

    IF ( iPass == 2 ) THEN

      usedSomewhere   = .FALSE.
      unusedSomewhere = .FALSE.

      DO iReg = LBOUND(regions,1),UBOUND(regions,1)
        nSpecies = regions(iReg)%specInput%nSpecies
        IF ( regions(iReg)%specInput%usedFlag ) THEN
! ------- if USED, then require number of sections to be consistent
          IF ( iSpecType /= nSpecies ) THEN
            CALL ErrorStop(global,ERR_SPEC_NTYPES,__LINE__)
          END IF ! iSpecType
          IF ( nSpecies == 0 ) regions(iReg)%specInput%usedFlag = .FALSE.
        ELSE
! ------- if not USED, then no consistency check, but check nSpecies set to 0
          IF ( nSpecies /= 0 )  THEN
            CALL ErrorStop(global,ERR_SPEC_NTYPES,__LINE__)
          END IF ! nSpecies
        END IF ! usedFlag

        usedSomewhere   = usedSomewhere   .OR. &
                                  regions(iReg)%specInput%usedFlag
        unusedSomewhere = unusedSomewhere .OR. &
                            .NOT. regions(iReg)%specInput%usedFlag
      END DO ! iReg

      IF (usedSomewhere.AND.unusedSomewhere) THEN
        CALL ErrorStop( global,ERR_MP_ALLORNONE,__LINE__ )
      END IF ! usedSomewhere.AND.unusedSomewhere

      global%specUsed = usedSomewhere

    END IF ! iPass

! ------------------------------------------------------------------------------
!   Check that a value of MATERIAL keyword in SPECIES_TYPE section is not used
!   by more than one species
! ------------------------------------------------------------------------------

    IF ( iPass == 2 ) THEN
      DO iReg = LBOUND(regions,1),UBOUND(regions,1)
        nSpecies = regions(iReg)%specInput%nSpecies

        IF ( nSpecies > 1 ) THEN
          DO iSpec = 1,nSpecies-1
            DO iSpec2 = iSpec+1,nSpecies
              IF (TRIM(regions(iReg)%specInput%specType(iSpec)%pMaterial%name) &
                  == &
                 TRIM(regions(iReg)%specInput%specType(iSpec2)%pMaterial%name) &
                 ) THEN
                CALL ErrorStop( global,ERR_SPEC_NUNIQUEMAT,__LINE__ )
              END IF ! 
            END DO ! iSpec2
          END DO ! iSpec
        END IF ! nSpecies
      END DO ! iReg
    END IF ! iPass

! ------------------------------------------------------------------------------
!   Close file
! ------------------------------------------------------------------------------

    CLOSE(IF_INPUT,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= 0 ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error
  END DO Outer

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_ReadInputFile

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_ReadInputFile.F90,v $
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
! Revision 1.4  2005/11/10 02:38:30  haselbac
! Clean-up
!
! Revision 1.3  2004/07/28 15:29:21  jferry
! created global variable for spec use
!
! Revision 1.2  2004/06/16 20:01:29  haselbac
! Added use of ModBuildFileNames, cosmetics
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
! ******************************************************************************

