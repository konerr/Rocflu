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
! *****************************************************************************
!
! Purpose: Read in user input related to physical materials.
!
! Description: None.
!
! Input: 
!   global 	Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! *****************************************************************************
!
! $Id: INRT_ReadMaterialInput.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! *****************************************************************************

SUBROUTINE INRT_ReadMaterialInput(global)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModMaterials
  USE ModError
  USE ModParameters

  USE ModBuildFileNames, ONLY: BuildFileNamePlain
    
  USE ModInterfaces, ONLY: ReadBothSection

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER, PARAMETER :: NSTRKEYS_MAX = 5
  INTEGER, PARAMETER :: NKEYS_MAX = 20
  LOGICAL :: strDefined(NSTRKEYS_MAX),defined(NKEYS_MAX)
  CHARACTER(20) :: strKeys(NSTRKEYS_MAX),keys(NKEYS_MAX)
  CHARACTER(256) :: line  
  CHARACTER(CHRLEN) :: RCSIdentString,strVals(NSTRKEYS_MAX)
  CHARACTER(CHRLEN) :: fname
  INTEGER :: errorFlag,nMat,iKeyMolw,iKeyDens,iKeyPr,iKeyRefVisc,iKeySpht, &
             iKeySurfTens,ikeySuthCoef,ikeySuthTemp,iKeyTboil,iKeyTmelt, &
             iKeyCond,iKeyDetonVel,iMat,iPass,iStrKeyName,iStrKeyPhase, &
             nStrKeys,nKeys
  REAL(RFREAL) :: vals(NKEYS_MAX)
  TYPE(t_material), POINTER :: material

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: INRT_ReadMaterialInput.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction( global,'INRT_ReadMaterialInput',__FILE__ )

! ******************************************************************************
! Define keys for string-valued quantities
! ******************************************************************************

  iStrKeyName  = 1
  iStrKeyPhase = 2
  nStrKeys     = 2

  IF ( nStrKeys > NSTRKEYS_MAX ) THEN
    CALL ErrorStop(global,ERR_EXCEEDS_DECL_MEM,__LINE__)
  END IF ! nStrKeys

  strKeys(iStrKeyName)  = 'NAME'
  strKeys(iStrKeyPhase) = 'PHASE'

! ******************************************************************************
! Define keys for real-valued quantities
! ******************************************************************************

  iKeyMolw     = 1
  iKeyDens     = 2
  iKeySpht     = 3
  iKeySurfTens = 4
  iKeyTboil    = 5
  iKeyTmelt    = 6
  iKeyRefVisc  = 7
  iKeySuthTemp = 8
  iKeySuthCoef = 9
  iKeyPr       = 10
  iKeyCond     = 11
  iKeyDetonVel = 12
  nKeys        = 12

  IF ( nKeys > NKEYS_MAX ) THEN 
    CALL ErrorStop(global,ERR_EXCEEDS_DECL_MEM,__LINE__)
  END IF ! nKeys

  keys(iKeyMolw)     = 'MOLW'
  keys(iKeyDens)     = 'DENS'
  keys(iKeySpht)     = 'SPHT'
  keys(iKeySurfTens) = 'SURFTENS'
  keys(iKeyTboil)    = 'TBOIL'
  keys(iKeyTmelt)    = 'TMELT'
  keys(iKeyRefVisc)  = 'VISCOSITY'
  keys(iKeySuthTemp) = 'SUTHTEMP'
  keys(iKeySuthCoef) = 'SUTHCOEF'
  keys(iKeyPr)       = 'PRANDTL'
  keys(iKeyCond)     = 'COND'
  keys(iKeyDetonVel) = 'DETONVEL'

! ******************************************************************************
! Search for MATERIAL sections
! ******************************************************************************

  CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.inp',fname)

  DO iPass = 1,2

! ==============================================================================
!   Open file
! ==============================================================================

    OPEN(IF_INPUT,FILE=TRIM(fname),FORM='FORMATTED',STATUS='OLD', &
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname))
    END IF ! global%error

! ==============================================================================
!   Read file looking for keywords
! ==============================================================================

    SELECT CASE ( iPass )

! -----------------------------------------------------------------------------
!     First pass: Count number of material sections and allocate materials
! -----------------------------------------------------------------------------

      CASE ( 1 ) 
        nMat = 0 ! initialize count of materials

        DO
          READ(IF_INPUT,'(A256)',ERR=10,END=86) line
          IF ( TRIM(line) == '# MATERIAL' ) THEN 
            nMat = nMat + 1
          END IF ! TRIM
        END DO ! <empty>

86      CONTINUE

        IF ( nMat > 0 ) THEN
          ALLOCATE(global%materials(nMat),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
          END IF ! global%error
        ELSE 
          NULLIFY(global%materials)
        END IF ! nMat

        global%nMaterials = nMat

! -----------------------------------------------------------------------------
!     Second pass: Search for real- and string-valued keys
! -----------------------------------------------------------------------------

      CASE ( 2 ) 
        iMat = 0 ! initialize index of materials

        DO
          READ(IF_INPUT,'(A256)',ERR=10,END=87) line
          
          IF ( TRIM(line) == '# MATERIAL' ) THEN
            iMat = iMat + 1
            material => global%materials(iMat)

! --------- Set values to indicate the status of undefined --------------------

            material%molw     = -1.0_RFREAL
            material%dens     = -1.0_RFREAL
            material%spht     = -1.0_RFREAL
            material%surftens = -1.0_RFREAL
            material%Tboil    = -1.0_RFREAL
            material%Tmelt    = -1.0_RFREAL
            material%Phase    = -1
            material%refVisc  = -1.0_RFREAL
            material%suthTemp = -1.0_RFREAL
            material%suthCoef = -1.0_RFREAL
            material%pr       = -1.0_RFREAL
            material%cond     = -1.0_RFREAL
            material%detonVel = REAL(CRAZY_VALUE_INT,KIND=RFREAL)

! --------- Read real- and string-valued keys ---------------------------------

            CALL ReadBothSection(global,IF_INPUT,nKeys,nStrKeys,keys, &
                                 strKeys,vals,strVals,defined,strDefined)

! --------- ensure that material is named -------------------------------------

            IF ( strDefined(iStrKeyName) .EQV. .TRUE. ) THEN
              material%name = TRIM(strVals(iStrKeyName))
            ELSE
              CALL ErrorStop(global,ERR_VAL_UNDEFINED,__LINE__)
            END IF ! strDefined(iStrKeyName)

! --------- Set material phase ------------------------------------------------

            IF ( strDefined(iStrKeyPhase) .EQV. .TRUE. ) THEN                                              
              SELECT CASE ( strVals(iStrKeyPhase)(1:1) )
                CASE ( 'G','g' )
                  material%phase = 1
                CASE ( 'L','l' )
                  material%phase = 2
                CASE ( 'S','s' )
                  material%phase = 3
                CASE DEFAULT
                  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
                END SELECT ! strVals
            END IF ! strDefined

! --------- Set material index to iMat (its index in global%materials(:)) -----

            material%index = iMat

! --------- Set real-valued quantities ----------------------------------------

            IF ( defined(iKeyMolw) .EQV. .TRUE. ) THEN 
              material%molw = vals(iKeyMolw)
            ELSE 
              CALL ErrorStop(global,ERR_VAL_UNDEFINED,__LINE__,'MOLW')
            END IF ! defined(iKeyMolw)

            IF ( defined(iKeyDens) .EQV. .TRUE. ) THEN 
              material%dens = vals(iKeyDens)
            ELSE 
              IF ( material%phase /= 1 ) THEN 
                CALL ErrorStop(global,ERR_VAL_UNDEFINED,__LINE__,'DENS')
              END IF ! material%phase
            END IF ! defined(iKeyDens)

            IF ( defined(iKeySpht) .EQV. .TRUE. ) THEN
              material%spht = vals(iKeySpht)
            ELSE
              CALL ErrorStop(global,ERR_VAL_UNDEFINED,__LINE__,'SPHT')
            END IF ! defined(iKeySpht)

            IF ( defined(iKeySurfTens) .EQV. .TRUE. ) THEN 
              material%surftens = vals(iKeySurfTens)
            END IF ! defined(iKeySurfTens)

            IF ( defined(iKeyTboil) .EQV. .TRUE. ) THEN
              material%Tboil = vals(iKeyTboil)
            END IF ! defined(iKeyTboil)

            IF ( defined(iKeyTmelt) .EQV. .TRUE. ) THEN
              material%Tmelt = vals(iKeyTmelt)
            END IF ! defined(iKeyTmelt)

            IF ( defined(iKeyRefVisc) .EQV. .TRUE. ) THEN
              material%refVisc = vals(iKeyRefVisc)
            END IF ! defined(iKeyRefVisc)

            IF ( defined(iKeySuthTemp) .EQV. .TRUE. ) THEN
              material%suthTemp = vals(iKeySuthTemp)
            END IF ! defined(iKeySuthTemp)

            IF ( defined(iKeySuthCoef) .EQV. .TRUE. ) THEN
              material%suthCoef = vals(iKeySuthCoef)
            END IF ! defined(iKeySuthCoef)

            IF ( defined(iKeyPr) .EQV. .TRUE. ) THEN
              material%pr = vals(iKeyPr)
            END IF ! defined(iKeyPr)

            IF ( defined(iKeyCond) .EQV. .TRUE. ) THEN
              material%cond = vals(iKeyCond)
            END IF ! defined(iKeyPr)

            IF ( defined(iKeyDetonVel) .EQV. .TRUE. ) THEN
              material%detonVel = vals(iKeyDetonVel)
            END IF ! defined(iKeyPr)
          END IF ! line
        END DO ! <empty>

87      CONTINUE

! -----------------------------------------------------------------------------
!     Default
! -----------------------------------------------------------------------------

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! iPass

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(IF_INPUT,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /=  ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname))
    END IF ! global%error
  END DO ! iPass

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)
  
  GOTO 999

10   CONTINUE
  CALL ErrorStop(global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname))

999  CONTINUE

END SUBROUTINE INRT_ReadMaterialInput

! *****************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_ReadMaterialInput.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2009/07/09 20:44:23  mparmar
! Added reading of viscosity and Prandtl number for material
!
! Revision 1.3  2008/12/06 08:43:50  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:02  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:50:12  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:15  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2005/11/14 17:01:44  haselbac
! Added check for DENS, clean-up
!
! Revision 1.2  2005/11/10 02:32:35  haselbac
! Clean-up
!
! Revision 1.1  2004/12/01 21:56:39  fnajjar
! Initial revision after changing case
!
! Revision 1.6  2004/07/27 21:28:32  jferry
! minor bug fix
!
! Revision 1.5  2004/07/23 22:43:17  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.4  2004/04/15 16:04:21  jferry
! minor formatting (removed trailing spaces)
!
! Revision 1.3  2004/03/02 21:49:46  jferry
! Added melting and boiling point to material definitions
!
! Revision 1.2  2003/09/13 20:17:31  fnajjar
! Added surface tension to Materials datastructure
!
! Revision 1.1  2003/03/24 23:23:25  jferry
! converted from libfloflu routine to rocinteract routine
!
! *****************************************************************************

