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
!******************************************************************************
!
! Purpose: Reads in default information for all interactions
!
! Description: none.
!
! Input: regions = data of all regions
!
! Output: fills user data into region%inrtInput
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_ReadDefaultSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_ReadDefaultSection( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInteract,   ONLY : t_inrt_input
  USE ModError
  USE ModParameters
  USE INRT_ModParameters

  USE ModInterfaces,      ONLY : ReadSection
  USE ModInterfaces,      ONLY : MakeNumberedKeys
  USE INRT_ModInterfaces, ONLY : INRT_SetActiveness,INRT_SetPermission
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg,iPlag,iPeul,iInrt

! ... local variables
  INTEGER, PARAMETER :: NKEYS_MAX = 30
  INTEGER, PARAMETER :: NPEUL_MAX = 10

  CHARACTER(CHRLEN)  :: RCSIdentString
  CHARACTER(CHRLEN)  :: keys(NKEYS_MAX)

  INTEGER :: brbeg,brend
  INTEGER :: nPlag,nPeul
  INTEGER :: nFixedKeys,nKeys
  INTEGER :: ind,indMixt,indPlag0,indPeul0
  INTEGER :: iKey,iKeyMixtActv,iKeyPlagActv,iKeyMixtPerm,iKeyPlagPerm
  INTEGER :: iKeyComputeAux,iKeyTwoDAverage
  INTEGER :: iKeyActv,iKeyPerm,iKeyPeulActv0,iKeyPeulPerm0

  LOGICAL :: defined(NKEYS_MAX)

  REAL(RFREAL) :: vals(NKEYS_MAX)

  TYPE(t_inrt_input), POINTER :: input
  TYPE(t_global),     POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: INRT_ReadDefaultSection.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'INRT_ReadDefaultSection',__FILE__ )

! begin -----------------------------------------------------------------------

! define Node keys and other fixed keys

  iKeyComputeAux  = 1
  iKeyTwoDAverage = 2
  iKeyMixtActv    = 3
  iKeyPlagActv    = 4
  iKeyMixtPerm    = 5
  iKeyPlagPerm    = 6
  nFixedKeys = 6

  keys(iKeyComputeAux)  = 'COMPUTE_AUX'
  keys(iKeyTwoDAverage) = '2D_AVERAGE'
  keys(iKeyMixtActv) = 'MIXT_ACTV'
  keys(iKeyPlagActv) = 'PLAG_ACTV'
  keys(iKeyMixtPerm) = 'MIXT_PERM'
  keys(iKeyPlagPerm) = 'PLAG_PERM'

  iKeyPeulActv0 = nFixedKeys
  iKeyPeulPerm0 = iKeyPeulActv0 + NPEUL_MAX

  CALL MakeNumberedKeys(keys,iKeyPeulActv0+1,'SPEC',1,NPEUL_MAX,1)
  CALL MakeNumberedKeys(keys,iKeyPeulPerm0+1,'SPEC',1,NPEUL_MAX,1)

  DO iPeul=1,NPEUL_MAX
    keys(iKeyPeulActv0+iPeul) = TRIM(keys(iKeyPeulActv0+iPeul))//'_ACTV'
    keys(iKeyPeulPerm0+iPeul) = TRIM(keys(iKeyPeulPerm0+iPeul))//'_PERM'
  END DO ! iPeul

  nKeys = iKeyPeulPerm0 + NPEUL_MAX

  IF (nKeys > NKEYS_MAX) CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

! Read default section from input file

  CALL ReadSection( global,IF_INPUT,nKeys,keys,vals,defined )
  brbeg = LBOUND(regions,1)
  brend = UBOUND(regions,1)

  DO iReg=brbeg,brend

    input => regions(iReg)%inrtInput

    IF (input%defaultRead) &
      CALL ErrorStop( global,ERR_INRT_DEFREAD,__LINE__ )

    input%defaultRead = .TRUE.

    nPlag = input%nPlag
    nPeul = input%nPeul
    IF (nPeul > NPEUL_MAX) &
      CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

    input%computeAux = .TRUE.

    IF (defined(iKeyComputeAux)) THEN

      SELECT CASE ( NINT(vals(iKeyComputeAux)) )

      CASE (0)
        input%computeAux = .FALSE.

      CASE (1)
        input%computeAux = .TRUE.

      CASE DEFAULT
        CALL ErrorStop( global,ERR_INRT_BADSWITCH,__LINE__ )

      END SELECT ! iKeyComputeAux

    END IF ! iKeyComputeAux

    input%twoDAverage = 0

    IF (defined(iKeyTwoDAverage)) THEN

      SELECT CASE ( NINT(vals(iKeyTwoDAverage)) )

      CASE (0:1)
        input%twoDAverage = NINT(vals(iKeyTwoDAverage))

      CASE DEFAULT
        CALL ErrorStop( global,ERR_INRT_BADSWITCH,__LINE__ )

      END SELECT ! iKeyTwoDAverage

    END IF ! iKeyTwoDAverage

! - Indices for Nodes

    indMixt  = input%indMixt
    indPlag0 = input%indPlag0
    indPeul0 = input%indPeul0

! - Check for Mixture controls

    IF (defined(iKeyMixtActv)) &
      CALL INRT_SetActiveness(global,vals(iKeyMixtActv), &
                              input%globActiveness(indMixt))

    IF (defined(iKeyMixtPerm)) &
      CALL INRT_SetPermission(global,vals(iKeyMixtPerm), &
                              input%globPermission(indMixt))

! - Check for Lagrangian particle controls

    DO iPlag = 1,nPlag+1

      ind = indPlag0 + iPlag

      IF (defined(iKeyPlagActv)) &
        CALL INRT_SetActiveness(global,vals(iKeyPlagActv), &
                                input%globActiveness(ind))

      IF (defined(iKeyPlagPerm)) &
        CALL INRT_SetPermission(global,vals(iKeyPlagPerm), &
                                input%globPermission(ind))

    END DO ! iPlag

! - Check for Smoke controls

    DO iPeul = 1,nPeul

      iKeyActv = iKeyPeulActv0 + iPeul
      iKeyPerm = iKeyPeulPerm0 + iPeul
      ind      = indPeul0      + iPeul

      IF (defined(iKeyActv)) &
        CALL INRT_SetActiveness(global,vals(iKeyActv), &
                                input%globActiveness(ind))

      IF (defined(iKeyPerm)) &
        CALL INRT_SetPermission(global,vals(iKeyPerm), &
                                input%globPermission(ind))

    END DO ! iPeul

! - Initialize local Activeness and Permission to global values

    DO iInrt = 1,INRT_TYPE_TOTAL

      input%inrts(iInrt)%activeness = input%globActiveness
      input%inrts(iInrt)%permission = input%globPermission

    END DO ! iInrt

  END DO ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_ReadDefaultSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_ReadDefaultSection.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
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
! Revision 1.1  2004/12/01 21:56:34  fnajjar
! Initial revision after changing case
!
! Revision 1.9  2004/07/28 15:42:12  jferry
! deleted defunct constructs: useDetangle, useSmokeDrag, useSmokeHeatTransfer
!
! Revision 1.8  2004/07/23 22:43:16  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.7  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.6  2004/03/02 21:47:29  jferry
! Added After Update interactions
!
! Revision 1.5  2003/09/25 15:48:43  jferry
! implemented Boiling Regulation interaction
!
! Revision 1.4  2003/04/02 22:32:04  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.3  2003/03/24 23:30:52  jferry
! overhauled rocinteract to allow interaction design to use user input
!
! Revision 1.2  2003/03/11 16:09:39  jferry
! Added comments
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!******************************************************************************

