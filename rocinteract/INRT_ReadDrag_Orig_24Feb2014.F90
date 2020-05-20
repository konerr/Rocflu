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
! Purpose: Reads in information related to the interaction Drag
!
! Description: none.
!
! Input: regions = data of all regions
!
! Output: fills user data into region%inrtInput%inrts
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_ReadDrag_Orig_24Feb2014.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_ReadDrag( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInteract,   ONLY : t_inrt_input,t_inrt_interact
  USE ModError
  USE ModParameters
  USE INRT_ModParameters

  USE ModInterfaces,      ONLY : ReadSection
  USE INRT_ModInterfaces, ONLY : INRT_SetActiveness,INRT_SetPermission, &
                                 INRT_DetermineTokens,INRT_DefineDrag
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg,iPlag

! ... local variables
  INTEGER, PARAMETER :: NKEYS_MAX = 20

  CHARACTER(CHRLEN)  :: RCSIdentString
  CHARACTER(CHRLEN)  :: keys(NKEYS_MAX)

  INTEGER :: brbeg,brend,nEdges
  INTEGER :: nPlag
  INTEGER :: nImplKeys,nNodeKeys,nKeys
  INTEGER :: ind,indMixt,indPlag0
  INTEGER :: iKeyUsed,iKeyModel,iKeyUnsteady,iKeyNode0
  INTEGER :: iKeyMixtActv,iKeyPlagActv,iKeyMixtPerm,iKeyPlagPerm

  LOGICAL :: defined(NKEYS_MAX)

  REAL(RFREAL) :: vals(NKEYS_MAX)

  TYPE(t_inrt_input),    POINTER :: input
  TYPE(t_inrt_interact), POINTER :: inrt
  TYPE(t_global),        POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_ReadDrag_Orig_24Feb2014.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'INRT_ReadDrag',__FILE__ )

! begin -----------------------------------------------------------------------

! define implementation-dependent keys

  iKeyUsed     = 1
  iKeyModel    = 2
  iKeyUnsteady = 3
  nImplKeys    = 3

  keys(iKeyUsed)     = 'USED'
  keys(iKeyModel)    = 'MODEL'
  keys(iKeyUnsteady) = 'UNSTEADY'

! define Node keys

  iKeyNode0 = nImplKeys
  iKeyMixtActv = iKeyNode0 + 1
  iKeyPlagActv = iKeyNode0 + 2
  iKeyMixtPerm = iKeyNode0 + 3
  iKeyPlagPerm = iKeyNode0 + 4
  nNodeKeys = 4

  keys(iKeyMixtActv) = 'MIXT_ACTV'
  keys(iKeyPlagActv) = 'PLAG_ACTV'
  keys(iKeyMixtPerm) = 'MIXT_PERM'
  keys(iKeyPlagPerm) = 'PLAG_PERM'

  nKeys = iKeyNode0 + nNodeKeys

  IF (nKeys > NKEYS_MAX) CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

! Read interaction section from input file

  CALL ReadSection( global,IF_INPUT,nKeys,keys,vals,defined )
  brbeg = LBOUND(regions,1)
  brend = UBOUND(regions,1)

  DO iReg=brbeg,brend

    input => regions(iReg)%inrtInput
    inrt  => input%inrts(INRT_TYPE_DRAG)

! - Check that INRT_DEFAULT section has been read, and that interaction has not

    IF (.NOT. input%defaultRead) &
      CALL ErrorStop( global,ERR_INRT_DEFUNREAD,__LINE__ )

    IF (inrt%used) CALL ErrorStop( global,ERR_INRT_READ,__LINE__ )

! - Use local variables for some useful quantities

    nPlag = input%nPlag

    indMixt  = input%indMixt
    indPlag0 = input%indPlag0

! - Check if interaction is used

    inrt%used = .TRUE. ! used by default when its section appears

    IF (defined(iKeyUsed)) THEN
      IF (NINT(vals(iKeyUsed)) == 0) inrt%used = .FALSE.
    END IF ! defined(iKeyUsed)

    IF (nPlag < 1) inrt%used = .FALSE. ! cannot occur without particles

    IF (.NOT. inrt%used) CYCLE ! do not bother with unused interactions

! - Define interaction (using any relevant information from input deck)

    CALL INRT_DefineDrag(regions(iReg))

! - Check for switches

    inrt%switches(INRT_SWI_DRAG_MODEL) = INRT_DRAG_MODEL_DEFAULT

    IF (defined(iKeyModel)) THEN

      SELECT CASE (NINT(vals(iKeyModel)))

      CASE (1)
        inrt%switches(INRT_SWI_DRAG_MODEL) = INRT_DRAG_MODEL_STOKES

      CASE (2)
        inrt%switches(INRT_SWI_DRAG_MODEL) = INRT_DRAG_MODEL_SN

      CASE (3)
        inrt%switches(INRT_SWI_DRAG_MODEL) = INRT_DRAG_MODEL_SMRFLD

      CASE (4)
        inrt%switches(INRT_SWI_DRAG_MODEL) = INRT_DRAG_MODEL_PARMAR

      CASE DEFAULT
        CALL ErrorStop( global,ERR_INRT_BADSWITCH,__LINE__ )

      END SELECT ! vals(iKeyModel)

    END IF ! defined(iKeyModel)

    inrt%switches(INRT_SWI_DRAG_UNSTEADY) = INRT_DRAG_UNSTEADY_DEFAULT

    IF (defined(iKeyUnsteady)) THEN
      
      SELECT CASE (NINT(vals(iKeyUnsteady)))

      CASE (0)
        inrt%switches(INRT_SWI_DRAG_UNSTEADY) = INRT_DRAG_UNSTEADY_NONE

      CASE (1)
        inrt%switches(INRT_SWI_DRAG_UNSTEADY) = INRT_DRAG_UNSTEADY_USE

      CASE DEFAULT
        CALL ErrorStop( global,ERR_INRT_BADSWITCH,__LINE__ )

      END SELECT ! vals(iKeyUnsteady)

    END IF ! defined(iKeyUnsteady)

! - Check for Mixture controls

    IF (defined(iKeyMixtActv)) &
      CALL INRT_SetActiveness(global,vals(iKeyMixtActv), &
                              inrt%activeness(indMixt))

    IF (defined(iKeyMixtPerm)) &
      CALL INRT_SetPermission(global,vals(iKeyMixtPerm), &
                              inrt%permission(indMixt))

! - Check for Lagrangian particle controls

    DO iPlag=1,nPlag+1

      ind = indPlag0 + iPlag

      IF (defined(iKeyPlagActv)) &
        CALL INRT_SetActiveness(global,vals(iKeyPlagActv), &
                                inrt%activeness(ind))

      IF (defined(iKeyPlagPerm)) &
        CALL INRT_SetPermission(global,vals(iKeyPlagPerm), &
                                inrt%permission(ind))

    END DO ! iPlag

! - Determine permission Tokens

    CALL INRT_DetermineTokens(regions(iReg),inrt)

  END DO ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_ReadDrag

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_ReadDrag_Orig_24Feb2014.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2009/07/09 20:44:22  mparmar
! Added reading of Parmar drag law
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
! Revision 1.2  2007/03/07 22:18:05  fnajjar
! Included Sommerfeld drag law
!
! Revision 1.1  2004/12/01 21:56:35  fnajjar
! Initial revision after changing case
!
! Revision 1.7  2004/07/23 22:43:16  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.6  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.5  2004/03/02 21:49:23  jferry
! Added inrtUsed flag to mixture data structure
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

