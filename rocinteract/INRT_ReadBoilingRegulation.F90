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
! Purpose: Reads in information related to the interaction Boiling Regulation
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
! $Id: INRT_ReadBoilingRegulation.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_ReadBoilingRegulation( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModMaterials,  ONLY : t_material
  USE ModInteract,   ONLY : t_inrt_input,t_inrt_interact
  USE ModError
  USE ModParameters
  USE INRT_ModParameters

  USE ModInterfaces,        ONLY : ReadBothSection
  USE ModInterfacesInteract,ONLY : INRT_SetMaterial
  USE INRT_ModInterfaces,   ONLY : INRT_DefineBoilingRegulation
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  INTEGER, PARAMETER :: NSTRKEYS_MAX = 5
  INTEGER, PARAMETER :: NKEYS_MAX    = 20

  CHARACTER(CHRLEN)  :: RCSIdentString
  CHARACTER(CHRLEN)  :: strKeys(NSTRKEYS_MAX),keys(NKEYS_MAX)
  CHARACTER(CHRLEN)  :: strVals(NSTRKEYS_MAX)

  INTEGER :: brbeg,brend
  INTEGER :: nStrKeys,nImplKeys,nKeys
  INTEGER :: iStrKeyMaterialLiq,iStrKeyMaterialGas
  INTEGER :: iKey,iKeyUsed,iKeyModel,iKeyBoilPt,iKeyEnPMs

  LOGICAL :: defined(NKEYS_MAX),strDefined(NSTRKEYS_MAX)

  REAL(RFREAL) :: boilPt,enPMs
  REAL(RFREAL) :: vals(NKEYS_MAX)

  TYPE(t_material),      POINTER :: matLiq,matGas
  TYPE(t_inrt_input),    POINTER :: input
  TYPE(t_inrt_interact), POINTER :: inrt
  TYPE(t_global),        POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_ReadBoilingRegulation.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'INRT_ReadBoilingRegulation',__FILE__ )

! begin -----------------------------------------------------------------------

! define string keys

  iStrKeyMaterialLiq = 1
  iStrKeyMaterialGas = 2
  nStrKeys = 2

  strKeys(iStrKeyMaterialLiq) = 'MATERIAL_LIQ'
  strKeys(iStrKeyMaterialGas) = 'MATERIAL_GAS'

  IF (nStrKeys > NSTRKEYS_MAX) &
    CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

! define implementation-dependent keys

  iKeyUsed   = 1
  iKeyModel  = 2
  iKeyBoilPt = 3
  iKeyEnPMs  = 4
  nImplKeys  = 4

  keys(iKeyUsed)   = 'USED'
  keys(iKeyModel)  = 'MODEL'
  keys(iKeyBoilPt) = 'BOILING_POINT'
  keys(iKeyEnPMs)  = 'ENERGY_PER_MASS'

  nKeys = nImplKeys

  IF (nKeys > NKEYS_MAX) CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

! Read interaction section from input file

  CALL ReadBothSection( global,IF_INPUT,nKeys,nStrKeys,keys,strKeys, &
                        vals,strVals,defined,strDefined )
  brbeg = LBOUND(regions,1)
  brend = UBOUND(regions,1)

  DO iReg=brbeg,brend

    input => regions(iReg)%inrtInput
    inrt  => input%inrts(INRT_TYPE_BOILRGN)

! - Check that INRT_DEFAULT section has been read, and that interaction has not

    IF (.NOT. input%defaultRead) &
      CALL ErrorStop( global,ERR_INRT_DEFUNREAD,__LINE__ )

    IF (inrt%used) CALL ErrorStop( global,ERR_INRT_READ,__LINE__ )

! - Check if interaction is used

    inrt%used = .TRUE. ! used by default when its section appears

    IF (defined(iKeyUsed)) THEN
      IF (NINT(vals(iKeyUsed)) == 0) inrt%used = .FALSE.
    END IF ! defined(iKeyUsed)

    IF (.NOT. inrt%used) CYCLE ! do not bother with unused interactions

! - Define interaction (using any relevant information from input deck)

    CALL INRT_DefineBoilingRegulation(regions(iReg))

! - Set pointers to liq and gas materials

    IF (strDefined(iStrKeyMaterialLiq)) THEN
      CALL INRT_SetMaterial(global,matLiq,strVals(iStrKeyMaterialLiq))
    ELSE
      CALL ErrorStop( global,ERR_INRT_MISSINGMAT,__LINE__ )
    END IF ! strDefined(iStrKeyMaterialLiq)

    IF (strDefined(iStrKeyMaterialGas)) THEN
      CALL INRT_SetMaterial(global,matGas,strVals(iStrKeyMaterialGas))
    ELSE
      CALL ErrorStop( global,ERR_INRT_MISSINGMAT,__LINE__ )
    END IF ! strDefined(iStrKeyMaterialGas)

! - Check that input and output materials have the same properties

    IF (matLiq%molw /= matGas%molw .OR. matLiq%dens /= matGas%dens .OR. &
        matLiq%spht /= matGas%spht) &
      CALL ErrorStop( global,ERR_INRT_BOIL_SAME,__LINE__ )

! - Check for switches

! - Material indices

    inrt%switches(INRT_SWI_BOILRGN_LIQIND) = matLiq%index
    inrt%switches(INRT_SWI_BOILRGN_GASIND) = matGas%index

! - Which model is used

    inrt%switches(INRT_SWI_BOILRGN_MODEL) = INRT_BOILRGN_MODEL_DEFAULT

    IF (defined(iKeyModel)) THEN

      SELECT CASE (NINT(vals(iKeyModel)))

      CASE (1)
        inrt%switches(INRT_SWI_BOILRGN_MODEL) = INRT_BOILRGN_MODEL_SHARP

      CASE DEFAULT
        CALL ErrorStop( global,ERR_INRT_BADSWITCH,__LINE__ )

      END SELECT ! vals(iKeyModel)

    END IF ! defined(iKeyModel)

! - Check for data

! - Boiling point of substance

    IF (defined(iKeyBoilPt)) THEN

      boilPt = vals(iKeyBoilPt)

      IF (boilPt <= 0._RFREAL) CALL ErrorStop( global,ERR_INRT_BADVAL,__LINE__ )

      inrt%data(INRT_DAT_BOILRGN_BOILPT) = boilPt

    ELSE

      CALL ErrorStop( global,ERR_INRT_MISSINGVAL,__LINE__ )

    END IF ! defined(iKeyBoilPt)

! - Energy released during condensation per unit mass of substance

    IF (defined(iKeyEnPMs)) THEN

      enPMs = vals(iKeyEnPMs)

      IF (enPMs <= 0._RFREAL) CALL ErrorStop( global,ERR_INRT_BADVAL,__LINE__ )

      inrt%data(INRT_DAT_BOILRGN_ENPMS) = enPMs

    ELSE

      CALL ErrorStop( global,ERR_INRT_MISSINGVAL,__LINE__ )

    END IF ! defined(iKeyEnPMs)

  END DO ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_ReadBoilingRegulation

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_ReadBoilingRegulation.F90,v $
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
! Revision 1.1  2007/04/09 18:50:11  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:15  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 21:56:32  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2004/07/23 22:43:16  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.4  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.3  2004/03/02 21:48:09  jferry
! First phase of replacing Detangle interaction
!
! Revision 1.2  2003/09/26 21:46:54  fnajjar
! Modified ModInterfaces call to ModInterfacesInteract
!
! Revision 1.1  2003/09/25 15:48:43  jferry
! implemented Boiling Regulation interaction
!
!******************************************************************************

