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
! Purpose: set explicit interfaces to subroutines and functions.
!
! Description: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_ModInterfaces.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

MODULE INRT_ModInterfaces

  IMPLICIT NONE

  INTERFACE

! These are PUBLIC (i.e., in ModInterfacesInteract)
!
! INRT_BurnStatusUpdate, INRT_PrintMaterialInput,
! INRT_PrintUserInput, INRT_ReadMaterialInput,
! INRT_SetMaterial, INRT_SetParticleTemp, INRT_SourceTerms,
! INRT_TwoDimAverage, INRT_UserInput, INRT_VaporEnergyConversion

! These are PRIVATE (i.e., not in ModInterfacesInteract)

  SUBROUTINE INRT_AllocateAuxillary( global,inrt,nEdges,nSwitches,nData )
    USE ModGlobal,   ONLY : t_global
    USE ModInteract, ONLY : t_inrt_interact
    TYPE(t_global),        POINTER    :: global
    TYPE(t_inrt_interact), POINTER    :: inrt
    INTEGER,               INTENT(IN) :: nEdges,nSwitches,nData
  END SUBROUTINE INRT_AllocateAuxillary

  SUBROUTINE INRT_AugmentConSources( region,iInrt )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
    INTEGER,        INTENT(IN)    :: iInrt
  END SUBROUTINE INRT_AugmentConSources

  SUBROUTINE INRT_AugmentDisSources( region,iInrt )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT), TARGET :: region
    INTEGER,        INTENT(IN)            :: iInrt
  END SUBROUTINE INRT_AugmentDisSources

  SUBROUTINE INRT_CalcBurning( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT), TARGET :: region
  END SUBROUTINE INRT_CalcBurning

  SUBROUTINE INRT_CalcDrag( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT), TARGET :: region
  END SUBROUTINE INRT_CalcDrag

  SUBROUTINE INRT_CalcDragUnsteady( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT), TARGET :: region
  END SUBROUTINE INRT_CalcDragUnsteady
 
  ! Subbu - Added added mass explicit & implicit
  SUBROUTINE INRT_CalcDragUnsteady_AMExplicit( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT), TARGET :: region
  END SUBROUTINE INRT_CalcDragUnsteady_AMExplicit

  SUBROUTINE INRT_CalcDragUnsteady_AMImplicit( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT), TARGET :: region
  END SUBROUTINE INRT_CalcDragUnsteady_AMImplicit
  ! Subbu - End added added mass explicit & implicit

  SUBROUTINE INRT_CalcHeatTransferNonBurn( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT), TARGET :: region
  END SUBROUTINE INRT_CalcHeatTransferNonBurn

  SUBROUTINE INRT_CalcScouring( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT), TARGET :: region
  END SUBROUTINE INRT_CalcScouring

  SUBROUTINE INRT_CheckUserInput( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE INRT_CheckUserInput

  SUBROUTINE INRT_ComputeMaxEdges( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE INRT_ComputeMaxEdges

  SUBROUTINE INRT_DefineBoilingRegulation( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE INRT_DefineBoilingRegulation

  SUBROUTINE INRT_DefineBurning( region,matIndIn,matIndOut,matIndOx, &
                                 oxUsed,plagOutExists )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
    INTEGER,        INTENT(IN)    :: matIndIn,matIndOut,matIndOx
    LOGICAL,        INTENT(INOUT) :: oxUsed
    LOGICAL,        INTENT(OUT)   :: plagOutExists
  END SUBROUTINE INRT_DefineBurning

  SUBROUTINE INRT_DefineDrag( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE INRT_DefineDrag

  SUBROUTINE INRT_DefineHeatTransferNonBurn( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE INRT_DefineHeatTransferNonBurn

  SUBROUTINE INRT_DefineScouring( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE INRT_DefineScouring

  SUBROUTINE INRT_DetermineTokens( region,inrt )
    USE ModDataStruct, ONLY : t_region
    USE ModInteract,   ONLY : t_inrt_interact
    TYPE(t_region),        INTENT(INOUT) :: region
    TYPE(t_inrt_interact), POINTER       :: inrt
  END SUBROUTINE INRT_DetermineTokens

  SUBROUTINE INRT_FinishMomentumEdge(global,inrt,iXEdge,iEnd)
    USE ModGlobal,   ONLY : t_global
    USE ModInteract, ONLY : t_inrt_interact
    TYPE(t_global),        POINTER    :: global
    TYPE(t_inrt_interact), POINTER    :: inrt
    INTEGER,               INTENT(IN) :: iXEdge,iEnd
  END SUBROUTINE INRT_FinishMomentumEdge

  SUBROUTINE INRT_Initialize( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE INRT_Initialize

  SUBROUTINE INRT_ReadBoilingRegulation( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE INRT_ReadBoilingRegulation

  SUBROUTINE INRT_ReadBurning( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE INRT_ReadBurning

  SUBROUTINE INRT_ReadDefaultSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE INRT_ReadDefaultSection

  SUBROUTINE INRT_ReadDrag( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE INRT_ReadDrag

  SUBROUTINE INRT_ReadHeatTransferNonBurn( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE INRT_ReadHeatTransferNonBurn

  SUBROUTINE INRT_ReadInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE INRT_ReadInputFile

  SUBROUTINE INRT_ReadScouring( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE INRT_ReadScouring

  SUBROUTINE INRT_SetActiveness( global,val,actv )
    USE ModGlobal, ONLY : t_global
    USE ModDataTypes
    TYPE(t_global), POINTER     :: global
    REAL(RFREAL),   INTENT(IN)  :: val
    INTEGER,        INTENT(OUT) :: actv
  END SUBROUTINE INRT_SetActiveness

  SUBROUTINE INRT_SetPermission( global,val,perm )
    USE ModGlobal, ONLY : t_global
    USE ModDataTypes
    TYPE(t_global), POINTER     :: global
    REAL(RFREAL),   INTENT(IN)  :: val
    INTEGER,        INTENT(OUT) :: perm
  END SUBROUTINE INRT_SetPermission

  END INTERFACE

END MODULE INRT_ModInterfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_ModInterfaces.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:49  mtcampbe
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
! Revision 1.13  2004/07/28 15:42:12  jferry
! deleted defunct constructs: useDetangle, useSmokeDrag, useSmokeHeatTransfer
!
! Revision 1.12  2004/07/27 21:27:13  jferry
! removed rocinteract allocation routines (moved to rocpart)
!
! Revision 1.11  2004/07/26 17:05:51  fnajjar
! moved allocation of inrtSources into Rocpart
!
! Revision 1.10  2004/04/15 16:04:21  jferry
! minor formatting (removed trailing spaces)
!
! Revision 1.9  2004/03/08 21:57:36  jferry
! better error checking for burning without smoke case
!
! Revision 1.8  2004/03/02 21:47:29  jferry
! Added After Update interactions
!
! Revision 1.7  2003/09/25 15:48:43  jferry
! implemented Boiling Regulation interaction
!
! Revision 1.6  2003/09/19 20:35:26  jferry
! Implemented oxidizer species for burning interaction
!
! Revision 1.5  2003/04/03 16:18:28  fnajjar
! Include routines for burning and scouring
!
! Revision 1.4  2003/04/02 22:32:03  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.3  2003/03/24 23:30:52  jferry
! overhauled rocinteract to allow interaction design to use user input
!
! Revision 1.2  2003/03/11 16:09:17  jferry
! Added new routines
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!******************************************************************************

