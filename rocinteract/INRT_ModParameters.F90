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
! Purpose: define parameters for interactions
!
! Description: none
!
! Notes: none
!
!******************************************************************************
!
! $Id: INRT_ModParameters.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

MODULE INRT_ModParameters

  IMPLICIT NONE

! Types of Interactions

  INTEGER, PARAMETER :: INRT_TYPE_DRAG      = 1, &
                        INRT_TYPE_HTRANSNB  = 2, &
                        INRT_TYPE_SCOURING  = 3, &
                        INRT_TYPE_BURNING   = 4, &
                        INRT_TYPE_BOILRGN   = 5, &
                        INRT_TYPE_TOTAL     = 5

! Switches and data

  INTEGER, PARAMETER :: INRT_SWI_DRAG_MODEL    = 1, &
                        INRT_SWI_DRAG_UNSTEADY = 2, &
                        !INRT_SWI_DRAG_TOTAL    = 2, &
                        ! Subbu - Add additional flags
                        INRT_SWI_DRAG_AMSCHEME = 3, &
                        INRT_SWI_DRAG_VU       = 4, &
                        INRT_SWI_DRAG_TOTAL    = 4, &
                        ! Subbu - End add additional flags
                        INRT_DAT_DRAG_TOTAL    = 0

  INTEGER, PARAMETER :: INRT_SWI_HTRANSNB_MODEL = 1, &
                        INRT_SWI_HTRANSNB_TOTAL = 1, &
                        INRT_DAT_HTRANSNB_TOTAL = 0

  INTEGER, PARAMETER :: INRT_SWI_SCOURING_TOTAL  = 0, &
                        INRT_DAT_SCOURING_COEF0  = 0, &
                        INRT_DAT_SCOURING_TOTAL0 = 0     ! TOTAL = this + nPeul

  INTEGER, PARAMETER :: INRT_SWI_BURNING_MODEL      = 1, &
                        INRT_SWI_BURNING_OXUSED     = 2, &
                        INRT_SWI_BURNING_VAPOR_METH = 3, &
                        INRT_SWI_BURNING_TOTAL      = 3, &
                        INRT_DAT_BURNING_VAPOR_TEMP = 1, &
                        INRT_DAT_BURNING_HEAT_COEF  = 2, &
                        INRT_DAT_BURNING_MFRC_PLAG  = 3, &
                        INRT_DAT_BURNING_MFRC_PEUL0 = 3, &
                        INRT_DAT_BURNING_TOTAL0     = 3  ! TOTAL = this + nPeul

  INTEGER, PARAMETER :: INRT_SWI_BOILRGN_MODEL  = 1, &
                        INRT_SWI_BOILRGN_LIQIND = 2, &
                        INRT_SWI_BOILRGN_GASIND = 3, &
                        INRT_SWI_BOILRGN_TOTAL  = 3, &
                        INRT_DAT_BOILRGN_BOILPT = 1, &
                        INRT_DAT_BOILRGN_ENPMS  = 2, &
                        INRT_DAT_BOILRGN_TOTAL  = 2

! Values of switches

  INTEGER, PARAMETER :: &
    INRT_DRAG_MODEL_STOKES  = 1, & ! Stokes drag
    INRT_DRAG_MODEL_SN      = 2, & ! Schiller-Naumann
    INRT_DRAG_MODEL_SMRFLD  = 3, & ! Sommerfeld
    INRT_DRAG_MODEL_PARMAR  = 4, & ! Parmar-Haselbacher-Bala
    INRT_DRAG_MODEL_DEFAULT = INRT_DRAG_MODEL_SN

  INTEGER, PARAMETER :: &
    INRT_DRAG_UNSTEADY_NONE    = 0, & ! No unsteady drag
    INRT_DRAG_UNSTEADY_USE     = 1, & ! Use all unsteady drag
    INRT_DRAG_UNSTEADY_DEFAULT = INRT_DRAG_UNSTEADY_NONE

  ! Subbu - Define values for  additional flags
  INTEGER, PARAMETER :: &
    INRT_DRAG_AMSCHEME_EXPLICIT = 1, & ! Explicit added-mass calc
    INRT_DRAG_AMSCHEME_IMPLICIT = 2, & ! Implicit added-mass cals
    INRT_DRAG_AMSCHEME_DEFAULT = INRT_DRAG_AMSCHEME_EXPLICIT

  INTEGER, PARAMETER :: &
    INRT_DRAG_VISCUNST_NONE    = 0, & ! No viscous unsteady drag
    INRT_DRAG_VISCUNST_USE     = 1, & ! Compute viscous unsteady drag
    INRT_DRAG_VISCUNST_DEFAULT = INRT_DRAG_VISCUNST_USE
  ! Subbu - End define values for  additional flags

  INTEGER, PARAMETER :: &
    INRT_HTRANSNB_MODEL_STOKES  = 1, & ! Stokes thermal drag
    INRT_HTRANSNB_MODEL_RM      = 2, & ! Ranz-Marshall
    INRT_HTRANSNB_MODEL_DEFAULT = INRT_HTRANSNB_MODEL_RM

  INTEGER, PARAMETER :: &
    INRT_BURNING_MODEL_BECKSTEAD = 1, & ! Beckstead model
    INRT_BURNING_MODEL_DEFAULT   = INRT_BURNING_MODEL_BECKSTEAD

  INTEGER, PARAMETER :: &
    INRT_BURNING_VAPOR_METH_NONE = 0, & ! Not used
    INRT_BURNING_VAPOR_METH_USED = 1    ! Used

  INTEGER, PARAMETER :: &
    INRT_BOILRGN_MODEL_SHARP   = 1, & ! Sharp cutoff
    INRT_BOILRGN_MODEL_DEFAULT = INRT_BOILRGN_MODEL_SHARP

! Miscellaneous

! Values of Burning Status for particles

  INTEGER, PARAMETER :: &
    INRT_BURNSTAT_OFF = 0, & ! Particle not burning
    INRT_BURNSTAT_ON  = 1    ! Particle burning

! Indices of Edges:
!
!   G = Gas Node
!   L = Lagrangian particle Node
!   S = Smoke Node
!   X = Internal Node

  INTEGER, PARAMETER :: INRT_DRAG_L_XMOM_G = 1, &
                        INRT_DRAG_L_YMOM_G = 2, &
                        INRT_DRAG_L_ZMOM_G = 3, &
                        INRT_DRAG_NEDGES   = 3

  INTEGER, PARAMETER :: INRT_HTRANSNB_L_ENER_G = 1, &
                        INRT_HTRANSNB_NEDGES   = 1

! no parameters of this type for Scouring:
! -> Edges are created dynamically by matching smoke and particle materials

  INTEGER, PARAMETER :: INRT_BURNING_G_MASS_X  = 1, &
                        INRT_BURNING_L_MASS_X  = 2, &
                        INRT_BURNING_S_MASS_X0 = 2, &
                        INRT_BURNING_X_ENER_G  = 3, &
                        INRT_BURNING_X_ENER_LV = 4, &
                        INRT_BURNING_X_MASS_G  = 5, &
                        INRT_BURNING_X_MASS_L  = 6, &
                        INRT_BURNING_X_MASS_S0 = 6, &
                        INRT_BURNING_NEDGES0   = 6

  INTEGER, PARAMETER :: INRT_BOILRGN_NEDGES = 0

!******************************************************************************

! Types of Edges

  INTEGER, PARAMETER :: INRT_EDGE_BAD      = 0, &  ! Bad Edge (uninitialized)
                        INRT_EDGE_MASS     = 1, &  ! Mass Edge
                        INRT_EDGE_MOME     = 2, &  ! Momentum Edge
                        INRT_EDGE_MOME_DUM = 3, &  ! Dummy momentum Edge
                        INRT_EDGE_ENER     = 5, &  ! Energy Edge
                        INRT_EDGE_MASS_GHO = 6     ! Ghost mass Edge

! Permission levels

  INTEGER, PARAMETER :: INRT_PERM_BAD   = -1, &  ! Bad value (uninitialized)
                        INRT_PERM_BLOCK =  0, &  ! Block
                        INRT_PERM_PMASS =  1, &  ! Permit Mass
                        INRT_PERM_PMOME =  2, &  ! Permit Mass and Momentum
                        INRT_PERM_PALL  =  3     ! Permit All

! Activeness (everything less than INRT_ACT_ACTIVE is passive)

  INTEGER, PARAMETER :: INRT_ACT_BAD    = 2, &  ! Bad value (uninitialized)
                        INRT_ACT_ACTIVE = 1     ! Node is active

END MODULE INRT_ModParameters

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_ModParameters.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2009/07/09 20:44:20  mparmar
! Added parameter Parmar drag law in parameter list
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
! Revision 1.10  2007/03/07 22:17:23  fnajjar
! Added Sommerfeld drag law in the parameter list
!
! Revision 1.9  2004/07/28 15:42:12  jferry
! deleted defunct constructs: useDetangle, useSmokeDrag, useSmokeHeatTransfer
!
! Revision 1.8  2004/04/15 16:04:21  jferry
! minor formatting (removed trailing spaces)
!
! Revision 1.7  2004/03/02 21:47:29  jferry
! Added After Update interactions
!
! Revision 1.6  2003/09/25 15:48:43  jferry
! implemented Boiling Regulation interaction
!
! Revision 1.5  2003/09/19 20:35:26  jferry
! Implemented oxidizer species for burning interaction
!
! Revision 1.4  2003/04/02 22:32:03  jferry
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

