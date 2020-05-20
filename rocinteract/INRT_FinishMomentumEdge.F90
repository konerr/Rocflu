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
! Purpose: Finishes definition for the special case of a momentum Edge.
!
! Description: none.
!
! Input:  iXEdge:  index of the momentum Edge (which is its x-component)
!         iEnd:    end of Edge to place an Insulate Token on (iEnd = 1 or 2)
!
! Output: inrt:    pointer to the updated interaction data structure
!
! Notes:
!
!   Defining an Edge requires that two things be specified:
!   the type of Edge (mass, momentum, or energy), and
!   the indices of the Nodes on both ends.
!
!   Momentum Edges require additional set-up, which this routine provides.
!   First, they require that two dummy momentum Edges be created (because
!   routines that store primary source terms can only store scalar quantities,
!   so there needs to be extra room to send the vector momentum).  Second,
!   they require that an Insulate Token be placed at one end as part of the
!   definition of the interaction.
!
!******************************************************************************
!
! $Id: INRT_FinishMomentumEdge.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_FinishMomentumEdge(global,inrt,iXEdge,iEnd)

  USE ModDataTypes
  USE ModGlobal,   ONLY : t_global
  USE ModInteract, ONLY : t_inrt_interact
  USE ModError
  USE INRT_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_global),        POINTER    :: global
  TYPE(t_inrt_interact), POINTER    :: inrt
  INTEGER,               INTENT(IN) :: iXEdge,iEnd

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: INRT_FinishMomentumEdge.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction( global,'INRT_FinishMomentumEdge',__FILE__ )

! begin -----------------------------------------------------------------------

! set permission Token (Permit Mass and Momemtum) on specified end

  inrt%edges(iXEdge)%token(iEnd) = INRT_PERM_PMOME

! fill in data for dummy Edges

  inrt%edges(iXEdge+1)%tEdge = INRT_EDGE_MOME_DUM
  inrt%edges(iXEdge+1)%iNode = inrt%edges(iXEdge)%iNode
  inrt%edges(iXEdge+1)%token = INRT_PERM_BLOCK

  inrt%edges(iXEdge+2)%tEdge = INRT_EDGE_MOME_DUM
  inrt%edges(iXEdge+2)%iNode = inrt%edges(iXEdge)%iNode
  inrt%edges(iXEdge+2)%token = INRT_PERM_BLOCK

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_FinishMomentumEdge

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_FinishMomentumEdge.F90,v $
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
! Revision 1.1  2004/12/01 21:56:27  fnajjar
! Initial revision after changing case
!
! Revision 1.3  2003/04/02 22:32:04  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.2  2003/03/11 16:09:39  jferry
! Added comments
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!******************************************************************************

