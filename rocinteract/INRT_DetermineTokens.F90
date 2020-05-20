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
! Purpose: puts permission Tokens on Edges
!
! Description: sets values of Tokens based in several criteria:
!
!   (a) the Permission level of Nodes
!   (b) the relative Activeness of the Nodes at either end
!   (c) if it is an upwind Node of a mass Edge
!   (d) if the Edge contains an internal Node
!
! Input: region = region data
!        inrt = interaction
!
! Output: modifies inrt
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_DetermineTokens.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_DetermineTokens( region,inrt )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInteract
  USE ModError
  USE INRT_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_region),        INTENT(INOUT) :: region
  TYPE(t_inrt_interact), POINTER       :: inrt

! ... loop variables
  INTEGER :: iEdge

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  TYPE(t_inrt_input), POINTER :: input
  TYPE(t_inrt_edge),  POINTER :: edge
  TYPE(t_global),     POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_DetermineTokens.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_DetermineTokens',__FILE__ )

! begin -----------------------------------------------------------------------

  input => region%inrtInput

  DO iEdge=1,inrt%nEdges

    edge => inrt%edges(iEdge)

! - Permission Tokens already placed on dummy Edges

    IF (edge%tEdge == INRT_EDGE_MOME_DUM) CYCLE

! - For Ghost mass Edge, set downwind permission Token to 0 (Block), and
! - Upwind Token to either 1 (Permit Mass) or 0 (Block).
! - Note that activeness plays no role for Ghost mass Edges.

    IF (edge%tEdge == INRT_EDGE_MASS_GHO) THEN

      edge%token(1) = MIN(INRT_PERM_PMASS,inrt%permission(edge%iNode(1)))
      edge%token(2) = INRT_PERM_BLOCK
      CYCLE

    END IF ! INRT_EDGE_MASS_GHO

! - Decrease permission Tokens to Permission level of corresponding Nodes

    edge%token(1) = MIN(edge%token(1),inrt%permission(edge%iNode(1)))
    edge%token(2) = MIN(edge%token(2),inrt%permission(edge%iNode(2)))

! - If Nodes on either end of Edge differ in Activeness, decrease the
! - permission Token on the more active end to 0 (Block)

    IF (inrt%activeness(edge%iNode(1)) > inrt%activeness(edge%iNode(2))) &
      edge%token(1) = MIN(edge%token(1),INRT_PERM_BLOCK)

    IF (inrt%activeness(edge%iNode(2)) > inrt%activeness(edge%iNode(1))) &
      edge%token(2) = MIN(edge%token(2),INRT_PERM_BLOCK)

! - Decrease permission Token of the upwind end of a mass Edge to 1
! - (Permit Mass only)

    IF (edge%tEdge == INRT_EDGE_MASS) &
      edge%token(1) = MIN(edge%token(1),INRT_PERM_PMASS)

! - Decrease permission Token of the upwind end of an Edge if there is
! - an Internal Node there

    IF (edge%iNode(1) == input%indIntl) &
      edge%token(1) = MIN(edge%token(1),INRT_PERM_BLOCK)

! - Decrease permission Token to 0 (Block) if it is equivalent to Block
! - for its Edge type

    SELECT CASE (edge%tEdge)

    CASE (INRT_EDGE_MOME)
      IF (edge%token(1) < INRT_PERM_PMOME) edge%token(1) = INRT_PERM_BLOCK
      IF (edge%token(2) < INRT_PERM_PMOME) edge%token(2) = INRT_PERM_BLOCK

    CASE (INRT_EDGE_ENER)
      IF (edge%token(1) < INRT_PERM_PALL ) edge%token(1) = INRT_PERM_BLOCK
      IF (edge%token(2) < INRT_PERM_PALL ) edge%token(2) = INRT_PERM_BLOCK

    END SELECT ! edge%tEdge

  END DO ! iEdge

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_DetermineTokens

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_DetermineTokens.F90,v $
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
! Revision 1.1  2004/12/01 21:56:26  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2003/09/19 20:35:26  jferry
! Implemented oxidizer species for burning interaction
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

