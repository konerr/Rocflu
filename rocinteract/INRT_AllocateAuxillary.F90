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
! Purpose: Allocates data structures for an individual interaction
!
! Description: none.
!
! Input: global    = global data
!        inrt      = interaction
!        nEdges    = number of Edges in interaction
!        nSwitches = number of integer parameters for interaction
!        nData     = number of real parameters for interaction
!
! Output: sets value of inrt%nEdges,
!         allocates and initializes inrt%edges, inrt%switches, and inrt%data
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_AllocateAuxillary.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_AllocateAuxillary( global,inrt,nEdges,nSwitches,nData )

  USE ModDataTypes
  USE ModGlobal,   ONLY : t_global
  USE ModInteract, ONLY : t_inrt_interact
  USE ModError
  USE INRT_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_global),        POINTER    :: global
  TYPE(t_inrt_interact), POINTER    :: inrt
  INTEGER,               INTENT(IN) :: nEdges,nSwitches,nData

! ... loop variables
  INTEGER :: iEdge

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER           :: errorFlag

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: INRT_AllocateAuxillary.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction( global,'INRT_AllocateAuxillary',__FILE__ )

! begin -----------------------------------------------------------------------

  inrt%nSwitches = nSwitches
  inrt%nData     = nData
  inrt%nEdges    = nEdges

  IF (nSwitches > 0) THEN

    ALLOCATE( inrt%switches(nSwitches),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    inrt%switches = -1

  END IF ! nSwitches

  IF (nData > 0) THEN

    ALLOCATE( inrt%data(nData),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    inrt%data = -1._RFREAL

  END IF ! nData

  IF (nEdges > 0) THEN

    ALLOCATE( inrt%edges(nEdges),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    DO iEdge = 1,nEdges
      inrt%edges(iEdge)%tEdge = INRT_EDGE_BAD
      inrt%edges(iEdge)%iNode = -1
      inrt%edges(iEdge)%token = INRT_PERM_PALL
    END DO ! iEdge

  END IF ! nEdges

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_AllocateAuxillary

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_AllocateAuxillary.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:49  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:01  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:50:11  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:14  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 21:56:08  fnajjar
! Initial revision after changing case
!
! Revision 1.4  2003/04/03 21:10:17  jferry
! implemented additional safety checks for rocinteract
!
! Revision 1.3  2003/04/02 22:32:03  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.2  2003/03/11 16:09:39  jferry
! Added comments
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!******************************************************************************

