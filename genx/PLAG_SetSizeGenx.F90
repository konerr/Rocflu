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
! Purpose: Set size of arrays in Genx.
!
! Description: none.
!
! Input: regions    = data of all regions,
!        wins, winp = GenX window registrations.
!
! Output: to Roccom.
!
! Notes: Surface registration for Tiles works only for External coupled bc.
!        Need to activate for both Internal and External bc.
!
!******************************************************************************
!
! $Id: PLAG_SetSizeGenx.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_SetSizeGenx( region )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  USE ModPartLag, ONLY    : t_plag
  IMPLICIT NONE
  INCLUDE 'roccomf90.h'

! ... parameters

  TYPE(t_region) :: region

! ... loop variables

! ... local variables

  INTEGER :: iLev, pid
  TYPE(t_global),     POINTER :: pGlobal
  TYPE(t_plag)  ,     POINTER :: pPlag

!******************************************************************************

  pGlobal => region%global

  CALL RegisterFunction( pGlobal,'PLAG_SetSizeGenx',__FILE__ )

!------------------------------------------------------------------------------
! Loop over all regions 
!------------------------------------------------------------------------------

  iLev   = region%currLevel
  pPlag => region%levels(iLev)%plag
  pid = region%iRegionGlobal*REGOFF+1

!------------------------------------------------------------------------------
! COM_set_size procedure must be called in RFLO_sendBoundaryValues
! if nPcls has changed.
!------------------------------------------------------------------------------
      
  CALL COM_set_size( TRIM(pGlobal%winp)//'.nc',pid,pPlag%nPcls)

!------------------------------------------------------------------------------
! finalize 
!------------------------------------------------------------------------------

  CALL DeregisterFunction( pGlobal )

END SUBROUTINE PLAG_SetSizeGenx

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_SetSizeGenx.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:30  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:47:51  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:58:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2004/11/16 08:51:13  jiao
! Fixed pane ID.
!
! Revision 1.2  2004/07/02 22:47:43  jiao
! Added definitions of iLev and pid.
!
! Revision 1.1  2004/07/02 22:08:04  fnajjar
! Initial import for Roccom3
!
!******************************************************************************

