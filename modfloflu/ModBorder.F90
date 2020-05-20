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
! ******************************************************************************
!
! Purpose: Definition of derived data type for borders.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ModBorder.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE ModBorder

  USE ModDataTypes
  
  IMPLICIT NONE

! ******************************************************************************
! Type definition
! ******************************************************************************

  TYPE t_border_data
    INTEGER :: sendRequest,sendRequestCount,sendRequestInt,tag,tagCount,tagInt
    REAL(RFREAL), DIMENSION(:), POINTER :: recvBuff1d,sendBuff1d
    INTEGER, DIMENSION(:,:), POINTER :: recvBuffInt,sendBuffInt
    REAL(RFREAL), DIMENSION(:,:), POINTER :: recvBuff,sendBuff
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: recvBuff3d,sendBuff3d
  END TYPE t_border_data

  TYPE t_border
    INTEGER :: iBorder,iProc,iRegionGlobal,iRegionLocal
    INTEGER :: nCellsRecv,nCellsSend,nVertRecv,nVertSend,nVertShared
#ifdef PLAG
    INTEGER :: nPclsRecv,nPclsSend,nPclsSendMax
#endif
    INTEGER, DIMENSION(:), POINTER :: icgRecv,icgRecvGlobalIds,icgSend, &
                                      ivgRecv,ivgSend,ivgShared
#ifdef PLAG
    INTEGER, DIMENSION(:,:), POINTER :: iPclSend
#endif
    TYPE(t_border_data) :: mixt,spec,plag
  END TYPE t_border

END MODULE ModBorder

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModBorder.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:37  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:16:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/11/28 23:05:04  mparmar
! Added recvBuff1d,sendBuff1d,recvBuff3d,sendBuff3d,icgRecvGlobalIds
!
! Revision 1.1  2007/04/09 18:49:10  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:16  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.7  2005/12/14 21:19:16  fnajjar
! Added nPclsSendMax for dynamic iPclSend
!
! Revision 1.6  2005/12/13 23:06:16  fnajjar
! Added tags and sendRequests pertinent for PLAG
!
! Revision 1.5  2005/05/18 22:05:26  fnajjar
! Added integer vars, fixed bug in declaration
!
! Revision 1.4  2005/04/30 13:47:51  haselbac
! Added arrays for parallelization of particle module
!
! Revision 1.3  2005/04/15 15:06:26  haselbac
! Added data arrays and variables
!
! Revision 1.2  2005/01/14 21:14:20  haselbac
! Added iProc
!
! Revision 1.1  2004/12/04 03:43:41  haselbac
! Initial revision
!
! ******************************************************************************

