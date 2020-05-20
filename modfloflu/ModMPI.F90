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
! Purpose: encapsulate MPI header file, define master process.
!
! Description: none
!
! Notes: used in order to avoid the INCLUDE statement in the source code
!        (not standard in Fortran 90)
!
! ******************************************************************************
!
! $Id: ModMPI.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001-2005 by the University of Illinois
!
! ******************************************************************************

MODULE ModMPI

  INTEGER, PARAMETER :: MASTERPROC   = 0, &   ! master process
                        MPI_PATCHOFF = 100    ! offset for patch numbers (tag)

  INCLUDE 'mpif.h'

  INTEGER, PARAMETER :: MPI_RFREAL = MPI_DOUBLE_PRECISION  ! goes with RFREAL

END MODULE ModMPI

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModMPI.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:37  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:10  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:17  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.10  2005/04/15 15:34:09  haselbac
! Bug fix
!
! Revision 1.9  2005/04/15 15:30:28  haselbac
! Changed so that mpif.h always included for RFLU
!
! Revision 1.8  2004/10/19 19:28:58  haselbac
! Cosmetics only
!
! Revision 1.7  2004/09/27 22:47:00  wasistho
! added RADI_TAG_SHIFT
!
! Revision 1.6  2004/03/06 02:33:04  wasistho
! moved mpi tag shifts from ModParameters to ModMPI
!
! Revision 1.5  2002/03/21 18:07:15  jblazek
! Added check of MPI_PATCHOFF (for tags).
!
! Revision 1.4  2002/03/18 23:07:19  jblazek
! Finished multiblock and MPI.
!
! Revision 1.3  2002/01/31 20:23:59  jblazek
! Added treatment of edge & corner cells.
!
! Revision 1.2  2001/12/11 21:59:29  jblazek
! memory allocation added.
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
! ******************************************************************************

