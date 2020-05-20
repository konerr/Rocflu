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
! Purpose: Compute information for proper use of DCUHRE
!
! Description: None.
!
! Input: 
!   NDIM        Number of dimensions
!   NF          Number of functions
!   KEY         Type of integration rule
!   MAXCLS      Maximum number of calls allowed (NOTE also output argument)
!
! Output:
!   MAXCLS      Maximum number of calls allowed (NOTE also input argument)
!   NW          Size of work array
!
! Notes:
!
!******************************************************************************
!
! $Id: RFLU_ComputeDCUHREInfo.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

  SUBROUTINE RFLU_ComputeDCUHREInfo(global,NDIM,NF,KEY,MAXCLS,NW)

    USE ModGlobal, ONLY: t_global
    USE ModError

    IMPLICIT NONE

! - parameters      

    INTEGER, INTENT(IN) :: KEY,NDIM,NF
    INTEGER, INTENT(OUT) :: NW
    INTEGER, INTENT(INOUT) :: MAXCLS    
    TYPE(t_global), POINTER :: global

! - locals

    INTEGER :: MAXSUB,NUM

!******************************************************************************

    CALL RegisterFunction(global,'RFLU_ComputeDCUHREInfo',__FILE__)

! - Compute information ------------------------------------------------------

    SELECT CASE (KEY) ! Select according to integration rule
      CASE (0)  
        NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) &
            + 4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
      CASE (4)
        NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! KEY

    MAXCLS = MAX(MAXCLS,4*NUM) ! documentation recommends MAXCLS >= 3*NUM      

    MAXSUB = (MAXCLS - NUM)/(2*NUM) + 1
    NW     = MAXSUB*(2*NDIM+2*NF+2) + 17*NF + 1

    CALL DeregisterFunction(global)

!******************************************************************************

  END SUBROUTINE RFLU_ComputeDCUHREInfo
  
!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeDCUHREInfo.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:35  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:52  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2002/09/09 14:15:01  haselbac
! global now under regions
!
! Revision 1.1  2002/07/25 14:34:59  haselbac
! Initial revision
!
!******************************************************************************
  

