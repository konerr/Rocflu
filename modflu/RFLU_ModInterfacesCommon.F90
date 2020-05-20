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
! Purpose: Set explicit interfaces to common subroutines and functions for
!   Rocflu.
!
! Description: None
!
! Notes: 
!   1. The subroutines contained in this interface file are common to at least
!      two Rocflu modules in the sense that these modules contain subroutines
!      with these names, but the actual source is NOT common. 
!
! ******************************************************************************
!
! $Id: RFLU_ModInterfacesCommon.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************
  
MODULE RFLU_ModInterfacesCommon

  IMPLICIT NONE

  INTERFACE

  SUBROUTINE RFLU_AllocateMemory(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_AllocateMemory

  SUBROUTINE RFLU_AllocateMemoryWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_AllocateMemoryWrapper
  
  SUBROUTINE RFLU_DeallocateMemory(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_DeallocateMemory

  SUBROUTINE RFLU_DeallocateMemoryWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_DeallocateMemoryWrapper
  
  SUBROUTINE RFLU_PrintHeader(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLU_PrintHeader

  END INTERFACE

END MODULE RFLU_ModInterfacesCommon

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModInterfacesCommon.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:41  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:41  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2004/10/19 19:28:04  haselbac
! Adapted interface of RFLU_AllocateMemoryWrapper
!
! Revision 1.3  2004/03/19 21:19:21  haselbac
! Removed interfaces for alloc/dealloc routines
!
! Revision 1.2  2004/02/26 21:02:01  haselbac
! Added/deleted entries for memory allocation due to PLAG integration
!
! Revision 1.1  2003/04/10 14:37:10  haselbac
! Initial revision
!
! ******************************************************************************

