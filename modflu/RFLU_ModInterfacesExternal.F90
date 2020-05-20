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
! Purpose: Set explicit interfaces to external subroutines and functions for 
!   Rocflu.
!
! Description: None
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModInterfacesExternal.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************
  
MODULE RFLU_ModInterfacesExternal

  IMPLICIT NONE

  INTERFACE

#ifdef GENX
  SUBROUTINE RFLU_CheckCouplingInput(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), DIMENSION(:) :: regions   
  END SUBROUTINE RFLU_CheckCouplingInput

  SUBROUTINE Fluid_compute_integrals(globalGenx,integ)
    USE ModGenx
    DOUBLE PRECISION, DIMENSION(MAN_INTEG_SIZE) :: integ
    TYPE(t_globalGenx), POINTER :: globalGenx
  END SUBROUTINE Fluid_compute_integrals
 
  SUBROUTINE Fluid_finalize(globalGenx)
    USE ModGenx, ONLY: t_globalGenx
    TYPE(t_globalGenx), POINTER :: globalGenx
  END SUBROUTINE Fluid_finalize

  SUBROUTINE Fluid_preHdfOutput( globalGenx )
    USE ModGenx, ONLY : t_globalGenx
    TYPE(t_globalGenx), POINTER :: globalGenx
  END SUBROUTINE Fluid_preHdfOutput

  SUBROUTINE Fluid_postHdfOutput( globalGenx )
    USE ModGenx, ONLY : t_globalGenx
    TYPE(t_globalGenx), POINTER :: globalGenx
  END SUBROUTINE Fluid_postHdfOutput
 
  SUBROUTINE RFLU_FlowSolverDummy(globalGenx,timeSystem,dTimeSystem, &
                                  genxHandleBc,genxHandleGm)
    USE ModDataTypes
    USE ModGenx, ONLY: t_globalGenx
    INTEGER, INTENT(in) :: genxHandleBc, genxHandleGm
    DOUBLE PRECISION, INTENT(IN) :: dTimeSystem,timeSystem
    TYPE(t_globalGenx), POINTER :: globalGenx
  END SUBROUTINE RFLU_FlowSolverDummy

  SUBROUTINE RFLU_GetBoundaryValues(region)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLU_GetBoundaryValues

  SUBROUTINE RFLU_GetDeformation(region)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLU_GetDeformation

  SUBROUTINE RFLU_PutBoundaryValues(region)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLU_PutBoundaryValues

  SUBROUTINE RFLU_PutBoundaryValuesAlpha(region)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLU_PutBoundaryValuesAlpha

  SUBROUTINE RFLU_UpdateInbuffGm(globalGenx,dAlpha)
    USE ModDataTypes
    USE ModGenx, ONLY: t_globalGenx
    DOUBLE PRECISION, INTENT(in) :: dAlpha
    TYPE(t_globalGenx), POINTER :: globalGenx
  END SUBROUTINE RFLU_UpdateInbuffGm
#endif

  END INTERFACE

END MODULE RFLU_ModInterfacesExternal

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModInterfacesExternal.F90,v $
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
! Revision 1.6  2006/01/07 10:18:16  wasistho
! added Fluid_pre/postHdfOutput
!
! Revision 1.5  2004/10/19 19:28:06  haselbac
! Removed interfaces of routines moved into GENX modules, cosmetics
!
! Revision 1.4  2003/10/03 20:45:08  haselbac
! Removed interface for RFLU_UpdateInbuffWrapper
!
! Revision 1.3  2003/05/01 14:10:18  haselbac
! Added RFLU_CheckCouplingInput
!
! Revision 1.2  2003/04/24 15:38:10  haselbac
! Deleted input argument in RFLU_PutBoundaryValues
!
! Revision 1.1  2003/04/10 14:37:10  haselbac
! Initial revision
!
! ******************************************************************************

