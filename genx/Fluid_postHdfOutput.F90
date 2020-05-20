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
! Purpose: treatments after solution dump to Hdf files.
!
! Description: none.
!
! Input: globalGenx = global data structure (contains name of the window)
!
! Output: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: Fluid_postHdfOutput.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE Fluid_postHdfOutput( globalGenx )

  USE ModDataTypes
  USE ModGenx, ONLY       : t_globalGenx
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE ModParameters

  IMPLICIT NONE
  INCLUDE 'roccomf90.h'

! ... parameters
  TYPE(t_globalGenx), POINTER  :: globalGenx
  TYPE(t_global), POINTER      :: global
  TYPE(t_region), POINTER      :: regions(:)

! ... loop variables
  INTEGER :: iReg, iStat

! ... local variables
  INTEGER :: iLev
  REAL(RFREAL) :: eps

!******************************************************************************

  global  => globalGenx%global

  regions => globalGenx%levels(1)%regions

  CALL RegisterFunction( global,'Fluid_preHdfOutput',__FILE__ )

! set tav from accumulated values ---------------------------------------------
! to be moved to a wrapper

  eps = 100._RFREAL*EPSILON( 1._RFREAL )

#ifdef STATS
  DO iReg=1,global%nRegionsLocal
      IF ((global%flowType==FLOW_UNSTEADY) .AND. (global%doStat==ACTIVE)) THEN
        IF (global%integrTime > eps) THEN
          IF (global%mixtNStat > 0) THEN
            DO iStat=1,global%mixtNStat
              regions(iReg)%mixt%tav(iStat,:) = &
              regions(iReg)%mixt%tav(iStat,:)*global%integrTime
            ENDDO
          ENDIF  ! mixtNstat
        ENDIF    ! integrTime
      ENDIF      ! unsteady and dostat
  ENDDO          ! iReg
#endif

  CALL DeregisterFunction( global )

END SUBROUTINE Fluid_postHdfOutput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: Fluid_postHdfOutput.F90,v $
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
! Revision 1.2  2006/01/07 10:19:56  wasistho
! added Rocflu treatment
!
! Revision 1.1  2005/12/08 19:58:26  wasistho
! initial import
!
!******************************************************************************

