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
! Purpose: register portable random number variables with GenX.
!
! Description: none.
!
! Input: regions = patches and region (volume) variables
!
! Output: to Roccom via randInitGenxInterface.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RandInitGenxInterface.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE randInitGenxInterface( regions, wins, winv )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModDataStruct, ONLY : t_region
#endif
  USE ModError
  USE ModParameters
  USE ModRandom
  IMPLICIT NONE
  INCLUDE 'roccomf90.h'

! ... parameters
  CHARACTER(CHRLEN) :: wins, winv
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  INTEGER :: pid 
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'randInitGenxInterface',__FILE__ )

! input data (currently none) -------------------------------------------------

! output surface data (currently none) ----------------------------------------
 
! output restart data ---------------------------------------------------------

  CALL COM_new_attribute( TRIM(winv)//'.rand', 'p',         &
                          COM_INTEGER,1, '')
  CALL COM_set_size( TRIM(winv)//'.rand', 0, RAND_TOTAL_SIZE)

! store pointers to variables, loop over all regions --------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE          ) THEN    ! on my processor

! --- volume data -------------------------------------------------------------

      pid  = iReg*REGOFF
      CALL COM_set_array( TRIM(winv)//'.rand',pid,regions(iReg)%randData%mt(0))

    ENDIF      ! region on this processor and active
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE randInitGenxInterface

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RandInitGenxInterface.F90,v $
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
! Revision 1.1  2004/12/01 21:23:47  haselbac
! Initial revision after changing case
!
! Revision 1.5  2004/06/30 04:07:48  wasistho
! moved Genx related parameter REGOFF to ModParameters
!
! Revision 1.4  2004/06/29 23:53:25  wasistho
! migrated to Roccom-3
!
! Revision 1.3  2004/03/01 23:47:35  jiao
! Changed the F90 implementation for COM_init_attribute and COM_init_mesh
! not to require registered scalars to be F90 pointers.
!
! Revision 1.2  2004/02/20 00:49:04  jiao
! Changed to use COM_INIT_ATTRIBUTE_SPECIAL in order to compile with older IBM compilers.
!
! Revision 1.1  2003/11/21 22:20:04  fnajjar
! Initial import
!
!******************************************************************************

