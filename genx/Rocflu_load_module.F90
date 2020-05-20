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
! Purpose: Open window to Rocflu and register interface functions.
!
! Description: None.
!
! Input: 
!   winName	Name of Rocflu window
!
! Output: None.
! 
! Notes: None.
!
! ******************************************************************************
!
! $Id: Rocflu_load_module.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE Rocflu_load_module(winName)

  USE ModGenx, ONLY: t_globalGenx, &
                     associate_pointer
  
  USE ModInterfaces, ONLY : RFLU_InitFlowSolver, &
                            RFLU_FlowSolverDummy, &
                            RFLU_FlowSolver, &
                            RFLU_UpdateInbuffGm, &
                            Fluid_compute_integrals, &
                            Fluid_finalize, &
                            Fluid_preHdfOutput, &
                            Fluid_postHdfOutput
  IMPLICIT NONE

  INCLUDE 'roccomf90.h'
  INCLUDE 'rocmanf90.h'

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Interface
! ==============================================================================

  INTERFACE
    SUBROUTINE COM_set_pointer(attr,ptr,asso)
      USE ModGenx, ONLY: t_globalGenx
      CHARACTER(*), INTENT(IN) :: attr
      TYPE(t_globalGenx), POINTER  :: ptr
      EXTERNAL asso
    END SUBROUTINE COM_set_pointer
  END INTERFACE

! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(*), INTENT(IN) :: winName

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: isDummy
  INTEGER :: types(7)
  TYPE(t_globalGenx), POINTER :: glb

! ******************************************************************************
! Start
! ******************************************************************************

  ALLOCATE(glb)
  ALLOCATE(glb%global)

  isDummy            = TRIM(winName) == "RocfluDummy"
  glb%isDummy        = isDummy
  glb%global%winName = winName

! ******************************************************************************
! Initialize fluids window
! ******************************************************************************

  CALL COM_new_window(winName)

! ******************************************************************************
! Create an attribute for global data
! ******************************************************************************

  CALL COM_new_attribute(winName//'.global','w',COM_F90POINTER,1,'')
  CALL COM_allocate_array(winName//'.global')

! ******************************************************************************
! Solver initialization
! ******************************************************************************

  types(1) = COM_F90POINTER
  types(2) = COM_DOUBLE_PRECISION
  types(3) = COM_MPI_COMM
  types(4) = COM_INTEGER
  types(5) = COM_STRING
  types(6) = COM_STRING
  types(7) = COM_INTEGER

  CALL COM_set_member_function(winName//'.initialize',RFLU_InitFlowSolver, &
                               winName//'.global','biiiiii',types)

! ******************************************************************************
! Solution update
! ******************************************************************************

  types(1) = COM_F90POINTER
  types(2) = COM_DOUBLE_PRECISION
  types(3) = COM_DOUBLE_PRECISION  
  types(4) = COM_INTEGER
  types(5) = COM_INTEGER

  IF ( isDummy .EQV. .TRUE. ) THEN
    CALL COM_set_member_function(winName//'.update_solution', &
                                RFLU_FlowSolverDummy,winName//'.global', & 
                                  'biiii',types)
  ELSE
    CALL COM_set_member_function(winName//'.update_solution', &
                                 RFLU_FlowSolver,winName//'.global','biiii', & 
                                 types)
  END IF ! isDummy

! ******************************************************************************
! Integral update
! ******************************************************************************

  types(1) = COM_F90POINTER
  types(2) = COM_DOUBLE_PRECISION

  CALL COM_set_member_function(winName//'.'//MAN_COMP_INTEG_NAME, &
                               Fluid_compute_integrals,winName//'.global', &
                               'bo',types)

! ******************************************************************************
! Standalone mode
! ******************************************************************************

  CALL COM_set_member_function(winName//'.update_inbuff_gm_fluid', &
                               RFLU_UpdateInbuffGm,winName//'.global','bi', & 
                               types)

! ******************************************************************************
! Shutdown
! ******************************************************************************

  CALL COM_set_member_function(winName//'.finalize',Fluid_finalize, &
                               winName//'.global','b',types)

! ******************************************************************************
! Pre and Post IO
! ******************************************************************************

  CALL COM_set_member_function(winName//'.pre_hdf_output', &
                               Fluid_preHdfOutput,winName//'.global','b', &
                               types)

  CALL COM_set_member_function(winName//'.post_hdf_output', &
                               Fluid_postHdfOutput,winName//'.global','b', &
                               types)

! ******************************************************************************
! Finish
! ******************************************************************************

  CALL COM_window_init_done(winName)

  CALL COM_set_pointer(winName//'.global',glb,associate_pointer )

! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE Rocflu_load_module


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: Rocflu_load_module.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:30  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:47:51  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:58:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2006/01/07 10:18:47  wasistho
! added Fluid_pre/postHdfOutput
!
! Revision 1.2  2005/02/03 23:06:58  jiao
! Changed third argument to RFLU_InitFlowSolver to COM_MPI_COMM.
!
! Revision 1.1  2004/12/01 21:23:58  haselbac
! Initial revision after changing case
!
! Revision 1.5  2004/10/19 19:23:17  haselbac
! Cosmetics only
!
! Revision 1.4  2003/12/07 04:06:26  jiao
! Changed the prototype of COM_set_pointer to be the same as COM_get_pointer,
! so that it takes ASSOCIATE_POINTER as the third argument. This change is
! needed to support PGI compilers.
!
! Revision 1.3  2002/11/15 21:23:01  haselbac
! Added integral update routine registration
!
! Revision 1.2  2002/10/17 19:57:15  haselbac
! Adapted argument lists of solver routines
!
! Revision 1.1  2002/10/05 18:28:58  haselbac
! Initial revision
!
! ******************************************************************************

