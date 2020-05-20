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
! Purpose: Write optimal LES statistics in ASCII ROCFLU format.
!
! Description: None.
!
! Input: 
!   region      region data
!
! Output: None.
!
! Notes: 
!   1. Works only for single regions right now
!
!******************************************************************************
!
! $Id: RFLU_WriteStatsFileOLES.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_WriteStatsFileOLES(global)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError    

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
! Local variables
! ==============================================================================
  
  CHARACTER(CHRLEN) :: RCSIdentString  
  
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global
  
! ******************************************************************************
! Start, set pointer
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_WriteStatsFileOLES.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction(global,'RFLU_WriteStatsFileOLES',__FILE__)

! ******************************************************************************
! Write data
! ******************************************************************************

  WRITE(IF_STATS_OLES,'(6(2X,E15.8))') & 
    global%currentTime,global%enerOLES,global%dissOLES, & 
    global%uVarOLES,global%vVarOLES,global%wVarOLES

  CALL DeregisterFunction(global)
     
! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE RFLU_WriteStatsFileOLES


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLU_WriteStatsFileOLES.F90,v $
!   Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
!   merged rocflu micro and macro
!
!   Revision 1.1.1.1  2014/07/15 14:31:38  brollin
!   New Stable version
!
!   Revision 1.3  2008/12/06 08:43:48  mtcampbe
!   Updated license.
!
!   Revision 1.2  2008/11/19 22:17:01  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.1  2007/04/09 18:49:57  haselbac
!   Initial revision after split from RocfloMP
!
!   Revision 1.1  2007/04/09 18:01:02  haselbac
!   Initial revision after split from RocfloMP
!
!   Revision 1.1  2002/09/09 16:28:02  haselbac
!   Initial revision
!
! ******************************************************************************

