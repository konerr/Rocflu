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
! Purpose: Get variable location for mixture conserved variables.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_GetCvLoc.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

INTEGER FUNCTION RFLU_GetCvLoc(global,fluidModel,var)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
      
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: fluidModel,var
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_GetCvLoc.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction(global,'RFLU_GetCvLoc',__FILE__)

! ******************************************************************************
! Set variable info depending on fluid model
! ******************************************************************************

  SELECT CASE ( fluidModel ) 

! ==============================================================================
!   Incompressible fluid
! ==============================================================================    

    CASE ( FLUID_MODEL_INCOMP ) 
      SELECT CASE ( var )
        CASE ( CV_MIXT_XVEL )  
          RFLU_GetCvLoc = 1
        CASE ( CV_MIXT_YVEL ) 
          RFLU_GetCvLoc = 2
        CASE ( CV_MIXT_ZVEL )
          RFLU_GetCvLoc = 3
        CASE ( CV_MIXT_PRES ) 
          RFLU_GetCvLoc = 4
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! var
      
! ==============================================================================
!   Compressible fluid
! ==============================================================================    
    
    CASE ( FLUID_MODEL_COMP ) 
      SELECT CASE ( var )
        CASE ( CV_MIXT_DENS ) 
          RFLU_GetCvLoc = 1
        CASE ( CV_MIXT_XMOM:CV_MIXT_XVEL )  
          RFLU_GetCvLoc = 2        
        CASE ( CV_MIXT_YMOM:CV_MIXT_YVEL ) 
          RFLU_GetCvLoc = 3     
        CASE ( CV_MIXT_ZMOM:CV_MIXT_ZVEL )
          RFLU_GetCvLoc = 4       
        CASE ( CV_MIXT_ENER:CV_MIXT_PRES ) 
          RFLU_GetCvLoc = 5    
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! var   

! ==============================================================================
!   Default
! ==============================================================================    

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! fluidModel

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END FUNCTION RFLU_GetCvLoc

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_GetCvLoc.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:35  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:52  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/03/26 20:21:33  haselbac
! Removed error trap for GL model
!
! Revision 1.1  2004/11/06 03:16:43  haselbac
! Initial revision
!
! ******************************************************************************

