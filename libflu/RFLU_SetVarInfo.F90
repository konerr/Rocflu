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
! Purpose: Set variable info for mixture.
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
! $Id: RFLU_SetVarInfo.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_SetVarInfo(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
      
  USE ModInterfaces, ONLY: RFLU_GetCvLoc    
      
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: cvMixtPres,cvMixtXVel,cvMixtYVel,cvMixtZVel
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_SetVarInfo.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SetVarInfo',__FILE__)

! ******************************************************************************
! Set variable info depending on fluid model
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%fluidModel ) 

! ==============================================================================
!   Incompressible fluid
! ==============================================================================    

    CASE ( FLUID_MODEL_INCOMP ) 
      cvMixtXVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_XVEL)
      cvMixtYVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_YVEL)
      cvMixtZVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_ZVEL)
      cvMixtPres = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_PRES)                  
      
      pRegion%mixt%cvInfo(cvMixtXVel) = VAR_INFO_POSNEG
      pRegion%mixt%cvInfo(cvMixtYVel) = VAR_INFO_POSNEG
      pRegion%mixt%cvInfo(cvMixtZVel) = VAR_INFO_POSNEG
      pRegion%mixt%cvInfo(cvMixtPres) = VAR_INFO_POS           
      
! ==============================================================================
!   Compressible fluid
! ==============================================================================    
    
    CASE ( FLUID_MODEL_COMP ) 
      pRegion%mixt%cvInfo(CV_MIXT_DENS) = VAR_INFO_POS
      pRegion%mixt%cvInfo(CV_MIXT_XVEL) = VAR_INFO_POSNEG
      pRegion%mixt%cvInfo(CV_MIXT_YVEL) = VAR_INFO_POSNEG
      pRegion%mixt%cvInfo(CV_MIXT_ZVEL) = VAR_INFO_POSNEG
      pRegion%mixt%cvInfo(CV_MIXT_PRES) = VAR_INFO_POS      

! ==============================================================================
!   Default
! ==============================================================================    

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixtInput%fluidModel

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SetVarInfo

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_SetVarInfo.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:36  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:53  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:50  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/03/26 20:21:47  haselbac
! Removed error trap for GL model
!
! Revision 1.3  2004/11/14 19:44:00  haselbac
! Generalized setting of variable info
!
! Revision 1.2  2004/11/06 03:17:34  haselbac
! Cosmetics only
!
! Revision 1.1  2004/11/02 02:26:50  haselbac
! Initial revision
!
! ******************************************************************************

