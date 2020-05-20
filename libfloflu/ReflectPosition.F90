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
! Purpose: Reflect position about plane.
!
! Description: None.
!
! Input: 
!   nx          x-component of normal vector of plane
!   ny          y-component of normal vector of plane
!   nz          z-component of normal vector of plane
!   xc          x-component of centroid of plane
!   yc          y-component of centroid of plane
!   zc          z-component of centroid of plane
!   xComp       x-component of position
!   yComp       y-component of position
!   zComp       z-component of position
!
! Output: 
!   xComp       x-component of reflected position
!   yComp       y-component of reflected position
!   zComp       z-component of reflected position
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ReflectPosition.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReflectPosition(nx,ny,nz,xc,yc,zc,xComp,yComp,zComp)

  USE ModDataTypes
  USE ModParameters
   
  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  REAL(RFREAL), INTENT(IN) :: nx,ny,nz,xc,yc,zc
  REAL(RFREAL), INTENT(INOUT) :: xComp,yComp,zComp

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  REAL(RFREAL) :: term 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: ReflectPosition.F90,v $ $Revision: 1.1.1.1 $'

! ******************************************************************************
! Reflect position
! ******************************************************************************

  term = 2.0_RFREAL*((xComp-xc)*nx + (yComp-yc)*ny + (zComp-zc)*nz)
  
  xComp = xComp - term*nx
  yComp = yComp - term*ny  
  zComp = zComp - term*nz
  
! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE ReflectPosition

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReflectPosition.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.1  2004/05/05 20:34:05  fnajjar
! Initial revision
!
! ******************************************************************************

