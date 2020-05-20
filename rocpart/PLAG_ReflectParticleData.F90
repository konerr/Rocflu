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
! Purpose: Reflect particle data about patch face. 
!
! Description: None.
!
! Input: 
!  pPatch       Pointer to patch
!  pPlag        Pointer to particle data
!  ifl          Local face index
!  iPcl         Particle index
!  xLocOld      x-component of old particle location
!  yLocOld      y-component of old particle location
!  zLocOld      z-component of old particle location    
!  xLoc         x-component of particle location
!  yLoc         y-component of particle location
!  zLoc         z-component of particle location
!  xTraj        x-component of particle trajectory
!  yTraj        y-component of particle trajectory
!  zTraj        z-component of particle trajectory
!
! Output: 
!  xLocOld      x-component of reflected old particle location
!  yLocOld      y-component of reflected old particle location
!  zLocOld      z-component of reflected old particle location  
!  xLoc         x-component of reflected particle location
!  yLoc         y-component of reflected particle location
!  zLoc         z-component of reflected particle location
!  xTraj        x-component of reflected particle trajectory
!  yTraj        y-component of reflected particle trajectory
!  zTraj        z-component of reflected particle trajectory
!
! Notes: 
!   1. Reflect particle data about centroid of face instead of about point of
!      intersection of particle path and plane defined by face. This is ok 
!      because vector between centroid and intersection lies in plane defined 
!      by face. 
!
! ******************************************************************************
!
! $Id: PLAG_ReflectParticleData.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_ReflectParticleData(pPatch,pPlag,ifl,iPcl,xLocOld,yLocOld, & 
                                    zLocOld,xLoc,yLoc,zLoc,xTraj,yTraj,zTraj)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModPartLag, ONLY: t_plag
  USE ModMPI
  
  USE PLAG_ModParameters    
   
  USE ModInterfaces, ONLY: ReflectPosition, & 
                           ReflectVector 
   
  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(IN) :: ifl,iPcl
  REAL(RFREAL), INTENT(INOUT) :: xLoc,xLocOld,yLoc,yLocOld,zLoc,zLocOld
  REAL(RFREAL), INTENT(INOUT), OPTIONAL :: xTraj,yTraj,zTraj                                 
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_plag), POINTER :: pPlag

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  REAL(RFREAL) :: nx,ny,nz,xc,xMom,xRhsSum,yc,yMom,yRhsSum,zc,zMom,zRhsSum

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ReflectParticleData.F90,v $ $Revision: 1.1.1.1 $'

! ******************************************************************************
! Reflect particle data
! ******************************************************************************

! ==============================================================================  
! Get face normal and centroid
! ==============================================================================  

  nx = pPatch%fn(XCOORD,ifl)
  ny = pPatch%fn(YCOORD,ifl)
  nz = pPatch%fn(ZCOORD,ifl)    

  xc = pPatch%fc(XCOORD,ifl)
  yc = pPatch%fc(YCOORD,ifl)
  zc = pPatch%fc(ZCOORD,ifl)    

! ==============================================================================  
! Reflect old position, new position, and trajectory
! ==============================================================================  

  CALL ReflectPosition(nx,ny,nz,xc,yc,zc,xLocOld,yLocOld,zLocOld)
  CALL ReflectPosition(nx,ny,nz,xc,yc,zc,xLoc,yLoc,zLoc)
  
  IF ( (PRESENT(xTraj) .EQV. .TRUE.) .AND. & 
       (PRESENT(yTraj) .EQV. .TRUE.) .AND. &
       (PRESENT(zTraj) .EQV. .TRUE.) ) THEN 
    CALL ReflectVector(nx,ny,nz,xTraj,yTraj,zTraj)
  END IF ! PRESENT
  
! ==============================================================================  
! Reflect particle momentum (and any other data)
! ==============================================================================  
  
  xMom = pPlag%cv(CV_PLAG_XMOM,iPcl)
  yMom = pPlag%cv(CV_PLAG_YMOM,iPcl) 
  zMom = pPlag%cv(CV_PLAG_ZMOM,iPcl)
  
  CALL ReflectVector(nx,ny,nz,xMom,yMom,zMom)  
  
  pPlag%cv(CV_PLAG_XMOM,iPcl) = xMom 
  pPlag%cv(CV_PLAG_YMOM,iPcl) = yMom
  pPlag%cv(CV_PLAG_ZMOM,iPcl) = zMom  

  CALL ReflectVector(nx,ny,nz,xMom,yMom,zMom)  
  xMom = pPlag%cvOld(CV_PLAG_XMOM,iPcl)
  yMom = pPlag%cvOld(CV_PLAG_YMOM,iPcl) 
  zMom = pPlag%cvOld(CV_PLAG_ZMOM,iPcl)
  
  CALL ReflectVector(nx,ny,nz,xMom,yMom,zMom)  
  
  pPlag%cvOld(CV_PLAG_XMOM,iPcl) = xMom 
  pPlag%cvOld(CV_PLAG_YMOM,iPcl) = yMom
  pPlag%cvOld(CV_PLAG_ZMOM,iPcl) = zMom  

! ==============================================================================  
! Reflect particle rhsSum, needed for RK4 and RK3
! ==============================================================================  
  
  xRhsSum = pPlag%rhsSum(CV_PLAG_XMOM,iPcl)
  yRhsSum = pPlag%rhsSum(CV_PLAG_YMOM,iPcl) 
  zRhsSum = pPlag%rhsSum(CV_PLAG_ZMOM,iPcl)
  
  CALL ReflectVector(nx,ny,nz,xRhsSum,yRhsSum,zRhsSum)  
  
  pPlag%rhsSum(CV_PLAG_XMOM,iPcl) = xRhsSum 
  pPlag%rhsSum(CV_PLAG_YMOM,iPcl) = yRhsSum
  pPlag%rhsSum(CV_PLAG_ZMOM,iPcl) = zRhsSum  
  
  xRhsSum = pPlag%rhsSum(CV_PLAG_XPOS,iPcl)
  yRhsSum = pPlag%rhsSum(CV_PLAG_YPOS,iPcl) 
  zRhsSum = pPlag%rhsSum(CV_PLAG_ZPOS,iPcl)
  
  CALL ReflectVector(nx,ny,nz,xRhsSum,yRhsSum,zRhsSum)  
  
  pPlag%rhsSum(CV_PLAG_XPOS,iPcl) = xRhsSum 
  pPlag%rhsSum(CV_PLAG_YPOS,iPcl) = yRhsSum
  pPlag%rhsSum(CV_PLAG_ZPOS,iPcl) = zRhsSum

! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE PLAG_ReflectParticleData

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ReflectParticleData.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:52  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:50:27  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:34  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.6  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.5  2005/01/01 21:34:29  haselbac
! Made traj arguments optional
!
! Revision 1.4  2004/11/30 23:34:10  fnajjar
! Included further description in comment for reflection of rhsSum
!
! Revision 1.3  2004/06/17 14:32:16  fnajjar
! Applied reflectivity to rhsSum vectors of positions and momenta
!
! Revision 1.2  2004/05/05 21:01:22  fnajjar
! Bug fixes: Reflect position properly and reflect old momentum
!
! Revision 1.1  2004/04/08 01:32:13  haselbac
! Initial revision
!
! ******************************************************************************

