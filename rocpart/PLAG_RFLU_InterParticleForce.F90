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
! Purpose: Correct the mixture properties using higher-order
!          interpolation schemed onto particle locations
!
! Description: none.
!
! Input: pRegion = current region.
!
! Output: region%levels%plag%dv 
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_RFLU_InterParticleForce.F90,v 1.4 2016/03/22 19:31:40 fred Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_InterParticleForce(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModMixture, ONLY    : t_mixt
  USE ModError  
  USE ModParameters
  USE PLAG_ModParameters

  USE ModInterfaces, ONLY: MixtPerf_R_M, & 
                           MixtPerf_T_DPR

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
  INTEGER :: icg,iPcl,nPcls
  INTEGER, POINTER, DIMENSION(:) :: pCvPlagMass
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv
  
  REAL(RFREAL) :: dx,dy,dz,massL
  REAL(RFREAL), POINTER, DIMENSION(:,:)    :: pCv,pPlagRhs,pVFracL
  REAL(RFREAL), POINTER, DIMENSION(:,:,:) :: pGradVFracL
  
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag),   POINTER :: pPlag 
  TYPE(t_mixt),   POINTER :: pMixt 

! TEMPORARY: Manoj: Adding interparticle force  
  REAL(RFREAL) :: beta,diamL,factor,phi,pi,phiCP,Ps,volL,C1,r
! END TEMPORARY

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_InterParticleForce.F90,v $'

  global => pRegion%global
  
  CALL RegisterFunction( global, 'PLAG_RFLU_InterParticleForce',__FILE__ )

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pMixt => pRegion%mixt 
  pGrid => pRegion%grid
  pPlag => pRegion%plag

  pi = global%pi

  nPcls = pRegion%plag%nPcls

! ******************************************************************************
! Set pointers for discrete particles 
! ******************************************************************************

  pAiv        => pPlag%aiv
  pCv         => pPlag%cv
  pVFracL     => pPlag%vFracL
  pPlagRhs    => pPlag%rhs
  pGradVFracL => pPlag%gradVFracL
  pCvPlagMass => pPlag%cvPlagMass

! ******************************************************************************
! Set parameters for interparticle force
! ******************************************************************************

  ! BBR - begin - getting rid unexplained collision model coding
  !Ps    = 1.0E5_RFREAL ! Value give by Stanley
  Ps    = pRegion%plagInput%collisionPs
 ! WRITE(999,*) 'This is my current Ps value = ',Ps
  C1 = 1.08D0
  r = 3.0D-1
  beta  = 3.0_RFREAL
  phiCP = 0.6_RFREAL

! ******************************************************************************        
! Correct discrete particle dv
! Force for a single particle NOT super-particle.
! -ve of force on the particle goes to rhs.
! ******************************************************************************
    
  DO iPcl = 1, nPcls
    icg = pAiv(AIV_PLAG_ICELLS,iPcl)

    diamL = pPlag%dv(DV_PLAG_DIAM,iPcl)
    volL  = (pi*(diamL**3.0_RFREAL)/6.0_RFREAL)

    massL = SUM( pCv(pCvPlagMass(:),iPcl) )

    phi = pVFracL(1,iPcl)
! TEMPORARY: Manoj
!    dx = pCv(CV_PLAG_XPOS,iPcl) - pGrid%cofg(XCOORD,icg)
!    dy = pCv(CV_PLAG_YPOS,iPcl) - pGrid%cofg(YCOORD,icg)
!    dz = pCv(CV_PLAG_ZPOS,iPcl) - pGrid%cofg(ZCOORD,icg)    
!    phi = pVFracL(1,iPcl) + pGradVFracL(XCOORD,1,iPcl)*dx &
!                          + pGradVFracL(YCOORD,1,iPcl)*dy &
!                          + pGradVFracL(ZCOORD,1,iPcl)*dz
! END TEMPORARY

! TEMPORARY: Manoj: 2010-03-24: Scaling up Ps according to volume fraction
!    ps = 10.0_RFREAL**(6.0_RFREAL+TANH(5.0_RFREAL*(phi-0.3_RFREAL))) ! Worked for 300 micro meter
!    ps = 10.0_RFREAL**(6.0_RFREAL+TANH(7.0_RFREAL*(phi-0.2_RFREAL))) ! Did not work for 80
!    ps = 10.0_RFREAL**(6.0_RFREAL+TANH(20.0_RFREAL*(phi-0.1_RFREAL))) ! Trying for 80
!    ps = 10.0_RFREAL**(8.0_RFREAL+TANH(20.0_RFREAL*(phi-0.1_RFREAL))) ! 2012-05-18
! BBR - begin - getting rid of experimental formulation of collision model
!    ps = 10.0_RFREAL**(LOG10(pRegion%plagInput%collisionPs) &
!                       +TANH(20.0_RFREAL*(phi-0.1_RFREAL))) ! 2012-05-18
! BBR - end 
!    ps = 10.0_RFREAL**( LOG10(pRegion%plagInput%collisionPs) &
!                        *0.5_RFREAL*( 1.0_RFREAL + TANH(20.0_RFREAL*(phi-0.1_RFREAL)) ) &
!                      ) ! 2012-05-23
!    WRITE(*,*) "LOG10(pRegion%plagInput%collisionPs)=",LOG10(pRegion%plagInput%collisionPs)
!    WRITE(*,*) "Ps=",Ps
!    STOP
! END TEMPORARY

!-- Manoj: 2012-05-29: Making denomenator +ve definite using ABS(phiCP-phi)
   IF (pRegion%mixtInput%prepIntVal19 == 0) THEN
    factor = volL*Ps*( beta*phi**(beta-2.0_RFREAL)/(ABS(phiCP-phi)) &
                          + phi**(beta-1.0_RFREAL)/(ABS(phiCP-phi))**2.0_RFREAL)
   ELSE IF (pRegion%mixtInput%prepIntVal19 == 1) THEN 
     factor = C1*(phi*phi*phi)*exp((r*phi)/(phiCP-phi))
   END IF

    pPlagRhs(CV_PLAG_XMOM,iPcl) = pPlagRhs(CV_PLAG_XMOM,iPcl) &
                                + factor*pGradVFracL(XCOORD,1,iPcl)
    pPlagRhs(CV_PLAG_YMOM,iPcl) = pPlagRhs(CV_PLAG_YMOM,iPcl) &
                                + factor*pGradVFracL(YCOORD,1,iPcl)
    pPlagRhs(CV_PLAG_ZMOM,iPcl) = pPlagRhs(CV_PLAG_ZMOM,iPcl) &
                                + factor*pGradVFracL(ZCOORD,1,iPcl)
! TEMPORARY: Manoj
! WRITE(*,*) "iPcl=",iPcl," =================================="
! WRITE(*,'(2(A,1X,E12.6,5X))') "x=",pCv(CV_PLAG_XPOS,iPcl),"    y=",pCv(CV_PLAG_YPOS,iPcl)
! WRITE(*,'(2(A,1X,E12.6,5X))') "acc-x=",-factor*pGradVFracL(XCOORD,1,iPcl)/massL, &
!                               "acc-y=",-factor*pGradVFracL(YCOORD,1,iPcl)/massL
!
!    pPlagRhs(CV_PLAG_XMOM,iPcl) = -8.0E7_RFREAL*massL 
!    pPlagRhs(CV_PLAG_YMOM,iPcl) = 0.0_RFREAL
!    pPlagRhs(CV_PLAG_ZMOM,iPcl) = 0.0_RFREAL
! END TEMPORARY    
  END DO ! iPcl
     
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLU_InterParticleForce

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_InterParticleForce.F90,v $
! Revision 1.4  2016/03/22 19:31:40  fred
! Adding Glasser's model as a choice for collision modeling
!
! Revision 1.3  2015/08/12 19:41:42  brollin
! Updating module declaration in rfluinit.F90
!
! Revision 1.2  2015/07/23 23:11:19  brollin
! 1) The pressure coefficient of the  collision model has been changed back to its original form
! 2) New options in the format of the grid and solutions have been added. Now the user can choose the endianness, and convert from one to the over in rfluconv.
! 3) The solutions are now stored in folders named by timestamp or iteration number.
! 4) The address enty in the hashtable has been changed to an integer(8) for cases when the grid becomes very large.
! 5) RFLU_WritePM can now compute PM2 on the fly for the Macroscale problem
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
!******************************************************************************

