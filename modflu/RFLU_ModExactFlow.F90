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
! Purpose: Suite of routines to compute exact flow solutions.
!
! Description: None.
!
! Notes: 
!   1. Collected routines in single module because computation of exact 
!      solutions is needed in at least three places: initialization of 
!      solution, setting of boundary profiles, and computation of errors. 
!
! ******************************************************************************
!
! $Id: RFLU_ModExactFlow.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModExactFlow

  USE ModParameters
  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
      
  IMPLICIT NONE
     
  PRIVATE
  PUBLIC :: RFLU_ComputeExactFlowAcoustic, &
            RFLU_ComputeExactFlowCAcoust, &
            RFLU_ComputeExactFlowCylPotential, &
            RFLU_ComputeExactFlowCulick, & 
            RFLU_ComputeExactFlowFreeVortex, &
            RFLU_ComputeExactFlowGaussianPulse, &
            RFLU_ComputeExactFlowPAcoust, &
            RFLU_ComputeExactFlowProudman, & 
            RFLU_ComputeExactFlowRadialPulse, &
            RFLU_ComputeExactFlowRayleighProblem, &
            RFLU_ComputeExactFlowRingleb, & 
            RFLU_ComputeExactFlowSphPotential, &
            RFLU_ComputeExactFlowSphStokes, &
            RFLU_ComputeExactFlowSsVortex, & 
            RFLU_ComputeExactFlowTaylorVortex, & 
            RFLU_SetExactFlowLinear, & 
            RFLU_SetExactFlowTrig

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModExactFlow.F90,v $ $Revision: 1.1.1.1 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS





! ******************************************************************************
!
! Purpose: Compute exact solution for Acoustic flow.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   x           x-coordinate
!   t           time 
!   ro          mean density
!   po          mean pressure
!   Mo          mean velocity
!   g           specific heat ratio
!   A1          amplitude
!   A2          amplitude
!
! Output:
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes:
!   1. Total pressure is not uniform!
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowAcoustic(global,x,t,ro,po,Mo,g,A1,A2,d,u,v, &
                                           w,p)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: A1,A2,g,Mo,po,ro,t,x
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    REAL(RFREAL) :: ao,theta1,theta2,uo

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

    ao = (g*po/ro)**0.5_RFREAL
    uo = Mo*ao

! Computing perturbation in mean flow
    theta1 = (x-(uo-ao)*t)
    theta2 = (x-(uo+ao)*t)
    theta1 = (x-(uo-ao)*(t+3.0E-2_RFREAL))
    theta2 = (x-(uo+ao)*(t+3.0E-2_RFREAL))

    u = A1*COS(theta1) + A2*COS(theta2)
    v = 0.0_RFREAL
    w = 0.0_RFREAL

    p = -u*ro*ao

    d = -u*ro/ao

! Adding perturbation to mean flow
    p = po + p
    u = uo + u
    d = ro + d 

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowAcoustic






! ******************************************************************************
!
! Purpose: Compute exact solution for channel acoustics.
!
! Description: None.
!
! Input: 
!   global      Pointer to global data
!   x           x-coordinate
!   y           y-coordinate
!   z           z-coordinate
!   t           time
!   Lx          Length of pipe in x-direction
!   Ly          Length of pipe in y-direction
!   Lz          Length of pipe in z-direction
!   n1          Index of mode in x-direction 
!   n2          Index of mode in y-direction
!   n3          Index of mode in z-direction
!   omega       Frequency
!   dTot        Total density
!   pTot        Total pressure
!   aTot        Total speed of sound
!   const       Constant is equal to amplitude
!
! Output: 
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes: 
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowCAcoust(global,x,y,z,t,Lx,Ly,Lz,n1,n2,n3, &
                                          omega,dTot,pTot,aTot,const,d,u,v,w,p)

    USE RFLU_ModBessel

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: n1,n2,n3 
    REAL(RFREAL), INTENT(IN) :: aTot,const,dTot,Lx,Ly,Lz,omega,pTot,t,x,y,z
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    REAL(RFREAL) :: cost,cosx,cosy,cosz,sint,sinx,siny,sinz,term

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

    cosx = COS(n1*global%pi*x/Lx)
    sinx = SIN(n1*global%pi*x/Lx)
    
    cosy = COS(n2*global%pi*y/Ly)
    siny = SIN(n2*global%pi*y/Ly)
    
    cosz = COS(n3*global%pi*z/Lz)
    sinz = SIN(n3*global%pi*z/Lz)
    
    cost  = COS(omega*t)
    sint  = SIN(omega*t)    
    
    term = const/(dTot*omega)
        
    p = const*cosx*cosy*cosz*cost
 
    d = dTot + p/(aTot*aTot)
    p = pTot + p

    u = (n1*global%pi/Lx)*term*sinx*cosy*cosz*sint
    v = (n2*global%pi/Ly)*term*cosx*siny*cosz*sint
    w = (n3*global%pi/Lz)*term*cosx*cosy*sinz*sint

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowCAcoust






! ******************************************************************************
!
! Purpose: Compute exact solution for potential flow around cylinder.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   x           x-coordinate
!   y           y-coordinate
!   ro          uniform density
!   uo          ambient flow velocity
!   po          ambient pressure
!   R           radius of cylinder
!
! Output:
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowCylPotential(global,x,y,ro,uo,po,R,d,u,v,w,p)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: po,R,ro,uo,x,y
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    REAL(RFREAL) :: rd2,theta

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

    theta = ATAN2(y,x)
    rd2   = R*R/(x*x+y*y)

    d = ro 
    u = uo-uo*rd2*COS(theta*2.0_RFREAL)
    v = -uo*rd2*SIN(theta*2.0_RFREAL)
    w = 0.0_RFREAL
    p = po+0.5_RFREAL*ro*uo*uo*rd2*(2.0_RFREAL*COS(theta*2.0_RFREAL) - rd2)

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowCylPotential








! ******************************************************************************
!
! Purpose: Compute exact solution for Culick flow.
!
! Description: None.
!
! Input: 
!   global      Pointer to global data
!   x           x-coordinate
!   y           y-coordinate
!   z		z-coordinate
!
! Output: 
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes: 
!   1. Assume cylindrical domain with axis along x-coordinate direction.
!   2. Assume radius equal to 0.01.
!   3. Assume injection velocity equal to unity.
!   4. Total pressure is not uniform but assumed so for now.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowCulick(global,x,y,z,d,u,v,w,p)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: x,y,z
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    REAL(RFREAL) :: r,ro,theta,vr

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

    r     = SQRT(y*y + z*z)
    ro    = 0.01_RFREAL
    theta = ATAN2(z,y)

    d  = 1.0_RFREAL
    
    vr = -SIN(0.5_RFREAL*global%pi*(r/ro)**2)/(r/ro)
    u  = global%pi*x/ro*COS(0.5_RFREAL*global%pi*(r/ro)**2)    
    v  = vr*COS(theta)
    w  = vr*SIN(theta)                               

!    p = 1.0E+5_RFREAL ! TEMPORARY
    p = 1.0E+5_RFREAL + 0.5_RFREAL*global%pi*(x/ro)**2.0_RFREAL

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowCulick






! ******************************************************************************
!
! Purpose: Compute exact solution for free vortex between two concentric 
!          cylinders.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   x           x-coordinate
!   y           y-coordinate
!   ro          uniform density
!   uo          ambient flow velocity
!   po          ambient pressure
!   R           radius of cylinder
!   g           Specific heat ratio
!
! Output:
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowFreeVortex(global,x,y,ro,uo,po,R,g,d,u,v,w,p)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: g,po,R,ro,uo,x,y
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    REAL(RFREAL) :: Mi2,radius,theta

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

    radius = SQRT(x*x+y*y)
    theta = ATAN2(y,x)

    Mi2 = uo*uo*ro/(g*po)

    d = ro*(1+0.5_RFREAL*(g-1.0_RFREAL)*Mi2*(1-(R/radius)**2.0_RFREAL) &
           )**(1.0_RFREAL/(g-1.0_RFREAL))
    u = (uo*R/radius)*SIN(theta)
    v = -(uo*R/radius)*COS(theta)
    w = 0.0_RFREAL
    p = po*(d/ro)**g

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowFreeVortex








! ******************************************************************************
!
! Purpose: Compute exact solution for Gaussian pulse flow.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   x           x-coordinate
!   ro          mean density
!   uo          mean velocity
!   po          mean pressure
!   c           speed of sound
!   g           specific heat ratio
!   L           length of domain
!   A           amplitude
!
! Output:
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes:
!   1. Total pressure is not uniform!
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowGaussianPulse(global,x,ro,uo,po,c,g,L,A,d, &
                                                u,v,w,p)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: A,c,g,L,po,ro,uo,x
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

    u = uo + A*(SIN(global%pi*x/L))**9.0_RFREAL
    v = 0.0_RFREAL
    w = 0.0_RFREAL

    p = po + (u - uo)*(ro*c)
! Comment : Both the following initialization are essentially same
!    d = ro*(p/po)**(1.0_RFREAL/g)
    d = ro + (u - uo)*(ro/c)

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowGaussianPulse






! ******************************************************************************
!
! Purpose: Compute exact solution for pipe acoustics.
!
! Description: None.
!
! Input: 
!   global      Pointer to global data
!   x           x-coordinate
!   y           y-coordinate
!   z           z-coordinate
!   t           time
!   L           Length of pipe
!   ro          Outer radius of pipe
!   iBc         Type of boundary condition in axial direction (0 = zero
!               pressure disturbance, 1 = zero pressure disturbance gradient)
!   im          Index of mode in circumferential direction
!   in          Index of mode in axial direction
!   iq          Index of mode in radial direction
!   etaqm       q-th root of m-th order Bessel function of first kind
!   omega       Frequency
!   dTot        Total density
!   pTot        Total pressure
!   aTot        Total speed of sound
!   const       Constant (appears in theoretical solution, but is not equal 
!               to amplitude)
!
! Output: 
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes: 
!   1. Exact solution assumes that z-axis corresponds to axis of symmetry of 
!      pipe.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowPAcoust(global,x,y,z,t,L,ro,iBc,im,in,iq, &
                                          etaqm,omega,dTot,pTot,aTot,const, & 
                                          d,u,v,w,p)

    USE RFLU_ModBessel

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iBc,im,in,iq
    REAL(RFREAL), INTENT(IN) :: aTot,const,dTot,etaqm,L,omega,pTot,ro,t,x,y,z
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    REAL(RFREAL) :: cimt,cinz,cot,dummyReal,Jm,Jmd,r,simt,sinz,sot,term, & 
                    theta,ur,ut,uz

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

    r = SQRT(x*x + y*y)
    theta = ATAN2(y,x)

    CALL RFLU_JYNDD(im,etaqm*r/ro,Jm,Jmd,dummyReal,dummyReal,dummyReal, & 
                    dummyReal)
    
    cimt = COS(im*theta)
    simt = SIN(im*theta)
    
    cinz = COS(in*global%pi*z/L)
    sinz = SIN(in*global%pi*z/L)    
    
    cot  = COS(omega*t)
    sot  = SIN(omega*t)    
    
    term = const/(dTot*omega)
        
    IF ( iBc == 0 ) THEN     
      p = const*Jm*cimt*sinz*cot
    ELSE IF ( iBc == 1 ) THEN 
      p = const*Jm*cimt*cinz*cot
    ELSE 
      p = 0.0_RFREAL
    END IF ! iBc
    
    d = dTot + p/aTot*aTot
    p = pTot + p

    IF ( iBc == 0 ) THEN   
      ur = -   term  *Jmd*(etaqm/ro)*cimt*sinz                 *sot
      ut =  im*term/r*Jm            *simt*sinz                 *sot
      uz = -   term  *Jm            *cimt*cinz*(in*global%pi/L)*sot    
    ELSE IF ( iBc == 1 ) THEN
      ur = -   term  *Jmd*(etaqm/ro)*cimt*cinz                 *sot
      ut =  im*term/r*Jm            *simt*cinz                 *sot
      uz =     term  *Jm            *cimt*sinz*(in*global%pi/L)*sot
    ELSE 
      ur = 0.0_RFREAL
      ut = 0.0_RFREAL
      uz = 0.0_RFREAL      
    END IF ! iBc

    u = ur*COS(theta) - ut*SIN(theta)
    v = ur*SIN(theta) + ut*COS(theta)
    w = uz

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowPAcoust






! ******************************************************************************
!
! Purpose: Compute exact solution for Proudman-Culick flow.
!
! Description: None.
!
! Input: 
!   global      Pointer to global data
!   x           x-coordinate
!   y           y-coordinate
!   height      Height of domain
!   dInc        Density
!   vInj        Injection velocity
!   pTot        Total pressure at (x,y) = (0,0)
!
! Output: 
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes: 
!   1. Total pressure is not uniform!
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowProudman(global,x,y,height,dInc,vInj,pTot, &
                                           d,u,v,w,p)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: dInc,height,pTot,vInj,x,y
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

    u = -0.5_RFREAL*global%pi*x/height*vInj*COS(0.5_RFREAL*global%pi*y/height)
    v =                                vInj*SIN(0.5_RFREAL*global%pi*y/height)

    w = 0.0_RFREAL

    d = dInc
    p = pTot - 0.25_RFREAL*vInj**2*(1.0_RFREAL &
                                  + 0.5_RFREAL*(global%pi*x/height)**2 & 
                                  -         COS(global%pi*y/height))

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowProudman









! ******************************************************************************
!
! Purpose: Compute exact solution for Radial pulse flow.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   radius      radial distance from origin
!   theta       angle w.r.t x-axis
!   ro          mean density
!   uo          mean velocity
!   po          mean pressure
!   c           speed of sound
!   g           specific heat ratio
!   L           length of domain
!   A           amplitude
!
! Output:
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes:
!   1. Total pressure is not uniform!
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowRadialPulse(global,radius,theta,ro,uo,po,c, &
                                              g,L,A,d,u,v,w,p)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: A,c,g,L,po,radius,ro,theta,uo
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

    u = uo + A*(COS(global%pi*radius/L))**9.0_RFREAL
    v = 0.0_RFREAL
    w = 0.0_RFREAL

    p = po + (u - uo)*(ro*c)
! Comment : Both the following initialization are essentially same
!    d = ro*(p/po)**(1.0_RFREAL/g)
    d = ro + (u - uo)*(ro/c)

    v = (u - uo)*SIN(theta)
    u = uo + (u - uo)*COS(theta)

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowRadialPulse








! ******************************************************************************
!
! Purpose: Compute exact solution for Rayleigh problem flow.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   t           time
!   y           vertical distance from flatplate
!   ro          uniform density
!   uo          max velocity
!   po          uniform pressure
!   muo         uniform viscosity
!
! Output:
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes:
!   1. Total pressure is not uniform!
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowRayleighProblem(global,t,y,ro,uo,po,muo, &
                                                  d,u,v,w,p)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: muo,po,ro,t,uo,y
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

    p = po
    d = ro 

    IF (t == 0.0_RFREAL) THEN
      u = uo
    ELSE
      u = uo*ERF(y/(2.0_RFREAL*SQRT(muo*t/ro)))
    END IF ! t

    v = 0.0_RFREAL
    w = 0.0_RFREAL

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowRayleighProblem








! ******************************************************************************
!
! Purpose: Compute exact solution for Ringleb flow.
!
! Description: The solution to the Ringleb flow is given in terms of qBar 
!   (normalized velocity magnitude) and the k (streamline constant). It is
!   not possible to determine these as an explicit function of x and y. 
!   Instead, given x and y, qBar and k are determined iteratively. Once 
!   qBar and k are known, the solution can be computed in terms of the 
!   primitive variables. 
!
! Input: 
!   x           x-coordinate
!   y           y-coordinate
!   rGas        Gas constant
!   pTot        Total pressure
!   tTot        Total temperature
!
! Output: 
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes: 
!   1. This routine assumes a perfect gas and that the ratio of specific heats
!      is equal to 1.4.
!   2. Depending on how the geometry is defined (which values of qBar and k 
!      were chosen), it may become necessary to alter the initial guess for 
!      qBar.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowRingleb(x,y,rGas,pTot,tTot,d,u,v,w,p)

    USE ModInterfaces, ONLY: MixtPerf_C_GRT,MixtPerf_D_PRT

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: pTot,rGas,tTot,x,y
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w

! ==============================================================================
! Locals
! ==============================================================================

    CHARACTER(CHRLEN) :: RCSIdentString
    INTEGER :: cntr,cntrMax
    REAL(RFREAL) :: a,aBar,aTot,alpha,daBardqBar,dFdqBar,dJdaBar,dJdqBar, &
                    dqBar,drhoBardqBar,F,J,k,qBar,qBarMax,rhoBar,t,tol

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

! ==============================================================================
!   Determine solution in terms of qBar and k with Newton-Raphson method
! ==============================================================================

    cntr    = 0
    cntrMax = 100 ! Maximum number of Newton-Raphson steps

    qBar    = 0.3_RFREAL       ! Initial guess for qBar
    qBarMax = SQRT(5.0_RFREAL) ! Physically possible maximum value for qBar 

    tol = 1.0E-16_RFREAL ! Convergence tolerance

! -------------------------------------------------------------------------------
!   Loop over iterations
! -------------------------------------------------------------------------------

    DO 
      cntr = cntr + 1

      aBar   = SQRT(1.0_RFREAL - 0.2_RFREAL*qBar**2)                        
      rhoBar = aBar**5.0_RFREAL

      J      = 1.0_RFREAL/aBar                 & 
             + 1.0_RFREAL/(3.0_RFREAL*aBar**3) & 
             + 1.0_RFREAL/(5.0_RFREAL*aBar**5) &
             - 0.5_RFREAL*LOG((1.0_RFREAL+aBar)/(1.0_RFREAL-aBar))              

! --- Compute function value ---------------------------------------------------                

      F = (x - 0.5_RFREAL*J)**2 + y**2 - 1.0_RFREAL/(4.0_RFREAL*rhoBar**2*qBar**4)                

! --- Compute derivative -------------------------------------------------------        

      drhoBardqBar = -qBar*aBar**3      
      dJdaBar = -1.0_RFREAL/aBar**2 & 
                -1.0_RFREAL/aBar**4 & 
                -1.0_RFREAL/aBar**6 &
                -1.0_RFREAL/((1.0_RFREAL+aBar)*(1.0_RFREAL-aBar))
      daBardqBar = -0.2_RFREAL*qBar/aBar
      dJdqBar = dJdaBar*daBardqBar        
      dFdqBar = -(x - 0.5_RFREAL*J)*dJdqBar & 
              + 0.5_RFREAL/(rhoBar**3*qBar**4)*drhoBardqBar & 
              + 1.0_RFREAL/(rhoBar**2*qBar**5)

! --- Update -------------------------------------------------------------------        

      dqBar = -F/dFdqBar         

      qBar = qBar + dqBar                                                                                 

! --- Limit to physically possible maximum -------------------------------------

      qBar = MIN(qBar,qBarMax)                   

! --- Convergence check and update k -------------------------------------------        

      IF ( (ABS(dqBar) < tol ) .OR. (cntr == cntrMax) ) THEN
        k = SQRT(2.0_RFREAL/(1.0_RFREAL/qBar**2 & 
                           - 2.0_RFREAL*rhoBar*(x - 0.5_RFREAL*J)))

        EXIT
      END IF ! ABS(dqBar)                   
    END DO ! <empty>

! ==============================================================================
!   Determine solution in terms of primitive variables
! ==============================================================================

    p    = pTot*rhoBar**1.4_RFREAL
    t    = tTot*rhoBar**0.4_RFREAL          
    d    = MixtPerf_D_PRT(p,rGas,t)      
    aTot = MixtPerf_C_GRT(1.4_RFREAL,rGas,tTot)
    a    = aTot*aBar

    alpha = ASIN(MIN(1.0_RFREAL,qBar/k))  
    u     = aTot*qBar*COS(alpha)*SIGN(1.0_RFREAL,y)
    v     = aTot*qBar*SIN(alpha) 
    w     = 0.0_RFREAL

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowRingleb








! ******************************************************************************
!
! Purpose: Compute exact solution for potential flow around sphere in x-y plane.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   x           x-coordinate
!   y           y-coordinate
!   ro          uniform density
!   uo          ambient flow velocity
!   po          ambient pressure
!   R           radius of cylinder
!
! Output:
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowSphPotential(global,x,y,ro,uo,po,R,d,u,v,w,p)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: po,R,ro,uo,x,y
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    REAL(RFREAL) :: rd3,rd6,theta

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

    theta = ATAN2(y,x)
    rd3   = R*R*R/((x*x+y*y)*SQRT(x*x+y*y))
    rd6   = rd3*rd3

    d = ro 
    u = uo*(1.0_RFREAL-0.25_RFREAL*rd3-0.75_RFREAL*rd3*COS(theta*2.0_RFREAL))
    v = -uo*0.75_RFREAL*rd3*SIN(theta*2.0_RFREAL)
    w = 0.0_RFREAL
    p = po+0.5_RFREAL*ro*uo*uo*((3.0_RFREAL*rd3-0.75_RFREAL*rd6)* &
                                (COS(theta))**2.0_RFREAL-(rd3+0.25*rd6))

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowSphPotential








! ******************************************************************************
!
! Purpose: Compute exact solution for stokes flow around sphere in x-y plane.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   x           x-coordinate
!   y           y-coordinate
!   ro          uniform density
!   uo          ambient flow velocity
!   po          ambient pressure
!   R           radius of cylinder
!   muo         coefficient of viscosity
!
! Output:
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowSphStokes(global,x,y,ro,uo,po,R,muo,d,u,v,w,p)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: muo,po,R,ro,uo,x,y
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    REAL(RFREAL) :: radius,rd,rd3,theta,ur,ut

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

    theta = ATAN2(y,x)
    radius= SQRT(x*x+y*y)
    rd    = R/SQRT(x*x+y*y)
    rd3   = R*R*R/((x*x+y*y)*SQRT(x*x+y*y))

    ur = uo*( 1.0_RFREAL+0.50_RFREAL*rd3-1.5_RFREAL*rd)*COS(theta)
    ut = uo*(-1.0_RFREAL+0.25_RFREAL*rd3+.75_RFREAL*rd)*SIN(theta)

    d = ro 
    u = ur*COS(theta) - ut*SIN(theta) 
    v = ur*SIN(theta) + ut*COS(theta) 
    w = 0.0_RFREAL
    p = po-1.5_RFREAL*muo*rd*uo*COS(theta)/radius

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowSphStokes








! ******************************************************************************
!
! Purpose: Compute exact solution for supersonic vortex.
!
! Description: None.
!
! Input: 
!   x           x-coordinate
!   y           y-coordinate
!   gGas        ratio of specific heats
!   rGas        Gas constant
!   ri          Inner radius
!   Mi          Mach number at inner radius
!   pTot        Total pressure
!   tTot        Total temperature
!
! Output: 
!   d           Density
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes: 
!   1. This routine assumes a perfect gas.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowSsVortex(x,y,gGas,rGas,ri,Mi,pTot,tTot, &
                                           d,u,v,w,p)

    USE ModInterfaces, ONLY: MixtPerf_C_GRT, & 
                             MixtPerf_D_DoGMa, &
                             MixtPerf_D_PRT, & 
                             MixtPerf_P_GMaPo, &
                             MixtPerf_P_DDoGPo, &
                             MixtPerf_T_DPR 

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: gGas,Mi,pTot,rGas,ri,tTot,x,y
    REAL(RFREAL), INTENT(OUT) :: d,p,u,v,w

! ==============================================================================
!   Locals
! ==============================================================================

    REAL(RFREAL) :: ai,alpha,di,dTot,pi,r,term,ti,Vi

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************

    dTot = MixtPerf_D_PRT(pTot,rGas,tTot)

    di = MixtPerf_D_DoGMa(dTot,gGas,Mi)      
    pi = MixtPerf_P_GMaPo(gGas,Mi,pTot)
    ti = MixtPerf_T_DPR(di,pi,rGas) 
    ai = MixtPerf_C_GRT(gGas,rGas,ti)
    Vi = Mi*ai

    r     = SQRT(x*x + y*y)
    alpha = ATAN(y/x) 

    term = 1.0_RFREAL - (ri/r)**2
    term = 1.0_RFREAL + 0.5_RFREAL*(gGas - 1.0_RFREAL)*Mi*Mi*term
    d    = di*term**(1.0_RFREAL/(gGas - 1.0_RFREAL))
    p    = MixtPerf_P_DDoGPo(d,dTot,gGas,pTot)

    u =  Vi*ri/r*SIN(alpha)
    v = -Vi*ri/r*COS(alpha)
    w =  0.0_RFREAL

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowSsVortex







! ******************************************************************************
!
! Purpose: Compute exact solution for Taylor vortex.
!
! Description: None.
!
! Input: 
!   t           time 
!   pi          value of pi 
!   x           x-coordinate
!   y           y-coordinate
!   L           Domain length/2
!   refL        reference length
!   refNu       reference coefficient of kinematic viscosity
!   refU        reference velocity
!   refD        reference density
!   refP        reference pressure
!
! Output: 
!   u           x-velocity component
!   v           y-velocity component
!   w           z-velocity component
!   p           Static pressure
!
! Notes: 
!   1. This routine assumes a perfect gas.
!
! ******************************************************************************

  SUBROUTINE RFLU_ComputeExactFlowTaylorVortex(t,pi,x,y,L,refL,refNu,refU,refD, &
                                               refP,u,v,w,p)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), INTENT(IN) :: L,pi,refD,refL,refNu,refP,refU,t,x,y
    REAL(RFREAL), INTENT(OUT) :: p,u,v,w

! ******************************************************************************
!   Start, compute exact solution
! ******************************************************************************
  
    u = -refU*COS(x/refL)*SIN(y/refL)*EXP(-2.0_RFREAL*refNu*t/(refL*refL))
    v =  refU*SIN(x/refL)*COS(y/refL)*EXP(-2.0_RFREAL*refNu*t/(refL*refL))
    w =  0.0_RFREAL
    p = -0.25_RFREAL*refD*refU*refU*(COS(2.0_RFREAL*x/refL) &
                                    +COS(2.0_RFREAL*y/refL)) &
                                   *EXP(-4.0_RFREAL*refNu*t/(refL*refL)) &
                                   + refP

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ComputeExactFlowTaylorVortex







! ******************************************************************************
!
! Purpose: Set linear variable behavior.
!
! Description: None.
!
! Input: 
!   x           x-coordinate
!   y           y-coordinate
!   z           z-coordinate
!   iVar        Variable index
!
! Output: 
!   var         Variable
!   gx          x-component of gradient
!   gy          y-component of gradient
!   gz          z-component of gradient
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_SetExactFlowLinear(x,y,z,iVar,var,gx,gy,gz)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iVar
    REAL(RFREAL), INTENT(IN) :: x,y,z
    REAL(RFREAL), INTENT(OUT) :: gx,gy,gz,var

! ==============================================================================
!   Locals
! ==============================================================================

    REAL(RFREAL) :: a,b,c

! ******************************************************************************
!   Start, set exact solution
! ******************************************************************************

    a = 1 + (iVar-1)*(ZCOORD-XCOORD)
    b = a + 1
    c = a + 2

    var = a*x + b*y + c*z
    
    gx = a
    gy = b
    gz = c 

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_SetExactFlowLinear








! ******************************************************************************
!
! Purpose: Set trigonometric variable behavior.
!
! Description: None.
!
! Input: 
!   global      Pointer to globa data
!   nx          Wave number for x-direction
!   ny          Wave number for y-direction
!   nz          Wave number for z-direction
!   x           x-coordinate
!   y           y-coordinate
!   z           z-coordinate
!   iVar        Variable index
!
! Output: 
!   var         Variable
!   gx          x-component of gradient
!   gy          y-component of gradient
!   gz          z-component of gradient
!
! Notes: 
!   1. At present, do not make use of iVar.
!
! ******************************************************************************

  SUBROUTINE RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,iVar,var,gx,gy,gz)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iVar
    REAL(RFREAL), INTENT(IN) :: nx,ny,nz,x,y,z
    REAL(RFREAL), INTENT(OUT) :: gx,gy,gz,var
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    REAL(RFREAL) :: a,b,c

! ******************************************************************************
!   Start, set exact solution
! ******************************************************************************

    a = nx*global%pi
    b = ny*global%pi
    c = nz*global%pi

    var = COS(a*x)*SIN(b*y)*COS(c*z)
    
    gx = -a*SIN(a*x)*SIN(b*y)*COS(c*z)
    gy =  b*COS(a*x)*COS(b*y)*COS(c*z)
    gz = -c*COS(a*x)*SIN(b*y)*SIN(c*z) 

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_SetExactFlowTrig






! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModExactFlow


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModExactFlow.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:36  brollin
! New Stable version
!
! Revision 1.9  2010/03/14 23:47:31  mparmar
! Added free-vortex and Stokes flow
!
! Revision 1.8  2009/09/28 14:21:40  mparmar
! Added rayleigh problem, potential flow over cylinder and sphere cases
!
! Revision 1.7  2009/08/13 01:21:21  mparmar
! Added pressure variation in RFLU_ComputeExactFlowCulick
!
! Revision 1.6  2009/07/08 19:11:47  mparmar
! Added RFLU_ComputeExactFlowRadialPulse and removed TEMPORARY
!
! Revision 1.5  2008/12/06 08:43:40  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:16:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/12/03 18:22:44  mparmar
! Changing exact solution to compute sin in gaussianpulse
!
! Revision 1.2  2007/11/28 23:05:22  mparmar
! Added acoustic, Gaussian pulse and Taylor-Green vortex cases
!
! Revision 1.1  2007/04/09 18:49:24  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:40  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.8  2006/04/13 18:07:34  haselbac
! Added routine for Culick flow
!
! Revision 1.7  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.6  2006/01/06 22:10:26  haselbac
! Added routines for linear and trig solution for grad testing
!
! Revision 1.5  2005/04/29 12:54:36  haselbac
! Adapted routine for exact pipe acoust solution to accept time argument
!
! Revision 1.4  2005/04/20 14:41:53  haselbac
! Bug fix and extensions for pipe acoustics case
!
! Revision 1.3  2005/03/31 17:02:44  haselbac
! Fixed bug in initialization of pressure for ONERA C0 case
!
! Revision 1.2  2005/03/15 20:44:44  haselbac
! Added routine to compute exact solution for pipe acoustics
!
! Revision 1.1  2004/07/06 15:14:27  haselbac
! Initial revision
!
! ******************************************************************************

