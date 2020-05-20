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
! Purpose: Collection of utility routines for program burn. 
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: SPEC_RFLU_ModPBAUtils.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2006-2007 by the University of Illinois
!
! ******************************************************************************

MODULE SPEC_RFLU_ModPBAUtils

  USE ModDataTypes  

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: SPEC_RFLU_ModPBAUtils.F90,v $'
  
! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: SPEC_RFLU_PBA_ComputeY, &
            SPEC_RFLU_PBA_Eo_GUd, &
            SPEC_RFLU_PBA_GetCJStates, &
            SPEC_RFLU_PBA_GetSolution

! ==============================================================================
! Private functions
! ==============================================================================
  
  !PRIVATE ::

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  




! ******************************************************************************
!
! Purpose: Compute Y fraction of products for program burn.
!
! Description: None.
!
! Input:
!   x         Centroid of the cell
!   xFront    Detonation front
!   xHalfW    Half of cell width
!   xWidth    Detonation width
!   xTail     Detonation tail
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  FUNCTION SPEC_RFLU_PBA_ComputeY(x,xFront,xHalfW,xWidth,xTail)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    REAL(RFREAL), INTENT(IN) :: x,xFront,xHalfW,xWidth,xTail
    REAL(RFREAL) :: SPEC_RFLU_PBA_ComputeY
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    REAL(RFREAL) :: xRange1,xRange2,xRange3,xx1,xx2,Y1,Y2,Y3

! ******************************************************************************
!   Start
! ******************************************************************************
   
! Case 1 ----------------------------------------------------------------------
    IF ( x <= xTail - xHalfW ) THEN
      SPEC_RFLU_PBA_ComputeY = 1.0_RFREAL

! Case 2 ----------------------------------------------------------------------
    ELSEIF ( x <= xTail ) THEN
      xRange1 = xTail - (x - xHalfW)
      Y1      = 1.0_RFREAL

      xRange2 = (x + xHalfW) - xTail
      xx1     = 0.0_RFREAL
      xx2     = ((x + xHalfW) - xTail)/xWidth
      Y2      = 1.0_RFREAL - (xx1 + xx2)/2.0_RFREAL ! Linear profile
!      Y2      = 1.0_RFREAL - (xx1*xx1 + xx2*xx2 + xx1*xx2)/3.0_RFREAL ! Quadratic profile

      xRange3 = 0.0_RFREAL
      Y3      = 0.0_RFREAL

      IF ( ABS(xWidth) < ABS(xRange2) ) THEN
        xRange2 = xWidth
        Y2      = 0.5_RFREAL ! Linear profile average
!        Y2      = 2.0_RFREAL/3.0_RFREAL ! Quadratic profile average

        xRange3 = (x + xHalfW) - xFront
        Y3      = 0.0_RFREAL
      END IF

      SPEC_RFLU_PBA_ComputeY = (xRange1*Y1+xRange2*Y2+xRange3*Y3)/(2.0_RFREAL*xHalfW) 

! Case 3 ----------------------------------------------------------------------
    ELSEIF ( x >= xFront + xHalfW ) THEN
      SPEC_RFLU_PBA_ComputeY = 0.0_RFREAL

! Case 4 ----------------------------------------------------------------------
    ELSEIF ( x >= xFront ) THEN
      xRange1 = 0.0_RFREAL
      Y1      = 1.0_RFREAL

      xRange2 = xFront - (x - xHalfW)
      xx1     = ((x - xHalfW) - xTail)/xWidth
      xx2     = 1.0_RFREAL
      Y2      = 1.0_RFREAL - (xx1 + xx2)/2.0_RFREAL ! Linear profile
!      Y2      = 1.0_RFREAL - (xx1*xx1 + xx2*xx2 + xx1*xx2)/3.0_RFREAL ! Quadratic profile

      xRange3 = (x + xHalfW) - xFront
      Y3      = 0.0_RFREAL

      IF ( ABS(xWidth) < ABS(xRange2) ) THEN
        xRange1 = xTail - (x - xHalfW)
        Y1      = 1.0_RFREAL

        xRange2 = xWidth
        Y2      = 0.5_RFREAL ! Linear profile average
!        Y2      = 2.0_RFREAL/3.0_RFREAL ! Quadratic profile average
      END IF

      SPEC_RFLU_PBA_ComputeY = (xRange1*Y1+xRange2*Y2+xRange3*Y3)/(2.0_RFREAL*xHalfW) 

! Case 5 ----------------------------------------------------------------------
    ELSEIF ( x <= xTail + xHalfW ) THEN
      xRange1 = xTail - (x - xHalfW)
      Y1      = 1.0_RFREAL

      xRange2 = (x + xHalfW) - xTail
      xx1     = 0.0_RFREAL
      xx2     = ((x + xHalfW) - xTail)/xWidth
      Y2      = 1.0_RFREAL - (xx1 + xx2)/2.0_RFREAL ! Linear profile
!      Y2      = 1.0_RFREAL - (xx1*xx1 + xx2*xx2 + xx1*xx2)/3.0_RFREAL ! Quadratic profile
      
      xRange3 = 0.0_RFREAL
      Y3      = 0.0_RFREAL

      IF ( ABS(xWidth) < ABS(xRange2) ) THEN
        xRange2 = xWidth
        Y2      = 0.5_RFREAL ! Linear profile average
!        Y2      = 2.0_RFREAL/3.0_RFREAL ! Quadratic profile average
        
        xRange3 = (x + xHalfW) - xFront
        Y3      = 0.0_RFREAL
      END IF

      SPEC_RFLU_PBA_ComputeY = (xRange1*Y1+xRange2*Y2+xRange3*Y3)/(2.0_RFREAL*xHalfW) 

! Case 6 ----------------------------------------------------------------------
    ELSEIF ( x >= xFront - xHalfW ) THEN
      xRange1 = 0.0_RFREAL
      Y1      = 1.0_RFREAL
      
      xRange2 = xFront - (x - xHalfW)
      xx1     = ((x - xHalfW) - xTail)/xWidth
      xx2     = 1.0_RFREAL
      Y2      = 1.0_RFREAL - (xx1 + xx2)/2.0_RFREAL ! Linear profile
!      Y2      = 1.0_RFREAL - (xx1*xx1 + xx2*xx2 + xx1*xx2)/3.0_RFREAL ! Quadratic profile
      
      xRange3 = (x + xHalfW) - xFront
      Y3      = 0.0_RFREAL

      IF ( ABS(xWidth) < ABS(xRange2) ) THEN
        xRange1 = xTail - (x - xHalfW)
        Y1      = 1.0_RFREAL
        
        xRange2 = xWidth
        Y2      = 0.5_RFREAL ! Linear profile average
!        Y2      = 2.0_RFREAL/3.0_RFREAL ! Quadratic profile average
      END IF

      SPEC_RFLU_PBA_ComputeY = (xRange1*Y1+xRange2*Y2+xRange3*Y3)/(2.0_RFREAL*xHalfW) 

! Case 7 ----------------------------------------------------------------------
    ELSE
      xx1 = ((x - xHalfW) - xTail)/xWidth
      xx2 = ((x + xHalfW) - xTail)/xWidth

      SPEC_RFLU_PBA_ComputeY = 1.0_RFREAL - (xx1 + xx2)/2.0_RFREAL ! Linear profile
!      SPEC_RFLU_PBA_ComputeY = 1.0_RFREAL - (xx1*xx1 + xx2*xx2 + xx1*xx2)/3.0_RFREAL ! Quadratic profile
    END IF

! ******************************************************************************
!   End  
! ******************************************************************************

  END FUNCTION SPEC_RFLU_PBA_ComputeY








! ******************************************************************************
!
! Purpose: Compute energy of explosive for program burn.
!
! Description: None.
!
! Input:
!   g             Ratio of specific heats
!   ud            Detonation velocity
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  FUNCTION SPEC_RFLU_PBA_Eo_GUd(g,ud)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    REAL(RFREAL), INTENT(IN) :: g,ud
    REAL(RFREAL) :: SPEC_RFLU_PBA_Eo_GUd
    
! ******************************************************************************
!   Start
! ******************************************************************************
   
    SPEC_RFLU_PBA_Eo_GUd = ((ud)**2.0_RFREAL)/(2.0_RFREAL*(g*g-1.0_RFREAL))

! ******************************************************************************
!   End  
! ******************************************************************************

  END FUNCTION SPEC_RFLU_PBA_Eo_GUd








! ******************************************************************************
!
! Purpose: Compute CJ state for program burn.
!
! Description: None.
!
! Input:
!   g             Ratio of specific heats
!   ro            Density of explosive
!   ud            Detonation velocity
!
! Output:
!   e             Internal energy of products
!   r             Density of products
!   u             Velocity of products
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE SPEC_RFLU_PBA_GetCJStates(g,ro,ud,e,r,u)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    REAL(RFREAL), INTENT(IN) :: g,ro,ud
    REAL(RFREAL), INTENT(OUT) :: e,r,u
    
! ******************************************************************************
!   Start
! ******************************************************************************
  
    e = g*ud*ud/((g+1.0_RFREAL)*(g*g-1.0_RFREAL))
    r = ro*(g+1.0_RFREAL)/g
    u = ud/(g+1.0_RFREAL)

! ******************************************************************************
!   End  
! ******************************************************************************

  END SUBROUTINE SPEC_RFLU_PBA_GetCJStates








! ******************************************************************************
!
! Purpose: Compute solution for program burn.
!
! Description: None.
!
! Input:
!   g             Ratio of specific heats
!   gc            Gas constant
!   ro            Density of explosive
!   po            Pressure in explosive
!   ud            Detonation velocity
!   Y             Mass fraction of products
!
! Output:
!   e             Internal energy
!   p             Pressure
!   r             Density
!   T             Temperature
!   u             Velocity
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE SPEC_RFLU_PBA_GetSolution(g,gc,ro,po,ud,Y,e,p,r,T,u)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    REAL(RFREAL), INTENT(IN) :: g,gc,ro,po,ud,Y
    REAL(RFREAL), INTENT(OUT) :: e,p,r,T,u
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    REAL(RFREAL) :: alpha,nTol
    
! ******************************************************************************
!   Start
! ******************************************************************************
  
    nTol = 1.0E-14_RFREAL

    alpha = po/(ro*ud*ud)

    u = (ud/(2.0_RFREAL+(g-1.0_RFREAL)*Y))*(1.0_RFREAL &
          -SQRT(1.0_RFREAL-((2.0_RFREAL+(g-1.0_RFREAL)*Y)*Y)/(g+1.0_RFREAL)))

    e = (u*(ud-u)/((g-1.0_RFREAL)*Y))*(1.0_RFREAL+alpha)*(1.0_RFREAL+alpha)
    p = ro*ud*u*(1.0_RFREAL+alpha)
    r = (ro*ud/(ud-u))/(1.0_RFREAL+alpha)

    u = (1.0_RFREAL+alpha)*(ud/(2.0_RFREAL+(g-1.0_RFREAL)*Y))*(1.0_RFREAL &
          -SQRT(1.0_RFREAL-((2.0_RFREAL+(g-1.0_RFREAL)*Y)*Y)/(g+1.0_RFREAL))) &
        - alpha*ud  

    IF ( ABS(Y) < nTol ) THEN
      u = 0.0_RFREAL
      e = ud*ud*(1.0_RFREAL+2.0_RFREAL*alpha+g*g*alpha*alpha) &
          /(2.0_RFREAL*(g*g-1.0_RFREAL))
      p = po
      r = ro
    END IF

    T = p/(gc*r)

! ******************************************************************************
!   End  
! ******************************************************************************

  END SUBROUTINE SPEC_RFLU_PBA_GetSolution








END MODULE SPEC_RFLU_ModPBAUtils

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_RFLU_ModPBAUtils.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
! ******************************************************************************

