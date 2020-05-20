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
! Purpose: Collection of routines for absorbing boundary condition.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModABC.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2009 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModABC

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModError
  USE ModMPI

  USE RFLU_ModExactFlow
  USE RFLU_ModFlowHardCode

  USE ModInterfaces, ONLY: MixtPerf_Eo_DGPUVW, &
                           MixtPerf_G_CpR, &
                           MixtPerf_R_M

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModABC.F90,v $'

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: RFLU_ABC_CreateSigma, &
            RFLU_ABC_DestroySigma, &
            RFLU_ABC_InitSigma, &
            RFLU_ABC_SetSigma, &
            RFLU_ABC_SetRefSoln, &
            RFLU_ABC_SourceTerms

! ==============================================================================
! Private functions
! ==============================================================================

  PRIVATE :: RFLU_ABC_NullifySigma

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS







! ******************************************************************************
!
! Purpose: Allocate memory for Sigma. 
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ABC_CreateSigma(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ABC_CreateSigma',__FILE__)

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  ALLOCATE(pRegion%mixt%sigma(pRegion%grid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%mixt%sigma')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ABC_CreateSigma









! ******************************************************************************
!
! Purpose: Destroy sigma variables.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ABC_DestroySigma(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ABC_DestroySigma',__FILE__)

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  DEALLOCATE(pRegion%mixt%sigma,STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pRegion%mixt%sigma')
  END IF ! global%error

  CALL RFLU_ABC_NullifySigma(pRegion)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ABC_DestroySigma









! ******************************************************************************
!
! Purpose: Initialize Sigma. 
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ABC_InitSigma(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: icg
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ABC_InitSigma',__FILE__)

! ******************************************************************************
! Initialize sigma
! ******************************************************************************

  DO icg = 1,pRegion%grid%nCellsTot
    pRegion%mixt%sigma(icg) = 0.0_RFREAL
  END DO ! icg

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ABC_InitSigma








! ******************************************************************************
!
! Purpose: Nullify sigma variable.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ABC_NullifySigma(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ABC_NullifySigma',__FILE__)

! ******************************************************************************
! Nullify memory
! ******************************************************************************

  NULLIFY(pRegion%mixt%sigma)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ABC_NullifySigma







! ******************************************************************************
!
! Purpose: Set Sigma. 
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ABC_SetSigma(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: icg
  REAL(RFREAL) :: C0,C1,C2,cx,cy,r,radius,rc,x1,x2,xc,y1,y2,yc
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ABC_SetSigma',__FILE__)

! ******************************************************************************
! Initialize sigma
! ******************************************************************************

  SELECT CASE( global%abcDomain )

! ==============================================================================
!   1D domain
! ==============================================================================

    CASE (0)
      x1 = global%abcLeft
      x2 = global%abcRight

      SELECT CASE( global%abcOrder )
        CASE (0)
         C0 = global%abcCoeff0

         DO icg = 1,pRegion%grid%nCellsTot
           xc = pRegion%grid%cofg(XCOORD,icg)

           IF ( (xc < x1) .OR. (xc > x2) ) THEN
             pRegion%mixt%sigma(icg) = global%abcSponge*(C0)
           END IF
         END DO ! icg

        CASE (1)
         C0 = global%abcCoeff0
         C1 = global%abcCoeff1

         DO icg = 1,pRegion%grid%nCellsTot
           xc = pRegion%grid%cofg(XCOORD,icg)

           IF ( xc < x1 ) THEN
             r = ABS(xc-x1) 
           ELSEIF (xc > x2 ) THEN
             r = ABS(xc-x2) 
           ELSE
             CYCLE
           END IF

           pRegion%mixt%sigma(icg) = global%abcSponge*(C0+C1*r)
         END DO ! icg

        CASE (2)
         C0 = global%abcCoeff0
         C1 = global%abcCoeff1
         C2 = global%abcCoeff2

         DO icg = 1,pRegion%grid%nCellsTot
           xc = pRegion%grid%cofg(XCOORD,icg)

           IF ( xc < x1 ) THEN
             r = ABS(xc-x1) 
           ELSEIF (xc > x2 ) THEN
             r = ABS(xc-x2) 
           ELSE
             CYCLE
           END IF

           pRegion%mixt%sigma(icg) = global%abcSponge*(C0+C1*r+C2*r*r)
         END DO ! icg

        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%abcOrder

! ==============================================================================
!   2D domain
! ==============================================================================

    CASE (1)
      x1 = global%abcLeft
      y1 = global%abcBottom
      x2 = global%abcRight
      y2 = global%abcTop

      SELECT CASE( global%abcOrder )
        CASE (0)
         C0 = global%abcCoeff0

         DO icg = 1,pRegion%grid%nCellsTot
           xc = pRegion%grid%cofg(XCOORD,icg)
           yc = pRegion%grid%cofg(YCOORD,icg)

           IF ( (xc < x1) .OR. (xc > x2) .OR. (yc < y1) .OR. (yc > y2) ) THEN
             pRegion%mixt%sigma(icg) = global%abcSponge*(C0)
           END IF
         END DO ! icg

        CASE (1)
         C0 = global%abcCoeff0
         C1 = global%abcCoeff1

         DO icg = 1,pRegion%grid%nCellsTot
           xc = pRegion%grid%cofg(XCOORD,icg)
           yc = pRegion%grid%cofg(YCOORD,icg)

           IF ( xc < x1 ) THEN
             IF ( yc < y1 ) THEN
               r = SQRT((xc-x1)*(xc-x1)+(yc-y1)*(yc-y1))
             ELSEIF (yc > y2 ) THEN
               r = SQRT((xc-x1)*(xc-x1)+(yc-y2)*(yc-y2))
             ELSE
               r = ABS(xc-x1)
             END IF
           ELSEIF (xc > x2 ) THEN
             IF ( yc < y1 ) THEN
               r = SQRT((xc-x2)*(xc-x2)+(yc-y1)*(yc-y1))
             ELSEIF (yc > y2 ) THEN
               r = SQRT((xc-x2)*(xc-x2)+(yc-y2)*(yc-y2))
             ELSE
               r = ABS(xc-x2)
             END IF
           ELSE
             IF ( yc < y1 ) THEN
               r = ABS(yc-y1)
             ELSEIF (yc > y2 ) THEN
               r = ABS(yc-y2)
             ELSE
               CYCLE
             END IF
           END IF

           pRegion%mixt%sigma(icg) = global%abcSponge*(C0+C1*r)
         END DO ! icg

        CASE (2)
         C0 = global%abcCoeff0
         C1 = global%abcCoeff1
         C2 = global%abcCoeff2

         DO icg = 1,pRegion%grid%nCellsTot
           xc = pRegion%grid%cofg(XCOORD,icg)
           yc = pRegion%grid%cofg(YCOORD,icg)

           IF ( xc < x1 ) THEN
             IF ( yc < y1 ) THEN
               r = SQRT((xc-x1)*(xc-x1)+(yc-y1)*(yc-y1))
             ELSEIF (yc > y2 ) THEN
               r = SQRT((xc-x1)*(xc-x1)+(yc-y2)*(yc-y2))
             ELSE
               r = ABS(xc-x1)
             END IF
           ELSEIF (xc > x2 ) THEN
             IF ( yc < y1 ) THEN
               r = SQRT((xc-x2)*(xc-x2)+(yc-y1)*(yc-y1))
             ELSEIF (yc > y2 ) THEN
               r = SQRT((xc-x2)*(xc-x2)+(yc-y2)*(yc-y2))
             ELSE
               r = ABS(xc-x2)
             END IF
           ELSE
             IF ( yc < y1 ) THEN
               r = ABS(yc-y1)
             ELSEIF (yc > y2 ) THEN
               r = ABS(yc-y2)
             ELSE
               CYCLE
             END IF
           END IF

           pRegion%mixt%sigma(icg) = global%abcSponge*(C0+C1*r+C2*r*r)
         END DO ! icg

        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%abcOrder

! ==============================================================================
!   Circular domain 
! ==============================================================================

    CASE (2)
      radius = global%abcRadius
      cx     = global%abcCenterX
      cy     = global%abcCenterY

      SELECT CASE( global%abcOrder )
        CASE (0)
         C0 = global%abcCoeff0

         DO icg = 1,pRegion%grid%nCellsTot
           xc = pRegion%grid%cofg(XCOORD,icg)
           yc = pRegion%grid%cofg(YCOORD,icg)
           rc = SQRT(xc*xc+yc*yc)

           IF ( rc > radius ) THEN
             pRegion%mixt%sigma(icg) = global%abcSponge*(C0)
           END IF
         END DO ! icg

        CASE (1)
         C0 = global%abcCoeff0
         C1 = global%abcCoeff1

         DO icg = 1,pRegion%grid%nCellsTot
           xc = pRegion%grid%cofg(XCOORD,icg)
           yc = pRegion%grid%cofg(YCOORD,icg)
           rc = SQRT(xc*xc+yc*yc)

           IF ( rc > radius ) THEN
             r = ABS(rc-radius) 

             pRegion%mixt%sigma(icg) = global%abcSponge*(C0+C1*r)
           END IF
         END DO ! icg

        CASE (2)
         C0 = global%abcCoeff0
         C1 = global%abcCoeff1
         C2 = global%abcCoeff2

         DO icg = 1,pRegion%grid%nCellsTot
           xc = pRegion%grid%cofg(XCOORD,icg)
           yc = pRegion%grid%cofg(YCOORD,icg)
           rc = SQRT(xc*xc+yc*yc)

           IF ( rc > radius ) THEN
             r = ABS(rc-radius) 

             pRegion%mixt%sigma(icg) = global%abcSponge*(C0+C1*r+C2*r*r)
           END IF
         END DO ! icg

        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%abcOrder

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pPatch%abcDomain

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ABC_SetSigma











! ******************************************************************************
!
! Purpose: Set reference solution for sponge layer. 
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ABC_SetRefSoln(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: errorString
  INTEGER :: iBc,icg,im,in,indCp,indMol,iq,n1,n2,n3
  REAL(RFREAL) :: A,A_c,A_cl,A1,A2,aTot,c,const,cp,cpRef,Cvm,cvg,cvl,cvv,d, &
                  dInc,dMin,dOffs,dTot,dummyReal,etaqm,g,gc,gcRef,gRef,gx,gy, &
                  gz,H,height,L,L_l,Lx,Ly,Lz,Mi,mInj,Mo,mw,nx,ny,nz,omega,p, &
                  pg,pi,pl,pMin,po,pOffs,psi,pTot,pv,r,r1,r2,Rc,Rc_l,radius, &
                  refD,refL,refNu,refP,refU,rg,rGas,ri,rl,ro,rv,rVap,t,tTot, &
                  theta,u,uo,uo_c,um,ur,ut,uz,v,vInj,Vm2,w,x,xMin,xo,xx,y,yp, &
                  yMin,z,zi,zMin
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvRef,pGv
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt_input), POINTER :: pMixtInput

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ABC_SetRefSoln',__FILE__)

  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Setting reference soln in sponge layer...'

    IF ( global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A,A)') SOLVER_NAME,'Case: ',TRIM(global%casename)
    END IF ! global%verbLevel                         
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid      => pRegion%grid
  pCvRef     => pRegion%mixt%cvRef
  pGv        => pRegion%mixt%gv
  pMixtInput => pRegion%mixtInput

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

! ******************************************************************************
! Set reference solution for sponge layer 
! ******************************************************************************

  IF ( global%abcDistrib == 0 ) THEN

! ******************************************************************************
!   Constant reference values in sponge layer
! ******************************************************************************

    IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
      pRegion%mixt%cvRef(CV_MIXT_DENS,0:1) = global%abcDens/global%refDensity
      pRegion%mixt%cvRef(CV_MIXT_XVEL,0:1) = global%abcVelX/global%refVelocity
      pRegion%mixt%cvRef(CV_MIXT_YVEL,0:1) = global%abcVelY/global%refVelocity
      pRegion%mixt%cvRef(CV_MIXT_ZVEL,0:1) = global%abcVelZ/global%refVelocity
      pRegion%mixt%cvRef(CV_MIXT_PRES,0:1) = (global%abcPress &
                                              - global%refPressure) &
                    /(global%refDensity*global%refVelocity*global%refVelocity)
    ELSE
      pRegion%mixt%cvRef(CV_MIXT_DENS,0:1) = global%abcDens
      pRegion%mixt%cvRef(CV_MIXT_XMOM,0:1) = global%abcVelX
      pRegion%mixt%cvRef(CV_MIXT_YMOM,0:1) = global%abcVelY
      pRegion%mixt%cvRef(CV_MIXT_ZMOM,0:1) = global%abcVelZ
      pRegion%mixt%cvRef(CV_MIXT_PRES,0:1) = global%abcPress
    END IF ! global%solverType

  ELSE

! ******************************************************************************
!   Variable reference values in sponge layer
! ******************************************************************************

    SELECT CASE ( pRegion%mixtInput%fluidModel )

! ==============================================================================
!     Compressible fluid model
! ==============================================================================  

      CASE ( FLUID_MODEL_COMP )

        SELECT CASE ( global%casename )

! ------------------------------------------------------------------------------
!       Proudman-Culick flow. NOTE this problem is two-dimensional and assumed 
!       to lie in the x-y plane, and that the injection boundary is located at 
!       y = -height.
! ------------------------------------------------------------------------------

        CASE ( "onera_c0", "onera_c0_2d_100x50" )
          CALL RFLU_GetParamsHardCodeProudman(dInc,mInj,vInj,pTot)

          height = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVertTot))
        
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            CALL RFLU_ComputeExactFlowProudman(global,x,y,height,dInc,vInj, &
                                               pTot,d,u,v,w,p)
          
            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              g  = global%refGamma
            ELSE
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)

              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)
            END IF ! solverType

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              pCvRef(CV_MIXT_DENS,icg) = d
              pCvRef(CV_MIXT_XVEL,icg) = u
              pCvRef(CV_MIXT_YVEL,icg) = v
              pCvRef(CV_MIXT_ZVEL,icg) = w
              pCvRef(CV_MIXT_PRES,icg) = p
            ELSE
              pCvRef(CV_MIXT_DENS,icg) = d
              pCvRef(CV_MIXT_XMOM,icg) = d*u
              pCvRef(CV_MIXT_YMOM,icg) = d*v
              pCvRef(CV_MIXT_ZMOM,icg) = d*w
              pCvRef(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
            END IF ! global%solverType
          END DO ! icg 

! ------------------------------------------------------------------------------
!       Culick flow 
! ------------------------------------------------------------------------------

        CASE ( "onera_c0_3d" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)

            CALL RFLU_ComputeExactFlowCulick(global,x,y,z,d,u,v,w,p)

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              g  = global%refGamma
            ELSE
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)

              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)
            END IF ! solverType

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              pCvRef(CV_MIXT_DENS,icg) = d
              pCvRef(CV_MIXT_XVEL,icg) = u
              pCvRef(CV_MIXT_YVEL,icg) = v
              pCvRef(CV_MIXT_ZVEL,icg) = w
              pCvRef(CV_MIXT_PRES,icg) = p
            ELSE
              pCvRef(CV_MIXT_DENS,icg) = d
              pCvRef(CV_MIXT_XMOM,icg) = d*u
              pCvRef(CV_MIXT_YMOM,icg) = d*v
              pCvRef(CV_MIXT_ZMOM,icg) = d*w
              pCvRef(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
            END IF ! global%solverType
          END DO ! icg          

! ------------------------------------------------------------------------------
!       Default - must be due to input error
! ------------------------------------------------------------------------------
          
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! global%casename

! ==============================================================================
!     Default   
! ==============================================================================  
              
      CASE DEFAULT 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%fluidModel

! ==============================================================================
! Non-dimensionalize reference solution for sponge layer for variable case 
! for Non-dissipative solver
! ==============================================================================

    IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
      DO icg = 1,pGrid%nCellsTot
        d = pCvRef(CV_MIXT_DENS,icg)
        u = pCvRef(CV_MIXT_XVEL,icg)
        v = pCvRef(CV_MIXT_YVEL,icg)
        w = pCvRef(CV_MIXT_ZVEL,icg)
        p = pCvRef(CV_MIXT_PRES,icg)

        x = pGrid%cofg(XCOORD,icg)
        zi = (x-0.4_RFREAL)/0.1_RFREAL
        IF ( zi > 1.0_RFREAL ) THEN
          zi = 1.0_RFREAL
        END IF ! zi

        d = d*(1.0_RFREAL-zi) + global%abcDens*(zi)
        u = u*(1.0_RFREAL-zi) + global%abcVelX*(zi)
        v = v*(1.0_RFREAL-zi) + global%abcVelY*(zi)
        w = w*(1.0_RFREAL-zi) + global%abcVelZ*(zi)
        p = p*(1.0_RFREAL-zi) + global%abcPress*(zi)

        pCvRef(CV_MIXT_DENS,icg) = d/global%refDensity
        pCvRef(CV_MIXT_XVEL,icg) = u/global%refVelocity
        pCvRef(CV_MIXT_YVEL,icg) = v/global%refVelocity
        pCvRef(CV_MIXT_ZVEL,icg) = w/global%refVelocity
        pCvRef(CV_MIXT_PRES,icg) = (p-global%refPressure) &
                     /(global%refDensity*global%refVelocity*global%refVelocity)
      END DO ! icg
    END IF ! global%solverType

  END IF ! global%abcDistrib 

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Setting reference soln in sponge layer done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ABC_SetRefSoln 









! *****************************************************************************
!
! Purpose: Add source terms due to absorbing boundary layer.
!
! Description: None.
!
! Input: 
!   region                 data of current region.
!
! Output: 
!   region%levels%mixt%rhs complete right-hand side (residual).
!
! Notes: None.
!
! *****************************************************************************

SUBROUTINE RFLU_ABC_SourceTerms(region)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModParameters

  USE ModInterfaces, ONLY: MixtPerf_Eo_DGPUVW, &
                           MixtPerf_G_CpR, &
                           MixtPerf_R_M

  IMPLICIT NONE

! *****************************************************************************
! Arguments
! *****************************************************************************

  TYPE(t_region), TARGET :: region

! *****************************************************************************
! Locals
! *****************************************************************************

  INTEGER :: distrib,errorFlag,icg,indCp,indMol
  REAL(RFREAL) :: cp,g,gc,mw,po,r,ro,rouo,rovo,rowo,roEo,ru,rv,rw,rE,sDens, &
                  sEner,sXmom,sYmom,sZmom,u,uo,v,vo,vol,w,wo
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pRhs
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_region), POINTER :: pRegion

! *****************************************************************************
! Start, set pointers
! *****************************************************************************

  pRegion => region
  global  => pRegion%global

  CALL RegisterFunction(global,'RFLU_ABC_SourceTerms',__FILE__)

  pGrid => pRegion%grid
  pCv   => pRegion%mixt%cv
  pRhs  => pRegion%mixt%rhs

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

  distrib = global%abcDistrib

! *****************************************************************************
! Loop over cells and compute source terms for absorbing boundary layer 
! *****************************************************************************

  DO icg = 1,pGrid%nCellsTot
    r  = pCv(CV_MIXT_DENS,icg)
    ru = pCv(CV_MIXT_XMOM,icg)
    rv = pCv(CV_MIXT_YMOM,icg)
    rw = pCv(CV_MIXT_ZMOM,icg)
    rE = pCv(CV_MIXT_ENER,icg)

    ro = pRegion%mixt%cvRef(CV_MIXT_DENS,distrib*icg)
    uo = pRegion%mixt%cvRef(CV_MIXT_XMOM,distrib*icg)
    vo = pRegion%mixt%cvRef(CV_MIXT_YMOM,distrib*icg)
    wo = pRegion%mixt%cvRef(CV_MIXT_ZMOM,distrib*icg)
    po = pRegion%mixt%cvRef(CV_MIXT_PRES,distrib*icg)

    mw = pRegion%mixt%gv(GV_MIXT_MOL,indMol*icg)
    cp = pRegion%mixt%gv(GV_MIXT_CP ,indCp *icg)
    gc = MixtPerf_R_M(mw)
    g  = MixtPerf_G_CpR(cp,gc)

    rouo = ro*uo
    rovo = ro*vo
    rowo = ro*wo
    roEo = ro*MixtPerf_Eo_DGPUVW(ro,g,po,uo,vo,wo) 

    vol = pGrid%vol(icg)

    sDens = -pRegion%mixt%sigma(icg)*vol*(r-ro)
    sXmom = -pRegion%mixt%sigma(icg)*vol*(ru-rouo)
    sYmom = -pRegion%mixt%sigma(icg)*vol*(rv-rovo)
    sZmom = -pRegion%mixt%sigma(icg)*vol*(rw-rowo)
    sEner = -pRegion%mixt%sigma(icg)*vol*(rE-roEo)

    pRhs(CV_MIXT_DENS,icg) = pRhs(CV_MIXT_DENS,icg) - sDens
    pRhs(CV_MIXT_XMOM,icg) = pRhs(CV_MIXT_XMOM,icg) - sXmom
    pRhs(CV_MIXT_YMOM,icg) = pRhs(CV_MIXT_YMOM,icg) - sYmom
    pRhs(CV_MIXT_ZMOM,icg) = pRhs(CV_MIXT_ZMOM,icg) - sZmom
    pRhs(CV_MIXT_ENER,icg) = pRhs(CV_MIXT_ENER,icg) - sEner
  END DO ! ic

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE RFLU_ABC_SourceTerms







END MODULE RFLU_ModABC

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModABC.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:36  brollin
! New Stable version
!
! Revision 1.1  2009/07/08 19:11:21  mparmar
! Initial revision
!
! ******************************************************************************
