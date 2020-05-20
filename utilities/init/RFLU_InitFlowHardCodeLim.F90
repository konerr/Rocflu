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
! Purpose: Initialize part of flow field in a region using hard-coded values.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine assumes a perfect gas.
!
! ******************************************************************************
!
! $Id: RFLU_InitFlowHardCodeLim.F90,v 1.4 2016/03/04 16:32:01 rahul Exp $
!
! Copyright: (c) 2005-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_InitFlowHardCodeLim(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModParameters

  USE ModInterfaces, ONLY: MixtPerf_Eo_DGPUVW, &
                           MixtPerf_G_CpR, &
                           MixtPerf_R_M
    
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
  INTEGER :: flag,icg,indCp,indMol,l,m,n
  REAL(RFREAL) :: cp,d,g,gaussAmp,gaussCenter,gaussWidth,gc,Lmn,p,perturb,mw,&
                  radius,the,u,v,w,x,y,z,zz,Eo
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: randNum
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pGv
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt_input), POINTER :: pMixtInput

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_InitFlowHardCodeLim.F90,v $ $Revision: 1.4 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_InitFlowHardCodeLim', &
                        __FILE__)
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Initializing flow field from limited hard code...'
                             
    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,A)') SOLVER_NAME,'Case: ',TRIM(global%casename)
    END IF ! global%verbLevel                         
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid      => pRegion%grid
  pCv        => pRegion%mixt%cv
  pGv        => pRegion%mixt%gv
  pMixtInput => pRegion%mixtInput
  
  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol    
  
! ******************************************************************************
! Initialize flow field based on user input and fluid model
! ******************************************************************************

  SELECT CASE ( pMixtInput%fluidModel )
  
! ==============================================================================
!   Incompressible fluid model
! ==============================================================================  
  
    CASE ( FLUID_MODEL_INCOMP ) 
      pRegion%mixt%cvState = CV_MIXT_STATE_PRIM
      
! TEMPORARY, to be replaced by proper code
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) 
! END TEMPORARY

! ==============================================================================
!   Compressible fluid model
! ==============================================================================  
    
    CASE ( FLUID_MODEL_COMP )
      pRegion%mixt%cvState = CV_MIXT_STATE_CONS
    
      SELECT CASE ( global%casename )


! Subbu - Dec 2011
! ------------------------------------------------------------------------------
!       Detonation wave in 3D - Cylindrical domain
! ------------------------------------------------------------------------------
        CASE ( "cyldet" )
          flag        = pMixtInput%prepRealVal2
          gaussAmp    = pMixtInput%prepRealVal6
          gaussCenter = pMixtInput%prepRealVal1
          gaussWidth  = pMixtInput%prepRealVal5*gaussCenter
            
          IF (flag == 1) THEN ! Random init
           ALLOCATE(randNum(pGrid%nCellsTot))
           DO icg = 1,pGrid%nCellsTot
             randNum(icg) = icg
           END DO
           CALL RANDOM_SEED()
           CALL RANDOM_NUMBER(randNum)

           DO icg = 1,pGrid%nCellsTot
             x = pGrid%cofg(XCOORD,icg)
             y = pGrid%cofg(YCOORD,icg)
             z = pGrid%cofg(ZCOORD,icg)
             radius  = (x**2.0_RFREAL+y**2.0_RFREAL)**0.5_RFREAL
             perturb = gaussAmp*EXP(-((radius-gaussCenter)**2.0_RFREAL) &
                      /(2.0_RFREAL*gaussWidth**2.0_RFREAL))*randNum(icg)
            
             pCv(CV_MIXT_DENS,icg) = pCv(CV_MIXT_DENS,icg)*(1.0_RFREAL+perturb)
             pCv(CV_MIXT_XMOM,icg) = pCv(CV_MIXT_XMOM,icg)*(1.0_RFREAL+perturb)
             pCv(CV_MIXT_YMOM,icg) = pCv(CV_MIXT_YMOM,icg)*(1.0_RFREAL+perturb)
             pCv(CV_MIXT_ZMOM,icg) = pCv(CV_MIXT_ZMOM,icg)*(1.0_RFREAL+perturb)
             pCv(CV_MIXT_ENER,icg) = pCv(CV_MIXT_ENER,icg)*(1.0_RFREAL+perturb)
           END DO ! icg

          ELSE
           m = pMixtInput%prepIntVal1 ! Wave no in theta
           n = pMixtInput%prepIntVal2 ! Wave number in z
           DO icg = 1,pGrid%nCellsTot
             x = pGrid%cofg(XCOORD,icg)
             y = pGrid%cofg(YCOORD,icg)
             z = pGrid%cofg(ZCOORD,icg)
             IF (ABS(m) == 0) THEN ! m=0; 2D-(r-z) 
              radius  = y
              zz      = x
              the     = DATAN2(z,y)
             ELSE                  ! 2D- (r-theta) OR 3D-(r-theta-z)
              radius  = (x**2.0_RFREAL+y**2.0_RFREAL)**0.5_RFREAL
              zz      = z
              the     = DATAN2(y,x)
             END IF
             perturb = gaussAmp*EXP(-((radius-gaussCenter)**2.0_RFREAL) &
                     /(2.0_RFREAL*gaussWidth**2.0_RFREAL))*DCOS(m*the)*DCOS(n*zz/gaussCenter) ! gaussCenter = r0
            
             pCv(CV_MIXT_DENS,icg) = pCv(CV_MIXT_DENS,icg)*(1.0_RFREAL+perturb)
             pCv(CV_MIXT_XMOM,icg) = pCv(CV_MIXT_XMOM,icg)*(1.0_RFREAL+perturb)
             pCv(CV_MIXT_YMOM,icg) = pCv(CV_MIXT_YMOM,icg)*(1.0_RFREAL+perturb)
             pCv(CV_MIXT_ZMOM,icg) = pCv(CV_MIXT_ZMOM,icg)*(1.0_RFREAL+perturb)
             pCv(CV_MIXT_ENER,icg) = pCv(CV_MIXT_ENER,icg)*(1.0_RFREAL+perturb)
           END DO ! icg
          END IF

! End Subbu - Jul 2012

! ------------------------------------------------------------------------------
!       Cylinder diffracting shock
! ------------------------------------------------------------------------------

        CASE ( "cylds" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN 
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3 
              v = 0.0_RFREAL 
              w = 0.0_RFREAL 
              p = pMixtInput%prepRealVal4 
            
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)
        
              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)            
            
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
            END IF ! x
          END DO ! icg           
  

! ------------------------------------------------------------------------------
!       Contact-Interface
! ------------------------------------------------------------------------------

        CASE ( "ctint" ) !Added by Chris Neal
          
          DO icg = 1,pGrid%nCellsTot
              x = pGrid%cofg(XCOORD,icg) !extract cell centroid x coordinate

              IF ( x < pMixtInput%prepRealVal1 ) THEN
                d = pMixtInput%prepRealVal2 !New density behind interface
                u = pCv(CV_MIXT_XMOM,icg)/pCv(CV_MIXT_DENS,icg)
                v = pCv(CV_MIXT_YMOM,icg)/pCv(CV_MIXT_DENS,icg)
                w = pCv(CV_MIXT_ZMOM,icg)/pCv(CV_MIXT_DENS,icg)


                mw = pGv(GV_MIXT_MOL,indMol*icg)
                cp = pGv(GV_MIXT_CP ,indCp *icg)

                gc = MixtPerf_R_M(mw)
                g  = MixtPerf_G_CpR(cp,gc)

                p = (g-1.0_RFREAL)*(pCv(CV_MIXT_ENER,icg) - 0.5_RFREAL*pCv(CV_MIXT_DENS,icg)*(u**2+v**2+w**2))
                Eo = p/(d*(g-1.0_RFREAL)) +0.5_RFREAL*(u**2+v**2+w**2)

                pCv(CV_MIXT_DENS,icg) = d
                pCv(CV_MIXT_XMOM,icg) = d*u
                pCv(CV_MIXT_YMOM,icg) = d*v
                pCv(CV_MIXT_ZMOM,icg) = d*w
                pCv(CV_MIXT_ENER,icg) = d*Eo

              END IF ! x

            END DO ! icg

! ------------------------------------------------------------------------------
!       Skews diffracting shock
! ------------------------------------------------------------------------------

! ----- Shock Mach number of 3.0 -----------------------------------------------  

        CASE ( "skews_ms2p0","skews_ms3p0","skews_ms4p0" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN 
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3 
              v = 0.0_RFREAL 
              w = 0.0_RFREAL 
              p = pMixtInput%prepRealVal4 

              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)
        
              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)  

              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
            END IF ! x
          END DO ! icg           
 
! ------------------------------------------------------------------------------
!       Shock-particle interaction
! ------------------------------------------------------------------------------

        CASE ( "shktb" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN 
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3 
              v = 0.0_RFREAL 
              w = 0.0_RFREAL 
              p = pMixtInput%prepRealVal4 

              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)
        
              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)  

              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
            END IF ! x
          END DO ! icg           
 


! Subbu - Sep 2011
! ------------------------------------------------------------------------------
!       Detonation wave in 3D
! ------------------------------------------------------------------------------
          CASE ( "sphdet", "detwav3D" )

          ALLOCATE(randNum(pGrid%nCellsTot))
          DO icg = 1,pGrid%nCellsTot
            randNum(icg) = icg
          END DO
          CALL RANDOM_SEED()
          CALL RANDOM_NUMBER(randNum)

          gaussAmp    = pMixtInput%prepRealVal6
          gaussCenter = pMixtInput%prepRealVal1
          gaussWidth  = pMixtInput%prepRealVal5*gaussCenter
          l           = pMixtInput%prepRealVal3 ! Spherical wave number
          m           = pMixtInput%prepRealVal4 ! Mode number

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)
              
            CALL AscLgndr(l,m,DCOS(z),Lmn)
            !radius  = (x**2.0_RFREAL+y**2.0_RFREAL+z**2.0_RFREAL)**0.5_RFREAL
            radius = x  
            !perturb = gaussAmp*EXP(-((radius-gaussCenter)**2.0_RFREAL) &
            !          /(2.0_RFREAL*gaussWidth**2.0_RFREAL))*randNum(icg)

            perturb = gaussAmp*EXP(-((radius-gaussCenter)**2.0_RFREAL) &
                      /(2.0_RFREAL*gaussWidth**2.0_RFREAL))*Lmn*DCOS(y)

            pCv(CV_MIXT_DENS,icg) = pCv(CV_MIXT_DENS,icg)*(1.0_RFREAL+perturb)
            pCv(CV_MIXT_XMOM,icg) = pCv(CV_MIXT_XMOM,icg)*(1.0_RFREAL+perturb)
            pCv(CV_MIXT_YMOM,icg) = pCv(CV_MIXT_YMOM,icg)*(1.0_RFREAL+perturb)
            pCv(CV_MIXT_ZMOM,icg) = pCv(CV_MIXT_ZMOM,icg)*(1.0_RFREAL+perturb)
            pCv(CV_MIXT_ENER,icg) = pCv(CV_MIXT_ENER,icg)*(1.0_RFREAL+perturb)
          END DO ! icg

! End Subbu - Sep 2011

 
! ------------------------------------------------------------------------------
!       Sphere diffracting shock
! ------------------------------------------------------------------------------

        CASE ( "sphds" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN 
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3 
              v = 0.0_RFREAL 
              w = 0.0_RFREAL 
              p = pMixtInput%prepRealVal4 
            
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)
        
              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)            
            
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
            END IF ! x
          END DO ! icg           
  
! ------------------------------------------------------------------------------
!       Default - must be due to input error
! ------------------------------------------------------------------------------

        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
      END SELECT ! global%casename
      
! ==============================================================================
!   Default
! ==============================================================================  
    
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) 
  END SELECT ! pMixtInput%fluidModel

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Initializing flow field from limited hard code done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InitFlowHardCodeLim



SUBROUTINE AscLgndr(l,m,x,plgndr)

 USE ModDataTypes

 IMPLICIT NONE

 INTEGER,INTENT(IN) :: l,m
 REAL(RFREAL),INTENT(IN) ::  x
 REAL(RFREAL),INTENT(OUT) ::  plgndr

 !Computes the associated Legendre polynomial Pm
 !l (x). Here m and l are integers satisfying
 !0 < m < l, while x lies in the range -1 <= x <= 1.

 INTEGER :: i,ll
 REAL(RFREAL) :: fact,pll,pmm,pmmp1,somx2,factDiff,factSum,norm
 REAL(RFREAL) , PARAMETER :: PI = 3.14159265358979

 IF(m.LT.0 .OR. m.GT.l .OR. ABS(x).GT.1.0_RFREAL) THEN !PAUSE 'bad arguments in plgndr'
  WRITE(*,*) 'Error in Associated Legendre computation'
  STOP
 END IF

 pmm=1.0_RFREAL !Compute Pmm
 IF(m .GT. 0) THEN
  somx2=SQRT((1-x)*(1+x))
  fact=1.0_RFREAL
  DO i=1,m
   pmm = -pmm*fact*somx2
   fact=fact+2.0_RFREAL
  END DO
 END IF

 IF(l.EQ.m) THEN
  plgndr=pmm
 ELSE
  pmmp1=x*(2.0_RFREAL*m+1.0_RFREAL)*pmm !Compute Pmm+1.
   IF(l.EQ.m+1) THEN
    plgndr=pmmp1
   ELSE !Compute Pml , l > m+ 1.
    DO ll=m+2,l
     pll=(x*(2.0_RFREAL*ll-1.0_RFREAL)*pmmp1-(ll+m-1.0_RFREAL)*pmm)/(ll-m)
     pmm=pmmp1
     pmmp1=pll
    END DO
    plgndr=pll
   END IF
 END IF

 IF (l-m .EQ. 0) THEN
  factDiff = 1.0_RFREAL
 ELSE
  factDiff = 1.0_RFREAL
  DO i=1,l-m
   factDiff = factDiff * i
  END DO
 END IF

 factSum = 1.0_RFREAL
 DO i=1,l+m
  factSum = factSum * i
 END DO

 norm = SQRT( ((2.0_RFREAL*l + 1.0_RFREAL)*factDiff)/ (4.0_RFREAL*PI*factSum ) )
 plgndr = norm*plgndr

END SUBROUTINE AscLgndr


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_InitFlowHardCodeLim.F90,v $
! Revision 1.4  2016/03/04 16:32:01  rahul
! Added shktb case.
!
! Revision 1.3  2015/02/16 23:00:33  neal
! Removed a useless allocation and reference to pDv data type for ctint initialization.
!
! Revision 1.2  2015/02/11 04:09:41  neal
! Added new initialization case, ctint for contact-interface simulations
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:55  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:07  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:55:41  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/03/26 20:22:23  haselbac
! Removed error trap for GL model
!
! Revision 1.1  2005/11/10 02:54:20  haselbac
! Renamed to shorten name
!
! Revision 1.6  2005/09/13 21:38:47  haselbac
! Removed hardcoded gamma value
!
! Revision 1.5  2005/07/05 19:49:04  mparmar
! Added init for diffracting shock over sphere
!
! Revision 1.4  2005/06/16 20:57:11  haselbac
! Now use gGas for Skews cases instead of hardcoded 1.4
!
! Revision 1.3  2005/06/16 20:55:54  haselbac
! Added new Skews cases
!
! Revision 1.2  2005/04/22 15:22:34  haselbac
! Changed message printed to screen
!
! Revision 1.1  2005/04/15 16:20:22  haselbac
! Initial revision
!
! ******************************************************************************

