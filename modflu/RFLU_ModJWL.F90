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
! Purpose: Collection of routines for JWL equation of state. 
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModJWL.F90,v 1.4 2016/02/08 17:00:26 fred Exp $
!
! Copyright: (c) 2006-2007 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModJWL

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input, &
                        t_mixt
  USE ModError
  USE ModMPI
  USE ModTools, ONLY: IsNan
  USE ModInterfaces, ONLY: MixtPerf_C_DGP, &
                           MixtPerf_Eo_DGPVm, &
                           MixtPerf_P_DEoGVm2, &
                           MixtPerf_T_DPR

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModJWL.F90,v $'
  REAL(RFREAL), PRIVATE :: ATNT,BTNT,CTNT,wTNT,R1TNT,R2TNT,rhoTNT,shcvTNT
    !ATNT      = 371.213E9_RFREAL, & ! Pa
    !BTNT      = 3.2306E9_RFREAL, &  ! Pa
    !CTNT      = 1.04527E9_RFREAL, & ! Pa
    !wTNT      = 0.30_RFREAL, &
    !R1TNT     = 4.15_RFREAL, &
    !R2TNT     = 0.95_RFREAL, &
    !rhoTNT    = 1.63E3_RFREAL, &            ! Kg/m3
    !shcvTNT   = 613.4969325153374_RFREAL, & ! 1e6/rhoTNT, J/(K-kg), Cv'=rhoTNT*Cv=1.0D6 (pa/K)
    !kTNT      = 0.297_RFREAL, &
    !muTNT     = 1.74E-6_RFREAL, &         ! NOTE: take air viscosity temporary
    !PrTNT     = 0.4672491789056_RFREAL, & ! shcvTNT*(wTNT+1.0_RFREAL)*muTNT/kTNT
    !MwTNT     = 24.6789_RFREAL, &
    !mlfN2TNT  = 0.2_RFREAL, & ! 1.5_RFREAL/7.5_RFREAL
    !mlfCOTNT  = 0.466666666666667_RFREAL, & ! 3.5_RFREAL/7.5_RFREAL
    !mlfH2OTNT = 0.333333333333333, & ! 2.5_RFREAL/7.5_RFREAL
    !MwAir     = 28.8501_RFREAL, &
    !mlfN2Air  = 0.79_RFREAL, &
    !mlfO2Air  = 0.21_RFREAL
  
! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: &
      RFLU_JWL_C_ER, &
      RFLU_JWL_ComputeEnergyMixt, &
      RFLU_JWL_ComputePressureMixt, &
      RFLU_JWL_E_PR, &
      RFLU_JWL_P_ER, &
      RFLU_JWL_T_PR

! ==============================================================================
! Private functions
! ==============================================================================
  
  PRIVATE :: RFLU_JWL_FindInverseMatrix

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  




! ******************************************************************************
!
! Purpose: Compute speed of sound from JWL equation of state.
!
! Description: None.
!
! Input:
!   e             Energy
!   r             Density
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  FUNCTION RFLU_JWL_C_ER(pRegion,e,r)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region),POINTER :: pRegion
    REAL(RFREAL), INTENT(IN) :: e,r
    REAL(RFREAL) :: RFLU_JWL_C_ER
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    REAL(RFREAL) :: a,CTNT,v
    
    ATNT   = pRegion%mixtInput%prepRealVal14
    BTNT   = pRegion%mixtInput%prepRealVal15
    wTNT   = pRegion%mixtInput%prepRealVal17
    R1TNT  = pRegion%mixtInput%prepRealVal18
    R2TNT  = pRegion%mixtInput%prepRealVal19
    rhoTNT = pRegion%mixtInput%prepRealVal4*pRegion%mixtInput%prepRealVal3

! ******************************************************************************
!   Start
! ******************************************************************************
    
    v    = 1.0_RFREAL/r
    CTNT = (rhoTNT*e - ATNT/R1TNT*EXP(-R1TNT*rhoTNT*v) &
                     - BTNT/R2TNT*EXP(-R2TNT*rhoTNT*v))*wTNT*(rhoTNT*v)**wTNT
    a    =  SQRT(rhoTNT*v**2.0_RFREAL*( ATNT*R1TNT*EXP(-R1TNT*rhoTNT*v)  &
                      + BTNT*R2TNT*EXP(-R2TNT*rhoTNT*v)  &
                      + CTNT*(1.0_RFREAL+wTNT)/(rhoTNT*v)**(2.0_RFREAL+wTNT)))
                                
    RFLU_JWL_C_ER = a

! ******************************************************************************
!   End  
! ******************************************************************************

  END FUNCTION RFLU_JWL_C_ER








! ******************************************************************************
!
! Purpose: Compute speed of sound, energy, and temperature
!          from JWL equation of state
!          for mixture of JWL and perfect gas.
!
! Description: None.
!
! Input:
!   icg           Cell Number 
!   g             Gamma, ratio of specific heats for mixture
!   gc            R, gas constant for mixture
!   p             Pressure
!   r             Density
!   Y             Mass fraction of explosive (JWL) gas
!
! Output:
!   a             Speed of sound
!   e             Energy
!   T             Temperature
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_JWL_ComputeEnergyMixt(pRegion,icg,g,gc,p,r,Y,a,e,T)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region),POINTER :: pRegion
    INTEGER :: icg
    REAL(RFREAL), INTENT(IN) :: g,gc,p,r,Y
    REAL(RFREAL), INTENT(OUT) :: a,e,T
    
! ==============================================================================
!   Locals
! ==============================================================================
    LOGICAL :: set_e_JWL
    INTEGER :: i,j,iter,iterMax
    DOUBLE PRECISION :: aa,ab,df,ea,eb,ebp,nTol,pa,pb, &
                    ra,rb,s,Ta,Tap,Tb,Tbp,term, &
                    tol,tolok,v,va,vb,Vm,T1,T2,T3,const,difff_n
    DOUBLE PRECISION  :: rJWL_old,rPerf_old
    DOUBLE PRECISION  :: f(2), f_old(2), sol(2), sol_old(2), Ja(2,2),InvJa(2,2), &
                   difff(2),temp_1(2),temp_2(2),&
                   res(2,2), diffsol(2),alt(2)

    ATNT   = pRegion%mixtInput%prepRealVal14
    BTNT   = pRegion%mixtInput%prepRealVal15
    wTNT   = pRegion%mixtInput%prepRealVal17
    R1TNT  = pRegion%mixtInput%prepRealVal18
    R2TNT  = pRegion%mixtInput%prepRealVal19
    rhoTNT = pRegion%mixtInput%prepRealVal4*pRegion%mixtInput%prepRealVal3
    shcvTNT  = pRegion%mixtInput%prepRealVal20

!******************************************************************************
!  Set Pointers
!******************************************************************************
   rJWL_old = pRegion%mixt%JWL_e_prev(icg,1)
   rPerf_old = pRegion%mixt%JWL_e_prev(icg,2)
   set_e_JWL = pRegion%mixt%set_e_JWL(icg) 
! ******************************************************************************
!   Start
! ******************************************************************************
    
    ! Subbu - Reduce tolerance for Broyden method
     nTol = 1.0D-3
    ! Subbu - End reduce tolerance for Broyden method

    Vm = 0.0D0

    v = 1.0D0/r

! ==============================================================================
!   Explosive gas (JWL) 
! ==============================================================================

    IF ( ABS(1.0_RFREAL - Y) < nTol .OR. Y > 1.0D0 ) THEN 
      e = RFLU_JWL_E_PR(pRegion,p,r)
      T = RFLU_JWL_T_PR(pRegion,p,r)
      a = RFLU_JWL_C_ER(pRegion,e,r)

      IF ( e < 0.0_RFREAL ) THEN
        WRITE(*,*) 'Energy is negative in RFLU_JWL_ComputeEnergyMixt JWL'
        WRITE(*,*) p,r,Y
        STOP
      END IF ! e

      RETURN
    END IF ! Y

! ==============================================================================
!   Perfect gas 
! ==============================================================================

    IF ( ABS( Y ) < nTol ) THEN                             
      e = MixtPerf_Eo_DGPVm(r,g,p,Vm) ! Since Vm=0, thus Eo is internal energy
      T = MixtPerf_T_DPR(r,p,gc)
      a = MixtPerf_C_DGP(r,g,p)

      IF ( e < 0.0D0 ) THEN
        WRITE(*,*) 'Energy is negative in RFLU_JWL_ComputeEnergyMixt Perf'
        STOP
      END IF

      RETURN
    END IF ! Y

! TO DO: Manoj-JWL, Clean up following code    
! ==============================================================================
! Explosive gas and perfect gas mixture (Secant Method)
! ==============================================================================
! FRED - Replacing old secant method with a Broyden method solver
! ------------------------------------------------------------------------------
! Initial Guesses
! ------------------------------------------------------------------------------

  IF (rJWL_old == 0.0D0 .AND. rPerf_old == 0.0D0) THEN
  sol(1) = v/(2.0D0*(1.0D0-Y)) !6.213065525944548D0 !v*8.0D0
  sol(2) = v/(2.0D0*Y)
  ELSE
  sol(1) = rJWL_old
  sol(2) = rPerf_old
  END IF  

! ------------------------------------------------------------------------------
! Secant Method Iteration
! ------------------------------------------------------------------------------

  iterMax = 100
  tol     = 1.0D-8
  tolok   = 1.0D-3
  
  va = sol(1)
  vb = sol(2)

  f(1) = va*Y + vb*(1.0D0-Y) - v
  f(2) = RFLU_JWL_T_PR(pRegion,p,1.0D0/va)-MixtPerf_T_DPR(1.0D0/vb,p,gc)

  Ja(1,1) = Y
  Ja(1,2) = 1.0D0-Y

  T1 = ATNT*exp(-R1TNT*rhoTNT*va)
  T2 = BTNT*exp(-R2TNT*rhoTNT*vb)
  T3 = 1.0D0/(wTNT*shcvTNT)

  Ja(2,1) = &
         T3*(p+((R1TNT*rhoTNT*va)-1.0D0)*T1+((R2TNT*rhoTNT*va)-1.0D0)*T2)
  Ja(2,2) = -p/gc

  CALL RFLU_JWL_FindInverseMatrix(Ja,InvJa,2,2)
  
  sol_old = sol

  DO i = 1,2
   alt(i) = 0.0D0
  END DO

  DO i = 1,2
   DO j = 1,2
      alt(i) = alt(i) + InvJa(i,j)*f(j)
   END DO
  END DO

  sol = sol_old - alt

  DO iter = 1,iterMax

   f_old = f
   va = sol(1)
   vb = sol(2)

   f(1) = va*Y + vb*(1.0D0-Y) - v
   f(2) = RFLU_JWL_T_PR(pRegion,p,1.0D0/va)-MixtPerf_T_DPR(1.0D0/vb,p,gc)

   difff = f - f_old
   diffsol = sol - sol_old

   difff_n = (difff(1)*difff(1))+(difff(2)*difff(2))
   
   IF (difff_n .LE. 1.0D-8) THEN
    EXIT
   END IF

!_______________________Step 1 - First Pseudo-inverse operation_____________________________

   DO i = 1,2
    temp_1(i) = 0.0D0
   END DO

   DO i = 1,2
    DO j = 1,2
     temp_1(i) = temp_1(i) + InvJa(i,j)*difff(j)
    END DO
   END DO

   temp_2 = diffsol - temp_1

   DO i = 1,2
    DO j = 1,2
     res(i,j) = temp_2(i)*difff(j)/difff_n
    END DO
   END DO

   !_______________________Step 2 - 2 by 2 Pseudo inverse operation______________________


   InvJa = InvJa + res 

   sol_old = sol

   DO i = 1,2
     alt(i) = 0.0D0
   END DO

   DO i = 1,2
    DO j = 1,2
     alt(i) = alt(i) + InvJa(i,j)*f(j)
    END DO
   END DO

   sol = sol_old - alt

   va = sol(1)
   vb = sol(2)

   s = sqrt((sol(1)-sol_old(1))**2.0D0+(sol(2)-sol_old(2))**2.0D0)

   IF(s < tol) THEN
     EXIT
   END IF

   END DO ! iter

  ra = 1.0D0/va
  rb = 1.0D0/vb
  T = RFLU_JWL_T_PR(pRegion,p,ra)
  Tb = MixtPerf_T_DPR(rb,p,gc)
  
  IF ( (IsNan(T) .EQV. .TRUE.) .OR. (IsNan(Tb) .EQV. .TRUE.)) THEN
    WRITE(*,*) 'Error: Temperature is NaN in energy iteration'
    WRITE(*,*) Y,r,p
    STOP
  END IF

  IF ( T < 0.0D0 .OR. Tb < 0.0D0 ) THEN
    WRITE(*,*) 'Error: negative temperature computed in energy iteration:'
    WRITE(*,*) Y,va,vb,v,T,Tb,p
    STOP
  END IF  
  
  ea = RFLU_JWL_E_PR(pRegion,p,ra)
  aa = RFLU_JWL_C_ER(pRegion,ea,ra)
  eb = MixtPerf_Eo_DGPVm(rb,g,p,Vm) ! Since Vm=0, thus Eo is internal energy
  ab = MixtPerf_C_DGP(rb,g,p)

  e  = ea*Y + eb*(1.0D0-Y)
  a  = aa*Y + ab*(1.0D0-Y)

  IF ( ABS(T-Tb) > ABS(tolok*T) ) THEN
    WRITE(*,*) 'Energy Broyden Method may not be converged! ra,rb,T,Tb,sNorm,Y:'
    WRITE(*,'(6(E15.8,1X))') ra,rb,T,Tb,s,Y
    WRITE(*,'(A,5(1XE15.8))') "Input values g,gc,p,r,Y ",g,gc,p,r,Y
    STOP
  END IF


  IF ( (IsNan(ea) .EQV. .TRUE.) .OR. (IsNan(e) .EQV. .TRUE.)) THEN
    WRITE(*,*) 'Error: Energy is NaN in energy iteration'
    WRITE(*,*) Y,r,p,ea,e
    STOP
  END IF

  IF ( ea < 0.0D0 .OR. e < 0.0D0 ) THEN
    WRITE(*,*) 'Error: negative energy computed in energy iteration:',Y,ra,rb,r,p,e,ea,T
    STOP
  END IF

  IF (set_e_JWL .EQV. .FALSE.) THEN
    set_e_JWL = .TRUE.
    rJWL_old = ra
    rPerf_old = rb
  END IF
    
! ******************************************************************************
!   End  
! ******************************************************************************

  END SUBROUTINE RFLU_JWL_ComputeEnergyMixt




! ******************************************************************************
!
! Purpose: Compute speed of sound, mass fractions of perfect gas and detonation
!          gas (JWL), pressure, and temperature from JWL equation of state
!          for mixture of JWL and perfect gas.
!
! Description: None.
!
! Input:
!   icg           Cell Number
!   g             Gamma, ratio of specific heats for mixture
!   gc            R, gas constant for mixture
!   e             Energy
!   r             Density
!   Y             Mass fraction of explosive (JWL) gas
!
! Output:
!   a             Speed of sound
!   eJWL          Mass fraction for detonation gas (JWL)
!   ePerf         Mass fraction for perfect gas
!   p             Pressure
!   T             Temperature
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_JWL_ComputePressureMixt(pRegion,icg,g,gc,e,r,Y,a,eJWL,ePerf,p,T)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    INTEGER :: icg
    REAL(RFREAL), INTENT(IN) :: g,gc,e,r,Y
    REAL(RFREAL), INTENT(OUT) :: a,eJWL,ePerf,p,T
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    LOGICAL :: set_p_JWL
    INTEGER :: i,iter,iterMax,j,k
    DOUBLE PRECISION :: aa,ab,dFdva,dGdva,ea,eb, &
                    fp,Gva,nTol,pa,pb,ra, &
                    rap,rb,rbp,shcpTNT,sNorm,sNormOld, &
                    sNormq,Ta,Tap,Tb,Tbp,term,tol,tolok,v,va,vap,vb,vbp,Vm2,&
                    difff_n,diffs_n
    DOUBLE PRECISION :: rJWL_old,rPerf_old,eJWL_old,ePerf_old
    DOUBLE PRECISION :: df(4),dJa(4,4),f(4),f_old(4),InvJa(4,4),Ja(4,4),&
                    Jadx(4),s(4),sol(4),sol_old(4),diffsol(4),temp_1(4),&
                    temp_2(4),difff(4),res(4,4),alt(4)
    ATNT   = pRegion%mixtInput%prepRealVal14
    BTNT   = pRegion%mixtInput%prepRealVal15
    wTNT   = pRegion%mixtInput%prepRealVal17
    R1TNT  = pRegion%mixtInput%prepRealVal18
    R2TNT  = pRegion%mixtInput%prepRealVal19
    rhoTNT = pRegion%mixtInput%prepRealVal4*pRegion%mixtInput%prepRealVal3
    shcvTNT  = pRegion%mixtInput%prepRealVal20
   

!**********************Pointers************************************************
    rJWL_old = pRegion%mixt%JWL_p_prev(icg,1)
    rPerf_old = pRegion%mixt%JWL_p_prev(icg,2)
    eJWL_old = pRegion%mixt%JWL_p_prev(icg,3)
    ePerf_old = pRegion%mixt%JWL_p_prev(icg,4)
    set_p_JWL = pRegion%mixt%set_p_JWL(icg)

! ******************************************************************************
!   Start
! ******************************************************************************
   
    ! Subbu - Reduce tolerance Broyden method
    nTol = 1.0D-3
    ! Subbu - End reduce tolerance Broyden method

    Vm2 = 0.0D0
    v = 1.0D0/r

! ==============================================================================
!   Explosive gas (JWL) 
! ==============================================================================

    IF ( ABS(1.0_RFREAL - Y) < nTol .OR. Y > 1.0_RFREAL ) THEN 
      p = RFLU_JWL_P_ER(pRegion,e,r)
      T = RFLU_JWL_T_PR(pRegion,p,r)
      a = RFLU_JWL_C_ER(pRegion,e,r)

      eJWL  = e
      ePerf = 0.0D0

      IF ( p < 0.0_RFREAL ) THEN
        WRITE(*,*) 'Pressure is negative in RFLU_JWL_ComputePressureMixt JWL'
        WRITE(*,*) Y,e,r,p,T
        STOP
      END IF ! e

      IF (T < 0.0_RFREAL) THEN
        WRITE(*,*) 'Temperature is negative in RFLU_JWL_ComputePressureMixt JWL'
        WRITE(*,*) Y,e,r,p,T
        STOP
      END IF !FRED - Adding kill statements for negative temperatures as well

  IF ( (IsNan(T) .EQV. .TRUE.) .OR. (IsNan(Tb) .EQV. .TRUE.)) THEN
    WRITE(*,*) 'Error: Temperature is NaN in JWL pressure'
    WRITE(*,*) Y,r,e
    STOP
  END IF

  IF ( (IsNan(p) .EQV. .TRUE.) .OR. (IsNan(pb) .EQV. .TRUE.)) THEN
    WRITE(*,*) 'Error: Pressure is NaN in JWL pressure'
    WRITE(*,*) Y,r,e
    STOP
  END IF


      RETURN
    END IF ! Y

! ==============================================================================
!   Perfect gas 
! ==============================================================================

    IF ( ABS( Y ) < nTol ) THEN                             
      p = MixtPerf_P_DEoGVm2(r,e,g,Vm2) ! Since Vm2=0, e is also total energy
      T = MixtPerf_T_DPR(r,p,gc)
      a = MixtPerf_C_DGP(r,g,p)

      eJWL  = 0.0_RFREAL
      ePerf = e

      IF ( p < 0.0_RFREAL ) THEN
        WRITE(*,*) 'Pressure is negative in RFLU_JWL_ComputePressureMixt Perf'
        WRITE(*,*) Y,e,r,p,T
        STOP
      END IF ! e

      IF (T < 0.0_RFREAL) THEN
        WRITE(*,*) 'Temperature negative in RFLU_JWL_ComputePressureMixt Perf'
        STOP
      END IF !FRED - Adding kill statements for negative temperatures as well

  IF ( (IsNan(T) .EQV. .TRUE.) .OR. (IsNan(Tb) .EQV. .TRUE.)) THEN
    WRITE(*,*) 'Error: Temperature is NaN in Ideal pressure'
    WRITE(*,*) Y,r,e
    STOP
  END IF

  IF ( (IsNan(p) .EQV. .TRUE.) .OR. (IsNan(pb) .EQV. .TRUE.)) THEN
    WRITE(*,*) 'Error: Pressure is NaN in Ideal pressure'
    WRITE(*,*) Y,r,e
    STOP
  END IF


      RETURN
    END IF ! Y

! TO DO: Manoj-JWL, Clean up following code    
! ==============================================================================
! Explosive gas and perfect gas mixture (Broyden's Method) 
! ==============================================================================

! ------------------------------------------------------------------------------
! Initial Guess
! ------------------------------------------------------------------------------
  
  IF (rJWL_old == 0.0D0 .AND. rPerf_old == 0.0D0 .AND. eJWL_old == 0.0D0 .AND. ePerf_old == 0.0D0) THEN
  sol(1) = v/(2.0D0*(1.0D0-Y)) !6.213065525944548D0 !v*8.0D0
  sol(2) = v/(2.0D0*Y) !0.857860966971079D0 !v
  sol(3) = e !9.323363026977179D+05 !e*4.0D0
  sol(4) = e !2.145518769638580D+05 !e
  ELSE
  sol(1) = rJWL_old
  sol(2) = rPerf_old
  sol(3) = eJWL_old
  sol(4) = ePerf_old
  END IF

  sol_old = sol !FRED - Modified initial guesses for specific volumes

  va = sol_old(1)
  vb = sol_old(2)
  ea = sol_old(3) 
  eb = sol_old(4)

! ------------------------------------------------------------------------------
! Broyden Method Iteration
! ------------------------------------------------------------------------------

  iterMax = 100
  tol     = 1.0D-8
  tolok   = 1.0D-3
  
    va = sol_old(1)
    vb = sol_old(2)
    ea = sol_old(3)
    eb = sol_old(4)

    f(1) = va*Y + vb*(1.0D0-Y) - v
    f(2) = ea*Y + eb*(1.0D0-Y) - e
    ra = 1.0D0/va
    rb = 1.0D0/vb
    pa = RFLU_JWL_P_ER(pRegion,ea,ra)
    Ta = RFLU_JWL_T_PR(pRegion,pa,ra)
    pb = MixtPerf_P_DEoGVm2(rb,eb,g,Vm2) ! Since Vm2=0, e is also total energy
    Tb = MixtPerf_T_DPR(rb,pb,gc)

    f(3) = pa - pb
    f(4) = Ta - Tb

!   Calculate Jacobian matrix approximately
    Ja(1,1)   = Y
    Ja(1,2)   = 1.0D0-Y
    Ja(1,3) = 0.0D0
    Ja(1,4) = 0.0D0

    Ja(2,1) = 0.0D0
    Ja(2,2) = 0.0D0
    Ja(2,3)   = Y
    Ja(2,4)   = 1.0D0-Y

    dFdva = ATNT*(wTNT/va-R1TNT*rhoTNT+wTNT/R1TNT/rhoTNT/va/va) &
           *exp(-R1TNT*rhoTNT*va)                                   &
          + BTNT*(wTNT/va-R2TNT*rhoTNT+wTNT/R2TNT/rhoTNT/va/va) &
           *exp(-R2TNT*rhoTNT*va)
    Ja(3,1)   = -wTNT*ea/va/va + dFdva
    Ja(3,2)   = (g-1.0D0)*eb/(vb*vb)
    Ja(3,3)   =  wTNT/va
    Ja(3,4)   = -(g-1.0D0)/vb

    dGdva = ATNT/shcvTNT*exp(-R1TNT*rhoTNT*va) &
          + BTNT/shcvTNT*exp(-R2TNT*rhoTNT*va)
    Ja(4,1) = dGdva
    Ja(4,2) = 0.0D0
    Ja(4,3) = 1.0D0/shcvTNT
    Ja(4,4) = (1.0D0-g)/gc

    CALL RFLU_JWL_FindInverseMatrix(Ja,InvJa,4,4)

    DO i=1,4
      alt(i)=0.0D0
    END DO

    DO i=1,4
     DO j=1,4
      alt(i) = alt(i) + InvJa(i,j)*f(j)
     END DO
    END DO

    sol = sol_old - alt

  DO iter = 1,iterMax

    f_old = f
    va = sol(1)
    vb = sol(2)
    ea = sol(3)
    eb = sol(4)

    f(1) = va*Y + vb*(1.0D0-Y) - v
    f(2) = ea*Y + eb*(1.0D0-Y) - e
    ra = 1.0D0/va
    rb = 1.0D0/vb
    pa = RFLU_JWL_P_ER(pRegion,ea,ra)
    Ta = RFLU_JWL_T_PR(pRegion,pa,ra)
    pb = MixtPerf_P_DEoGVm2(rb,eb,g,Vm2) ! Since Vm2=0, e is also total energy
    Tb = MixtPerf_T_DPR(rb,pb,gc)

    IF ( (IsNan(T) .EQV. .TRUE.) .OR. (IsNan(Tb) .EQV. .TRUE.)) THEN
    WRITE(*,*) 'Temperature is negative in a Broyden step - pressure'
    EXIT
    END IF


    f(3) = pa - pb
    f(4) = Ta - Tb

    difff = f-f_old
    diffsol = sol-sol_old
!______________________________Step 1 - Pseudo-inverse Calculation_______________


    difff_n = (difff(1)*difff(1))+(difff(2)*difff(2))&
              +(difff(3)*difff(3))+(difff(4)*difff(4))


    IF (difff_n .LE. 1.0D-8) THEN
    va = sol(1)
    vb = sol(2)
    ea = sol(3)
    eb = sol(4)
    EXIT
    ENDIF

    DO i=1,4
      temp_1(i)=0.0D0
    END DO

    DO i=1,4
     DO j=1,4
      temp_1(i) = temp_1(i) + InvJa(i,j)*difff(j)
     END DO
    END DO

    temp_2 = diffsol - temp_1

   !____________________________Step 2 - 4 by 4 Pseudo-Inverse Calculation_______



    DO j = 1,4
      DO k = 1,4
        res(j,k) = temp_2(j)*difff(k)/difff_n
      END DO
    END DO



    InvJa = InvJa + res

    sol_old = sol

    DO i=1,4
      alt(i)=0.0D0
    END DO

    DO i=1,4
     DO j=1,4
      alt(i) = alt(i) + InvJa(i,j)*f(j)
     END DO
    END DO

    sol = sol_old - alt

    va = sol(1)
    vb = sol(2)
    ea = sol(3)
    eb = sol(4)

    sNorm = sqrt((sol(1)-sol_old(1))*(sol(1)-sol_old(1)) + &
         (sol(2)-sol_old(2))*(sol(2)-sol_old(2)) + &
         (sol(3)-sol_old(3))*(sol(3)-sol_old(3)) + &
         (sol(4)-sol_old(4))*(sol(4)-sol_old(4)))

    IF (sNorm < tol .OR. IsNan(ra) .EQV. .TRUE. .OR. IsNan(rb) .EQV. .TRUE.) THEN
      EXIT
    END IF

  END DO ! iter

! ------------------------------------------------------------------------------
! Calculate dependent variables
! ------------------------------------------------------------------------------

  ra = 1.0D0/va
  rb = 1.0D0/vb
  p  = RFLU_JWL_P_ER(pRegion,ea,ra)
  T  = RFLU_JWL_T_PR(pRegion,p,ra)
  aa = RFLU_JWL_C_ER(pRegion,ea,ra)
  pb = MixtPerf_P_DEoGVm2(rb,eb,g,Vm2) ! Since Vm2=0, e is also total energy
  Tb = MixtPerf_T_DPR(rb,pb,gc)
  ab = MixtPerf_C_DGP(rb,g,p)
   
  IF ( (IsNan(T) .EQV. .TRUE.) .OR. (IsNan(Tb) .EQV. .TRUE.)) THEN
    WRITE(*,*) 'Error: Temperature is NaN in pressure iteration'
    WRITE(*,*) Y,r,e
    STOP
  END IF

  IF ( (IsNan(p) .EQV. .TRUE.) .OR. (IsNan(pb) .EQV. .TRUE.)) THEN
    WRITE(*,*) 'Error: Pressure is NaN in pressure iteration'
    WRITE(*,*) Y,r,e
    STOP
  END IF

 IF (p < 1.0D-1 .OR. pb < 1.0D-1) THEN
    WRITE(*,*) 'Computed pressure converged, but is less than 1'
    WRITE(*,*) Y,r,e,p
    STOP
 END IF

  IF (T < 0.0D0 .OR. p < 0.0D0) THEN
  WRITE(*,*) 'Negative quantities computed in Pressure Iteration'
    IF (T < 0.0D0) THEN
       WRITE(*,*) 'Negative Temperature computed:',T,Tb

    ELSE
       WRITE(*,*) 'Negative Pressure computed:',p,pb
    END IF

  STOP
  END IF 

  IF ( ABS(p-pb) > ABS(tolok*p) .OR. ABS(T-Tb) > ABS(tolok*T) ) THEN
    WRITE(*,*) 'Broyden Method may not be converged! p,pb,T,Tb,sNorm,Y:'
    WRITE(*,'(6(E15.8,1X))') p,pb, T,Tb, sNorm,Y
    WRITE(*,'(A,5(1XE15.8))') "Input values g,gc,e,r,Y ",g,gc,e,r,Y
    T = Y*T + (1.0_RFREAL-Y)*Tb
    p = Y*p + (1.0_RFREAL-Y)*pb
    !STOP
  END IF


  a     = aa*Y + ab*(1.0D0 - Y)
  eJWL  = ea
  ePerf = eb

  IF (set_p_JWL .EQV. .FALSE.) THEN
  set_p_JWL = .TRUE.
  rJWL_old = ra
  rPerf_old = rb
  eJWL_old = ea
  ePerf_old = eb
  END IF
    
! ******************************************************************************
!   End  
! ******************************************************************************

  END SUBROUTINE RFLU_JWL_ComputePressureMixt








! ******************************************************************************
!
! Purpose: Compute internal energy from JWL equation of state.
!
! Description: None.
!
! Input:
!   p             Pressure
!   r             Density
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  FUNCTION RFLU_JWL_E_PR(pRegion,p,r)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region),POINTER :: pRegion
    REAL(RFREAL), INTENT(IN) :: p,r
    REAL(RFREAL) :: RFLU_JWL_E_PR
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    REAL(RFREAL) :: e,v
    ATNT   = pRegion%mixtInput%prepRealVal14
    BTNT   = pRegion%mixtInput%prepRealVal15
    wTNT   = pRegion%mixtInput%prepRealVal17
    R1TNT  = pRegion%mixtInput%prepRealVal18
    R2TNT  = pRegion%mixtInput%prepRealVal19
    rhoTNT = pRegion%mixtInput%prepRealVal4*pRegion%mixtInput%prepRealVal3
    
! ******************************************************************************
!   Start
! ******************************************************************************
    
    v = 1.0_RFREAL/r
    e = p*v/wTNT - ATNT*(v/wTNT-1.0_RFREAL/rhoTNT/R1TNT)*EXP(-R1TNT*rhoTNT*v) &
                 - BTNT*(v/wTNT-1.0_RFREAL/rhoTNT/R2TNT)*EXP(-R2TNT*rhoTNT*v) 

    RFLU_JWL_E_PR = e

! ******************************************************************************
!   End  
! ******************************************************************************

  END FUNCTION RFLU_JWL_E_PR








! ###############################################
! TO DO: Manoj-JWL, Clean up following code    
!Subroutine to find the inverse of a square matrix
!Author : Yoshifumi Nozaki

  SUBROUTINE RFLU_JWL_FindInverseMatrix(matrix,inverse,n,nmax)
  implicit none
        
!---Declarations
        INTEGER, INTENT(IN) :: n,nmax
        double precision, INTENT(IN), DIMENSION(n,n) :: matrix  !Input A matrix
        double precision, INTENT(OUT), DIMENSION(nmax,nmax) :: inverse !Inverted matrix
        
        INTEGER :: i, j, k, l
        double precision :: m
        double precision, DIMENSION(n,2*n) :: augmatrix !augmented matrix

    
        !Augment input matrix with an identity matrix
        DO i = 1,n
          DO j = 1,2*n
            IF (j <= n ) THEN
              augmatrix(i,j) = matrix(i,j)
            ELSE IF ((i+n) == j) THEN
              augmatrix(i,j) = 1.0d0
            Else
              augmatrix(i,j) = 0.0d0
            ENDIF
          END DO
        END DO   
		
        !Ensure diagonal elements are non-zero
        DO k = 1,n-1
          DO j = k+1,n
            IF (augmatrix(k,k) == 0) THEN
!			  write(*,*) 'There exists diagonal elements'
               DO i = k+1, n
                 IF (augmatrix(i,k) /= 0) THEN
                   DO  l = 1, 2*n
                     augmatrix(k,l) = augmatrix(k,l)+augmatrix(i,l)
                   END DO
                 ENDIF
               END DO
            ENDIF
          END DO
        END DO   
		
        !Reduce augmented matrix to upper traingular form
        DO k =1, n-1
          DO j = k+1, n   
            m = augmatrix(j,k)/augmatrix(k,k)
			!write(*,*) k, j, m
            DO i = k, 2*n
              augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
            END DO
          END DO
        END DO
        
        !Test for invertibility
        DO i = 1, n
          IF (augmatrix(i,i) == 0) THEN
            write(*,*) "ERROR-Matrix is non-invertible"
            inverse = 0.d0
            WRITE(11,'(4(E23.16,1X))') matrix(1,1:4)
            WRITE(11,'(4(E23.16,1X))') Matrix(2,1:4)
            WRITE(11,'(4(E23.16,1X))') Matrix(3,1:4)
            WRITE(11,'(4(E23.16,1X))') Matrix(4,1:4)
            STOP 
            return
          ENDIF
        END DO
       		
        !Make diagonal elements as 1
        DO i = 1 , n
          m = augmatrix(i,i)
          DO j = i , (2 * n)                                
            augmatrix(i,j) = (augmatrix(i,j) / m)
          END DO
        END DO
		
        !Reduced right side half of augmented matrix to identity matrix
        DO k = n-1, 1, -1
          DO i =1, k
            m = augmatrix(i,k+1)
            DO j = k, (2*n)
              augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
            END DO
          END DO
        END DO        
		
		!store answer
        DO i =1, n
          DO j = 1, n
            inverse(i,j) = augmatrix(i,j+n)
          END DO
        END DO

        !do i=1,n
	    !write(26,99) augmatrix (i,1:2*n)
	    !end do
		
!99  format(100(1X,E30.15)) 

  END SUBROUTINE RFLU_JWL_FindInverseMatrix 
! ###############################################








! ******************************************************************************
!
! Purpose: Compute pressure from JWL equation of state.
!
! Description: None.
!
! Input:
!   e             Energy
!   r             Density
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  FUNCTION RFLU_JWL_P_ER(pRegion,e,r)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    REAL(RFREAL), INTENT(IN) :: e,r
    REAL(RFREAL) :: RFLU_JWL_P_ER
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    REAL(RFREAL) :: p,v
    ATNT   = pRegion%mixtInput%prepRealVal14
    BTNT   = pRegion%mixtInput%prepRealVal15
    wTNT   = pRegion%mixtInput%prepRealVal17
    R1TNT  = pRegion%mixtInput%prepRealVal18
    R2TNT  = pRegion%mixtInput%prepRealVal19
    rhoTNT = pRegion%mixtInput%prepRealVal4*pRegion%mixtInput%prepRealVal3
    
! ******************************************************************************
!   Start
! ******************************************************************************
    
    v = 1.0_RFREAL/r
    p = ATNT*(1.0_RFREAL-wTNT/R1TNT/rhoTNT/v)*EXP(-R1TNT*rhoTNT*v) &
      + BTNT*(1.0_RFREAL-wTNT/R2TNT/rhoTNT/v)*EXP(-R2TNT*rhoTNT*v) + wTNT*e/v

    RFLU_JWL_P_ER = p

! ******************************************************************************
!   End  
! ******************************************************************************

  END FUNCTION RFLU_JWL_P_ER








! ******************************************************************************
!
! Purpose: Compute temperature from JWL equation of state.
!
! Description: None.
!
! Input:
!   p             Pressure
!   r             Density
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  FUNCTION RFLU_JWL_T_PR(pRegion,p,r)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region),POINTER :: pRegion
    REAL(RFREAL), INTENT(IN) :: p,r
    REAL(RFREAL) :: RFLU_JWL_T_PR
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    REAL(RFREAL) :: T,v
    ATNT   = pRegion%mixtInput%prepRealVal14
    BTNT   = pRegion%mixtInput%prepRealVal15
    wTNT   = pRegion%mixtInput%prepRealVal17
    R1TNT  = pRegion%mixtInput%prepRealVal18
    R2TNT  = pRegion%mixtInput%prepRealVal19
    rhoTNT = pRegion%mixtInput%prepRealVal4*pRegion%mixtInput%prepRealVal3
    shcvTNT  = pRegion%mixtInput%prepRealVal20
    
! ******************************************************************************
!   Start
! ******************************************************************************
    
    v = 1.0_RFREAL/r
    T = (p - ATNT*EXP(-R1TNT*rhoTNT*v)  &
           - BTNT*EXP(-R2TNT*rhoTNT*v))*v/wTNT/shcvTNT    

    RFLU_JWL_T_PR = T

! ******************************************************************************
!   End  
! ******************************************************************************

  END FUNCTION RFLU_JWL_T_PR








END MODULE RFLU_ModJWL

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModJWL.F90,v $
! Revision 1.4  2016/02/08 17:00:26  fred
! Fixing Vulcan compiling issue
!
! Revision 1.3  2016/02/06 17:23:43  fred
! Adding error statements in case P,T are Nan in non-iterative parts of functions
!
! Revision 1.2  2016/02/04 19:56:53  fred
! Adding iterative JWL EOS capabilities for the cylindrical detonation case
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
!
! ******************************************************************************

