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
! Purpose: Compute errors of computed solution relative to exact solution.
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
! $Id: RFLU_ComputeExactFlowError.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ComputeExactFlowError(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModParameters
    
  USE RFLU_ModBessel  
  USE RFLU_ModBoundaryTests, ONLY: RFLU_TestIsBoundaryCell
  USE RFLU_ModDifferentiationCells, ONLY: RFLU_ComputeGradCellsWrapper
  USE RFLU_ModExactFlow
  USE RFLU_ModFlowHardCode      
    
  USE ModInterfaces, ONLY: MixtPerf_D_CGP, &  
                           MixtPerf_Eo_DGPUVW, &
                           MixtPerf_G_CPR, & 
                           MixtPerf_P_DRT, & 
                           MixtPerf_R_CpG, &
                           MixtPerf_R_M, &
                           RFLU_PrintLocInfo 
         
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

  LOGICAL :: printErrorNorms
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,iBc,icg,im,in,indCp,indMol,iq,iVar,iVarBeg,iVarEnd, &
             l8normLoc,nCells
  INTEGER, DIMENSION(:), ALLOCATABLE :: varInfo
  INTEGER, DIMENSION(1,MIN_VAL:MAX_VAL) :: locDummy
  REAL(RFREAL) :: aTot,A,A1,A2,A_c,c,const,cp,cpGas,dc,de,dInc,dTot, &
                  dummyReal,etaqm,g,gc,gGas,gmc,gme,gxc,gxe,gyc,gye,gzc,gze, &
                  height,idc,L,L1,l1norm,l2norm,l8norm,L_l,Mi,mInj,Mo,muo,mw, &
                  nx,ny,nz,omega,pc,pe,pi,po,pTot,R,radius,refD,refL,refNu, &
                  refP,refU,rGas,ri,ro,t,term,tTot,uc,ue,uo,uo_c,var,vc,ve, &
                  vInj,we,x,y,z
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt_input), POINTER :: pMixtInput

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = &
    '$RCSfile: RFLU_ComputeExactFlowError.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeExactFlowError', &
                        __FILE__)
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Computing errors in flow solution...'
                             
    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,A)') SOLVER_NAME,'Case: ',TRIM(global%casename)
    END IF ! global%verbLevel                         
  END IF ! global%verbLevel

! ==============================================================================
! Set pointers
! ==============================================================================

  pGrid      => pRegion%grid
  pCv        => pRegion%mixt%cv
  pDv        => pRegion%mixt%dv
  pGv        => pRegion%mixt%gv
  pMixtInput => pRegion%mixtInput

! ==============================================================================
! Set constants and initialize variables
! ==============================================================================

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

  cpGas = global%refCp
  gGas  = global%refGamma  
  rGas  = MixtPerf_R_CpG(cpGas,gGas)
  
  nCells = 0
  
  l1norm =  0.0_RFREAL
  l2norm =  0.0_RFREAL
  l8norm = -HUGE(1.0_RFREAL)  
  
  printErrorNorms = .TRUE.  
  
! ******************************************************************************
! Initialize flow field based on user input
! ******************************************************************************

  SELECT CASE ( global%casename )

! ==============================================================================
!   Acoustic flow case
! ==============================================================================

    CASE ( "acoustic" )
      A1 = pMixtInput%prepRealVal1
      A2 = pMixtInput%prepRealVal2
      Mo = pMixtInput%prepRealVal3
      ro = pMixtInput%prepRealVal4
      po = pMixtInput%prepRealVal5

      DO icg = 1,pGrid%nCellsTot
        nCells = nCells + 1      
                           
        x = pGrid%cofg(XCOORD,icg)
                           
        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          g  = global%refGamma
        ELSE
          mw = pGv(GV_MIXT_MOL,indMol*icg)
          cp = pGv(GV_MIXT_CP ,indCp *icg)
                       
          gc = MixtPerf_R_M(mw)
          g  = MixtPerf_G_CpR(cp,gc)
        END IF ! solverType

        ue = uo
        pe = po
        de = ro
        ve = 0.0_RFREAL
        we = 0.0_RFREAL

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          t = global%currentTime + global%dtImposed/2.0_RFREAL
          CALL RFLU_ComputeExactFlowAcoustic(global,x,t,ro,po,Mo,g, &
                                             A1,A2,de,ue,ve,we,pe)

          dc = pCv(CV_MIXT_DENS,icg)
          uc = pCv(CV_MIXT_XVEL,icg)
        ELSE
          t = global%currentTime 
          CALL RFLU_ComputeExactFlowAcoustic(global,x,t,ro,po,Mo,g, &
                                             A1,A2,de,ue,ve,we,pe)

          dc  = pCv(CV_MIXT_DENS,icg)
          idc = 1.0_RFREAL/dc

          uc = idc*pCv(CV_MIXT_XMOM,icg)
        END IF ! solverType

        term = ABS(SQRT(dc*dc)/SQRT(de*de) - 1.0_RFREAL)
        l1norm = l1norm + term
        l2norm = l2norm + term**2

        IF ( term > l8norm ) THEN
          l8norm    = term
          l8normLoc = icg
        END IF ! term
      END DO ! icg

! ==============================================================================
!   Potential flow around cylinder
! ==============================================================================  
  
    CASE ( "cylpotential" )
      ro  = pMixtInput%prepRealVal1
      uo  = pMixtInput%prepRealVal2
      po  = pMixtInput%prepRealVal3
      R   = pMixtInput%prepRealVal4

      DO icg = 1,pGrid%nCellsTot
        nCells = nCells + 1      
                           
        x = pGrid%cofg(XCOORD,icg)
        y = pGrid%cofg(YCOORD,icg)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          g  = global%refGamma
        ELSE
          mw = pGv(GV_MIXT_MOL,indMol*icg)
          cp = pGv(GV_MIXT_CP ,indCp *icg)
                       
          gc = MixtPerf_R_M(mw)
          g  = MixtPerf_G_CpR(cp,gc)
        END IF ! solverType

        CALL RFLU_ComputeExactFlowCylPotential(global,x,y,ro,uo,po,R, &
                                               de,ue,ve,we,pe)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          dc = pCv(CV_MIXT_DENS,icg)
          uc = pCv(CV_MIXT_XVEL,icg)
          pc = pCv(CV_MIXT_PRES,icg)
        ELSE
          dc  = pCv(CV_MIXT_DENS,icg)
          idc = 1.0_RFREAL/dc

          uc = idc*pCv(CV_MIXT_XMOM,icg)
          pc = pDv(DV_MIXT_PRES,icg)
        END IF ! solverType

        term = ABS(SQRT(pc*pc)/SQRT(pe*pe) - 1.0_RFREAL)
        term = ABS(SQRT(uc*uc)/SQRT(ue*ue) - 1.0_RFREAL)
        l1norm = l1norm + term
        l2norm = l2norm + term**2

        IF ( term > l8norm ) THEN
          l8norm    = term
          l8normLoc = icg
        END IF ! term
      END DO ! icg
 
! ==============================================================================
!   Free vortex flow between two concentric cylinders
! ==============================================================================  
  
    CASE ( "freevortex" )
      ro  = pMixtInput%prepRealVal1
      uo  = pMixtInput%prepRealVal2
      po  = pMixtInput%prepRealVal3
      R   = pMixtInput%prepRealVal4

      DO icg = 1,pGrid%nCellsTot
        nCells = nCells + 1      

        x = pGrid%cofg(XCOORD,icg)
        y = pGrid%cofg(YCOORD,icg)
                           
        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          g  = global%refGamma
        ELSE
          mw = pGv(GV_MIXT_MOL,indMol*icg)
          cp = pGv(GV_MIXT_CP ,indCp *icg)
                       
          gc = MixtPerf_R_M(mw)
          g  = MixtPerf_G_CpR(cp,gc)
        END IF ! solverType

        CALL RFLU_ComputeExactFlowFreeVortex(global,x,y,ro,uo,po,R,g, &
                                             de,ue,ve,we,pe)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          dc = pCv(CV_MIXT_DENS,icg)
          uc = pCv(CV_MIXT_XVEL,icg)
          pc = pCv(CV_MIXT_PRES,icg)
        ELSE
          dc  = pCv(CV_MIXT_DENS,icg)
          idc = 1.0_RFREAL/dc

          uc = idc*pCv(CV_MIXT_XMOM,icg)
          pc = pDv(DV_MIXT_PRES,icg)
        END IF ! solverType

        term = ABS(SQRT(pc*pc)/SQRT(pe*pe) - 1.0_RFREAL)
        term = ABS(SQRT(uc*uc)/SQRT(ue*ue) - 1.0_RFREAL)
        l1norm = l1norm + term
        l2norm = l2norm + term**2

        IF ( term > l8norm ) THEN
          l8norm    = term
          l8normLoc = icg
        END IF ! term
      END DO ! icg
 
! ==============================================================================
!   Gradient test
! ==============================================================================  
  
    CASE ( "gtlin","gttri" )
    
! ------------------------------------------------------------------------------
!     Set variables and allocate memory
! ------------------------------------------------------------------------------    
    
      iVarBeg = CV_MIXT_DENS ! NOTE at present works only for 1 var
      iVarEnd = iVarBeg       

      ALLOCATE(grad(XCOORD:ZCOORD,iVarBeg:iVarEnd,pGrid%nCellsTot), &
               STAT=errorFlag) 
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'grad')
      END IF ! global%error  

! ------------------------------------------------------------------------------
!     Compute gradients
! ------------------------------------------------------------------------------    
    
      ALLOCATE(varInfo(iVarBeg:iVarEnd),STAT=errorFlag) 
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varInfo')
      END IF ! global%error    
    
      varInfo(iVarBeg:iVarEnd) = iVarBeg ! NOTE at present works only for 1 var
    
      CALL RFLU_ComputeGradCellsWrapper(pRegion,iVarBeg,iVarEnd,iVarBeg, &
                                        iVarEnd,varInfo,pCv,grad)

      DEALLOCATE(varInfo,STAT=errorFlag) 
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varInfo')
      END IF ! global%error   
      
! ------------------------------------------------------------------------------
!     Compute error norms
! ------------------------------------------------------------------------------    
        
      SELECT CASE ( global%casename) 
      
! ----- Linear function --------------------------------------------------------      
      
        CASE ( "gtlin" )  
          SELECT CASE ( pRegion%mixtInput%dimens ) 
            CASE ( 1 )   
              DO icg = 1,pGrid%nCellsTot
                nCells = nCells + 1      

                x = pGrid%cofg(XCOORD,icg)
                y = pGrid%cofg(YCOORD,icg)
                z = pGrid%cofg(ZCOORD,icg)

                DO iVar = iVarBeg,iVarEnd
                  CALL RFLU_SetExactFlowLinear(x,y,z,iVar,var,gxe,gye,gze)

                  gxc = grad(XCOORD,iVar,icg)

                  gmc = ABS(gxc)
                  gme = ABS(gxe)

                  term = ABS(gmc/gme-1.0_RFREAL)

                  l1norm = l1norm + term 
                  l2norm = l2norm + term**2

                  IF ( term > l8norm ) THEN 
                    l8norm    = term
                    l8normLoc = icg
                  END IF ! term    
                END DO ! iVar 
              END DO ! icg 
            CASE ( 2 )   
              DO icg = 1,pGrid%nCellsTot
                nCells = nCells + 1      

                x = pGrid%cofg(XCOORD,icg)
                y = pGrid%cofg(YCOORD,icg)
                z = pGrid%cofg(ZCOORD,icg)

                DO iVar = iVarBeg,iVarEnd
                  CALL RFLU_SetExactFlowLinear(x,y,z,iVar,var,gxe,gye,gze)

                  gxc = grad(XCOORD,iVar,icg)
                  gyc = grad(YCOORD,iVar,icg)                    

                  gmc = SQRT(gxc*gxc + gyc*gyc)
                  gme = SQRT(gxe*gxe + gye*gye)          

                  term = ABS(gmc/gme-1.0_RFREAL)

                  l1norm = l1norm + term 
                  l2norm = l2norm + term**2

                  IF ( term > l8norm ) THEN 
                    l8norm    = term
                    l8normLoc = icg
                  END IF ! term    
                END DO ! iVar 
              END DO ! icg 
            CASE ( 3 ) 
              DO icg = 1,pGrid%nCellsTot
                nCells = nCells + 1      

                x = pGrid%cofg(XCOORD,icg)
                y = pGrid%cofg(YCOORD,icg)
                z = pGrid%cofg(ZCOORD,icg)

                DO iVar = iVarBeg,iVarEnd
                  CALL RFLU_SetExactFlowLinear(x,y,z,iVar,var,gxe,gye,gze)

                  gxc = grad(XCOORD,iVar,icg)
                  gyc = grad(YCOORD,iVar,icg)
                  gzc = grad(ZCOORD,iVar,icg)                    

                  gmc = SQRT(gxc*gxc + gyc*gyc + gzc*gzc)
                  gme = SQRT(gxe*gxe + gye*gye + gze*gze)          

                  term = ABS(gmc/gme-1.0_RFREAL)

                  l1norm = l1norm + term 
                  l2norm = l2norm + term**2

                  IF ( term > l8norm ) THEN 
                    l8norm    = term
                    l8normLoc = icg
                  END IF ! term    
                END DO ! iVar 
              END DO ! icg         
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pRegion%mixtInput%dimens    

! ----- Trigonometric function -------------------------------------------------      

        CASE ( "gttri" ) 
          nx = pRegion%mixtInput%prepRealVal1
          ny = pRegion%mixtInput%prepRealVal2
          nz = pRegion%mixtInput%prepRealVal3            
        
          SELECT CASE ( pRegion%mixtInput%dimens )
            CASE ( 1 )   
              DO icg = 1,pGrid%nCellsTot
                nCells = nCells + 1      

                x = pGrid%cofg(XCOORD,icg)
                y = pGrid%cofg(YCOORD,icg)
                z = pGrid%cofg(ZCOORD,icg)

                DO iVar = iVarBeg,iVarEnd
                  CALL RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,iVar,var, &
                                             gxe,gye,gze)

                  gxc = grad(XCOORD,iVar,icg)
                  gyc = grad(YCOORD,iVar,icg)                    

                  gmc = ABS(gxc)
                  gme = ABS(gxe)

                  term = ABS(gmc/gme-1.0_RFREAL)

                  l1norm = l1norm + term 
                  l2norm = l2norm + term**2

                  IF ( term > l8norm ) THEN 
                    l8norm    = term
                    l8normLoc = icg
                  END IF ! term    
                END DO ! iVar 
              END DO ! icg 
            CASE ( 2 )   
              DO icg = 1,pGrid%nCellsTot
                nCells = nCells + 1      

                x = pGrid%cofg(XCOORD,icg)
                y = pGrid%cofg(YCOORD,icg)
                z = pGrid%cofg(ZCOORD,icg)

                DO iVar = iVarBeg,iVarEnd
                  CALL RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,iVar,var, &
                                             gxe,gye,gze)

                  gxc = grad(XCOORD,iVar,icg)
                  gyc = grad(YCOORD,iVar,icg)                    

                  gmc = SQRT(gxc*gxc + gyc*gyc)
                  gme = SQRT(gxe*gxe + gye*gye)          

                  term = ABS(gmc/gme-1.0_RFREAL)

                  l1norm = l1norm + term 
                  l2norm = l2norm + term**2

                  IF ( term > l8norm ) THEN 
                    l8norm    = term
                    l8normLoc = icg
                  END IF ! term    
                END DO ! iVar 
              END DO ! icg 
            CASE ( 3 ) 
              DO icg = 1,pGrid%nCellsTot
                nCells = nCells + 1      

                x = pGrid%cofg(XCOORD,icg)
                y = pGrid%cofg(YCOORD,icg)
                z = pGrid%cofg(ZCOORD,icg)

                DO iVar = iVarBeg,iVarEnd
                  CALL RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,iVar,var, &
                                             gxe,gye,gze)

                  gxc = grad(XCOORD,iVar,icg)
                  gyc = grad(YCOORD,iVar,icg)
                  gzc = grad(ZCOORD,iVar,icg)                    

                  gmc = SQRT(gxc*gxc + gyc*gyc + gzc*gzc)
                  gme = SQRT(gxe*gxe + gye*gye + gze*gze)          

                  term = ABS(gmc/gme-1.0_RFREAL)

                  l1norm = l1norm + term 
                  l2norm = l2norm + term**2

                  IF ( term > l8norm ) THEN 
                    l8norm    = term
                    l8normLoc = icg
                  END IF ! term    
                END DO ! iVar 
              END DO ! icg         
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pRegion%mixtInput%dimens            
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! global%casename
  
! ------------------------------------------------------------------------------
!     Deallocate memory
! ------------------------------------------------------------------------------    
    
      DEALLOCATE(grad,STAT=errorFlag) 
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'grad')
      END IF ! global%error  
   
! ==============================================================================
!   Gaussian pulse flow case
! ==============================================================================

    CASE ( "gaussianpulse" )
      A_c  = 0.00001_RFREAL 
! TEMPORARY : Have to be consistent with initialization code
!      uo_c = 0.1_RFREAL 
      uo_c = 0.0_RFREAL
  
      L1 = pMixtInput%prepRealVal1
      c  = pMixtInput%prepRealVal2
      A  = A_c*c
      uo = uo_c*c
  
      ro = pMixtInput%prepRealVal3
      po = pMixtInput%prepRealVal4

      DO icg = 1,pGrid%nCellsTot
        nCells = nCells + 1      
                           
        x = pGrid%cofg(XCOORD,icg)
                           
        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          g  = global%refGamma
        ELSE
          mw = pGv(GV_MIXT_MOL,indMol*icg)
          cp = pGv(GV_MIXT_CP ,indCp *icg)
                       
          gc = MixtPerf_R_M(mw)
          g  = MixtPerf_G_CpR(cp,gc)
        END IF ! solverType

        ue = uo
        pe = po
        de = ro
        ve = 0.0_RFREAL
        we = 0.0_RFREAL

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN

          x = x - (uo+c)*global%dtImposed/2.0_RFREAL

          IF ( x > -L1/2.0_RFREAL .AND. x < L1/2.0_RFREAL ) THEN
            x = x + L1/2.0_RFREAL
            CALL RFLU_ComputeExactFlowGaussianPulse(global,x,ro,uo,po,c,g,L1, &
                                                    A,de,ue,ve,we,pe)
            x = x - L1/2.0_RFREAL
          END IF ! x

          dc = pCv(CV_MIXT_DENS,icg)
          uc = pCv(CV_MIXT_XVEL,icg)

        ELSE
          IF ( x > -L1/2.0_RFREAL .AND. x < L1/2.0_RFREAL ) THEN
            x = x + L1/2.0_RFREAL
            CALL RFLU_ComputeExactFlowGaussianPulse(global,x,ro,uo,po,c,g,L1, &
                                                    A,de,ue,ve,we,pe)
            x = x - L1/2.0_RFREAL
          END IF ! x

          dc  = pCv(CV_MIXT_DENS,icg)
          idc = 1.0_RFREAL/dc

          uc = idc*pCv(CV_MIXT_XMOM,icg)
        END IF ! solverType

        term = ABS(SQRT(dc*dc)/SQRT(de*de) - 1.0_RFREAL)
        l1norm = l1norm + term
        l2norm = l2norm + term**2

        IF ( term > l8norm ) THEN
          l8norm    = term
          l8normLoc = icg
        END IF ! term
      END DO ! icg
 
! ==============================================================================
!   Rayleigh problem case
! ==============================================================================

    CASE ( "rayleighproblem" )
      ro  = pMixtInput%prepRealVal1
      uo  = pMixtInput%prepRealVal2
      po  = pMixtInput%prepRealVal3
      muo = global%refVisc
      t   = global%currentTime

      DO icg = 1,pGrid%nCells
        nCells = nCells + 1      
                           
        y = pGrid%cofg(YCOORD,icg)
                           
        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          g  = global%refGamma
        ELSE
          mw = pGv(GV_MIXT_MOL,indMol*icg)
          cp = pGv(GV_MIXT_CP ,indCp *icg)
                       
          gc = MixtPerf_R_M(mw)
          g  = MixtPerf_G_CpR(cp,gc)
        END IF ! solverType

        CALL RFLU_ComputeExactFlowRayleighProblem(global,t,y,ro,uo, &
                                                  po,muo,de,ue,ve,we,pe)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          uc = pCv(CV_MIXT_XVEL,icg)
        ELSE
          dc  = pCv(CV_MIXT_DENS,icg)
          idc = 1.0_RFREAL/dc

          uc = idc*pCv(CV_MIXT_XMOM,icg)
        END IF ! solverType

        IF ( ue == 0.0_RFREAL ) THEN
          term = ABS(SQRT(uc*uc))
        ELSE
          term = ABS(SQRT(uc*uc)/SQRT(ue*ue) - 1.0_RFREAL)
        END IF ! ue == 0
        l1norm = l1norm + term
        l2norm = l2norm + term**2

        IF ( term > l8norm ) THEN
          l8norm    = term
          l8normLoc = icg
        END IF ! term
      END DO ! icg
 
! ==============================================================================
!   Proudman-Culick flow. NOTE this problem is two-dimensional and assumed to 
!   lie in the x-y plane, and that the injection boundary is located at 
!   y = -height.
! ==============================================================================  
  
    CASE ( "onera_c0", "onera_c0_2d_100x50" )
      CALL RFLU_GetParamsHardCodeProudman(dInc,mInj,vInj,pTot)

      height = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVertTot))
    
      DO icg = 1,pGrid%nCellsTot
        nCells = nCells + 1      
      
        x = pGrid%cofg(XCOORD,icg)
        y = pGrid%cofg(YCOORD,icg)

        dc  = pCv(CV_MIXT_DENS,icg)
        idc = 1.0_RFREAL/dc
        
        uc = idc*pCv(CV_MIXT_XMOM,icg)
        vc = idc*pCv(CV_MIXT_YMOM,icg)  

        CALL RFLU_ComputeExactFlowProudman(global,x,y,height,dInc,vInj, &
                                           pTot,de,ue,ve,we,pe) 

        term = ABS(SQRT(uc*uc + vc*vc)/SQRT(ue*ue + ve*ve) - 1.0_RFREAL)              
        l1norm = l1norm + term 
        l2norm = l2norm + term**2

        IF ( term > l8norm ) THEN 
          l8norm    = term
          l8normLoc = icg
        END IF ! term  
      END DO ! icg     
  
! ==============================================================================
!   Pipe acoustics. NOTE the pipe is assumed to have the x-coordinate 
!   running down the axis. 
! ==============================================================================

    CASE ( "pipeacoust" )
      CALL RFLU_GetParamsHardCodePAcoust(pTot,aTot)

      dTot = MixtPerf_D_CGP(aTot,gGas,pTot)        

      L  = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVert)) 
      ro = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))        

      im  = MAX(pMixtInput%prepIntVal1,1)
      in  = MAX(pMixtInput%prepIntVal2,1)
      iq  = MAX(pMixtInput%prepIntVal3,1)
      iBc = MAX(MIN(pMixtInput%prepIntVal4,1),0)
                  
      const = MAX(pMixtInput%prepRealVal1,0.0_RFREAL)

      CALL RFLU_JYZOM(im,iq,dummyReal,etaqm,dummyReal,dummyReal)       

      omega = aTot*SQRT((in*global%pi/L)**2 + (etaqm/ro)**2)         

      IF ( global%verbLevel > VERBOSE_LOW ) THEN           
        WRITE(STDOUT,'(A,5X,A,1X,I2)'   ) SOLVER_NAME, &
              'Boundary condition:',iBc          
        WRITE(STDOUT,'(A,5X,A,3(1X,I2))') SOLVER_NAME, &
              'Mode:',im,in,iq                   
        WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
              'Total density (kg/m^3):   ',dTot
        WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
              'Total pressure (N/m^2):   ',pTot            
        WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
              'Angular frequency (rad/s):',omega
        WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
              'Constant (-):             ',const                  
      END IF ! global%verbLevel 
          
      DO icg = 1,pGrid%nCellsTot
        nCells = nCells + 1         
      
        x = pGrid%cofg(XCOORD,icg)
        y = pGrid%cofg(YCOORD,icg)
        z = pGrid%cofg(ZCOORD,icg)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          pc = pCv(CV_MIXT_PRES,icg)
        ELSE
          pc = pDv(DV_MIXT_PRES,icg)
        END IF ! global%solverType 

        CALL RFLU_ComputeExactFlowPAcoust(global,z,y,x,global%currentTime, & 
                                          L,ro,iBc,im,in,iq,etaqm,omega,dTot, &
                                          pTot,aTot,const,de,ue,ve,we,pe) 

        term = ABS(pc/pe - 1.0_RFREAL)              
        l1norm = l1norm + term 
        l2norm = l2norm + term**2

        IF ( term > l8norm ) THEN 
          l8norm    = term
          l8normLoc = icg
        END IF ! term               
      END DO ! icg    
    
! ==============================================================================
!   Ringleb flow. NOTE this problem is two-dimensional and assumed to lie in 
!   the x-y plane and that the exact solution is restricted to gamma = 1.4.
! ==============================================================================  
  
    CASE ( "ringleb" )
      CALL RFLU_GetParamsHardCodeRingleb(pTot,tTot)
    
      DO icg = 1,pGrid%nCellsTot
        nCells = nCells + 1
      
        x = pGrid%cofg(XCOORD,icg)
        y = pGrid%cofg(YCOORD,icg)

        dc = pCv(CV_MIXT_DENS,icg)

        CALL RFLU_ComputeExactFlowRingleb(x,y,rGas,pTot,tTot,de,ue,ve,we,pe)  

        term = ABS(dc/de-1.0_RFREAL)
        l1norm = l1norm + term 
        l2norm = l2norm + term**2

        IF ( term > l8norm ) THEN 
          l8norm    = term
          l8normLoc = icg
        END IF ! term                
      END DO ! icg    
    
! ==============================================================================
!   Supersonic vortex flow. NOTE this problem is two-dimensional and assumed 
!   to lie in the x-y plane. 
! ==============================================================================  
  
    CASE ( "ssvorth20x5l1"   ,"ssvortp20x5l1"   , &
           "ssvorth20x5l3"   ,"ssvortp20x5l3"   , &
           "ssvorth40x10l1"  ,"ssvortp40x10l1"  , & 
           "ssvorth40x10l3"  ,"ssvortp40x10l3"  , & 
           "ssvorth80x20l1"  ,"ssvortp80x20l1"  , & 
           "ssvorth80x20l3"  ,"ssvortp80x20l3"  , & 
           "ssvorth160x40l1" ,"ssvortp160x40l1" , &
           "ssvorth160x40l3" ,"ssvortp160x40l3" , &
           "ssvorth320x80l1" ,"ssvortp320x80l1" , &
           "ssvorth320x80l3" ,"ssvortp320x80l3" , &
           "ssvorth640x160l1","ssvortp640x160l1", &
           "ssvorth640x160l3","ssvortp640x160l3" )
      CALL RFLU_GetParamsHardCodeSsVortex(ri,Mi,pTot,tTot)     
      
      DO icg = 1,pGrid%nCellsTot
      
! TEMPORARY
!        IF ( RFLU_TestIsBoundaryCell(pRegion,icg) .EQV. .FALSE. ) THEN
! END TEMPORARY             
      
        nCells = nCells + 1

        x = pGrid%cofg(XCOORD,icg)
        y = pGrid%cofg(YCOORD,icg)

        dc = pCv(CV_MIXT_DENS,icg)

        CALL RFLU_ComputeExactFlowSsVortex(x,y,gGas,rGas,ri,Mi,pTot,tTot, & 
                                           de,ue,ve,we,pe)
   
        term = ABS(dc/de-1.0_RFREAL)
        l1norm = l1norm + term 
        l2norm = l2norm + term**2

        IF ( term > l8norm ) THEN 
          l8norm    = term
          l8normLoc = icg
        END IF ! term 

! TEMPORARY
!        END IF ! RFLU_TestIsBoundaryCell
! END TEMPORARY

      END DO ! icg
        
! ==============================================================================
!   Supersonic vortex flow. NOTE this problem is two-dimensional and assumed 
!   to lie in the x-y plane. 
! ==============================================================================  
  
    CASE ( "taylorvortex" ) 
      refL  = global%refLength
      refNu = global%refVisc/global%refDensity
      refU  = global%refVelocity   
      refD  = global%refDensity
      refP  = global%refPressure
 
      pi = global%pi
      L  = pMixtInput%prepRealVal1
 
      DO icg = 1,pGrid%nCellsTot
        nCells = nCells + 1      
                           
        x = pGrid%cofg(XCOORD,icg)
        y = pGrid%cofg(YCOORD,icg)

        IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
          t = global%currentTime + global%dtImposed/2.0_RFREAL
          CALL RFLU_ComputeExactFlowTaylorVortex(t,pi,x,y,L,refL,refNu,refU, &
                                                 refD,refP,ue,ve,we,pe)
          pc = pCv(CV_MIXT_PRES,icg)
        ELSE
          t = global%currentTime
          CALL RFLU_ComputeExactFlowTaylorVortex(t,pi,x,y,L,refL,refNu,refU, &
                                                 refD,refP,ue,ve,we,pe)
          pc = pDv(DV_MIXT_PRES,icg)
        END IF ! solverType

        term = ABS(SQRT(pc*pc)/SQRT(pe*pe) - 1.0_RFREAL)
        l1norm = l1norm + term
        l2norm = l2norm + term**2

        IF ( term > l8norm ) THEN
          l8norm    = term
          l8normLoc = icg
        END IF ! term
      END DO ! icg
 
! ==============================================================================
!   Default - due to input error or missing CALL in this routine
! ==============================================================================  
            
    CASE DEFAULT 
      printErrorNorms = .FALSE. 
    
      global%warnCounter = global%warnCounter + 1
    
      IF ( global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'*** WARNING ***', & 
              'Exact solution not available. Returning to calling procedure.'                                   
      END IF ! global%verbLevel
  END SELECT ! global%casename

! ******************************************************************************
! Print error norms and location
! ******************************************************************************

  IF ( printErrorNorms .EQV. .TRUE. ) THEN 
    l1norm = l1norm/REAL(nCells)
    l2norm = SQRT(l2norm/REAL(nCells))

    WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME,'Number of cells:',nCells 

    WRITE(STDOUT,'(A,3X,A)')                SOLVER_NAME,'Error norms:'
    WRITE(STDOUT,'(A,5X,A,1X,E13.6)')       SOLVER_NAME,'L1 norm:',l1norm
    WRITE(STDOUT,'(A,5X,A,1X,E13.6)')       SOLVER_NAME,'L2 norm:',l2norm
    WRITE(STDOUT,'(A,5X,A,1X,E13.6,1X,I6)') SOLVER_NAME,'L8 norm:',l8norm, &
                                            l8normLoc

    locDummy(1,MIN_VAL) = l8normLoc
    locDummy(1,MAX_VAL) = l8normLoc
  
    CALL RFLU_PrintLocInfo(pRegion,locDummy,1,LOCINFO_MODE_VERBOSE, & 
                           OUTPUT_MODE_MASTER_ONLY)
  END IF ! printErrorNorms

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Computing errors in flow solution done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ComputeExactFlowError

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeExactFlowError.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.8  2010/03/15 01:03:19  mparmar
! Added cylpotential and freevortex cases
!
! Revision 1.7  2009/09/28 14:22:05  mparmar
! Added error computation for rayleighproblem
!
! Revision 1.6  2009/07/08 19:12:26  mparmar
! Minor changes in logic for exact solution computation for gaussianpulse case
!
! Revision 1.5  2008/12/06 08:43:57  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:11  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/12/03 18:33:43  mparmar
! Changed computation of exact solution of gaussian pulse to get sinusoidal pulse
!
! Revision 1.2  2007/11/28 23:05:45  mparmar
! Added acoustic, gaussianpulse and taylorvortex
!
! Revision 1.1  2007/04/09 18:58:08  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.17  2007/02/27 13:19:53  haselbac
! Enabled 1d computations for gradient checking cases
!
! Revision 1.16  2006/08/04 03:07:35  haselbac
! Adapted to changes, now use module RFLU_ModBoundaryTests
!
! Revision 1.15  2006/01/06 22:19:43  haselbac
! Added comp of error for linear and trig cases
!
! Revision 1.14  2005/10/09 15:12:17  haselbac
! Added 2d C0 case
!
! Revision 1.13  2005/04/29 12:52:17  haselbac
! Adapted to changes in interface of RFLU_ComputeExactFlowPAcoust
!
! Revision 1.12  2005/04/20 14:48:10  haselbac
! Adapted pipeacoust section to changes in init
!
! Revision 1.11  2005/03/23 01:53:15  haselbac
! Added setting of modes for pipeacoust case
!
! Revision 1.10  2005/03/15 20:48:04  haselbac
! Added error computation for pipe acoustics
!
! Revision 1.9  2005/03/09 15:09:05  haselbac
! Added 1-cell wide ssvortex cases
!
! Revision 1.8  2004/10/19 19:30:11  haselbac
! Location info now always written
!
! Revision 1.7  2004/07/06 15:15:07  haselbac
! Adapted to changes in libflu and modflu, cosmetics
!
! Revision 1.6  2004/03/08 22:50:14  haselbac
! Bug fix: Adding missing interfaces
!
! Revision 1.5  2004/03/03 01:53:35  haselbac
! Added printing of location info for cell with largest error
!
! Revision 1.4  2004/02/23 23:05:42  haselbac
! Added Proudman solution for ONERA C0 case
!
! Revision 1.3  2004/02/13 03:04:35  haselbac
! Added more casenames, nCells var to allow interior error comp
!
! Revision 1.2  2004/02/02 01:12:57  haselbac
! Changed initialization
!
! Revision 1.1  2004/01/29 22:58:38  haselbac
! Initial revision
!
! ******************************************************************************

