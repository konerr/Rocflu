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
! Purpose: Collection of routines to help with imposition of constraints.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModConstraintUtils.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModConstraintUtils

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModGrid, ONLY: t_grid
  USE ModMPI
  USE ModTools

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_GetConstrType, & 
            RFLU_GetConstrValue
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModConstraintUtils.F90,v $ $Revision: 1.1.1.1 $' 
              
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS





! ******************************************************************************
!
! Purpose: Get constraint type and value for given variable and patch.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   pPatch      Pointer to patch data
!   var         Variable index
!   ifl         Face index
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

INTEGER FUNCTION RFLU_GetConstrType(pRegion,pPatch,var,ifl)

  USE RFLU_ModGridSpeedUtils, ONLY: RFLU_DescaleGridSpeed

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: ifl,var
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: dist,indCp,indMol,gasModel
  REAL(RFREAL) :: eqTol,minj,nx,ny,nz,rl
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_GetConstrType',__FILE__)

  eqTol = 1.0E-6_RFREAL

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  SELECT CASE ( var ) 
  
! ==============================================================================  
!   Mixture
! ==============================================================================  

! ------------------------------------------------------------------------------
!   Density
! ------------------------------------------------------------------------------

    CASE ( V_MIXT_DENS ) 
      SELECT CASE ( pPatch%bcType )             
        CASE ( BC_INJECTION )
          dist = pPatch%mixt%distrib

          minj = pPatch%mixt%vals(BCDAT_INJECT_MFRATE,dist*ifl)           

          IF ( minj > 0.0_RFREAL ) THEN
            RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET 
          ELSE 
            RFLU_GetConstrType = CONSTR_TYPE_NONE
          END IF ! minj                   
        CASE DEFAULT ! defensive coding
          RFLU_GetConstrType = CONSTR_TYPE_NONE
      END SELECT ! pPatch%bcType   

! ------------------------------------------------------------------------------
!   x-velocity
! ------------------------------------------------------------------------------
  
    CASE ( V_MIXT_XVEL ) 
      SELECT CASE ( pPatch%bcType )
        CASE ( BC_SLIPWALL )
          nx = pPatch%fn(XCOORD,ifl)
        
          IF ( FloatEqual(ABS(nx),1.0_RFREAL,eqTol) .EQV. .TRUE. ) THEN 
            RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET
          ELSE 
            RFLU_GetConstrType = CONSTR_TYPE_NONE
          END IF ! FloatEqual      
        CASE ( BC_NOSLIPWALL_HFLUX ) 
          RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET
        CASE ( BC_NOSLIPWALL_TEMP ) 
          RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET 
        CASE ( BC_INJECTION ) 
          nx = pPatch%fn(XCOORD,ifl)
        
          IF ( FloatEqual(ABS(nx),1.0_RFREAL,eqTol) .EQV. .TRUE. ) THEN 
            RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET
          ELSE 
            RFLU_GetConstrType = CONSTR_TYPE_NONE
          END IF ! FloatEqual                             
        CASE DEFAULT ! defensive coding
          RFLU_GetConstrType = CONSTR_TYPE_NONE
      END SELECT ! pPatch%bcType 

! ------------------------------------------------------------------------------
!   y-velocity
! ------------------------------------------------------------------------------

    CASE ( V_MIXT_YVEL ) 
      SELECT CASE ( pPatch%bcType )
        CASE ( BC_SLIPWALL )
          ny = pPatch%fn(YCOORD,ifl)
        
          IF ( FloatEqual(ABS(ny),1.0_RFREAL,eqTol) .EQV. .TRUE. ) THEN 
            RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET
          ELSE 
            RFLU_GetConstrType = CONSTR_TYPE_NONE
          END IF ! FloatEqual      
        CASE ( BC_NOSLIPWALL_HFLUX ) 
          RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET
        CASE ( BC_NOSLIPWALL_TEMP ) 
          RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET
        CASE ( BC_INJECTION ) 
          ny = pPatch%fn(YCOORD,ifl)
        
          IF ( FloatEqual(ABS(ny),1.0_RFREAL,eqTol) .EQV. .TRUE. ) THEN 
            RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET
          ELSE 
            RFLU_GetConstrType = CONSTR_TYPE_NONE
          END IF ! FloatEqual                                     
        CASE DEFAULT ! defensive coding
          RFLU_GetConstrType = CONSTR_TYPE_NONE
      END SELECT ! pPatch%bcType   

! ------------------------------------------------------------------------------
!   z-velocity
! ------------------------------------------------------------------------------

    CASE ( V_MIXT_ZVEL )
      SELECT CASE ( pPatch%bcType )
        CASE ( BC_SLIPWALL )
          nz = pPatch%fn(ZCOORD,ifl)
        
          IF ( FloatEqual(ABS(nz),1.0_RFREAL,eqTol) .EQV. .TRUE. ) THEN 
            RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET
          ELSE 
            RFLU_GetConstrType = CONSTR_TYPE_NONE
          END IF ! FloatEqual      
        CASE ( BC_NOSLIPWALL_HFLUX ) 
          RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET
        CASE ( BC_NOSLIPWALL_TEMP ) 
          RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET 
        CASE ( BC_INJECTION ) 
          nz = pPatch%fn(ZCOORD,ifl)
        
          IF ( FloatEqual(ABS(nz),1.0_RFREAL,eqTol) .EQV. .TRUE. ) THEN 
            RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET
          ELSE 
            RFLU_GetConstrType = CONSTR_TYPE_NONE
          END IF ! FloatEqual                                        
        CASE DEFAULT ! defensive coding
          RFLU_GetConstrType = CONSTR_TYPE_NONE
      END SELECT ! pPatch%bcType 
      
! ------------------------------------------------------------------------------
!   Temperature
! ------------------------------------------------------------------------------

    CASE ( V_MIXT_TEMP ) 
      SELECT CASE ( pPatch%bcType )      
        CASE ( BC_NOSLIPWALL_HFLUX ) 
          RFLU_GetConstrType = CONSTR_TYPE_VONNEUMANN
        CASE ( BC_NOSLIPWALL_TEMP ) 
          RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET        
        CASE ( BC_INJECTION )
          dist = pPatch%mixt%distrib

          minj = pPatch%mixt%vals(BCDAT_INJECT_MFRATE,dist*ifl)           

          IF ( minj > 0.0_RFREAL ) THEN
            RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET 
          ELSE
#ifndef GENX       
            RFLU_GetConstrType = CONSTR_TYPE_NONE
#else
            IF ( pRegion%mixtInput%flowModel == FLOW_EULER ) THEN 
              RFLU_GetConstrType = CONSTR_TYPE_NONE 
            ELSE 
              RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET
            END IF ! pRegion%mixtInput%flowModel
#endif      
          END IF ! minj                    
        CASE DEFAULT ! defensive coding
          RFLU_GetConstrType = CONSTR_TYPE_NONE
      END SELECT ! pPatch%bcType   

#ifdef SPEC
! ==============================================================================  
!   Species
! ==============================================================================  

    CASE ( V_SPEC_VAR1:V_SPEC_VAR9 ) 
      SELECT CASE ( pPatch%bcType )   
        CASE ( BC_INFLOW_TOTANG )
          RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET         
        CASE ( BC_INFLOW_VELTEMP )
          RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET         
        CASE ( BC_INJECTION )
          RFLU_GetConstrType = CONSTR_TYPE_DIRICHLET         
        CASE DEFAULT        
          RFLU_GetConstrType = CONSTR_TYPE_NONE                      
      END SELECT ! pPatch%bcType                
#endif 
 
! ==============================================================================  
!   Default - always unconstrained
! ==============================================================================  

    CASE DEFAULT
      RFLU_GetConstrType = CONSTR_TYPE_NONE
  END SELECT ! var

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END FUNCTION RFLU_GetConstrType











! ******************************************************************************
!
! Purpose: Get constraint value for given variable and patch.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   pPatch      Pointer to patch data
!   var         Variable index
!   ifl         Face index
!
! Output: 
!   constrType  Constraint type
!   constrVal   Constraint value
!
! Notes: None.
!
! ******************************************************************************

FUNCTION RFLU_GetConstrValue(pRegion,pPatch,var,ifl)

  USE RFLU_ModGridSpeedUtils, ONLY: RFLU_DescaleGridSpeed
  USE RFLU_ModRindStates, ONLY: RFLU_SetRindStateInjectPerf

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

  REAL(RFREAL) :: RFLU_GetConstrValue

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: ifl,var
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: dist,icg,indCp,indGs,indMol,gasModel
#ifdef SPEC
  INTEGER :: iSpec
#endif
  REAL(RFREAL) :: cp,eqTol,fs,fsu,hflux,Hl,minj,mm,nx,ny,nz,pl,rl,tinj,ul,vl,wl
#ifdef GENX
  REAL(RFREAL) :: tb
#endif
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_GetConstrValue',__FILE__)

  eqTol = 1.0E-6_RFREAL

  indCp    = pRegion%mixtInput%indCp
  indMol   = pRegion%mixtInput%indMol
  gasModel = pRegion%mixtInput%gasModel
  
  indGs = pRegion%grid%indGs

! ******************************************************************************
! Set dimensions and pointers
! ******************************************************************************

  SELECT CASE ( var ) 
  
! ==============================================================================  
!   Mixture
! ==============================================================================  

! ------------------------------------------------------------------------------
!   Density
! ------------------------------------------------------------------------------

    CASE ( V_MIXT_DENS ) 
      SELECT CASE ( pPatch%bcType )             
        CASE ( BC_INJECTION )
          nx = pPatch%fn(XCOORD,ifl)
          ny = pPatch%fn(YCOORD,ifl)
          nz = pPatch%fn(ZCOORD,ifl)          
        
          dist = pPatch%mixt%distrib

          IF ( gasModel == GAS_MODEL_TCPERF ) THEN
            icg = pPatch%bf2c(ifl)

            fs  = pPatch%gs(indGs*ifl)
            fsu = RFLU_DescaleGridSpeed(pRegion,fs)

            cp = pRegion%mixt%gv(GV_MIXT_CP ,indCp *icg)
            mm = pRegion%mixt%gv(GV_MIXT_MOL,indMol*icg)

            minj = pPatch%mixt%vals(BCDAT_INJECT_MFRATE,dist*ifl)           

            IF ( minj > 0.0_RFREAL ) THEN 
              tinj = pPatch%mixt%vals(BCDAT_INJECT_TEMP,dist*ifl)

              pl = pRegion%mixt%dv(DV_MIXT_PRES,icg)

              CALL RFLU_SetRindStateInjectPerf(cp,mm,nx,ny,nz,minj,tinj, & 
                                               pl,fsu,rl,ul,vl,wl,Hl)

              RFLU_GetConstrValue = rl
            ELSE ! defensive coding
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END IF ! minj
          ELSE 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END IF ! gasModel                                
        CASE DEFAULT ! defensive coding
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%bcType   

! ------------------------------------------------------------------------------
!   x-velocity
! ------------------------------------------------------------------------------
  
    CASE ( V_MIXT_XVEL ) 
      SELECT CASE ( pPatch%bcType )
        CASE ( BC_SLIPWALL )
          nx = pPatch%fn(XCOORD,ifl)
        
          IF ( FloatEqual(ABS(nx),1.0_RFREAL,eqTol) .EQV. .TRUE. ) THEN 
            RFLU_GetConstrValue = 0.0_RFREAL ! TEMPORARY - should be rel vel
          ELSE ! defensive coding
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END IF ! FloatEqual      
        CASE ( BC_NOSLIPWALL_HFLUX ) 
          RFLU_GetConstrValue = 0.0_RFREAL ! TEMPORARY - should be rel vel
        CASE ( BC_NOSLIPWALL_TEMP ) 
          RFLU_GetConstrValue = 0.0_RFREAL ! TEMPORARY - should be rel vel                             
        CASE ( BC_INJECTION ) 
          nx = pPatch%fn(XCOORD,ifl)
          ny = pPatch%fn(YCOORD,ifl)
          nz = pPatch%fn(ZCOORD,ifl)          
        
          IF ( FloatEqual(ABS(nx),1.0_RFREAL,eqTol) .EQV. .TRUE. ) THEN
            dist = pPatch%mixt%distrib
           
            IF ( gasModel == GAS_MODEL_TCPERF ) THEN
              icg = pPatch%bf2c(ifl)

              fs  = pPatch%gs(indGs*ifl)
              fsu = RFLU_DescaleGridSpeed(pRegion,fs)
              
              cp = pRegion%mixt%gv(GV_MIXT_CP ,indCp *icg)
              mm = pRegion%mixt%gv(GV_MIXT_MOL,indMol*icg)
                           
              minj = pPatch%mixt%vals(BCDAT_INJECT_MFRATE,dist*ifl)           

              IF ( minj > 0.0_RFREAL ) THEN 
                tinj = pPatch%mixt%vals(BCDAT_INJECT_TEMP,dist*ifl)

                pl = pRegion%mixt%dv(DV_MIXT_PRES,icg)
                
                CALL RFLU_SetRindStateInjectPerf(cp,mm,nx,ny,nz,minj,tinj, & 
                                                 pl,fsu,rl,ul,vl,wl,Hl)

                RFLU_GetConstrValue = ul
              ELSE 
                RFLU_GetConstrValue = fsu*nx
              END IF ! minj
            ELSE 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END IF ! gasModel
          ELSE ! defensive coding
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END IF ! FloatEqual               
        CASE DEFAULT ! defensive coding
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%bcType 

! ------------------------------------------------------------------------------
!   y-velocity
! ------------------------------------------------------------------------------

    CASE ( V_MIXT_YVEL ) 
      SELECT CASE ( pPatch%bcType )
        CASE ( BC_SLIPWALL )
          ny = pPatch%fn(YCOORD,ifl)
        
          IF ( FloatEqual(ABS(ny),1.0_RFREAL,eqTol) .EQV. .TRUE. ) THEN 
            RFLU_GetConstrValue = 0.0_RFREAL ! TEMPORARY - should be rel vel
          ELSE ! defensive coding
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END IF ! FloatEqual      
        CASE ( BC_NOSLIPWALL_HFLUX ) 
          RFLU_GetConstrValue = 0.0_RFREAL ! TEMPORARY - should be rel vel
        CASE ( BC_NOSLIPWALL_TEMP ) 
          RFLU_GetConstrValue = 0.0_RFREAL ! TEMPORARY - should be rel vel                             
        CASE ( BC_INJECTION ) 
          nx = pPatch%fn(XCOORD,ifl)
          ny = pPatch%fn(YCOORD,ifl)
          nz = pPatch%fn(ZCOORD,ifl)          
        
          IF ( FloatEqual(ABS(ny),1.0_RFREAL,eqTol) .EQV. .TRUE. ) THEN
            dist = pPatch%mixt%distrib
           
            IF ( gasModel == GAS_MODEL_TCPERF ) THEN
              icg = pPatch%bf2c(ifl)

              fs  = pPatch%gs(indGs*ifl)
              fsu = RFLU_DescaleGridSpeed(pRegion,fs)
              
              cp = pRegion%mixt%gv(GV_MIXT_CP ,indCp *icg)
              mm = pRegion%mixt%gv(GV_MIXT_MOL,indMol*icg)
                           
              minj = pPatch%mixt%vals(BCDAT_INJECT_MFRATE,dist*ifl)           

              IF ( minj > 0.0_RFREAL ) THEN 
                tinj = pPatch%mixt%vals(BCDAT_INJECT_TEMP,dist*ifl)

                pl = pRegion%mixt%dv(DV_MIXT_PRES,icg)
                
                CALL RFLU_SetRindStateInjectPerf(cp,mm,nx,ny,nz,minj,tinj, & 
                                                 pl,fsu,rl,ul,vl,wl,Hl)

                RFLU_GetConstrValue = vl
              ELSE 
                RFLU_GetConstrValue = fsu*ny
              END IF ! minj
            ELSE 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END IF ! gasModel
          ELSE ! defensive coding
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END IF ! FloatEqual            
        CASE DEFAULT ! defensive coding
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%bcType 

! ------------------------------------------------------------------------------
!   z-velocity
! ------------------------------------------------------------------------------

    CASE ( V_MIXT_ZVEL )
      SELECT CASE ( pPatch%bcType )
        CASE ( BC_SLIPWALL )
          nz = pPatch%fn(ZCOORD,ifl)
        
          IF ( FloatEqual(ABS(nz),1.0_RFREAL,eqTol) .EQV. .TRUE. ) THEN 
            RFLU_GetConstrValue = 0.0_RFREAL ! TEMPORARY - should be rel vel
          ELSE ! defensive coding
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END IF ! FloatEqual      
        CASE ( BC_NOSLIPWALL_HFLUX ) 
          RFLU_GetConstrValue = 0.0_RFREAL ! TEMPORARY - should be rel vel
        CASE ( BC_NOSLIPWALL_TEMP ) 
          RFLU_GetConstrValue = 0.0_RFREAL ! TEMPORARY - should be rel vel                             
        CASE ( BC_INJECTION ) 
          nx = pPatch%fn(XCOORD,ifl)
          ny = pPatch%fn(YCOORD,ifl)
          nz = pPatch%fn(ZCOORD,ifl)          
        
          IF ( FloatEqual(ABS(nz),1.0_RFREAL,eqTol) .EQV. .TRUE. ) THEN
            dist = pPatch%mixt%distrib
           
            IF ( gasModel == GAS_MODEL_TCPERF ) THEN
              icg = pPatch%bf2c(ifl)

              fs  = pPatch%gs(indGs*ifl)
              fsu = RFLU_DescaleGridSpeed(pRegion,fs)
              
              cp = pRegion%mixt%gv(GV_MIXT_CP ,indCp *icg)
              mm = pRegion%mixt%gv(GV_MIXT_MOL,indMol*icg)
                           
              minj = pPatch%mixt%vals(BCDAT_INJECT_MFRATE,dist*ifl)           

              IF ( minj > 0.0_RFREAL ) THEN 
                tinj = pPatch%mixt%vals(BCDAT_INJECT_TEMP,dist*ifl)

                pl = pRegion%mixt%dv(DV_MIXT_PRES,icg)
                
                CALL RFLU_SetRindStateInjectPerf(cp,mm,nx,ny,nz,minj,tinj, & 
                                                 pl,fsu,rl,ul,vl,wl,Hl)

                RFLU_GetConstrValue = wl
              ELSE 
                RFLU_GetConstrValue = fsu*nz
              END IF ! minj
            ELSE 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END IF ! gasModel
          ELSE ! defensive coding
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END IF ! FloatEqual            
        CASE DEFAULT ! defensive coding
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%bcType
      
! ------------------------------------------------------------------------------
!   Temperature
! ------------------------------------------------------------------------------

    CASE ( V_MIXT_TEMP ) 
      SELECT CASE ( pPatch%bcType )      
        CASE ( BC_NOSLIPWALL_HFLUX ) 
          SELECT CASE ( pPatch%mixt%distrib ) 
            CASE ( BCDAT_CONSTANT ) 
              icg   = pPatch%bf2c(ifl)              
              hflux = pPatch%mixt%vals(BCDAT_NOSLIP_Q,0)              
              
              RFLU_GetConstrValue = -hflux/pRegion%mixt%tv(TV_MIXT_TCOL,icg)
            CASE ( BCDAT_DISTRIB ) 
              icg   = pPatch%bf2c(ifl)              
              hflux = pPatch%mixt%vals(BCDAT_NOSLIP_Q,ifl)
              
              RFLU_GetConstrValue = -hflux/pRegion%mixt%tv(TV_MIXT_TCOL,icg)             
            CASE DEFAULT ! defensive coding
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pPatch%mixt%distrib
        CASE ( BC_NOSLIPWALL_TEMP ) 
          SELECT CASE ( pPatch%mixt%distrib ) 
            CASE ( BCDAT_CONSTANT ) 
              RFLU_GetConstrValue = pPatch%mixt%vals(BCDAT_NOSLIP_T,0)
            CASE ( BCDAT_DISTRIB )          
              RFLU_GetConstrValue = pPatch%mixt%vals(BCDAT_NOSLIP_T,ifl)
            CASE DEFAULT ! defensive coding
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pPatch%mixt%distrib         
        CASE ( BC_INJECTION ) 
          dist = pPatch%mixt%distrib
        
          minj = pPatch%mixt%vals(BCDAT_INJECT_MFRATE,dist*ifl)           

#ifndef GENX    
          IF ( minj > 0.0_RFREAL ) THEN
            tinj = pPatch%mixt%vals(BCDAT_INJECT_TEMP,dist*ifl)
             
            RFLU_GetConstrValue = tinj
          ELSE ! defensive coding
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END IF ! minj
#else
          IF ( pRegion%mixtInput%flowModel == FLOW_EULER ) THEN 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
	  ELSE
	    tb = pPatch%mixt%vals(BCDAT_INJECT_TEMP,dist*ifl)
	    
	    RFLU_GetConstrValue = tinj
	  END IF ! pRegion%mixtInput%flowModel
#endif	  	  
        CASE DEFAULT ! defensive coding
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pPatch%bcType   

#ifdef SPEC
! ==============================================================================  
!   Species
! ==============================================================================  

    CASE ( V_SPEC_VAR1:V_SPEC_VAR9 )
      iSpec = var - V_SPEC_VAR1 + 1
    
      SELECT CASE ( pPatch%bcType ) 
        CASE ( BC_INFLOW_TOTANG )
          SELECT CASE ( pPatch%spec%distrib )           
            CASE ( BCDAT_CONSTANT ) 
              RFLU_GetConstrValue = pPatch%spec%vals(iSpec,0)
            CASE ( BCDAT_DISTRIB )          
              RFLU_GetConstrValue = pPatch%spec%vals(iSpec,ifl)
            CASE DEFAULT ! defensive coding
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pPatch%spec%distrib                         
        CASE ( BC_INFLOW_VELTEMP )
          SELECT CASE ( pPatch%spec%distrib )           
            CASE ( BCDAT_CONSTANT ) 
              RFLU_GetConstrValue = pPatch%spec%vals(iSpec,0)
            CASE ( BCDAT_DISTRIB )          
              RFLU_GetConstrValue = pPatch%spec%vals(iSpec,ifl)
            CASE DEFAULT ! defensive coding
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pPatch%spec%distrib        
        CASE ( BC_INJECTION )
          SELECT CASE ( pPatch%spec%distrib )           
            CASE ( BCDAT_CONSTANT ) 
              RFLU_GetConstrValue = pPatch%spec%vals(iSpec,0)
            CASE ( BCDAT_DISTRIB )          
              RFLU_GetConstrValue = pPatch%spec%vals(iSpec,ifl)
            CASE DEFAULT ! defensive coding
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pPatch%spec%distrib          
        CASE DEFAULT ! defensive coding       
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)    
      END SELECT ! pPatch%bcType     
#endif

! ==============================================================================  
!   Default 
! ==============================================================================  

    CASE DEFAULT ! defensive coding
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! var

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END FUNCTION RFLU_GetConstrValue








! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModConstraintUtils


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModConstraintUtils.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:36  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:24  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:39  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.5  2006/08/19 15:39:01  mparmar
! Renamed patch variables
!
! Revision 1.4  2006/04/07 15:19:18  haselbac
! Removed tabs
!
! Revision 1.3  2005/10/31 21:09:36  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.2  2005/10/14 14:05:38  haselbac
! Added setting of boundary temperature for GENX simulations
!
! Revision 1.1  2005/10/05 14:33:44  haselbac
! Initial revision
!
! ******************************************************************************
  
  





