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
! Purpose: Set gas variables.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!  icgBeg       Beginning cell index
!  icgEnd       Ending cell index
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_SetGasVars.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_SetGasVars(pRegion,icgBeg,icgEnd)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid

#ifdef SPEC
  USE ModSpecies, ONLY: t_spec_type
#endif

  USE ModInterfaces, ONLY: MixtPerf_M_R, & 
                           MixtPerf_R_CpG, &
                           MixtPerf_R_M
                                 
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(IN) :: icgBeg,icgEnd
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: icg,indCp,indMol
  REAL(RFREAL) :: cpGasRef,gc,gGas,gGasRef
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pGv
  TYPE(t_global), POINTER :: global

#ifdef SPEC
  INTEGER :: iSpec
  REAL(RFREAL) :: cp,gcg,gcm,imm,immg,mmg,phip,r,Yg,Yi
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
  TYPE(t_spec_type), POINTER :: pSpecType
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_SetGasVars.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SetGasVars',__FILE__)

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pCv => pRegion%mixt%cv 
  pGv => pRegion%mixt%gv 

#ifdef SPEC
  pCvSpec => pRegion%spec%cv
#endif

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol
  
  cpGasRef = global%refCp
  gGasRef  = global%refGamma
  
! ******************************************************************************  
! Compute dependent variables
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%fluidModel ) 

! ==============================================================================  
!   Incompressible fluid model
! ==============================================================================  

    CASE ( FLUID_MODEL_INCOMP ) 

! ==============================================================================  
!   Compressible fluid model
! ==============================================================================  
    
    CASE ( FLUID_MODEL_COMP ) 
      SELECT CASE ( pRegion%mixtInput%gasModel )
      
! ------------------------------------------------------------------------------
!       Thermally and calorically perfect gas
! ------------------------------------------------------------------------------     
       
        CASE ( GAS_MODEL_TCPERF ) 
          gc = MixtPerf_R_CpG(cpGasRef,gGasRef)
  
          DO icg = icgBeg,icgEnd
            pGv(GV_MIXT_MOL,icg*indMol) = MixtPerf_M_R(gc)
            pGv(GV_MIXT_CP ,icg*indCp ) = cpGasRef    
          END DO ! icg

! ------------------------------------------------------------------------------
!       Mixture of thermally and calorically perfect gases
! ------------------------------------------------------------------------------     

        CASE ( GAS_MODEL_MIXT_TCPERF )         
#ifdef SPEC
          IF ( pRegion%spec%cvState /= CV_MIXT_STATE_PRIM ) THEN 
            CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
          END IF ! pRegion%spec%cvState  

          DO icg = icgBeg,icgEnd            
            imm = 0.0_RFREAL
            cp  = 0.0_RFREAL
            
            DO iSpec = 1,pRegion%specInput%nSpecies
              pSpecType => pRegion%specInput%specType(iSpec)
                        
              imm = imm + pCvSpec(iSpec,icg)/pSpecType%pMaterial%molw
              cp  = cp  + pCvSpec(iSpec,icg)*pSpecType%pMaterial%spht
            END DO ! iSpec
          
            pGv(GV_MIXT_MOL,icg) = 1.0_RFREAL/imm
            pGv(GV_MIXT_CP ,icg) = cp
          END DO ! icg
#else
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) ! Defensive coding
#endif              

! ------------------------------------------------------------------------------
!       Mixture of thermally and calorically perfect gaseous species and 
!       particulate species, treated as pseudo-gas
! ------------------------------------------------------------------------------     

        CASE ( GAS_MODEL_MIXT_PSEUDO )         
#ifdef SPEC
          IF ( pRegion%spec%cvState /= CV_MIXT_STATE_PRIM ) THEN 
            CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
          END IF ! pRegion%spec%cvState  

          DO icg = icgBeg,icgEnd
            r = pCv(CV_MIXT_DENS,icg)

            immg = 0.0_RFREAL
            cp   = 0.0_RFREAL
            Yg   = 0.0_RFREAL
            phip = 0.0_RFREAL

            DO iSpec = 1,pRegion%specInput%nSpecies
              pSpecType => pRegion%specInput%specType(iSpec)

              Yi = pCvSpec(iSpec,icg)

              cp = cp + Yi*pSpecType%pMaterial%spht

              IF ( pSpecType%discreteFlag .EQV. .TRUE. ) THEN 
                phip = phip + r*Yi/pSpecType%pMaterial%dens
              ELSE 
                immg = immg + Yi/pSpecType%pMaterial%molw
                Yg   = Yg   + Yi                                                              
              END IF ! pSpecType%discreteFlag
            END DO ! iSpec

            mmg = Yg/immg
            gcg = MixtPerf_R_M(mmg)
            gcm = Yg*gcg

            pGv(GV_MIXT_MOL,icg) = MixtPerf_M_R(gcm)
            pGv(GV_MIXT_CP ,icg) = cp
          END DO ! icg
#else
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) ! Defensive coding
#endif                   

! ------------------------------------------------------------------------------     
!       Gas-liquid mixture fluid model
! ------------------------------------------------------------------------------     

        CASE ( GAS_MODEL_MIXT_GASLIQ ) 
#ifdef SPEC
          IF ( pRegion%spec%cvState /= CV_MIXT_STATE_PRIM ) THEN 
            CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
          END IF ! pRegion%spec%cvState  

          DO icg = icgBeg,icgEnd
            pGv(GV_MIXT_MOL,icg*indMol) = &
              pRegion%specInput%specType(1)%pMaterial%molw
            pGv(GV_MIXT_CP ,icg*indCp ) = &
              pRegion%specInput%specType(1)%pMaterial%spht 
          END DO ! icg

#else
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) ! Defensive coding
#endif 

! ------------------------------------------------------------------------------
!       Mixture of thermally and calorically perfect gas and detonation products
! ------------------------------------------------------------------------------     

        CASE ( GAS_MODEL_MIXT_JWL )         
#ifdef SPEC
          IF ( pRegion%spec%cvState /= CV_MIXT_STATE_PRIM ) THEN 
            CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
          END IF ! pRegion%spec%cvState  

          DO icg = icgBeg,icgEnd            
            imm = 0.0_RFREAL
            cp  = 0.0_RFREAL
            
            DO iSpec = 1,pRegion%specInput%nSpecies
              pSpecType => pRegion%specInput%specType(iSpec)
                        
              imm = imm + pCvSpec(iSpec,icg)/pSpecType%pMaterial%molw
              cp  = cp  + pCvSpec(iSpec,icg)*pSpecType%pMaterial%spht
            END DO ! iSpec

! Subbu - See RFLU_SetTransportVars for cp & molwt calculation
! ************* Currently the Gas variables cp & molwt assigned here will be overridden
! in the function RFLU_SetTransportVars as RFLU_SetTransportVars is called after
! RFLU_SetGasVars subroutine ********************************************************* !
 
            pGv(GV_MIXT_MOL,icg) = 1.0_RFREAL/imm
            pGv(GV_MIXT_CP ,icg) = cp
          END DO ! icg
#else
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) ! Defensive coding
#endif   

! ------------------------------------------------------------------------------
!       Default
! ------------------------------------------------------------------------------     

        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pRegion%mixtInput%gasModel
      
! ==============================================================================  
!   Default
! ==============================================================================  

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)    
  END SELECT ! pRegion%mixtInput%fluidModel


! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SetGasVars

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_SetGasVars.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:36  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:53  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.8  2006/08/15 15:23:54  haselbac
! Cosmetic changes only
!
! Revision 1.7  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.6  2006/03/26 20:21:44  haselbac
! Added support for GL model
!
! Revision 1.5  2005/11/14 16:55:17  haselbac
! Added support for pseudo-gas model
!
! Revision 1.4  2005/11/10 02:15:33  haselbac
! Added support for different gas models
!
! Revision 1.3  2005/04/15 15:06:21  haselbac
! Added range arguments, removed pGrid declaration
!
! Revision 1.2  2004/11/14 19:43:09  haselbac
! Added setting of gv for incompressible fluid model
!
! Revision 1.1  2004/10/19 19:23:55  haselbac
! Initial revision
!
! ******************************************************************************

