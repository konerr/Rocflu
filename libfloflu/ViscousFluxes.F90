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
! Purpose: Compute viscous fluxes and add them to dissipation.
!
! Description: Distinction made between laminar and turbulent (total)
!              viscous fluxes; laminar flux is computed by generic routine
!              ViscousFluxEddy if TurbModel is inactive, 
!              else turbulent flux routine is called; upon output 
!              viscous fluxes added to numerical dissipation fluxes
!
! Input: 
!   region      Data of current region.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ViscousFluxes.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ViscousFluxes(region)

  USE ModError
  USE ModParameters
  USE ModDataTypes
  USE ModBndPatch, ONLY: t_patch
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch
  
  USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, & 
                               RFLU_ConvertCvPrim2Cons  
  USE RFLU_ModDifferentiationBFaces
  USE RFLU_ModDifferentiationFaces
  USE RFLU_ModViscousFlux  
  USE ModInterfaces, ONLY : RFLU_DecideNeedBGradFace 
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region) :: region

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: iLev, turbModel
  INTEGER :: iPatch
  INTEGER :: varInfo(CV_MIXT_XVEL:CV_MIXT_TEMP)
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => region%global

  CALL RegisterFunction(global,'ViscousFluxes',__FILE__)

! ******************************************************************************
! Define input parameters and gridlevel
! ******************************************************************************
 
  turbModel = region%mixtInput%turbModel 
 
! ******************************************************************************
! Add viscous residual from laminar or laminar+turbulent stress contributions
! ******************************************************************************

  pRegion => region%pRegion

  CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWT)
  
  varInfo(CV_MIXT_XVEL) = V_MIXT_XVEL
  varInfo(CV_MIXT_YVEL) = V_MIXT_YVEL
  varInfo(CV_MIXT_ZVEL) = V_MIXT_ZVEL
  varInfo(CV_MIXT_TEMP) = V_MIXT_TEMP
  
  CALL RFLU_ComputeGradFacesWrapper(pRegion,CV_MIXT_XVEL,CV_MIXT_TEMP, & 
                                    GRF_MIXT_XVEL,GRF_MIXT_TEMP, &
                                    pRegion%mixt%cv,pRegion%mixt%gradFace)
                              
  IF ( pRegion%grid%nFacesConstr > 0 ) THEN                              
    CALL RFLU_ComputeGradFacesConstr(pRegion,CV_MIXT_XVEL,CV_MIXT_TEMP, & 
                                     GRF_MIXT_XVEL,GRF_MIXT_TEMP,varInfo, &
                                     pRegion%mixt%cv,pRegion%mixt%gradFace)                              
  END IF ! pRegion%grid%nFacesConstr      

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
    
      CALL RFLU_ComputeGradBFacesWrapper(pRegion,pPatch,CV_MIXT_XVEL, &
                                         CV_MIXT_TEMP,GRBF_MIXT_XVEL, &
                                         GRBF_MIXT_TEMP,pRegion%mixt%cv, &
                                         pPatch%mixt%gradFace)
    
      IF ( pPatch%cReconst /= CONSTR_NONE ) THEN                                        
        CALL RFLU_ComputeBFGradConstrWrapper(pRegion,pPatch,CV_MIXT_XVEL, &
                                             CV_MIXT_TEMP,GRBF_MIXT_XVEL, &
                                             GRBF_MIXT_TEMP,varInfo, &
                                             pRegion%mixt%cv, &
                                             pPatch%mixt%gradFace)                           
      END IF ! pPatch%cReconst                                           
    END IF ! RFLU_DecideNeedBGradFace
  END DO ! iPatch
                                                                   
  IF ( turbModel == TURB_MODEL_NONE ) THEN ! Only laminar                                    
    CALL RFLU_EnforceHeatFlux(pRegion,pRegion%mixt%tv,TV_MIXT_TCOL)  
    
    IF ( pRegion%mixtInput%axiFlag .EQV. .TRUE. ) THEN
      CALL RFLU_AXI_ViscousFluxes(pRegion,pRegion%mixt%tv,TV_MIXT_MUEL,TV_MIXT_TCOL)
      CALL RFLU_AXI_ViscousFluxesPatches(pRegion,pRegion%mixt%tv,TV_MIXT_MUEL, &
                                     TV_MIXT_TCOL)
    ELSE
      CALL RFLU_ViscousFluxes(pRegion,pRegion%mixt%tv,TV_MIXT_MUEL,TV_MIXT_TCOL)
      CALL RFLU_ViscousFluxesPatches(pRegion,pRegion%mixt%tv,TV_MIXT_MUEL, &
                                     TV_MIXT_TCOL)
    END IF ! pRegion%mixtInput%axiFlag
  ELSE 
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END IF ! turbModel                                                                
                                                              
  CALL RFLU_ConvertCvPrim2Cons(pRegion,CV_MIXT_STATE_CONS)

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(region%global)

END SUBROUTINE ViscousFluxes

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ViscousFluxes.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:33  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:16:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2008/05/29 01:35:13  mparmar
! Added calls to viscous flux routines for axi-symm mode
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.12  2006/08/19 15:38:37  mparmar
! Moved region%mixt%bGradFace to patch%mixt%gradFace
!
! Revision 1.11  2006/04/15 16:53:35  haselbac
! Added IF on cReconst flag on patch
!
! Revision 1.10  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.9  2006/04/07 14:39:22  haselbac
! Changed calls to bface grad routines, now inside patch loop
!
! Revision 1.8  2005/10/27 18:53:54  haselbac
! Adapted to new cell grad routines
!
! Revision 1.7  2005/10/14 14:02:47  haselbac
! Added call to RFLU_EnforceHeatFlux
!
! Revision 1.6  2005/10/05 13:49:57  haselbac
! Adapted to new face grads, cosmetics
!
! Revision 1.5  2005/05/16 20:39:51  haselbac
! Now USE RFLU_ModViscousFlux
!
! Revision 1.4  2005/03/09 14:52:04  haselbac
! Cosmetics only
!
! Revision 1.3  2005/01/05 01:45:44  haselbac
! Added better comments
!
! Revision 1.2  2005/01/05 01:43:31  haselbac
! Added init of global coeffs, now also defined for all regions
!
! Revision 1.1  2004/12/01 16:52:17  haselbac
! Initial revision after changing case
!
! Revision 1.30  2004/08/02 23:13:14  wasistho
! mv libfloflu/viscousFluxEddy(Patch) to rocflo/RFLO_viscousFlux(Patch)
!
! Revision 1.29  2004/07/03 01:26:59  wasistho
! activate turbulent viscous fluxes in Rocflu
!
! Revision 1.28  2004/03/13 03:09:52  wasistho
! TURB_CoViscousFluxesFlo to TURB_CoViscousFluxes
!
! Revision 1.27  2004/03/08 23:32:41  wasistho
! changed turb file/routine name
!
! Revision 1.26  2004/01/29 22:52:50  haselbac
! Added ONLY to USE module statement
!
! Revision 1.25  2003/12/04 03:23:10  haselbac
! Added RFLU specific routines
!
! Revision 1.24  2003/11/20 16:40:36  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.21  2003/10/03 20:49:45  haselbac
! Renamed ip variable to iPatch
!
! Revision 1.20  2003/08/28 20:32:33  wasistho
! excluced ModInterfacesTurbulence,Radiation,Periodic from ModInterfaces
!
! Revision 1.19  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.18  2003/04/10 01:22:41  jblazek
! Got rid of pRegion in ViscousFluxesMP.
!
! Revision 1.17  2002/12/12 02:54:40  wasistho
! Added TURB_ViscousFluxes interface
!
! Revision 1.16  2002/09/09 14:09:08  haselbac
! mixtInput now under regions, added gradient routines
!
! Revision 1.15  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.13  2002/08/27 23:58:20  wasistho
! TURB_TotalVisFluxes is set back to TURB_ViscousFluxes
!
! Revision 1.12  2002/08/16 16:51:49  wasistho
! TURB_ViscousFluxes renamed to TURB_TotalVisFluxes
!
! Revision 1.11  2002/08/15 19:48:05  jblazek
! Implemented grid deformation capability.
!
! Revision 1.10  2002/08/01 01:34:04  wasistho
! Included RFLU capability and changed CalcGradVelT
!
! Revision 1.9  2002/07/29 17:13:08  jblazek
! Clean up after RFLU and TURB.
!
! Revision 1.8  2002/07/27 08:13:16  wasistho
! prepared for rocturb implementation
!
! Revision 1.7  2002/07/23 01:08:33  wasistho
! Turbulent part of viscous flux moved to rocturb
!
! Revision 1.6  2002/07/22 22:59:11  jblazek
! Some more clean up.
!
! Revision 1.5  2002/07/19 23:38:24  wasistho
! made compliant with CODING RULE
!
! Revision 1.4  2002/06/14 21:17:01  wasistho
! Added #ifdef RFLO until RFLU ready
!
! Revision 1.3  2002/05/21 01:51:56  wasistho
! add viscous terms
!
! Revision 1.1  2002/01/23 03:51:24  jblazek
! Added low-level time-stepping routines.
!
! ******************************************************************************

