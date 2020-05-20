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
! Purpose: Collection of routines for Eulerian-based particle infrastructure.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: PLAG_ModEulerian.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE PLAG_ModEulerian

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY: t_plag 
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModMPI
  USE PLAG_ModParameters

  IMPLICIT NONE
  
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: PLAG_ModEulerian.F90,v $ $Revision: 1.1.1.1 $'        

! ==============================================================================
! Public functions
! ==============================================================================
  
  PUBLIC :: PLAG_CreateEulerianField,  & 
            PLAG_DestroyEulerianField, & 
            PLAG_InitEulerianField,    &
            PLAG_RFLU_CalcEulerianField

! ==============================================================================
! Private functions
! ==============================================================================

!  PRIVATE :: 

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  
  
  
  
  
  
  
  

  





! *******************************************************************************
!
! Purpose: Create eulerian field.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPlag               Pointer to plag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_CreateEulerianField(pRegion,pPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion

    TYPE(t_plag), POINTER :: pPlag

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,ibc,iec,nEv
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_CreateEulerianField',__FILE__)

    nEv = pPlag%nEv        

    pGrid => pRegion%grid
    ibc = 1
    iec = pGrid%nCellsTot
  
! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    ALLOCATE(pPlag%ev(nEv,ibc:iec),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%ev')
    END IF ! global%error                       

! ******************************************************************************
!   Initialize
! ******************************************************************************

    CALL PLAG_InitEulerianField(pRegion,pPlag)

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_CreateEulerianField






! *******************************************************************************
!
! Purpose: Destroy Eulerian field.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPlag               Pointer to plag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_DestroyEulerianField(pRegion,pPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion

    TYPE(t_plag), POINTER :: pPlag
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag
    TYPE(t_global), POINTER :: global    

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_DestroyEulerianField',__FILE__)

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************
        
    DEALLOCATE(pPlag%ev,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPlag%ev')
    END IF ! global%error

! ******************************************************************************
!   Nullify memory
! ******************************************************************************
 
    NULLIFY(pPlag%ev)

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_DestroyEulerianField

  



! *******************************************************************************
!
! Purpose: Initialize Eulerian field.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPlag               Pointer to plag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_InitEulerianField(pRegion,pPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
     
    TYPE(t_region), POINTER :: pRegion

    TYPE(t_plag), POINTER :: pPlag
  
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,ibc,iec,iCell,iVar,nVars
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'PLAG_InitEulerianField',__FILE__)

    nVars = pPlag%nEv       

    pGrid => pRegion%grid
    ibc = 1
    iec = pGrid%nCellsTot

! ******************************************************************************
!   Initialize memory
! ******************************************************************************
        
    DO iCell = ibc,iec
      DO iVar = 1, nVars
        pPlag%ev(iVar,iCell)  = 0.0_RFREAL
      END DO ! iVar                  
    END DO ! iCell                      

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_InitEulerianField







! *******************************************************************************
!
! Purpose: Compute Eulerian-based field from Lagrangian datastructure.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   ev                  array containing Eulerian field
!
! Output: None.
!
! Notes: Routine relevant to write TECPLOT files in RFLU.
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLU_CalcEulerianField(pRegion,ev)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), DIMENSION(:,:) :: ev     
    TYPE(t_region), POINTER :: pRegion
  
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,ibc,icg,iCont,iec,iPcl,iVar,nCells1D,nCont,nPcls, & 
               nPlagCell,nVars,nVarsPlot
    INTEGER, POINTER, DIMENSION(:)   :: cvMass    

    REAL(RFREAL) :: diamMicron,massL,massSqrL,nPlagCellInv
    REAL(RFREAL) :: massFrac,densMixt,spLoadL,volL,volMixt
    REAL(RFREAL) :: densLTot,eta,etaAvg,etaLTot,massFracLTot,massFracSum, &
                    massGasTot,massLTot,massMixt,massPart,numDensTot, &
                    relVel(3),relVelMagL,reyL, spLoadLTot,volFrac, &
                    volFracSum,volFracLTot,volLTot,volGasTot
                    
    REAL(RFREAL), POINTER, DIMENSION(:)   :: dens
    REAL(RFREAL), POINTER, DIMENSION(:,:) :: cv,dv,tv

    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'PLAG_RFLU_CalcEulerianField',__FILE__)

    pGrid => pRegion%grid
    pPlag => pRegion%plag

    ibc = 1
    iec = pGrid%nCellsTot
   
    nCont = pRegion%plagInput%nCont
    nPcls = pPlag%nPcls
    nVars = pPlag%nEv

    cvMass => pPlag%cvPlagMass
    cv => pPlag%cv
    dv => pPlag%dv
    tv => pPlag%tv
    dens => pRegion%plagInput%dens

    IF ( pPlag%nPcls > 0 ) THEN

! ******************************************************************************
!     Check if size of array sent is aligned with the expected size
! ******************************************************************************

! CHECK
    nVarsPlot = pPlag%nEv -nCont
    
    IF ( nVarsPlot /= SIZE(ev,1) ) THEN 
      WRITE(*,*) 'nVarsPlot   = ',nVarsPlot
      WRITE(*,*) 'SIZE(ev,1) = ',SIZE(ev,1)
      CALL ErrorStop(global,ERR_DATADIM_MISMATCH,__LINE__)
    END IF ! nVars 

    IF ( iec /= SIZE(ev,2) ) THEN 
      WRITE(*,*) 'iec        = ',iec
      WRITE(*,*) 'SIZE(ev,2) = ',SIZE(ev,2)
      CALL ErrorStop(global,ERR_DATADIM_MISMATCH,__LINE__)
    END IF ! nVars 
        
! END CHECK

! ******************************************************************************
!     Reinitialize Eulerian field to compute instantaneous variables
! ******************************************************************************

      DO icg = ibc,iec
      DO iVar = 1, nVarsPlot
        ev(iVar,icg) = 0.0_RFREAL
      END DO ! iVar
      END DO ! icg

! ******************************************************************************
!     Compute instantaneous Eulerian field from Lagrangian datastructure
! ******************************************************************************

      DO iPcl = 1,nPcls
        icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcl)

        diamMicron = dv(DV_PLAG_DIAM,iPcl)*1.0E+06_RFREAL      
        massL      = SUM(cv(cvMass(:),iPcl))
        spLoadL    = pPlag%arv(ARV_PLAG_SPLOAD,iPcl)
        
        volL       = 0.0_RFREAL 
        DO iCont = 1, nCont
          volL = volL + cv(cvMass(iCont),iPcl)/dens(iCont)
        END DO ! iCont   

        ev(EV_PLAG_DIA3,icg) = ev(EV_PLAG_DIA3,icg) +diamMicron**3.0_RFREAL
        ev(EV_PLAG_DIA4,icg) = ev(EV_PLAG_DIA4,icg) +diamMicron**4.0_RFREAL

        ev(EV_PLAG_NDNS,icg) = ev(EV_PLAG_NDNS,icg) +1.0_RFREAL
 
        ev(EV_PLAG_UVEL,icg) = ev(EV_PLAG_UVEL,icg) +dv(DV_PLAG_UVEL,iPcl)
        ev(EV_PLAG_VVEL,icg) = ev(EV_PLAG_VVEL,icg) +dv(DV_PLAG_VVEL,iPcl)
        ev(EV_PLAG_WVEL,icg) = ev(EV_PLAG_WVEL,icg) +dv(DV_PLAG_WVEL,iPcl)
        
        ev(EV_PLAG_TEMP,icg) = ev(EV_PLAG_TEMP,icg) +dv(DV_PLAG_TEMP,iPcl)

        ev(EV_PLAG_MFRC,icg) = ev(EV_PLAG_MFRC,icg) +massL *spLoadL
        ev(EV_PLAG_VFRC,icg) = ev(EV_PLAG_VFRC,icg) +volL  *spLoadL

        relVel(1)   = dv(DV_PLAG_UVELMIXT,iPcl)-dv(DV_PLAG_UVEL,iPcl)
        relVel(2)   = dv(DV_PLAG_VVELMIXT,iPcl)-dv(DV_PLAG_VVEL,iPcl)
        relVel(3)   = dv(DV_PLAG_WVELMIXT,iPcl)-dv(DV_PLAG_WVEL,iPcl)
        relVelMagL  = SQRT( relVel(1)*relVel(1)+ &
                            relVel(2)*relVel(2)+ &
                            relVel(3)*relVel(3)  )

        reyL = dv(DV_PLAG_DIAM,iPcl) *relVelMagL *dv(DV_PLAG_DENSMIXT,iPcl) / &
                                                  tv(TV_PLAG_MUELMIXT,iPcl)

        ev(EV_PLAG_REYN,icg) = ev(EV_PLAG_REYN,icg) +reyL
      END DO ! iPcl

! ******************************************************************************
!     Scale field by cell-based number of particles
!       Also compute mass and volume fractions
! ******************************************************************************

      massFracSum = 0.0_RFREAL
      volFracSum  = 0.0_RFREAL

      DO icg = ibc,iec
        IF ( ev(EV_PLAG_NDNS,icg) > 0 ) THEN
          nPlagCellInv = 1.0_RFREAL/ev(EV_PLAG_NDNS,icg)

          ev(EV_PLAG_DIA3,icg) = ev(EV_PLAG_DIA3,icg) *nPlagCellInv
          ev(EV_PLAG_DIA4,icg) = ev(EV_PLAG_DIA4,icg) *nPlagCellInv
          ev(EV_PLAG_UVEL,icg) = ev(EV_PLAG_UVEL,icg) *nPlagCellInv
          ev(EV_PLAG_VVEL,icg) = ev(EV_PLAG_VVEL,icg) *nPlagCellInv
          ev(EV_PLAG_WVEL,icg) = ev(EV_PLAG_WVEL,icg) *nPlagCellInv
          ev(EV_PLAG_TEMP,icg) = ev(EV_PLAG_TEMP,icg) *nPlagCellInv
          ev(EV_PLAG_REYN,icg) = ev(EV_PLAG_REYN,icg) *nPlagCellInv
        ENDIF ! nPlagCell

        densMixt = pRegion%mixt%cv(CV_MIXT_DENS,icg)
        volMixt  = pRegion%grid%vol(icg)

        volFrac  = ev(EV_PLAG_VFRC,icg)/volMixt
        massFrac = ev(EV_PLAG_MFRC,icg)

        massFracSum = massFracSum +massFrac 
        volFracSum  = volFracSum  +densMixt *volMixt *(1.0_RFREAL-volFrac)

        ev(EV_PLAG_VFRC,icg) = volFrac 
      END DO ! icg  

!      WRITE(*,*) '>>> massFrac <<<',massFracSum/(massFracSum+volFracSum)

! ******************************************************************************
!     Compute total mass, volume in gas for whole volume in each region
!      Consider all cells irrespective of particles being present
! ******************************************************************************

      massLTot   = 0.0_RFREAL
      volLTot    = 0.0_RFREAL
      spLoadLTot = 0.0_RFREAL

      massGasTot = 0.0_RFREAL
      numDensTot = 0.0_RFREAL
      volGasTot  = 0.0_RFREAL

      DO icg = ibc,iec
       densMixt   = pRegion%mixt%cv(CV_MIXT_DENS,icg)
        volMixt    = pRegion%grid%vol(icg)
        massGasTot = massGasTot +densMixt *volMixt
        volGasTot  = volGasTot  +volMixt

        volLTot  = volLTot  +ev(EV_PLAG_VFRC,icg) *volMixt 
        massLTot = massLTot +ev(EV_PLAG_MFRC,icg)

        numDensTot = numDensTot +ev(EV_PLAG_NDNS,icg)
      END DO ! icg

      volFracLTot  = volLTot/volGasTot
      massFracLTot = massLTot/(massLTot +massGasTot *(1.0-volFracLTot))
      etaLTot      = massLTot/massGasTot

!      WRITE(*,*) 'massLTot, volLTot = ',massLTot,volLTot
!      WRITE(*,*) 'massGasTot , volGasTot = ',massGasTot,volGasTot
!      WRITE(*,*) 'Mass and Volume Fractions in Region for all cells'   
!      WRITE(*,*) 'Region           = ',pRegion%iRegionGlobal
!      WRITE(*,*) '  numDensTot/Vol   = ',numDensTot/volGasTot
!      WRITE(*,*) '  volFracTot       = ',volFracLTot
!      WRITE(*,*) '  massFracTot      = ',massFracLTot
!      WRITE(*,*) '  etaTot           = ',etaLTot

! ******************************************************************************
!     Compute total mass, volume in gas for whole volume in each region
!      Consider cells ONLY where particles are present
! ******************************************************************************

      massLTot   = 0.0_RFREAL
      volLTot    = 0.0_RFREAL

      massGasTot = 0.0_RFREAL
      numDensTot = 0.0_RFREAL
      volGasTot  = 0.0_RFREAL     
     
      DO icg = ibc,iec
        IF ( ev(EV_PLAG_NDNS,icg) > 0 ) THEN
          densMixt   = pRegion%mixt%cv(CV_MIXT_DENS,icg)
          volMixt    = pRegion%grid%vol(icg)
          massGasTot = massGasTot +densMixt *volMixt
          volGasTot  = volGasTot  +volMixt

          volLTot  = volLTot  +ev(EV_PLAG_VFRC,icg) *volMixt 
          massLTot = massLTot +ev(EV_PLAG_MFRC,icg)

          numDensTot = numDensTot + ev(EV_PLAG_NDNS,icg)
        END IF ! ev
      END DO ! icg

      volFracLTot  = volLTot/volGasTot
      massFracLTot = massLTot/(massLTot +massGasTot *(1.0-volFracLTot))
      etaLTot      = massLTot/massGasTot

!      WRITE(*,*) 'Mass and Volume Fractions in Region for cells where Particles are present'     
!      WRITE(*,*) 'Region           = ',pRegion%iRegionGlobal
!      WRITE(*,*) '  numDensTot/Vol   = ',numDensTot/volGasTot
!      WRITE(*,*) '  volFracTot       = ',volFracLTot
!      WRITE(*,*) '  massFracTot      = ',massFracLTot  
!      WRITE(*,*) '  etaTot           = ',etaLTot   

! ******************************************************************************
!    Compute Mass Fraction in each cell
! ******************************************************************************
     
      DO icg = ibc,iec
        ev(EV_PLAG_MFRC,icg) = 0.0_RFREAL
      END DO ! icg

      DO iPcl = 1,nPcls
        icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcl)
        
        massL      = SUM(cv(cvMass(:),iPcl))
        spLoadL    = pPlag%arv(ARV_PLAG_SPLOAD,iPcl)

        ev(EV_PLAG_MFRC,icg) = ev(EV_PLAG_MFRC,icg) +massL *spLoadL
      END DO ! icg  

      DO icg = ibc,iec
        densMixt = pRegion%mixt%cv(CV_MIXT_DENS,icg)
        volMixt  = pRegion%grid%vol(icg)
        
        massPart = ev(EV_PLAG_MFRC,icg)
        massMixt = densMixt *volMixt
        ev(EV_PLAG_MFRC,icg) = massPart/(massPart +massMixt) 
      END DO ! icg  

! ******************************************************************************
!    Scale number density by cell volume
! ******************************************************************************
     
      DO icg = ibc,iec
        volMixt  = pRegion%grid%vol(icg)

        ev(EV_PLAG_NDNS,icg) = ev(EV_PLAG_NDNS,icg)/volMixt
      END DO ! icg  
          
    END IF ! pPlag%nPcls


! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_RFLU_CalcEulerianField





END MODULE PLAG_ModEulerian

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModEulerian.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.6  2008/12/06 08:43:51  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2007/04/26 20:45:34  fnajjar
! Removed dv-based calculations that are irrelevant and modified the needed ones
!
! Revision 1.3  2007/04/24 14:07:05  fnajjar
! Cleaned up and added mass loading
!
! Revision 1.2  2007/04/16 23:21:41  fnajjar
! Deleted all RFLO-pertinent calls and statements
!
! Revision 1.1  2007/04/09 18:50:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.6  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.5  2006/02/16 14:51:58  fnajjar
! Bug fix for missing ev pointer array in PLAG_CalcEulerian
!
! Revision 1.4  2006/01/31 17:20:21  fnajjar
! Bug fix to embed PLAG_RFLU_CalcEulerian in ifdef RFLU
!
! Revision 1.3  2005/11/30 22:20:21  fnajjar
! Added PLAG_RFLU_CalcEulerianField
!
! Revision 1.2  2005/02/16 14:46:48  fnajjar
! Bug fix for w-velocity component of ev and cosmetics cleanup
!
! Revision 1.1  2005/01/08 20:44:32  fnajjar
! Initial import for PLAG statistics
!
! ******************************************************************************

