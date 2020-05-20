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
!******************************************************************************
!
! Purpose: calculate convective fluxes for RocfluidMP framework.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: regions%levels%mixt = new solution after one time step.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ConvectiveFluxesMP.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ConvectiveFluxesMP( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  USE ModInterfaces, ONLY: ConvectiveFluxes

#ifdef SPEC
  USE RFLU_ModConvertCv, ONLY: RFLU_ScalarConvertCvCons2Prim, &
                               RFLU_ScalarConvertCvPrim2Cons
  USE ModInterfaces, ONLY: RFLU_ScalarFirst, &
                           RFLU_ScalarSecond, &
                           RFLU_ScalarFirstPatch, &
                           RFLU_ScalarSecondPatch, &
                           RFLU_ScalarInitRhs
#endif

  IMPLICIT NONE

! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region) :: region

! =============================================================================
! Locals
! =============================================================================

  INTEGER :: iPatch,spaceDiscr,spaceOrder
  TYPE(t_global), POINTER :: global

#ifdef RFLU
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion
#endif

! *****************************************************************************
! Start
! *****************************************************************************

  global => region%global

  CALL RegisterFunction( global,'ConvectiveFluxesMP',__FILE__ )

! *****************************************************************************
! Get parameters
! *****************************************************************************

  spaceDiscr = region%mixtInput%spaceDiscr
  spaceOrder = region%mixtInput%spaceOrder

! *****************************************************************************
! Compute fluxes
! *****************************************************************************

! =============================================================================
! Mixture
! =============================================================================

  CALL ConvectiveFluxes( region )

#ifdef SPEC
! =============================================================================
! Species
! =============================================================================

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    pRegion => region%pRegion

    CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv, &
                                       pRegion%spec%cvState)

    CALL RFLU_ScalarInitRhs(pRegion,pRegion%specInput%nSpecies, &
                            pRegion%spec%diss,pRegion%spec%rhs)

    SELECT CASE ( spaceOrder ) 
      CASE ( 1 ) 
        CALL RFLU_ScalarFirst(pRegion,pRegion%specInput%nSpecies, &
                              pRegion%spec%cv,pRegion%spec%rhs)

        DO iPatch = 1,pRegion%grid%nPatches
          pPatch => pRegion%patches(iPatch)

          CALL RFLU_ScalarFirstPatch(pRegion,pPatch, &
                                     pRegion%specInput%nSpecies, &
                                     pRegion%spec%cv,pPatch%spec, &
                                     pRegion%spec%rhs)
        END DO ! iPatch                                     
      CASE ( 2 )     
        CALL RFLU_ScalarSecond(pRegion,pRegion%specInput%nSpecies, &
                               pRegion%spec%cv,pRegion%spec%gradCell, &
                               pRegion%spec%rhs)

        DO iPatch = 1,pRegion%grid%nPatches
          pPatch => pRegion%patches(iPatch)
          
          SELECT CASE ( pPatch%spaceOrder ) 
            CASE ( 1 ) 
              CALL RFLU_ScalarFirstPatch(pRegion,pPatch, &
                                         pRegion%specInput%nSpecies, &
                                         pRegion%spec%cv,pPatch%spec, &
                                         pRegion%spec%rhs)            
            CASE ( 2 )           
              CALL RFLU_ScalarSecondPatch(pRegion,pPatch, &
                                          pRegion%specInput%nSpecies, &
                                          pRegion%spec%cv, &
                                          pRegion%spec%gradCell, &
                                          pPatch%spec,pRegion%spec%rhs)
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pPatch%spaceOrder
        END DO ! iPatch        
      CASE DEFAULT 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! spaceOrder

    CALL RFLU_ScalarConvertCvPrim2Cons(pRegion,pRegion%spec%cv, &
                                       pRegion%spec%cvState)
  END IF ! global%specUsed
#endif

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE ConvectiveFluxesMP

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ConvectiveFluxesMP.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:24  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2006/08/19 15:38:28  mparmar
! Renamed patch variables
!
! Revision 1.2  2006/04/15 16:53:04  haselbac
! Added option of 1st order bfluxes with 2nd order scheme
!
! Revision 1.1  2004/12/01 16:48:27  haselbac
! Initial revision after changing case
!
! Revision 1.12  2004/07/28 15:29:18  jferry
! created global variable for spec use
!
! Revision 1.11  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.10  2004/01/29 22:52:40  haselbac
! Added second-order for species
!
! Revision 1.9  2003/11/25 21:01:38  haselbac
! Added rocspecies support
!
! Revision 1.8  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/10/03 20:11:55  wasistho
! initial installation of turbModel SA and DES
!
! Revision 1.4  2003/10/01 23:52:09  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.3  2003/09/26 21:45:34  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.2  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.1  2003/03/28 19:45:46  fnajjar
! Initial import for RocfluidMP
!
!******************************************************************************

