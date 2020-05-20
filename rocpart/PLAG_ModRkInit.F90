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
!*******************************************************************************
!
! Purpose: Suite of routines to initialize Runge-Kutta schemes.
!
! Description: None.
!
! Notes: None.
!
!*******************************************************************************
!
! $Id: PLAG_ModRkInit.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
!*******************************************************************************

MODULE PLAG_ModRkInit
  
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModDataTypes
  USE ModGlobal,     ONLY: t_global
  USE ModGrid,       ONLY: t_grid
  USE ModPartLag,    ONLY: t_plag, t_plag_input, t_tile_plag
  USE ModBndPatch,   ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModError

  USE PLAG_ModParameters
  USE INRT_ModParameters

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: PLAG_InjcTileRkInit, &
            PLAG_RkInitPrimary
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: PLAG_ModRkInit.F90,v $ $Revision: 1.1.1.1 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS




!******************************************************************************
!
! Purpose: Initialize Runge-Kutta scheme for Integer variables.
!
! Description: None.
!
! Input: 
!   region      Region data
!   iStage      Runge-Kutta stage
!   icBeg       Beginning index for cell update
!   icEnd       Ending index for cell update
!   ivBeg       Beginning index for variable update
!   ivEnd       Ending index for variable update
!   aiv         Integer variables
!   aivOld      Old integer variables
!
! Output: 
!   ivOld       Old integer variables
!
! Notes: None.
!
!******************************************************************************

SUBROUTINE PLAG_RkInitGenericInt(region,iStage,icBeg,icEnd,ivBeg,ivEnd,&
                                 aiv,aivOld)

! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: icBeg,icEnd,iStage,ivBeg,ivEnd
  INTEGER, DIMENSION(:,:), POINTER :: aiv,aivOld
  TYPE(t_region) :: region

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: ic,iv
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ModRkInit.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction(global,'PLAG_RkInitGenericInt',__FILE__)

! *****************************************************************************
! Initialize Runge-Kutta scheme
! *****************************************************************************

  SELECT CASE ( global%rkScheme ) 
    CASE ( RK_SCHEME_4_CLASSICAL ) 
      IF ( iStage == 1 ) THEN
        DO ic = icBeg,icEnd
          DO iv = ivBeg,ivEnd
            aivOld(iv,ic) = aiv(iv,ic)
          END DO ! iv
        END DO ! ic
      END IF ! iStage
    CASE ( RK_SCHEME_3_WRAY ) 
      DO ic = icBeg,icEnd
        DO iv = ivBeg,ivEnd
          aivOld(iv,ic) = aiv(iv,ic)
        END DO ! iv
      END DO ! ic   
    CASE ( RK_SCHEME_1_EULER ) 
      DO ic = icBeg,icEnd
        DO iv = ivBeg,ivEnd
          aivOld(iv,ic) = aiv(iv,ic)
        END DO ! iv
      END DO ! ic   
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%rkScheme

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RkInitGenericInt





!******************************************************************************
!
! Purpose: Initialize Runge-Kutta scheme for Real variables.
!
! Description: None.
!
! Input: 
!   region      Region data
!   iStage      Runge-Kutta stage
!   icBeg       Beginning index for cell update
!   icEnd       Ending index for cell update
!   ivBeg       Beginning index for variable update
!   ivEnd       Ending index for variable update
!   rv          Real variables
!   rvOld       Old real variables
!
! Output: 
!   rvOld       Old real variables
!
! Notes: None.
!
!******************************************************************************

SUBROUTINE PLAG_RkInitGenericReal(region,iStage,icBeg,icEnd,ivBeg,ivEnd,&
                                  rv,rvOld)

! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: icBeg,icEnd,iStage,ivBeg,ivEnd
  REAL(RFREAL), DIMENSION(:,:), POINTER :: rv,rvOld
  TYPE(t_region) :: region

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: ic,iv
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ModRkInit.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction(global,'PLAG_RkInitGenericReal',__FILE__)

! *****************************************************************************
! Initialize Runge-Kutta scheme
! *****************************************************************************

  SELECT CASE ( global%rkScheme ) 
    CASE ( RK_SCHEME_4_CLASSICAL ) 
      IF ( iStage == 1 ) THEN
        DO ic = icBeg,icEnd
          DO iv = ivBeg,ivEnd
            rvOld(iv,ic) = rv(iv,ic)
          END DO ! iv
        END DO ! ic
      END IF ! iStage
    CASE ( RK_SCHEME_3_WRAY ) 
      DO ic = icBeg,icEnd
        DO iv = ivBeg,ivEnd
          rvOld(iv,ic) = rv(iv,ic)
        END DO ! iv
      END DO ! ic   
    CASE ( RK_SCHEME_1_EULER ) 
      DO ic = icBeg,icEnd
        DO iv = ivBeg,ivEnd
          rvOld(iv,ic) = rv(iv,ic)
        END DO ! iv
      END DO ! ic   
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%rkScheme

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RkInitGenericReal





!******************************************************************************
!
! Purpose: Driver that initializes primary variables.
!
! Description: none.
!
! Input: istage = RK stage
!        region = data of current region.
!
! Output: region%levels%plag%cvOld 
!         region%levels%plag%aivOld
!         region%levels%plag%arvOld
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE PLAG_RkInitPrimary( region, iStage )

! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: iStage
  TYPE(t_region), INTENT(INOUT), TARGET :: region

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: nAiv,nArv,nCv,nPcls
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv, pAivOld

  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv, pArvOld, pCv, pCvOld
  
  TYPE(t_plag),   POINTER :: pPlag  
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ModRkInit.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction(global,'PLAG_RkInitPrimary',__FILE__)

!******************************************************************************
! Set pointers
!******************************************************************************

  pPlag   => region%plag

  pAiv    => pPlag%aiv
  pArv    => pPlag%arv   
  pCv     => pPlag%cv
  
  pAivOld => pPlag%aivOld
  pArvOld => pPlag%arvOld
  pCvOld  => pPlag%cvOld
  
!******************************************************************************
! Get dimensions
!******************************************************************************

  nPcls = region%plag%nPcls
  
  nAiv = pPlag%nAiv
  nArv = pPlag%nArv           
  nCv  = pPlag%nCv

!******************************************************************************
! Initialize previous solution
!******************************************************************************
  
  CALL PLAG_RkInitGenericInt(  region,iStage,1,nPcls,1,nAiv,pAiv,pAivOld )
  CALL PLAG_RkInitGenericReal( region,iStage,1,nPcls,1,nArv,pArv,pArvOld )
  CALL PLAG_RkInitGenericReal( region,iStage,1,nPcls,1,nCv ,pCv ,pCvOld  )

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RkInitPrimary







!******************************************************************************
!
! Purpose: Driver that initializes variables for tiles.
!
! Description: none.
!
! Input: region = current region.
!        iStage = current RK stage.
!
! Output: region%levels%tilePlag%cvOld 
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE PLAG_InjcTileRkInit( region, iStage )


! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: iStage
  TYPE(t_region), INTENT(INOUT), TARGET :: region

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: bcType,iPatch,nCv,nPatches,nTiles
  
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv, pCvOld

  TYPE(t_patch),     POINTER :: pPatch
  TYPE(t_tile_plag), POINTER :: pTilePlag  
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ModRkInit.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction(global,'PLAG_InjcTileRkInit',__FILE__)
  
!****************************************************************************** 
! Get dimensions
!****************************************************************************** 

  nPatches = region%grid%nPatches

!******************************************************************************     
! Loop over patches
!****************************************************************************** 

  DO iPatch=1,nPatches

    pPatch => region%patches(iPatch)

    bcType = pPatch%bcType

! =============================================================================
!   Select injection boundary condition
! =============================================================================

    IF ( bcType == BC_INJECTION .OR. &
         bcType == BC_INFLOW    .OR. &
	 bcType == BC_INFLOW_VELTEMP ) THEN 

!  ----------------------------------------------------------------------------   
!     Get tile dimensions and set pointers
!  ----------------------------------------------------------------------------   

      nTiles  = pPatch%nBFaces
      
      pTilePlag => pPatch%tilePlag
      
      nCv  = pTilePlag%nCv

      pCv     => pTilePlag%cv
      pCvOld  => pTilePlag%cvOld

!  ----------------------------------------------------------------------------  
!     Initialize previous solution
!  ---------------------------------------------------------------------------- 
      
      CALL PLAG_RkInitGenericReal( region,iStage,1,nTiles,1,nCv,pCv,pCvOld  )

    END IF !bcType
          
  END DO ! iPatch

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_InjcTileRkInit

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE PLAG_ModRkInit

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModRkInit.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.5  2008/12/06 08:43:51  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/08/07 21:54:01  fnajjar
! Removed statements with BC_RANGE since obsolete for RocfluMP
!
! Revision 1.2  2007/04/16 23:21:41  fnajjar
! Deleted all RFLO-pertinent calls and statements
!
! Revision 1.1  2007/04/09 18:50:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.5  2006/09/18 20:30:18  fnajjar
! Activated tile datastructure for inflow bc
!
! Revision 1.4  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.3  2005/05/26 23:02:05  fnajjar
! Bug fix in setting iLev before defining pPlag pointer
!
! Revision 1.2  2005/05/23 18:42:15  fnajjar
! Bug fix to define pointers before setting dimensions
!
! Revision 1.1  2005/05/19 16:02:43  fnajjar
! Initial import
!
!******************************************************************************

