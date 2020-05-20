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
! Purpose: Update the particle datastructure by removing particles flagged for
!   deletion.
!
! Description: None.
!
! Input: 
!  pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: PLAG_UpdateDataStruct.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_UpdateDataStruct(pRegion)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModPartLag, ONLY: t_plag
  USE ModMPI
 
  USE PLAG_ModParameters    
   
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
! =============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iReg

  INTEGER :: iAiv,iArv,iCv,iDv,iGap,iPcl,iTv,nAiv,nArv,nCv,nDv,&
             nPclsBeg,nPclsEnd,nPclsPrev,nTv
  INTEGER, DIMENSION(:,:), POINTER :: pAiv,pAivOld
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pArv,pArvOld,pCv,pCvOld,pDv,&
                                           pRhs,pRhsSum,pTv
  TYPE(t_global), POINTER :: global
  TYPE(t_plag),   POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_UpdateDataStruct.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_UpdateDataStruct',__FILE__)
  
! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pPlag => pRegion%plag

  pAiv  => pPlag%aiv
  pArv  => pPlag%arv
  pCv   => pPlag%cv
  pDv   => pPlag%dv
  pRhs  => pPlag%rhs
  pRhsSum => pPlag%rhsSum
  pTv   => pPlag%tv

  pAivOld => pPlag%aivOld
  pArvOld => pPlag%arvOld
  pCvOld  => pPlag%cvOld

! ******************************************************************************
! Get dimensions
! ******************************************************************************

  nAiv = pPlag%nAiv
  nArv = pPlag%nArv
  nCv  = pPlag%nCv
  nDv  = pPlag%nDv
  nTv  = pPlag%nTv

  iReg  = pRegion%iRegionGlobal

! ******************************************************************************
! Loop over particles
! ******************************************************************************

  iGap = 1
  iPcl = pPlag%nPcls 

  searchLoop: DO

! ==============================================================================  
!   Copy data if DELETE status in iGap and KEEP status in iPcl
! ==============================================================================  

    IF ( ( pAiv(AIV_PLAG_STATUS,iGap) /= PLAG_STATUS_KEEP   .AND. &
           pAiv(AIV_PLAG_STATUS,iGap) /= PLAG_STATUS_COMM ) .AND. &
         ( pAiv(AIV_PLAG_STATUS,iPcl) == PLAG_STATUS_KEEP   .OR.  &
           pAiv(AIV_PLAG_STATUS,iPcl) == PLAG_STATUS_COMM )       ) THEN

! ------------------------------------------------------------------------------
!      Copy data and set status to keep
! ------------------------------------------------------------------------------       
   
       DO iAiv = 1, nAiv
         pAiv(iAiv,iGap)    = pAiv(iAiv,iPcl)
         pAivOld(iAiv,iGap) = pAivOld(iAiv,iPcl)
       END DO ! iAiv
    
       DO iArv = 1, nArv
         pArv(iArv,iGap)    = pArv(iArv,iPcl)
         pArvOld(iArv,iGap) = pArvOld(iArv,iPcl)
       END DO ! iArv
  
       DO iCv = 1, nCv
         pCv(iCv,iGap)     = pCv(iCv,iPcl)
         pCvOld( iCv,iGap) = pCvOld(iCv,iPcl)
         pRhs(iCv,iGap)    = pRhs(iCv,iPcl)
         pRhsSum(iCv,iGap) = pRhsSum(iCv,iPcl)   
       END DO ! iCv

       DO iDv = 1, nDv
         pDv(iDv,iGap)     = pDv(iDv,iPcl)
       END DO ! iDv
  
       DO iTv = 1, nTv
         pTv(iTv,iGap)     = pTv(iTv,iPcl)
       END DO ! iDv
  
       pAiv(AIV_PLAG_STATUS,iGap) = pAiv(AIV_PLAG_STATUS,iPcl)
   
! ------------------------------------------------------------------------------
!      Increment/decrement counters
! ------------------------------------------------------------------------------    
    
       iGap = iGap + 1
       iPcl = iPcl - 1
       
! ==============================================================================
!   Update counters if KEEP or COMM status in iGap and DELETE status in iPcl
! ==============================================================================       
       
    ELSE
      IF ( pAiv(AIV_PLAG_STATUS,iGap) == PLAG_STATUS_KEEP .OR. &
           pAiv(AIV_PLAG_STATUS,iGap) == PLAG_STATUS_COMM      ) THEN        
        iGap = iGap + 1 
      END IF ! pAiv
      
      IF ( pAiv(AIV_PLAG_STATUS,iPcl) == PLAG_STATUS_DELETE ) THEN 
        iPcl = iPcl - 1
      END IF ! pAiv 
    END IF ! pAiv(AIV_PLAG_STATUS,iGap)

! ==============================================================================  
!   Exit update
! ==============================================================================  

    IF ( iGap > iPcl ) THEN
      EXIT searchLoop
    END IF ! iGap
  END DO searchLoop

! ******************************************************************************
! Update number of particles
! ******************************************************************************
 
  nPclsPrev   = pPlag%nPcls  
  pPlag%nPcls = iPcl

  nPclsBeg = MAX(1,pPlag%nPcls+1)
  nPclsEnd = nPclsPrev

! ******************************************************************************
! Initialize datastructure to crazy for particles between nPclsBeg and nPclsEnd
! ******************************************************************************

  DO iPcl = nPclsBeg, nPclsEnd
    DO iAiv = 1, nAiv
      pAiv(iAiv,iPcl)    = CRAZY_VALUE_INT 
      pAivOld(iAiv,iPcl) = CRAZY_VALUE_INT 
    END DO ! iAiv

    DO iArv = 1, nArv
      pArv(iArv,iPcl)    = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      pArvOld(iArv,iPcl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    END DO ! iArv

    DO iCv = 1, nCv
      pCv(iCv,iPcl)     = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      pCvOld(iCv,iPcl)  = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      pRhs(iCv,iPcl)    = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      pRhsSum(iCv,iPcl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    END DO ! iCv

    DO iDv = 1, nDv
      pDv(iDv,iPcl)     = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    END DO ! iDv

    DO iTv = 1, nTv
      pTv(iTv,iPcl)     = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    END DO ! iTv
  END DO ! iPcl 

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_UpdateDataStruct

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_UpdateDataStruct.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:52  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/04/16 23:22:44  fnajjar
! Deleted all RFLO-pertinent calls and statements
!
! Revision 1.1  2007/04/09 18:50:27  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:34  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.10  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.9  2005/12/30 16:29:59  fnajjar
! Added dv and tv in data reshuffle for consistency in GatherSurfStats
!
! Revision 1.8  2005/05/18 22:22:33  fnajjar
! Adapted to parallelization within RFLU
!
! Revision 1.7  2004/07/01 23:05:11  fnajjar
! Added aivOld to data shuffle kernel to inline with particle trajectory motion
!
! Revision 1.6  2004/05/22 00:15:16  fnajjar
! Bug fixes for proper definition of nPclsBeg and nPclsEnd
!
! Revision 1.5  2004/05/05 21:47:44  fnajjar
! Bug Fix: updating cvOld and arvOld in datastructure
!
! Revision 1.4  2004/05/05 21:20:38  fnajjar
! Clean up and simplification
!
! Revision 1.3  2004/04/09 22:58:59  fnajjar
! Made routine RFLO-aware with ifdef and included bug fix for datastructure reinitialization
!
! Revision 1.2  2004/04/07 20:26:07  fnajjar
! Redesigned kernel to account for degeneracy case of last particle to be removed
!
! Revision 1.1  2004/03/26 21:33:48  fnajjar
! Initial import for plag datastructure update
!
! ******************************************************************************

