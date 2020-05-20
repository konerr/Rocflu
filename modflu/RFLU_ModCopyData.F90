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
! Purpose: Suite of routines to copy data from and to partitioned regions.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModCopyData.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModCopyData

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGrid, ONLY: t_grid  
  USE ModMPI

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: RFLU_COPY_CellDataP2S_R2D, & 
            RFLU_COPY_CellDataP2S_R3D, &
            RFLU_COPY_CellDataS2P_I1D, &  
            RFLU_COPY_CellDataS2P_R2D

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModCopyData.F90,v $ $Revision: 1.1.1.1 $'


! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS







! ******************************************************************************
!
! Purpose: Copy real cell data from a partitioned region to serial region for 
!   2d arrays.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pGrid       Pointer to grid of partitioned region
!   var         Data on partitioned region
!   varSerial   Data on serial region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_COPY_CellDataP2S_R2D(global,pGrid,var,varSerial)

    IMPLICIT NONE
      
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), DIMENSION(:,:) :: var,varSerial
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: icg,icg2,iVar,iVarBeg,iVarEnd
      
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_COPY_CellDataP2S_R2D',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    iVarBeg = LBOUND(var,1)
    iVarEnd = UBOUND(var,1)    

    IF ( (iVarBeg /= LBOUND(varSerial,1)) .OR. & 
         (iVarEnd /= UBOUND(varSerial,1)) ) THEN 
      CALL ErrorStop(global,ERR_LUBOUND_MISMATCH,__LINE__)
    END IF ! iVarBeg
  
! TEMPORARY: Manoj: 2012-06-04: Looping over only real cells as virtual cell
!                               data is not correct because of no virtual
!                               particles.
!    DO icg = 1,pGrid%nCellsTot
    DO icg = 1,pGrid%nCells
      icg2 = pGrid%pc2sc(icg)
        
      DO iVar = iVarBeg,iVarEnd
        varSerial(iVar,icg2) = var(iVar,icg)
      END DO ! iVar
    END DO ! icg

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_COPY_CellDataP2S_R2D








! ******************************************************************************
!
! Purpose: Copy real cell data from a partitioned region to serial region for 
!   2d arrays.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pGrid       Pointer to grid of partitioned region
!   var         Data on partitioned region
!   varSerial   Data on serial region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_COPY_CellDataP2S_R3D(global,pGrid,var,varSerial)

    IMPLICIT NONE
      
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), DIMENSION(:,:,:) :: var,varSerial
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: icg,icg2,iCmp,iCmpBeg,iCmpEnd,iVar,iVarBeg,iVarEnd
      
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_COPY_CellDataP2S_R3D',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    iCmpBeg = LBOUND(var,1)
    iCmpEnd = UBOUND(var,1)

    iVarBeg = LBOUND(var,2)
    iVarEnd = UBOUND(var,2)    

    IF ( (iCmpBeg /= LBOUND(varSerial,1)) .OR. & 
         (iCmpEnd /= UBOUND(varSerial,1)) ) THEN 
      CALL ErrorStop(global,ERR_LUBOUND_MISMATCH,__LINE__)
    END IF ! iCmpBeg
    
    IF ( (iVarBeg /= LBOUND(varSerial,2)) .OR. & 
         (iVarEnd /= UBOUND(varSerial,2)) ) THEN 
      CALL ErrorStop(global,ERR_LUBOUND_MISMATCH,__LINE__)
    END IF ! iVarBeg
  
! TEMPORARY: Manoj: 2012-06-04: Looping over only real cells as virtual cell
!                               data is not correct because of no virtual
!                               particles.
!    DO icg = 1,pGrid%nCellsTot
    DO icg = 1,pGrid%nCells
      icg2 = pGrid%pc2sc(icg)
        
      DO iCmp = iCmpBeg,iCmpEnd  
        DO iVar = iVarBeg,iVarEnd
          varSerial(iCmp,iVar,icg2) = var(iCmp,iVar,icg)
        END DO ! iVar
      END DO ! iCmp
    END DO ! icg

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_COPY_CellDataP2S_R3D






! ******************************************************************************
!
! Purpose: Copy integer cell data from serial region to a partitioned region 
!   for 1d arrays.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pGrid       Pointer to grid of partitioned region
!   var         Data on partitioned region
!   varSerial   Data on serial region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_COPY_CellDataS2P_I1D(global,pGrid,var,varSerial)

    IMPLICIT NONE
      
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, DIMENSION(:) :: var,varSerial
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: icg,icg2
      
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_COPY_CellDataS2P_I1D',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************
  
    DO icg = 1,pGrid%nCellsTot
      icg2 = pGrid%pc2sc(icg)
        
      var(icg) = varSerial(icg2)
    END DO ! icg

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_COPY_CellDataS2P_I1D








! ******************************************************************************
!
! Purpose: Copy real cell data from serial region to a partitioned region for 
!   2d arrays.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pGrid       Pointer to grid of partitioned region
!   var         Data on partitioned region
!   varSerial   Data on serial region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_COPY_CellDataS2P_R2D(global,pGrid,var,varSerial)

    IMPLICIT NONE
      
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), DIMENSION(:,:) :: var,varSerial
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: icg,icg2,iVar,iVarBeg,iVarEnd
      
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_COPY_CellDataS2P_R2D',__FILE__)

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    iVarBeg = LBOUND(var,1)
    iVarEnd = UBOUND(var,1)    

    IF ( (iVarBeg /= LBOUND(varSerial,1)) .OR. & 
         (iVarEnd /= UBOUND(varSerial,1)) ) THEN 
      CALL ErrorStop(global,ERR_LUBOUND_MISMATCH,__LINE__)
    END IF ! iVarBeg
  
    DO icg = 1,pGrid%nCellsTot
      icg2 = pGrid%pc2sc(icg)
        
      DO iVar = iVarBeg,iVarEnd
        var(iVar,icg) = varSerial(iVar,icg2)
      END DO ! iVar
    END DO ! icg

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_COPY_CellDataS2P_R2D







! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModCopyData


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModCopyData.F90,v $
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
! Revision 1.1  2007/04/09 18:00:40  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.5  2006/12/18 16:18:44  haselbac
! Bug fix: Merging of cases with sype patches needs nCellsTot copied
!
! Revision 1.4  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.3  2005/11/27 01:51:00  haselbac
! Added routine RFLU_COPY_CellDataP2S_R3D
!
! Revision 1.2  2005/08/19 02:33:18  haselbac
! Renamed routines and added new routine
!
! Revision 1.1  2005/04/15 15:06:42  haselbac
! Initial revision
!
! ******************************************************************************

