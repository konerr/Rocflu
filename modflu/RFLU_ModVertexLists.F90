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
! Purpose: Suite of routines for vertex lists.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModVertexLists.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModVertexLists

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY: t_grid
  USE ModStencil, ONLY: t_stencil
  USE ModSortSearch
  USE ModMPI

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_BuildVert2CellList, &
            RFLU_CreateVert2CellList, &   
            RFLU_DestroyVert2CellList, & 
            RFLU_NullifyVert2CellList
            
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModVertexLists.F90,v $ $Revision: 1.1.1.1 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  




! *******************************************************************************
!
! Purpose: Build vertex-to-cell list.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildVert2CellList(pRegion)

    USE RFLU_ModCellFaceEdgeInfo

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg,icl,ict,ivg,ivl,v2cBeg,v2cDegrSum,v2cEnd,v2cIndx
    INTEGER, DIMENSION(:), ALLOCATABLE :: v2cDegr
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildVert2CellList',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Building vertex-to-cell list...'            
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set grid pointer and set constants
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate temporary memory
! ******************************************************************************

    ALLOCATE(v2cDegr(pGrid%nVertTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'v2cDegr')
    END IF ! global%error    

    DO ivg = 1,pGrid%nVertTot
      v2cDegr(ivg) = 0
    END DO ! ivg    

! ******************************************************************************
!   Loop over cells and determine degree of each vertex
! ******************************************************************************

    DO icg = 1,pGrid%nCellsTot
      ict = RFLU_GetGlobalCellType(global,pGrid,icg)
      icl = pGrid%cellGlob2Loc(2,icg)

      SELECT CASE ( ict ) 
        CASE ( CELL_TYPE_TET ) 
          DO ivl = 1,4
            ivg = pGrid%tet2v(ivl,icl)              
            v2cDegr(ivg) = v2cDegr(ivg) + 1 
          END DO ! ivl
        CASE ( CELL_TYPE_HEX ) 
          DO ivl = 1,8
            ivg = pGrid%hex2v(ivl,icl)
            v2cDegr(ivg) = v2cDegr(ivg) + 1              
          END DO ! ivl
        CASE ( CELL_TYPE_PRI )           
          DO ivl = 1,6
            ivg = pGrid%pri2v(ivl,icl)
            v2cDegr(ivg) = v2cDegr(ivg) + 1              
          END DO ! ivl                  
        CASE ( CELL_TYPE_PYR ) 
          DO ivl = 1,5
            ivg = pGrid%pyr2v(ivl,icl)
            v2cDegr(ivg) = v2cDegr(ivg) + 1              
          END DO ! ivl            
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)                
      END SELECT ! ict
    END DO ! icg      

! ******************************************************************************
!   Build information array
! ******************************************************************************

    v2cDegrSum = 0      

    DO ivg = 1,pGrid%nVertTot
      IF ( ivg > 1 ) THEN 
        pGrid%v2cInfo(V2C_BEG,ivg) = pGrid%v2cInfo(V2C_END,ivg-1) + 1
      ELSE 
        pGrid%v2cInfo(V2C_BEG,ivg) = 1
      END IF ! ivg

      v2cDegrSum = v2cDegrSum + v2cDegr(ivg)

      pGrid%v2cInfo(V2C_END,ivg) = v2cDegrSum                                                                           
    END DO ! ivg          

! ******************************************************************************
!   Allocate actual v2c array
! ******************************************************************************

    ALLOCATE(pGrid%v2c(v2cDegrSum),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%v2c')
    END IF ! global%error              

! ******************************************************************************
!   Build actual array
! ******************************************************************************

    DO ivg = 1,pGrid%nVertTot ! Reinitialize degree counter     
      v2cDegr(ivg) = 0 
    END DO ! ivg

    DO icg = 1,pGrid%nCellsTot
      ict = RFLU_GetGlobalCellType(global,pGrid,icg)
      icl = pGrid%cellGlob2Loc(2,icg)

      SELECT CASE ( ict ) 
        CASE ( CELL_TYPE_TET ) 
          DO ivl = 1,4
            ivg = pGrid%tet2v(ivl,icl)
            v2cIndx = pGrid%v2cInfo(V2C_BEG,ivg) + v2cDegr(ivg)              
            pGrid%v2c(v2cIndx) = icg 
            v2cDegr(ivg) = v2cDegr(ivg) + 1  
          END DO ! ivl
        CASE ( CELL_TYPE_HEX ) 
          DO ivl = 1,8
            ivg = pGrid%hex2v(ivl,icl)
            v2cIndx = pGrid%v2cInfo(V2C_BEG,ivg) + v2cDegr(ivg)              
            pGrid%v2c(v2cIndx) = icg 
            v2cDegr(ivg) = v2cDegr(ivg) + 1  
          END DO ! ivl
        CASE ( CELL_TYPE_PRI )           
          DO ivl = 1,6
            ivg = pGrid%pri2v(ivl,icl)
            v2cIndx = pGrid%v2cInfo(V2C_BEG,ivg) + v2cDegr(ivg)              
            pGrid%v2c(v2cIndx) = icg 
            v2cDegr(ivg) = v2cDegr(ivg) + 1  
          END DO ! ivl                  
        CASE ( CELL_TYPE_PYR ) 
          DO ivl = 1,5
            ivg = pGrid%pyr2v(ivl,icl)
            v2cIndx = pGrid%v2cInfo(V2C_BEG,ivg) + v2cDegr(ivg)              
            pGrid%v2c(v2cIndx) = icg 
            v2cDegr(ivg) = v2cDegr(ivg) + 1  
          END DO ! ivl            
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)                
      END SELECT ! ict
    END DO ! icg                     

! ******************************************************************************
!   Deallocate temporary memory
! ******************************************************************************

    DEALLOCATE(v2cDegr,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'v2cDegr')
    END IF ! global%error    


#ifdef CHECK_DATASTRUCT
! ******************************************************************************
!   Data structure output for checking
! ******************************************************************************

    WRITE(STDOUT,'(A)') SOLVER_NAME
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Vertex-to-cell list'      

    DO ivg = 1,pGrid%nVertTot
      v2cBeg = pGrid%v2cInfo(V2C_BEG,ivg)
      v2cEnd = pGrid%v2cInfo(V2C_END,ivg)

      WRITE(STDOUT,'(A,1X,I6,1X,I3,3X,20(1X,I6))') SOLVER_NAME,ivg, & 
            v2cEnd-v2cBeg+1,pGrid%v2c(v2cBeg:v2cEnd)        
    END DO ! ivg          
#endif                  

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Building vertex-to-cell list done.'
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_BuildVert2CellList  








! *******************************************************************************
!
! Purpose: Create vertex-to-cell list.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_CreateVert2CellList(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,ic
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateVert2CellList',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Creating vertex-to-cell list...'            
    END IF ! global%verbLevel

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyVert2CellList(pRegion)

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    ALLOCATE(pGrid%v2cInfo(V2C_BEG:V2C_END,pGrid%nVertTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'v2cInfo')
    END IF ! global%error                

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Creating vertex-to-cell list done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_CreateVert2CellList 
 
 
 
 
 


! *******************************************************************************
!
! Purpose: Destroy vertex-to-cell list.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_DestroyVert2CellList(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyVert2CellList',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Destroying vertex-to-cell list...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    DEALLOCATE(pGrid%v2c,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%v2c')
    END IF ! global%error                

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyVert2CellList(pRegion)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Destroying vertex-to-cell list done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyVert2CellList 







! *******************************************************************************
!
! Purpose: Nullify vertex-to-cell list.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_NullifyVert2CellList(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================


    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_NullifyVert2CellList',__FILE__)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Nullifying vertex-to-cell list...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    NULLIFY(pGrid%v2c)
    NULLIFY(pGrid%v2cInfo)   

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Nullifying vertex-to-cell list done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_NullifyVert2CellList 







! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModVertexLists


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModVertexLists.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:36  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:46  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:58  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:27  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:42  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/04/07 15:19:21  haselbac
! Removed tabs
!
! Revision 1.1  2005/10/05 14:33:44  haselbac
! Initial revision
!
! ******************************************************************************

