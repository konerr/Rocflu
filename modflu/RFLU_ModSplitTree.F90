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
! Purpose: Suite of routines to carry out split tree operations.
!
! Description: None.
!
! Notes: To create and use a split tree, one has to take the following steps:
!   1. 
!   2. 
!   3. Deallocate the table by calling DestroyHashTable.
!
!*******************************************************************************
!
! $Id: RFLU_ModSplitTree.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!*******************************************************************************

MODULE RFLU_ModSplitTree

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModSortSearch

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: RFLU_CreateSplitTree, & 
            RFLU_BuildSplitTree, & 
            RFLU_DestroySplitTree

  SAVE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
  
  INTEGER, PARAMETER :: NLEVELS_MAX  = 50,    & ! Must be greater than 1
                        NBUCKETS_MAX = 10000, & ! Must be greater than 2
                        NPOINTS_MAX  = 5       
  
  INTEGER, PRIVATE :: nPoints
  
  INTEGER, DIMENSION(:), ALLOCATABLE :: pointList     
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: tree     
       
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: pointXyz     
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
! ==============================================================================
!   Create split tree
! ==============================================================================  

    SUBROUTINE RFLU_CreateSplitTree(global,nDataPoints)
    
      IMPLICIT NONE   
        
      INTEGER, INTENT(IN) :: nDataPoints 
      TYPE(t_global), POINTER :: global  
        
      INTEGER :: errorFlag,ip         
        
      CALL RegisterFunction(global,'RFLU_CreateSplitTree',__FILE__)

! ------------------------------------------------------------------------------
!     Copy argument into nPoints variable
! ------------------------------------------------------------------------------ 

      nPoints = nDataPoints

! ------------------------------------------------------------------------------
!     Allocate memory
! ------------------------------------------------------------------------------ 

      ALLOCATE(tree(7,NBUCKETS_MAX),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
      END IF ! global%error)

      ALLOCATE(pointList(nPoints),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
      END IF ! global%error

      ALLOCATE(pointXyz(nPoints),STAT=errorFlag)
      global%error = errorFlag  
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
      END IF ! global%error

! ------------------------------------------------------------------------------
!     Initialize
! ------------------------------------------------------------------------------ 

      tree(:,:) = 0

      DO ip = 1,nPoints
        pointList(ip) = ip
      END DO ! ip
    
      CALL DeregisterFunction(global)    

    END SUBROUTINE RFLU_CreateSplitTree

! ==============================================================================
!   Build split tree
! ==============================================================================  

    SUBROUTINE RFLU_BuildSplitTree(global,xyz)

      IMPLICIT NONE
      
      INTEGER :: errorFlag,iBranch1,iBranch2,iBucket,iBucketLast,ipl,ipm,ipg, & 
                 iSplitDir,nBuckets,nLevels
      
      REAL(RFREAL), INTENT(IN) :: xyz(3,nPoints)
      TYPE(t_global), POINTER :: global
      
! ------------------------------------------------------------------------------
!     Start
! ------------------------------------------------------------------------------      
      
      CALL RegisterFunction(global,'RFLU_BuildSplitTree',__FILE__)      

      nLevels  = 1
      nBuckets = 1
      
      iBucket     = 1
      iBucketLast = 1
      iSplitDir   = 1
      
      tree(1,iBucket) = 1
      tree(2,iBucket) = nPoints
      tree(3,iBucket) = nPoints
      tree(4,iBucket) = iSplitDir
      
! ------------------------------------------------------------------------------
!     Loop over levels
! ------------------------------------------------------------------------------      
      
      OUTER: DO
        iSplitDir = iSplitDir + 1
        IF ( iSplitDir > 3 ) THEN 
          iSplitDir = 1
        END IF ! iSplitDir
      
! ----- Loop over buckets      
      
        INNER: DO 
     
! ------- Determine whether bucket should be split     
           
          IF ( tree(3,iBucket) > NPOINTS_MAX ) THEN ! Split bucket

! --------- Check whether more buckets can be generated

            IF ( nBuckets > NBUCKETS_MAX-2 ) THEN 
              EXIT OUTER
            END IF ! nBuckets
          
! --------- Determine where bucket should be split - use median 
          
            ipl = 0
            DO ipg = tree(1,iBucket),tree(2,iBucket)
              ipl = ipl + 1
              pointXyz(ipl) = xyz(iSplitDir,pointList(ipg))
            END DO ! ipg
            
            CALL QuickSortRFREALInteger(pointXyz(1:tree(3,iBucket)), &
                          pointList(tree(1,iBucket):tree(2,iBucket)), & 
                          tree(3,iBucket))
          
            ipm = tree(3,iBucket)/2 + 1 ! NOTE integer division
          
! --------- Update information for current bucket          
          
            tree(5,iBucket) = pointList(ipm)          
          
            iBranch1 = iBucket + 1
            iBranch2 = iBucket + 2
            
            tree(6,iBucket) = iBranch1
            tree(7,iBucket) = iBranch2          
          
! --------- Create new buckets

            nBuckets = nBuckets + 2
            
            tree(1,iBranch1) = tree(1,iBucket)
            tree(2,iBranch1) = ipm - 1
            tree(3,iBranch1) = tree(2,iBranch1) - tree(1,iBranch1) + 1        
            tree(4,iBranch1) = iSplitDir
          
            tree(1,iBranch2) = ipm
            tree(2,iBranch2) = tree(3,iBucket)
            tree(3,iBranch2) = tree(2,iBranch2) - tree(1,iBranch2) + 1
            tree(4,iBranch2) = tree(4,iBranch1)           
          END IF ! tree(3,iBucket)  
          
! ------- Check whether more buckets to be split          
        
          IF ( iBucket < nBuckets ) THEN 
            IF ( iBucket /= iBucketLast ) THEN
              iBucket = iBucket + 1             
            ELSE 
              iBucketLast = nBuckets
              EXIT INNER
            END IF ! iBucket
          ELSE
            EXIT OUTER  
          END IF ! iBucket
        END DO INNER 
      
! ----- Check whether more levels can be inserted      
      
        IF ( nLevels < NLEVELS_MAX ) THEN 
          nLevels   = nLevels   + 1
          iBucket   = iBucket   + 1
        ELSE  
          EXIT OUTER
        END IF ! nLevels
      END DO OUTER

      WRITE(*,*) nLevels,nBuckets

      DO ipl = 1,nBuckets
        WRITE(*,*) ipl,tree(1:7,ipl)
      END DO ! ipl

      CALL DeregisterFunction(global)  

    END SUBROUTINE RFLU_BuildSplitTree

! ==============================================================================
!   Destroy split tree
! ==============================================================================  
   
    SUBROUTINE RFLU_DestroySplitTree(global)
    
      IMPLICIT NONE

      TYPE(t_global), POINTER :: global

      INTEGER :: errorFlag

      CALL RegisterFunction(global,'RFLU_DestroySplitTree',__FILE__)

      DEALLOCATE(tree,STAT=errorFlag)
      global%error = errorFlag  
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__)
      END IF ! global%error)

      DEALLOCATE(pointList,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__)
      END IF ! global%error

      DEALLOCATE(pointXyz,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__)
      END IF ! global%error

      CALL DeregisterFunction(global)

    END SUBROUTINE RFLU_DestroySplitTree


  

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModSplitTree


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLU_ModSplitTree.F90,v $
!   Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
!   merged rocflu micro and macro
!
!   Revision 1.1.1.1  2014/07/15 14:31:37  brollin
!   New Stable version
!
!   Revision 1.3  2008/12/06 08:43:45  mtcampbe
!   Updated license.
!
!   Revision 1.2  2008/11/19 22:16:58  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.1  2007/04/09 18:49:26  haselbac
!   Initial revision after split from RocfloMP
!
!   Revision 1.1  2007/04/09 18:00:42  haselbac
!   Initial revision after split from RocfloMP
!
!   Revision 1.4  2004/01/22 16:03:59  haselbac
!   Made contents of modules PRIVATE, only procs PUBLIC, to avoid errors on ALC and titan
!
!   Revision 1.3  2002/10/08 15:49:21  haselbac
!   {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
!   Revision 1.2  2002/09/09 15:12:12  haselbac
!   global now under regions
!
!   Revision 1.1  2002/04/11 18:48:48  haselbac
!   Initial revision
!
!
! ******************************************************************************

