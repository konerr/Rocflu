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
! Purpose: Read optimal LES integrals in ASCII ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. No checking of section strings.
!
!******************************************************************************
!
! $Id: RFLU_ReadIntegrals1245OLES.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_ReadIntegrals1245OLES(pRegion)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region    
  USE ModMPI

  USE RFLU_ModOLES

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
! Local variables
! ==============================================================================

  INTEGER :: errorFlag,hLoc,i,ic1l,ic2l,ic3l,ic4l,icmp,ifcp,iFile,iFun,j,k,l, & 
             m,nCells,vLoc
  CHARACTER(CHRLEN) :: iFileName,sectionString,RCSIdentString
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global
  
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion  
  
! ******************************************************************************
! Start, set pointer
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ReadIntegrals1245OLES.F90,v $ $Revision: 1.1.1.1 $'

  pGrid => pRegion%grid

! ******************************************************************************
! Open file
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ReadIntegrals1245OLES',__FILE__)
  
  nCells = SIZE(pGrid%fsOLES,1)

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A)') SOLVER_NAME  
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading optimal LES integrals...'
  END IF ! global%verbLevel

  WRITE(iFileName,'(A,I1)') TRIM(global%inDir)// & 
                            TRIM(global%casename)//'.int_',nCells

  iFile = IF_INTEG_OLES
  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD", & 
       IOSTAT=errorFlag)
  global%error = errorFlag        
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
  END IF ! global%error  
      
! ==============================================================================
! Integral 1
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Integral 1...'
  END IF ! global%verbLevel

  READ(iFile,'(A)') sectionString
  
!  DO icmp = 1,3
!    DO ifcp = 1,3
!      DO i = 1,3*nCells
!        READ(iFile,'(E15.6)') pGrid%int1OLES(icmp,ifcp,i)
!      END DO ! i
!    END DO ! ifcp
!  END DO ! icmp
  
  DO ifcp = 1,3
    DO ic1l = 1,nCells
      DO iFun = 1,9
        CALL RFLU_MapK2IJ(iFun,l,i)
        vLoc = RFLU_GetI1PosOLES(l,ic1l)
        
        READ(iFile,'(E15.6)') pGrid%int1OLES(i,ifcp,vLoc)
      END DO ! iFun
    END DO ! ic1l
  END DO ! ifcp  
  
! ==============================================================================
! Integral 2
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Integral 2...'
  END IF ! global%verbLevel

  READ(iFile,'(A)') sectionString
    
!  DO ifcp = 1,3
!    DO i = 1,3*nCells
!      DO j = 1,3*nCells
!        READ(iFile,'(2(E15.6,1X))') pGrid%int20OLES(ifcp,i,j), & 
!                                    pGrid%int21OLES(ifcp,i,j)
!      END DO ! j
!    END DO ! i
!  END DO ! ifcp
  
  DO ifcp = 1,3
    DO ic1l = 1,nCells
      DO ic2l = 1,nCells
        DO iFun = 1,9
          CALL RFLU_MapK2IJ(iFun,l,j)   
          vLoc = RFLU_GetI1PosOLES(l,ic1l)
          hLoc = RFLU_GetLPosOLES(j,ic2l)   
          
          READ(iFile,'(2(E15.6,1X))') pGrid%int20OLES(ifcp,vLoc,hLoc), & 
                                      pGrid%int21OLES(ifcp,vLoc,hLoc)                    
        END DO ! iFun
      END DO ! ic2l    
    END DO ! ic1l    
  END DO ! ifcp  
  
! ==============================================================================
! Integral 4
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Integral 4...'
  END IF ! global%verbLevel

  READ(iFile,'(A)') sectionString
  
!  DO icmp = 1,3
!    DO ifcp = 1,3
!      DO i = 1,9*nCells*nCells
!        READ(iFile,'(3(E15.6,1X))') pGrid%int40OLES(icmp,ifcp,i), & 
!                                    pGrid%int41OLES(icmp,ifcp,i), &
!                                    pGrid%int42OLES(icmp,ifcp,i)
!      END DO ! i
!    END DO ! ifcp
!  END DO ! icmp  
  
  DO ifcp = 1,3
    DO ic1l = 1,nCells
      DO ic2l = 1,nCells
        DO iFun = 1,27
          CALL RFLU_MapL2IJK(iFun,l,m,i)    
          vLoc = RFLU_GetI4PosOLES(l,m,ic1l,ic2l,nCells)   
          
          READ(iFile,'(3(E15.6,1X))') pGrid%int40OLES(i,ifcp,vLoc), & 
                                      pGrid%int41OLES(i,ifcp,vLoc), & 
                                      pGrid%int42OLES(i,ifcp,vLoc)                             
        END DO ! iFun
      END DO ! ic2l    
    END DO ! ic1l    
  END DO ! ifcp    
  
! ==============================================================================
! Integral 5
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Integral 5...'
  END IF ! global%verbLevel

  READ(iFile,'(A)') sectionString
  
!  DO ifcp = 1,3
!    DO i = 1,9*nCells*nCells
!      DO j = 1,9*nCells*nCells
!       READ(iFile,'(3(E15.6,1X))') pGrid%int50OLES(ifcp,i,j), & 
!                                   pGrid%int51OLES(ifcp,i,j), &
!                                   pGrid%int52OLES(ifcp,i,j)
!      END DO ! j
!    END DO ! i
!  END DO ! ifcp    

  DO ifcp = 1,3
    DO ic1l = 1,nCells
      DO ic2l = 1,nCells
        DO ic3l = 1,nCells
          DO ic4l = 1,nCells
            DO iFun = 1,81
              CALL RFLU_MapM2IJKL(iFun,l,m,j,k)
              vLoc = RFLU_GetQPosOLES(j,k,ic3l,ic4l,nCells)
              hLoc = RFLU_GetI4PosOLES(l,m,ic1l,ic2l,nCells) 

              READ(iFile,'(3(E15.6,1X))') pGrid%int50OLES(ifcp,vLoc,hLoc), & 
                                          pGrid%int51OLES(ifcp,vLoc,hLoc), & 
                                          pGrid%int52OLES(ifcp,vLoc,hLoc)                                    
            END DO ! iFun
          END DO ! ic4l
        END DO ! ic3l
      END DO ! ic2l    
    END DO ! ic1l    
  END DO ! ifcp

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(iFile,IOSTAT=errorFlag)
  global%error = errorFlag   
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
  END IF ! global%error
   
  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading optimal LES integrals done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)
    
! ******************************************************************************
! End
! ******************************************************************************
 
END SUBROUTINE RFLU_ReadIntegrals1245OLES


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLU_ReadIntegrals1245OLES.F90,v $
!   Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
!   merged rocflu micro and macro
!
!   Revision 1.1.1.1  2014/07/15 14:31:38  brollin
!   New Stable version
!
!   Revision 1.3  2008/12/06 08:43:48  mtcampbe
!   Updated license.
!
!   Revision 1.2  2008/11/19 22:17:00  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.1  2007/04/09 18:49:57  haselbac
!   Initial revision after split from RocfloMP
!
!   Revision 1.1  2007/04/09 18:01:01  haselbac
!   Initial revision after split from RocfloMP
!
!   Revision 1.4  2004/01/22 16:04:33  haselbac
!   Changed declaration to eliminate warning on ALC
!
!   Revision 1.3  2002/10/08 15:49:29  haselbac
!   {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
!   Revision 1.2  2002/10/05 19:30:30  haselbac
!   GENX integration: added outDir to file name
!
!   Revision 1.1  2002/09/09 16:28:02  haselbac
!   Initial revision
!
! ******************************************************************************

