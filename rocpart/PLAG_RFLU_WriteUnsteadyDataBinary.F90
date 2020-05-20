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
! Purpose: Write unsteady data file for particles in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: PLAG_RFLU_WriteUnsteadyDataBinary.F90,v 1.2 2015/07/23 23:11:19 brollin Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_WriteUnsteadyDataBinary(pRegion)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region 
  USE ModPartLag, ONLY: t_plag,t_tile_plag     
  USE ModMPI

  USE PLAG_ModParameters

  USE ModBuildFileNames, ONLY: BuildFileNameUnsteady

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName,sectionString,RCSIdentString
  INTEGER :: errorFlag,i,iCont,iFile,ifl,iMass,iPatch,iVar,j,nCont,nData,nVars 
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_plag), POINTER :: pPlag
  TYPE(t_tile_plag), POINTER :: pTilePlag  
  
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion  
  
! ******************************************************************************
! Start, open file
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_WriteUnsteadyDataBinary.F90,v $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_WriteUnsteadyDataBinary',__FILE__)

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary particle unsteady data file...'
  END IF ! global%verbLevel

  CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.plag_unsd', & 
                             pRegion%iRegionGlobal,global%currentTime, & 
                             iFileName)

  iFile = IF_SOLUT
!  OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN",IOSTAT=errorFlag)
! BBR - begin
    IF( global%solutFormat .EQ. FORMAT_BINARY )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
         IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_L )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
         ACCESS="SEQUENTIAL",CONVERT="LITTLE_ENDIAN",IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_B )THEN
    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
         ACCESS="SEQUENTIAL",CONVERT="BIG_ENDIAN",IOSTAT=errorFlag)
    END IF
! BBR - end 
  global%error = errorFlag          
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
  END IF ! global%error  

! ==============================================================================
! Header and general information
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
  END IF ! global%verbLevel

  sectionString = '# ROCFLU particle unsteady data file'
  WRITE(iFile) sectionString
  
  sectionString = '# Precision and range'
  WRITE(iFile) sectionString
  WRITE(iFile) PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)
  
  sectionString = '# Physical time'
  WRITE(iFile) sectionString
  WRITE(iFile) global%currentTime 

! ==============================================================================
! Dimensions
! ==============================================================================
  
  pGrid => pRegion%grid  
  pPlag => pRegion%plag
   
  nCont = pRegion%plagInput%nCont   
  
  nVars = 7
  nData = pRegion%plagInput%nTimeBH
   
  sectionString = '# Dimensions'
  WRITE(iFile) sectionString
  WRITE(iFile) pPlag%nPcls,nVars,nData
  
! ==============================================================================
! State vector
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Conserved variables...'
  END IF ! global%verbLevel

  sectionString = '# Relevant time history'
  WRITE(iFile) sectionString
  WRITE(iFile) (pPlag%timeBH(i),i=1,nData)

  sectionString = '# Particle x-velocity substantial derivative'
  WRITE(iFile) sectionString
  WRITE(iFile) ((pPlag%dudtPlag(XCOORD,i,j),i=1,nData),j=1,pPlag%nPcls)

  sectionString = '# Particle y-velocity substantial derivative'
  WRITE(iFile) sectionString
  WRITE(iFile) ((pPlag%dudtPlag(YCOORD,i,j),i=1,nData),j=1,pPlag%nPcls)

  sectionString = '# Particle z-velocity substantial derivative'
  WRITE(iFile) sectionString
  WRITE(iFile) ((pPlag%dudtPlag(ZCOORD,i,j),i=1,nData),j=1,pPlag%nPcls)

  sectionString = '# Fluid x-velocity substantial derivative'
  WRITE(iFile) sectionString
  WRITE(iFile) ((pPlag%dudtMixt(XCOORD,i,j),i=1,nData),j=1,pPlag%nPcls)

  sectionString = '# Fluid y-velocity substantial derivative'
  WRITE(iFile) sectionString
  WRITE(iFile) ((pPlag%dudtMixt(YCOORD,i,j),i=1,nData),j=1,pPlag%nPcls)

  sectionString = '# Fluid z-velocity substantial derivative'
  WRITE(iFile) sectionString
  WRITE(iFile) ((pPlag%dudtMixt(ZCOORD,i,j),i=1,nData),j=1,pPlag%nPcls)

! ==============================================================================
! End marker
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
  END IF ! global%verbLevel

  sectionString = '# End'
  WRITE(iFile) sectionString  

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
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary particle unsteady data file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)
    
! ******************************************************************************
! End
! ******************************************************************************
 
END SUBROUTINE PLAG_RFLU_WriteUnsteadyDataBinary


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_WriteUnsteadyDataBinary.F90,v $
! Revision 1.2  2015/07/23 23:11:19  brollin
! 1) The pressure coefficient of the  collision model has been changed back to its original form
! 2) New options in the format of the grid and solutions have been added. Now the user can choose the endianness, and convert from one to the over in rfluconv.
! 3) The solutions are now stored in folders named by timestamp or iteration number.
! 4) The address enty in the hashtable has been changed to an integer(8) for cases when the grid becomes very large.
! 5) RFLU_WritePM can now compute PM2 on the fly for the Macroscale problem
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
!
! ******************************************************************************

