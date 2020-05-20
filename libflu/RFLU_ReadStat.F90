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
! Purpose: read in time averaged solution
!
! Description: file to be read in contains statistics solution in 
!              interior cells; only binary format is supported
!
! Input: region      = region data (dimensions)
!        integrTime  = integrated time from file
!        mixttav     = time averaged mixture data from file
!        turbtav     = time averaged TURB data from file
!
! Output: region%mixt%tav     = time averaged mixture variables
!         region%turb%tav     = time averaged TURB variables
!         global%integrTime   = integrated time
!
! Notes: each region has statistics solution file
!
!******************************************************************************
!
! $Id: RFLU_ReadStat.F90,v 1.2 2015/07/23 23:11:18 brollin Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_ReadStat( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: ijk, ind, l

! ... local variables
  CHARACTER(CHRLEN+23) :: fname
  CHARACTER(CHRLEN)    :: RCSIdentString,msg,timeString1,timeString2

  INTEGER :: errorFlag,iReg, ijkbeg, ijkend, nCells
  INTEGER :: nTavgVar
  INTEGER :: mixtVarId(2,region%global%mixtNStat)
  INTEGER :: turbVarId(2,region%global%turbNStat)

  REAL(RFREAL)          :: currentTime
  REAL(RFREAL), POINTER :: mixttav(:,:), turbtav(:,:) 

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ReadStat.F90,v $'

  global => region%global

  CALL RegisterFunction(global,'RFLU_ReadStat',__FILE__)

! open file -----------------------------------------------------------------
  
  iReg = region%iRegionGlobal
  WRITE(fname,'(A,I5.5,A,1PE11.5)') TRIM(global%outDir)// & 
                                    TRIM(global%casename)//'.statb_',iReg, &
                                    '_',global%currentTime

!  OPEN(IF_STAT,file=fname,form='unformatted',status='old',iostat=errorFlag)
! BBR - begin
    IF( global%solutFormat .EQ. FORMAT_BINARY )THEN
    OPEN(IF_STAT,file=fname,form='unformatted',status='old',iostat=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_L )THEN
    OPEN(IF_STAT,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
         ACCESS="SEQUENTIAL",CONVERT="LITTLE_ENDIAN",IOSTAT=errorFlag)
    ELSEIF( global%solutFormat .EQ. FORMAT_BINARY_B )THEN
    OPEN(IF_STAT,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
         ACCESS="SEQUENTIAL",CONVERT="BIG_ENDIAN",IOSTAT=errorFlag)
    END IF
! BBR - end 
  global%error = errorFlag   
  IF (global%error /= 0) &
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname))

! read physical time and integrated time ------------------------------------

  READ(IF_STAT,err=10,end=10) currentTime
  READ(IF_STAT,err=10,end=10) global%integrTime
  WRITE(timeString1,'(1PE11.5)') global%currentTime
  WRITE(timeString2,'(1PE11.5)') currentTime          

  IF ( TRIM(timeString1) /= TRIM(timeString2) ) THEN
    WRITE(msg,1010) iReg,currentTime,global%currentTime
    CALL ErrorStop(global,ERR_TIME_SOLUTION,__LINE__, & 
                   msg//' File: '//TRIM(fname) )
  ENDIF

! read dimensions & check them ----------------------------------------------

  READ(IF_STAT,err=10,end=10) nCells

  IF (nCells /= region%grid%nCells) THEN
    WRITE(msg,1000) iReg, nCells, region%grid%nCells
    CALL ErrorStop(global,ERR_GRID_DIMENSIONS,__LINE__,msg)
  ENDIF

! mixture NSTAT and ID

  IF (global%mixtNStat > 0) THEN
    READ(IF_STAT,err=10,end=10) nTavgVar,mixtVarId(1,:)
    mixtVarId(2,:) = mod(mixtVarId(1,:),10)
    mixtVarId(1,:) = (mixtVarId(1,:)-mixtVarId(2,:))/10

    IF (nTavgVar /= global%mixtNStat) THEN
      CALL ErrorStop(global,ERR_STATS_RESTART,__LINE__)
    ENDIF

    DO ind=1,2
      DO l=1,global%mixtNStat
        IF (mixtVarId(ind,l)  /= global%mixtStatId(ind,l)) THEN
          CALL ErrorStop(global,ERR_STATS_RESTART,__LINE__)
        ENDIF
      ENDDO
    ENDDO
  ENDIF

! turbulence NSTAT and ID

#ifdef TURB
  IF (global%turbNStat > 0) THEN
    READ(IF_STAT,err=10,end=10) nTavgVar,turbVarId(1,:)
    turbVarId(2,:) = mod(turbVarId(1,:),10)
    turbVarId(1,:) = (turbVarId(1,:)-turbVarId(2,:))/10

    IF (nTavgVar /= global%turbNStat) THEN
      CALL ErrorStop(global,ERR_STATS_RESTART,__LINE__)
    ENDIF

    DO ind=1,2
      DO l=1,global%turbNStat
        IF (turbVarId(ind,l)  /= global%turbStatId(ind,l)) THEN
          CALL ErrorStop(global,ERR_STATS_RESTART,__LINE__)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
#endif

! read time averaged variables ----------------------------------------------

  IF (global%myProcid==MASTERPROC .AND. global%verbLevel==VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME,'  - read statistics'

  ijkbeg = 1
  ijkend = region%grid%nCells

  IF (global%mixtNStat > 0) THEN
    mixttav => region%mixt%tav
    READ(IF_STAT,err=10,end=10) ((mixttav(l,ijk), ijk=ijkbeg,ijkend), &
                                  l=1,global%mixtNStat)
  ENDIF
#ifdef TURB
  IF (global%turbNStat > 0) THEN
    turbtav => region%turb%tav
    READ(IF_STAT,err=10,end=10) ((turbtav(l,ijk), ijk=ijkbeg,ijkend), &
                                  l=1,global%turbNStat)
  ENDIF
#endif


! close file ----------------------------------------------------------------

  CLOSE(IF_STAT,iostat=errorFlag)
  global%error = errorFlag   
  IF (global%error /= 0) &
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname))

! finalization & error handling ---------------------------------------------

  CALL DeregisterFunction(global)
  GOTO 999

10   CONTINUE
  CALL ErrorStop(global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname))

999  CONTINUE
1000 FORMAT('Region ',I5,', nCells= ',I10,', nCellExpected= ',I10)
1010 FORMAT('Region ',I5,', time in file is= ',1PE12.5, &
            ' but it should be= ',E12.5,'.')

END SUBROUTINE RFLU_ReadStat

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ReadStat.F90,v $
! Revision 1.2  2015/07/23 23:11:18  brollin
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
! Revision 1.3  2008/12/06 08:43:36  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:53  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.9  2006/01/02 11:22:09  wasistho
! changed timeStamp to currentTime
!
! Revision 1.8  2002/11/02 01:51:40  wasistho
! Added TURB statistics
!
! Revision 1.7  2002/10/08 15:48:57  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.6  2002/10/05 18:52:34  haselbac
! Cosmetic changes only
!
! Revision 1.5  2002/09/09 14:15:01  haselbac
! global now under regions
!
! Revision 1.4  2002/07/22 15:45:50  wasistho
! Cleaned-up conforming Coding Rule
!
! Revision 1.3  2002/06/18 00:37:07  wasistho
! Added prefix SOLVER NAME to satistics STDOutput
!
! Revision 1.2  2002/06/17 18:33:34  wasistho
! modified times matching check
!
! Revision 1.1  2002/06/14 21:23:56  wasistho
! Added time avg statistics
!
!
!******************************************************************************

