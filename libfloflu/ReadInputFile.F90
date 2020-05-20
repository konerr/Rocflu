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
! Purpose: Read in user input (done on all processors).
!
! Description: None.
!
! Input: 
!   regions     Region data
!
! Output: 
!   None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ReadInputFile.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReadInputFile(regions)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  
  USE ModBuildFileNames, ONLY: BuildFileNamePlain  
  
  USE ModInterfaces, ONLY: ReadFormatsSection, &
                           ReadReferenceSection, &
                           ReadFlowSection, & 
                           ReadProbeSection, & 
                           ReadForcesSection, &
                           ReadTimestepSection, &
                           ReadMiscSection, &
                           ReadMixtureSection, & 
                           ReadMultigridSection, &
                           ReadABCSection, &
                           ReadGFMSection, &
                           ReadMvFrameSection, &
                           ReadNumericsSection, & 
                           ReadTransformSection, &
                           ReadViscositySection, &
                           ReadThrustSection, &
                           ReadInitFlowSection, &
                           ReadGridMotionSection, &
                           ReadAccelerationSection, &
                           ReadRandomSection, &
                           ReadPostSection, &
                           ReadPrepSection, &
                           ReadTimeZoomingSection, &
                           ReadRocketSection

#ifdef INRT
  USE ModInterfacesInteract, ONLY: INRT_ReadMaterialInput
#endif
#ifdef STATS
  USE ModInterfacesStatistics, ONLY: ReadStatisticSection
#endif
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: fname
  CHARACTER(256) :: line
  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction(global,'ReadInputFile',__FILE__)

! ******************************************************************************
! Open file
! ******************************************************************************

  CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.inp',fname)
  
  OPEN(IF_INPUT,FILE=TRIM(fname),FORM='FORMATTED',STATUS='OLD',IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname))
  END IF ! global%error

! ******************************************************************************
! Read file looking for keywords
! ******************************************************************************

  DO
    READ(IF_INPUT,'(A256)',ERR=10,END=86) line

    SELECT CASE( TRIM(line) )
      CASE ('# FORMATS')
        CALL ReadFormatsSection(global)
      CASE ('# REFERENCE')
        CALL ReadReferenceSection(global)
      CASE ('# RANDOM')
        CALL ReadRandomSection(global)
      CASE ('# FLOWMODEL')
        CALL ReadFlowSection(regions)
      CASE ('# MIXTURE')
        CALL ReadMixtureSection(regions)
      CASE ('# VISCMODEL')
        CALL ReadViscositySection(regions)
      CASE ('# ACCELERATION')
        CALL ReadAccelerationSection(global)
      CASE ('# PROBE')
        CALL ReadProbeSection(global)
      CASE ('# FORCES')
        CALL ReadForcesSection(global)
      CASE ('# TIMESTEP')
        CALL ReadTimestepSection(global)
      CASE ('# MULTIGRID')
        CALL ReadMultigridSection(global)
      CASE ('# ABC')
        CALL ReadABCSection(global)
      CASE ('# GFM')
        CALL ReadGFMSection(global)
      CASE ('# MVFRAME')
        CALL ReadMvFrameSection(global)
      CASE ('# NUMERICS')
        CALL ReadNumericsSection(regions)
      CASE ('# THRUST')
        CALL ReadThrustSection(global)
      CASE ('# INITFLOW')
        CALL ReadInitFlowSection(regions) 
      CASE ('# POST') 
        CALL ReadPostSection(global)
      CASE ('# TRANSFORM')
        CALL ReadTransformSection(global)      
      CASE ('# GRIDMOTION')
        CALL ReadGridMotionSection(regions)
      CASE ('# PREP')
        CALL ReadPrepSection(global)
      CASE ('# MISC')
        CALL ReadMiscSection(global)
      CASE ('# ROCKET')
        CALL ReadRocketSection(global)
      CASE ('# TIMEZOOMING')
        CALL ReadTimeZoomingSection(global)
#ifdef STATS
      CASE ('# STATISTICS')
        CALL ReadStatisticSection(global)
#endif
    END SELECT ! TRIM(line)
  END DO ! <empty>

86   CONTINUE

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_INPUT,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname))
  END IF ! global%error

! ******************************************************************************
! Read material sections
! ******************************************************************************

#ifdef INRT
  CALL INRT_ReadMaterialInput(global)
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)
  GOTO 999

10   CONTINUE
  CALL ErrorStop(global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname))

999  CONTINUE

END SUBROUTINE ReadInputFile

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadInputFile.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.7  2009/08/28 18:29:39  mtcampbe
! RocfluMP integration with Rocstar and some makefile tweaks.  To build
! Rocstar with new Rocflu:
! make ROCFLU=RocfluMP
! To build Rocstar with the new RocfluND:
! make ROCFLU=RocfluMP HYPRE=/the/hypre/install/path
!
! Revision 1.6  2009/07/08 19:11:23  mparmar
! Added call to ReadABCSection
!
! Revision 1.5  2008/12/06 08:43:32  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:16:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/11/28 23:17:44  mparmar
! Added reading of # MISC section
!
! Revision 1.2  2007/06/18 17:34:15  mparmar
! Added call to read moving frame section
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.4  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.3  2005/11/10 01:57:00  haselbac
! Added reading of MIXTURE section for Rocflu, clean-up
!
! Revision 1.2  2005/10/31 19:24:37  haselbac
! Added call to ReadMixtureSection, clean-up
!
! Revision 1.1  2004/12/01 16:50:28  haselbac
! Initial revision after changing case
!
! Revision 1.36  2004/11/29 17:13:57  wasistho
! use ModInterfacesStatistics
!
! Revision 1.35  2004/07/24 03:45:28  wasistho
! release readPostSection to flo and flu
!
! Revision 1.34  2004/04/08 03:15:13  wasistho
! nDummyCells in Rocflo read from INITFLOW section
!
! Revision 1.33  2003/11/21 22:33:10  fnajjar
! Updated Random Number Generator
!
! Revision 1.32  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.29  2003/09/26 21:45:34  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.28  2003/08/28 20:05:38  jblazek
! Added acceleration terms.
!
! Revision 1.27  2003/08/11 21:51:17  jblazek
! Added basic global grid smoothing scheme.
!
! Revision 1.26  2003/06/02 17:12:00  jblazek
! Added computation of thrust.
!
! Revision 1.25  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.24  2003/04/28 22:39:12  haselbac
! Added readPostSection and readPrepSection for RFLU
!
! Revision 1.23  2003/04/10 23:25:53  fnajjar
! Added readViscositySection
!
! Revision 1.22  2003/03/24 23:25:48  jferry
! moved some routines from libfloflu to rocinteract
!
! Revision 1.21  2003/03/17 19:37:02  jblazek
! Added inDir to path to the input file.
!
! Revision 1.20  2003/03/11 16:04:19  jferry
! Enclosed USE statements for multi-physics routines within ifdefs
!
! Revision 1.19  2003/02/26 23:38:30  jferry
! eliminated end=999 convention to ensure that files get closed
!
! Revision 1.18  2003/01/28 16:16:32  haselbac
! Read transform section only for RFLU
!
! Revision 1.17  2003/01/28 16:08:31  haselbac
! Added new calls transform, initial solution (RFLU), and grid motion (RFLU)
!
! Revision 1.16  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.15  2002/10/08 15:48:35  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.14  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.13  2002/09/09 14:02:21  haselbac
! mixtInput now under regions
!
! Revision 1.12  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.10  2002/08/18 02:27:47  wasistho
! Added ReadTurbulenceSection
!
! Revision 1.9  2002/06/14 21:17:01  wasistho
! Added time avg statistics
!
! Revision 1.8  2002/06/07 16:40:37  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.7  2002/05/04 16:20:55  haselbac
! Cosmetic changes
!
! Revision 1.6  2002/03/26 19:03:37  haselbac
! Added ROCFLU functionality
!
! Revision 1.5  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.4  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.3  2002/01/11 17:18:31  jblazek
! Updated description of I/O variables.
!
! Revision 1.2  2002/01/02 16:00:03  jblazek
! Added input for multigrid parameters.
!
! Revision 1.1  2001/12/07 16:54:31  jblazek
! Added files to read user input.
!
! ******************************************************************************

