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
! Purpose: Initialize global information.
!
! Description: None.
!
! Input:
!   casename            Case name
!   verbLevel           Verbosity level
!   communicator        Communicator
!   global              Pointer to type global
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_InitGlobal.F90,v 1.3 2015/02/11 16:58:26 brollin Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_InitGlobal(casename,verbLevel,communicator,global)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModMPI
  USE ModParameters

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(*), INTENT(IN) :: casename
  INTEGER, INTENT(IN) :: communicator,verbLevel
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString

  !BBR - begin
  INTEGER :: clock_rate,clock_max,dummy
  !BBR -end

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_InitGlobal.F90,v $ $Revision: 1.3 $'

! ******************************************************************************
! Initialize global variables
! ******************************************************************************

! ==============================================================================
! Main variables
! ==============================================================================

  global%nFunTree = 0

  global%casename = casename

  global%verbLevel = verbLevel
#ifdef GENX
  global%verbLevelCOM = 0
#endif

  global%checkLevel = CHECK_HIGH ! Default

  global%nLevels  = 1
  global%nRegions = 1

  global%nRegionsLocal = 1
  global%nRegionsProc  = 1

! ==============================================================================
! Grid and solution file formats
! ==============================================================================

  global%gridFormat  = FORMAT_ASCII
  global%solutFormat = FORMAT_ASCII

! ==============================================================================
! Parallel stuff
! ==============================================================================

  global%nProcAlloc = 1
  global%myProcid   = MASTERPROC ! Default process number (if not MPI)
  global%mpiComm    = communicator
  global%mpierr     = ERR_NONE

! ==============================================================================
!  ABC Absorbing Boundary Condition 
! ==============================================================================

  global%abcFlag    = .FALSE.
  global%abcKind    = 0 ! 0=sponge, 1=PML 
  global%abcDomain  = 0 ! 0=1D, 1=2D, 2=Circular
  global%abcOrder   = 0 ! 0=Constant, 1=Linear, 2=Quadratic
  global%abcDistrib = 0 ! 0=Constant, 1=Variable reference values

  global%abcSponge  = REAL(CRAZY_VALUE_INT,KIND=RFREAL) 
  global%abcDens    = REAL(CRAZY_VALUE_INT,KIND=RFREAL) 
  global%abcVelX    = REAL(CRAZY_VALUE_INT,KIND=RFREAL) 
  global%abcVelY    = REAL(CRAZY_VALUE_INT,KIND=RFREAL) 
  global%abcVelZ    = REAL(CRAZY_VALUE_INT,KIND=RFREAL) 
  global%abcPress   = REAL(CRAZY_VALUE_INT,KIND=RFREAL) 
  global%abcCenterX = REAL(CRAZY_VALUE_INT,KIND=RFREAL) 
  global%abcCenterY = REAL(CRAZY_VALUE_INT,KIND=RFREAL) 
  global%abcRadius  = REAL(CRAZY_VALUE_INT,KIND=RFREAL) 
  global%abcLeft    = REAL(CRAZY_VALUE_INT,KIND=RFREAL) 
  global%abcRight   = REAL(CRAZY_VALUE_INT,KIND=RFREAL) 
  global%abcBottom  = REAL(CRAZY_VALUE_INT,KIND=RFREAL) 
  global%abcTop     = REAL(CRAZY_VALUE_INT,KIND=RFREAL) 
  global%abcCoeff0  = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%abcCoeff1  = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%abcCoeff2  = REAL(CRAZY_VALUE_INT,KIND=RFREAL)

! ==============================================================================
!  GFM Ghost Fluid Method
! ==============================================================================

  global%gfmFlag    = .FALSE.

  global%gfmNRigids    = 0
  global%gfmNFluids    = 0
  global%gfmNParticles = 0

  global%gfmGridX1 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%gfmGridX2 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%gfmGridY1 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%gfmGridY2 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%gfmGridZ1 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%gfmGridZ2 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%gfmGridNx = CRAZY_VALUE_INT
  global%gfmGridNy = CRAZY_VALUE_INT
  global%gfmGridNz = CRAZY_VALUE_INT

! ==============================================================================
!  Moving Frame
! ==============================================================================

  global%mvFrameFlag  = .FALSE.
  global%mvfAccFlag   = .FALSE.
  global%mvfAccType   = ACC_CONSTANT 
 
  global%iPatchGlobalMvFrame = -1 ! Cannot be CRAZY to avoid format error
  
  global%mvfOmega     = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%mvfPhase     = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%mvfAccTe     = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%mvfAccTs     = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%mvfAccX      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%mvfAccY      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%mvfAccZ      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%mvfLocInitX  = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%mvfLocInitY  = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%mvfLocInitZ  = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%mvfMass      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%mvfVelInitX  = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%mvfVelInitY  = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  global%mvfVelInitZ  = REAL(CRAZY_VALUE_INT,KIND=RFREAL)

! ==============================================================================
! Time stepping
! ==============================================================================

  global%flowType    = FLOW_STEADY
  global%solverType  = SOLV_EXPLICIT
  global%hypreSolver = HYPRE_SOLV_BAMG

  global%tolPC        = 1.0E-12_RFREAL
  global%tolHypre     = 1.0E-16_RFREAL
  global%maxIterPC    = 1000
  global%maxIterHypre = 1000

  global%currentIter       = 0
  global%iterSinceRestart  = 0
  global%maxIter           = 10000
  global%printIter         = 1
  global%restartIter       = 0
  global%resTol            = 1.0E-5_RFREAL
  global%writeIter         = 1000  
  global%fomWriteIter      = 1000  

  global%currentTime       = 0.0_RFREAL ! No physical time set
  global%currentTimeRK     = 0.0_RFREAL
  global%dtImposed         = 1.0E-5_RFREAL
  global%dtMinLimit        = 1.0E-10_RFREAL
  global%fomWriteTime      = 5.0E-5_RFREAL
  global%maxTime           = 1.0E-4_RFREAL
  global%printTime         = 1.0E-12_RFREAL
  global%restartTime       = 0.0_RFREAL
  global%runTime           = global%maxTime
  global%timeSinceRestart  = 0.0_RFREAL
  global%timeSincePrint    = 0.0_RFREAL
  global%timeSinceProbe    = 0.0_RFREAL
  global%timeSinceWrite    = 0.0_RFREAL

  !begin BBR
  global%halftimeSinceWrite = 0.0_RFREAL
  global%wallTime           = 99999999_RFREAL ! Walltime Larger than what any 
                                              ! cluster allows (3years) except
                                              ! maybe Hipergator... (hehehe...)
  global%safeWriteTime      = 99999999_RFREAL ! Same as above
  global%previousTime      = -1.0_RFREAL

  CALL SYSTEM_CLOCK ( dummy, clock_rate, clock_max )
  global%timeruninit = dummy

  global%TimeCycle = .FALSE.
  global%Test_BBR = .FALSE.
  global%Test_Subbu = .FALSE.
  global%Test_Chris = .FALSE.
  global%Test_Yash = .FALSE.
  global%Test_Fred = .FALSE.
  global%Test_Sappy = .FALSE.
  global%Test_Rahul = .FALSE.
  global%Test_Brad = .FALSE.

  global%SafeWrite = 0

  !end BBR

  global%timeSinceFomWrite = 0.0_RFREAL
  global%writeTime         = 5.0E-5_RFREAL  

  global%rkScheme = RK_SCHEME_4_CLASSICAL

#ifndef GENX
  global%timeStamp = 0.0E+0_RFREAL
#endif

  global%restartFromScratch = .FALSE.

  global%dualTstSource = .FALSE.

! ==============================================================================
! Optimal LES stuff
! ==============================================================================

  global%dissOLES = 0.0_RFREAL
  global%enerOLES = 0.0_RFREAL

  global%progressCounter = 1

#ifdef STATS
! ==============================================================================
! Time averaged statistics
! ==============================================================================

  global%doStat     = 0
  global%reStat     = 0
  global%mixtNStat  = 0
  global%turbNStat  = 0
  global%integrTime = 0.0_RFREAL
#endif

! ==============================================================================
! Error variable
! ==============================================================================

  global%error = ERR_NONE

! ==============================================================================
! Mathematical constants (these should really be parameters...)
! ==============================================================================

  global%pi      = 4.0_RFREAL*ATAN(1.0_RFREAL)
  global%rad2deg = 180.0_RFREAL/global%pi
  global%deg2rad = 1.0_RFREAL/global%rad2deg
  global%rad     = global%deg2rad ! Should not be used in Rocflu

! ==============================================================================
! Multigrid
! ==============================================================================

  global%startLevel = 1
  global%cycleType  = MGCYCLE_NO
  global%refineIter = 9999999

! ==============================================================================
! Reference values
! ==============================================================================

  global%refVelocity = 100.0_RFREAL
  global%refPressure = 1.0E+5_RFREAL
  global%refDensity  = 1.2_RFREAL
  global%refCp       = 1004.5_RFREAL
  global%refGamma    = 1.4_RFREAL
  global%refLength   = 1.0_RFREAL
  global%refReNum    = 100.0_RFREAL
  global%prLam       = 0.72_RFREAL
  global%prTurb      = 0.9_RFREAL
  global%scnLam      = 0.22_RFREAL
  global%scnTurb     = 0.9_RFREAL

! ------------------------------------------------------------------------------
! Reference values specific to gas-liquid mixture computations
! ------------------------------------------------------------------------------
  
  global%refBetaPLiq    =  4.5E-7_RFREAL  
  global%refBetaTLiq    = -4.5E-7_RFREAL 
  global%refCvLiq       =  4179.0_RFREAL
  global%refDensityLiq  =  1000.0_RFREAL
  global%refPressLiq    =  1.0E+5_RFREAL
  global%refTempLiq     =  273.0_RFREAL

! ==============================================================================
! Geometric transformations
! ==============================================================================

  global%transformFlag = .FALSE. ! No transformation
  global%distortFlag   = .FALSE. ! No distortion
  global%enforceFlag   = .FALSE.
  
  global%angleX   = 0.0_RFREAL
  global%angleY   = 0.0_RFREAL
  global%angleZ   = 0.0_RFREAL
  global%scaleX   = 1.0_RFREAL
  global%scaleY   = 1.0_RFREAL
  global%scaleZ   = 1.0_RFREAL
  global%distortX = 1.0_RFREAL
  global%distortY = 1.0_RFREAL
  global%distortZ = 1.0_RFREAL  

! ==============================================================================
! Program burn
! ==============================================================================

  global%pbaFlag      = .FALSE.
  global%pbaBurnFlag  = .FALSE.

! ==============================================================================
! Probe
! ==============================================================================

  global%nProbes = 0

  global%probeSaveIter = 1
  global%probeSaveTime = 0.0_RFREAL

! ==============================================================================
! Forces and acceleration
! ==============================================================================

  global%accelOn        = .FALSE.
  global%forceFlag      = .FALSE.
  global%patchCoeffFlag = .TRUE.
  
  global%forceWriteCntr  = 0
  global%thrustWriteCntr = 0

  global%accelX = 0.0_RFREAL
  global%accelY = 0.0_RFREAL
  global%accelZ = 0.0_RFREAL    

  global%gravity = 9.81_RFREAL    

  global%forceRefArea   = 1.0_RFREAL
  global%forceRefLength = 1.0_RFREAL
  global%forceRefXCoord = 0.0_RFREAL
  global%forceRefYCoord = 0.0_RFREAL
  global%forceRefZCoord = 0.0_RFREAL    

! ==============================================================================
! Random number generator
! ==============================================================================

  global%randSeedOffset = 0
  global%randSeedType   = RAND_SEED_TYPE_FIXED

! ==============================================================================
! Multiphysics modules
! ==============================================================================

  global%inrtUsed = .FALSE.
  global%peulUsed = .FALSE.
  global%plagUsed = .FALSE.
  global%specUsed = .FALSE. 

  global%plagUsedSave = global%plagUsed

#ifdef SPEC
  global%nSpecies = 0
#endif

! ==============================================================================
! Miscellaneous
! ==============================================================================

  global%resInit = 0.0_RFREAL

  global%warnCounter = 0
  
  global%moduleType = MODULE_TYPE_NONE 

! ==============================================================================
! Misc section related variables
! ==============================================================================

  global%nBVertEstCorrFactor = 1.0_RFREAL

! ==============================================================================
! Flow initialization
! ==============================================================================

  global%initFlowFlag = CRAZY_VALUE_INT ! To force setting of value
#ifdef PLAG
  global%initPlagFlag = CRAZY_VALUE_INT ! To force setting of value
#endif

! ==============================================================================
! Preprocessor
! ==============================================================================

  global%prepPartMode    = PARTITION_MODE_PROPER
  global%syPePatchesFlag = .FALSE.

! ==============================================================================
! Postprocessor
! ==============================================================================

  global%postOutputFormat = POST_OUTPUT_FORMAT_TECPLOT

  global%postInterpOrder = 1
  global%postInterpType  = INTERP_TYPE_PROPER
  global%postPlotType    = PLOT_GRID_FLOW

  global%postCompErrFlag    = .FALSE.
  global%postMergeFlag      = .TRUE.
  global%postPlotMixtCvFlag = .TRUE.
  global%postPlotMixtDvFlag = .TRUE.
  global%postPlotMixtGvFlag = .TRUE.
  global%postPlotPatchFlag  = .TRUE.
  global%postPlotVolFlag    = .TRUE.
  global%postWriteMergeFlag = .FALSE.
  global%postZoneVirtFlag   = .FALSE.

  global%postDiscFlag     = .FALSE.
  global%postExtractFlag  = .FALSE.
  global%postLag2EulFlag  = .FALSE.
  global%postSpecFlag     = .FALSE.
  global%postVortFlag     = .FALSE.
  global%postVortCoreFlag = .FALSE.

  global%postSchType  =  0
  global%postSchExp   = 10.0_RFREAL
  global%postNFringes = 32

  global%postPlagFrac = 1.0_RFREAL

  global%postNServers = 1

! ==============================================================================
! Picking
! ==============================================================================

  global%pickCoordFlag = .FALSE.

  global%pickXCoordLow = -HUGE(1.0_RFREAL)
  global%pickXCoordUpp =  HUGE(1.0_RFREAL)
  global%pickYCoordLow = -HUGE(1.0_RFREAL)
  global%pickYCoordUpp =  HUGE(1.0_RFREAL)
  global%pickZCoordLow = -HUGE(1.0_RFREAL)
  global%pickZCoordUpp =  HUGE(1.0_RFREAL)

! ==============================================================================
! GENX; if GENX is defined, winName is set in rocflu_load_module.F90
! ==============================================================================

#ifndef GENX
  global%surfWinName      = ''
  global%surfWinNameInput = ''
  global%volWinName       = ''
  global%volWinNameInput  = ''
#endif

! The following will be overwritten by RFLU_ReadGENXControlFile

  global%inDir   = './'
  global%outDir  = './'

! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE RFLU_InitGlobal

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_InitGlobal.F90,v $
! Revision 1.3  2015/02/11 16:58:26  brollin
! New fix to Safe Restart Feature. Now Safe restart communicated to all ranks
! once any rank hit the Safe restart condition.
!
! Revision 1.2  2015/02/10 18:24:22  brollin
! adding missing Metis library files
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.13  2010/03/14 23:22:39  mparmar
! Added initialization of hypreSolver
!
! Revision 1.12  2009/07/08 19:11:30  mparmar
! Added initialization of absorbing layer data
!
! Revision 1.11  2008/12/06 08:43:35  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:16:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2008/05/29 01:35:17  mparmar
! Added mvfAccType, mvfOmega, mvfPhase for moving reference frame
!
! Revision 1.8  2008/01/07 14:03:38  haselbac
! Added init of postZoneVirtFlag
!
! Revision 1.7  2007/11/28 23:04:58  mparmar
! Added initialization for SOLV_IMPLICIT_HM solver and nBVertEstCorrFactor
!
! Revision 1.6  2007/09/04 13:15:40  haselbac
! Added init of new plotting variables
!
! Revision 1.5  2007/08/07 17:16:47  haselbac
! Minor bug fix
!
! Revision 1.4  2007/06/18 17:45:39  mparmar
! Added initialization of moving reference frame data
!
! Revision 1.3  2007/06/14 01:48:09  haselbac
! Added init of global%runTime
!
! Revision 1.2  2007/04/12 17:54:00  haselbac
! Added init of postPlagFrac
!
! Revision 1.1  2007/04/09 18:48:52  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.50  2006/10/20 21:23:29  mparmar
! Added initialization of gravity,thrustWriteCntr
!
! Revision 1.49  2006/08/04 03:03:23  haselbac
! Added init of grid distortion parameters, fixed typo
!
! Revision 1.48  2006/05/06 15:20:09  haselbac
! Bug fix: Missing ifdef PLAG
!
! Revision 1.47  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.46  2006/03/26 20:21:35  haselbac
! Added ref values for GL model
!
! Revision 1.45  2006/03/25 21:42:47  haselbac
! Added init of syPePatchesFlag
!
! Revision 1.44  2006/02/06 23:56:37  haselbac
! Added communicator to arg list to fix bug in GENX init with Panda
!
! Revision 1.43  2005/12/10 23:28:05  haselbac
! Renamed geom post variables to pick variables
!
! Revision 1.42  2005/12/10 16:54:30  haselbac
! Added init for postLag2EulFlag
!
! Revision 1.41  2005/12/01 17:09:53  fnajjar
! Added appropriate initialization, checking and printing of random seed type
!
! Revision 1.40  2005/11/10 02:01:55  haselbac
! Added init of acceleration components
!
! Revision 1.39  2005/10/28 19:17:07  haselbac
! Added setting of postPlotPatchFlag
!
! Revision 1.38  2005/10/05 20:02:28  haselbac
! Added init for post output format and nservers
!
! Revision 1.37  2005/08/10 00:30:55  haselbac
! Added init of postVortCoreFlag
!
! Revision 1.36  2005/08/09 00:54:26  haselbac
! Added init of patchCoeffFlag, postCompErrFlag, postWriteMergeFlag
!
! Revision 1.35  2005/07/25 12:20:22  haselbac
! Added init of postVortFlag, changed postSchExp to float
!
! Revision 1.34  2005/07/05 19:26:38  haselbac
! Added init for postSchType and postSchExp
!
! Revision 1.33  2005/07/01 15:13:12  haselbac
! Added init of verbLevelCOM
!
! Revision 1.32  2005/05/01 14:18:28  haselbac
! Added init of postDiscFlag and postNFringes
!
! Revision 1.31  2005/04/15 15:06:17  haselbac
! Removed Charm/FEM routines, added init of mpiComm and plagUsedSave
!
! Revision 1.30  2005/03/31 16:49:55  haselbac
! Added initialization of currentTimeRK
!
! Revision 1.29  2004/11/17 16:27:54  haselbac
! Added setting of rkScheme
!
! Revision 1.28  2004/10/26 15:16:16  haselbac
! Added initialization of postExtractFlag
!
! Revision 1.27  2004/10/19 19:24:22  haselbac
! Added and removed several variables
!
! Revision 1.26  2004/10/09 16:36:15  fnajjar
! Initialized initPlagFlag
!
! Revision 1.25  2004/07/21 14:54:31  haselbac
! Added initialization of postInterpType
!
! Revision 1.24  2004/06/16 20:00:19  haselbac
! Added force variables, clean-up, cosmetics
!
! Revision 1.23  2004/03/15 21:02:57  haselbac
! Changed init value of checkLevel to use parameter
!
! Revision 1.22  2004/03/05 22:09:01  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.21  2004/01/29 22:56:33  haselbac
! Added initialization of rad2deg and deg2rad
!
! Revision 1.20  2003/12/04 03:23:54  haselbac
! Collect initialization of ALL global variables, clean-up
!
! Revision 1.19  2003/11/25 21:02:47  haselbac
! Added initialization of initFlowFlag
!
! Revision 1.18  2003/10/29 21:38:14  haselbac
! Added init for new global variables
!
! Revision 1.17  2003/10/15 02:40:07  haselbac
! Added init for dtMinLimit
!
! Revision 1.16  2003/08/07 15:30:32  haselbac
! Added and changed var names
!
! Revision 1.15  2003/07/22 01:55:54  haselbac
! Added init of postInterpOrder and warnCounter
!
! Revision 1.14  2003/05/05 18:39:03  haselbac
! Added init of plotType
!
! Revision 1.13  2003/04/29 21:48:05  haselbac
! Added initialization of surfDiverFlag
!
! Revision 1.12  2003/04/28 22:40:01  haselbac
! Added init of post and prep vars, cosmetics
!
! Revision 1.11  2003/04/12 21:36:22  haselbac
! Added init of FEMRocfluGrid
!
! Revision 1.10  2003/04/07 14:21:07  haselbac
! Initialize new probe parameters, cosmetic changes
!
! Revision 1.9  2003/03/15 17:00:03  haselbac
! Removed splitFlag
!
! Revision 1.8  2003/02/01 00:27:17  haselbac
! Added initialization of transform variables and enforceFlag
!
! Revision 1.7  2002/10/27 18:52:16  haselbac
! Removed tabs
!
! Revision 1.6  2002/10/16 21:11:40  haselbac
! Added initialization of resInit
!
! Revision 1.5  2002/10/13 21:41:40  jiao
! Compiled Rocflu on IBM SP.
!
! Revision 1.4  2002/10/05 18:46:28  haselbac
! Integration into GENX, some clean-up
!
! Revision 1.3  2002/09/09 14:15:01  haselbac
! global now under regions
!
! Revision 1.2  2002/07/25 14:42:07  haselbac
! Added arguments, status of routine has changed, now called from main
!
! Revision 1.1  2002/06/14 20:05:48  haselbac
! Former RFLU_InitLocalGlobal.F90, renamed as ModLocal.F90 no longer used
!
! Revision 1.3  2002/05/28 13:49:10  haselbac
! Changed initial value for currentTime
!
! Revision 1.2  2002/05/07 18:51:01  haselbac
! Added initialization of transform flag
!
! Revision 1.1  2002/05/04 16:09:00  haselbac
! Initial revision
!
! ******************************************************************************

