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
! Purpose: Read in user input related to flow initialization.
!
! Description: None.
!
! Input: 
!   regions     Data for regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ReadInitFlowSection.F90,v 1.4 2016/03/22 19:32:14 fred Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReadInitFlowSection(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters
  USE ModMixture, ONLY: t_mixt_input
  
  USE ModInterfaces, ONLY: ReadSection, & 
                           ReadRegionSection  
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Local variables
! ==============================================================================

  INTEGER :: i,iReg,nVals
!  INTEGER, PARAMETER :: NVALS_MAX = 43

! Rahul- I/P for grid dx, dy and dz
! Fred - Variable JWL density
  INTEGER, PARAMETER :: NVALS_MAX = 50
! Rahul - end

  CHARACTER(10) :: keys(NVALS_MAX)
  LOGICAL :: defined(NVALS_MAX)
  REAL(RFREAL) :: vals(NVALS_MAX)
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction(global,'ReadInitFlowSection',__FILE__)

! specify keywords and search for them

  nVals = NVALS_MAX

  keys( 1) = 'FLAG'
  keys( 2) = 'VELX'
  keys( 3) = 'VELY'
  keys( 4) = 'VELZ'
  keys( 5) = 'PRESS'
  keys( 6) = 'DENS'
  keys( 7) = 'TEMP'
  keys( 8) = 'IVAL1'
  keys( 9) = 'IVAL2'  
  keys(10) = 'IVAL3'
  keys(11) = 'IVAL4'
  keys(12) = 'IVAL5'
  keys(13) = 'IVAL6' 
  keys(14) = 'IVAL7' 
  keys(15) = 'IVAL8' 
  keys(16) = 'IVAL9' 
  keys(17) = 'IVAL10'
  keys(18) = 'IVAL11'
  keys(19) = 'IVAL12'
  keys(20) = 'IVAL13'
  keys(21) = 'IVAL14'
  keys(22) = 'IVAL15'
  keys(23) = 'IVAL16'
  keys(24) = 'RVAL1'
  keys(25) = 'RVAL2'  
  keys(26) = 'RVAL3' 
  keys(27) = 'RVAL4'
  keys(28) = 'RVAL5'  
  keys(29) = 'RVAL6'
  keys(30) = 'RVAL7'
  keys(31) = 'RVAL8'         
  keys(32) = 'RVAL9'
  keys(33) = 'RVAL10'
  keys(34) = 'RVAL11' 
  keys(35) = 'RVAL12'
  keys(36) = 'RVAL13'  
  keys(37) = 'RVAL14'
  keys(38) = 'RVAL15'
  keys(39) = 'RVAL16'         
  keys(40) = 'RVAL17'         
  keys(41) = 'RVAL18'         
  keys(42) = 'RVAL19'         
  keys(43) = 'RVAL20'         

! Rahul - grid spacing of shktb case
  keys(44) = 'DX'         
  keys(45) = 'DY'         
  keys(46) = 'DZ'         
! Rahul - end

! Fred - Variable JWL PETN density  
  keys(47) = 'RVAL24'
! Fred - end

! Fred - Choose collision model/multimode
  keys(48) = 'IVAL17'
  keys(49) = 'IVAL18'
  keys(50) = 'IVAL19'
! Fred - end

  CALL ReadSection(global,IF_INPUT,nVals,keys,vals,defined )

  IF ( defined(1) .EQV. .FALSE. ) THEN
    CALL ErrorStop(global,ERR_VAL_UNDEFINED,__LINE__,'INITFLOW-FLAG')
  ELSE
    IF ( NINT(vals(1)) == INITFLOW_FROMSCRATCH ) THEN 
      global%initFlowFlag = INITFLOW_FROMSCRATCH
    
      DO iReg = LBOUND(regions,1),UBOUND(regions,1) 
        SELECT CASE ( regions(iReg)%mixtInput%fluidModel ) 
          CASE ( FLUID_MODEL_INCOMP ) 
            DO i = 2,5
              IF ( defined(i) .EQV. .FALSE. ) THEN 
                CALL ErrorStop(global,ERR_VAL_UNDEFINED,__LINE__,keys(i))
              END IF ! defined
            END DO ! i          
          
            regions(iReg)%mixtInput%iniVelX  = vals(2)
            regions(iReg)%mixtInput%iniVelY  = vals(3)   
            regions(iReg)%mixtInput%iniVelZ  = vals(4)  
            regions(iReg)%mixtInput%iniPress = vals(5)            
          CASE ( FLUID_MODEL_COMP ) 
            SELECT CASE ( regions(iReg)%mixtInput%gasModel ) 
              CASE ( GAS_MODEL_TCPERF, & 
                     GAS_MODEL_MIXT_TCPERF, & 
                     GAS_MODEL_MIXT_TPERF, & 
                     GAS_MODEL_MIXT_PSEUDO ) 
                DO i = 2,6
                  IF ( defined(i) .EQV. .FALSE. ) THEN 
                    CALL ErrorStop(global,ERR_VAL_UNDEFINED,__LINE__,keys(i))
                  END IF ! defined
                END DO ! i          

                regions(iReg)%mixtInput%iniVelX  = vals(2)
                regions(iReg)%mixtInput%iniVelY  = vals(3)   
                regions(iReg)%mixtInput%iniVelZ  = vals(4)  
                regions(iReg)%mixtInput%iniPress = vals(5)
                regions(iReg)%mixtInput%iniDens  = vals(6) 
              CASE ( GAS_MODEL_MIXT_GASLIQ )                         
                DO i = 2,4
                  IF ( defined(i) .EQV. .FALSE. ) THEN 
                    CALL ErrorStop(global,ERR_VAL_UNDEFINED,__LINE__,keys(i))
                  END IF ! defined
                END DO ! i

                DO i = 6,7
                  IF ( defined(i) .EQV. .FALSE. ) THEN 
                    CALL ErrorStop(global,ERR_VAL_UNDEFINED,__LINE__,keys(i))
                  END IF ! defined
                END DO ! i

                regions(iReg)%mixtInput%iniVelX  = vals(2)
                regions(iReg)%mixtInput%iniVelY  = vals(3)   
                regions(iReg)%mixtInput%iniVelZ  = vals(4)  
                regions(iReg)%mixtInput%iniPress = vals(5)
                regions(iReg)%mixtInput%iniTemp  = vals(7)
              CASE DEFAULT
                CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END SELECT ! regions(iReg)%mixtInput%gasModel
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! regions(iReg)%mixtInput%fluidModel       
      END DO ! iReg
    ELSE IF ( NINT(vals(1)) == INITFLOW_FROMFILE ) THEN 
      global%initFlowFlag = INITFLOW_FROMFILE
    ELSE IF ( NINT(vals(1)) == INITFLOW_FROMHARDCODE ) THEN  
      global%initFlowFlag = INITFLOW_FROMHARDCODE      
    ELSE IF ( NINT(vals(1)) == INITFLOW_FROMCOMBO_SERIAL ) THEN
      global%initFlowFlag = INITFLOW_FROMCOMBO_SERIAL
    ELSE IF ( NINT(vals(1)) == INITFLOW_FROMCOMBO_PARALLEL ) THEN
      global%initFlowFlag = INITFLOW_FROMCOMBO_PARALLEL     
    ELSE IF ( NINT(vals(1)) == INITFLOW_FROMFILE_DIFFNTOT ) THEN
      global%initFlowFlag = INITFLOW_FROMFILE_DIFFNTOT     
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! NINT(vals(1))     
  END IF ! defined

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)  
    IF ( defined(8) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepIntVal1 = CRAZY_VALUE_INT
    ELSE
      regions(iReg)%mixtInput%prepIntVal1 = vals(8)
    END IF ! defined  

    IF ( defined(9) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepIntVal2 = CRAZY_VALUE_INT
    ELSE
      regions(iReg)%mixtInput%prepIntVal2 = vals(9)
    END IF ! defined  

    IF ( defined(10) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepIntVal3 = CRAZY_VALUE_INT
    ELSE
      regions(iReg)%mixtInput%prepIntVal3 = vals(10)
    END IF ! defined
    
    IF ( defined(11) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepIntVal4 = CRAZY_VALUE_INT
    ELSE
      regions(iReg)%mixtInput%prepIntVal4 = vals(11)
    END IF ! defined  
    
    IF ( defined(12) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepIntVal5 = CRAZY_VALUE_INT
    ELSE
      regions(iReg)%mixtInput%prepIntVal5 = vals(12)
    END IF ! defined  
    
    IF ( defined(13) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepIntVal6 = CRAZY_VALUE_INT
    ELSE
      regions(iReg)%mixtInput%prepIntVal6 = vals(13)
    END IF ! defined                

    IF ( defined(14) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepIntVal7 = CRAZY_VALUE_INT
    ELSE
      regions(iReg)%mixtInput%prepIntVal7 = vals(14)
    END IF ! defined                

    IF ( defined(15) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepIntVal8 = CRAZY_VALUE_INT
    ELSE
      regions(iReg)%mixtInput%prepIntVal8 = vals(15)
    END IF ! defined                

    IF ( defined(16) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepIntVal9 = CRAZY_VALUE_INT
    ELSE
      regions(iReg)%mixtInput%prepIntVal9 = vals(16)
    END IF ! defined                

    IF ( defined(17) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepIntVal10 = CRAZY_VALUE_INT
    ELSE
      regions(iReg)%mixtInput%prepIntVal10 = vals(17)
    END IF ! defined                

    IF ( defined(18) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepIntVal11 = CRAZY_VALUE_INT
    ELSE
      regions(iReg)%mixtInput%prepIntVal11 = vals(18)
    END IF ! defined                

    IF ( defined(19) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepIntVal12 = CRAZY_VALUE_INT
    ELSE
      regions(iReg)%mixtInput%prepIntVal12 = vals(19)
    END IF ! defined                

    IF ( defined(20) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepIntVal13 = CRAZY_VALUE_INT
    ELSE
      regions(iReg)%mixtInput%prepIntVal13 = vals(20)
    END IF ! defined                

    IF ( defined(21) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepIntVal14 = CRAZY_VALUE_INT
    ELSE
      regions(iReg)%mixtInput%prepIntVal14 = vals(21)
    END IF ! defined                

    IF ( defined(22) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepIntVal15 = CRAZY_VALUE_INT
    ELSE
      regions(iReg)%mixtInput%prepIntVal15 = vals(22)
    END IF ! defined                

    IF ( defined(23) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepIntVal16 = CRAZY_VALUE_INT
    ELSE
      regions(iReg)%mixtInput%prepIntVal16 = vals(23)
    END IF ! defined                

    IF ( defined(24) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal1 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal1 = vals(24)
    END IF ! defined  

    IF ( defined(25) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal2 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal2 = vals(25)
    END IF ! defined  

    IF ( defined(26) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal3 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal3 = vals(26)
    END IF ! defined 

    IF ( defined(27) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal4 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal4 = vals(27)
    END IF ! defined 

    IF ( defined(28) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal5 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal5 = vals(28)
    END IF ! defined  

    IF ( defined(29) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal6 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal6 = vals(29)
    END IF ! defined
    
    IF ( defined(30) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal7 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal7 = vals(30)
    END IF ! defined
    
    IF ( defined(31) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal8 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal8 = vals(31)
    END IF ! defined        

    IF ( defined(32) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal9 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal9 = vals(32)
    END IF ! defined        

    IF ( defined(33) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal10 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal10 = vals(33)
    END IF ! defined        

    IF ( defined(34) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal11 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal11 = vals(34)
    END IF ! defined        

    IF ( defined(35) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal12 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal12 = vals(35)
    END IF ! defined        

    IF ( defined(36) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal13 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal13 = vals(36)
    END IF ! defined        

    IF ( defined(37) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal14 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal14 = vals(37)
    END IF ! defined        

    IF ( defined(38) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal15 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal15 = vals(38)
    END IF ! defined        

    IF ( defined(39) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal16 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal16 = vals(39)
    END IF ! defined        

    IF ( defined(40) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal17 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal17 = vals(40)
    END IF ! defined        

    IF ( defined(41) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal18 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal18 = vals(41)
    END IF ! defined        

    IF ( defined(42) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal19 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal19 = vals(42)
    END IF ! defined        

    IF ( defined(43) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal20 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal20 = vals(43)
    END IF ! defined        
! Rahul - dx, dy and dz for shktb case
    IF ( defined(44) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal21 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal21 = vals(44)
    END IF ! defined        

    IF ( defined(45) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal22 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal22 = vals(45)
    END IF ! defined        

    IF ( defined(46) .EQV. .FALSE. ) THEN 
      regions(iReg)%mixtInput%prepRealVal23 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal23 = vals(46)
    END IF ! defined        
! Rahul - end

! Fred - Variable PETN density for JWL 
    IF ( defined(47) .EQV. .FALSE. ) THEN
      regions(iReg)%mixtInput%prepRealVal24 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepRealVal24 = vals(47)
    END IF ! defined
! Fred - end

    IF ( defined(48) .EQV. .FALSE. ) THEN
      regions(iReg)%mixtInput%prepIntVal17 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepIntVal17 = vals(48)
    END IF ! defined

    IF ( defined(49) .EQV. .FALSE. ) THEN
      regions(iReg)%mixtInput%prepIntVal18 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepIntVal18 = vals(49)
    END IF ! defined
 
    IF ( defined(50) .EQV. .FALSE. ) THEN
      regions(iReg)%mixtInput%prepIntVal19 = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    ELSE
      regions(iReg)%mixtInput%prepIntVal19 = vals(50)
    END IF ! defined

  END DO ! iReg           

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE ReadInitFlowSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadInitFlowSection.F90,v $
! Revision 1.4  2016/03/22 19:32:14  fred
! Adding Glasser's model as a choice for collision modelling
!
! Revision 1.3  2016/02/06 17:24:25  fred
! Adding Variable JWL density RVAL into read list
!
! Revision 1.2  2016/02/03 20:38:34  rahul
! Added 3 additional keys to read dx, dy and dz of the shktb grid from the
! .inp file. This is a temporary fix.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2009/07/08 19:11:22  mparmar
! Adapted to new init option
!
! Revision 1.3  2008/12/06 08:43:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.9  2007/04/05 00:56:57  haselbac
! Added additional RVALxy params for 2p shocktube problems
!
! Revision 1.8  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.7  2006/03/26 20:21:20  haselbac
! Added TEMP input argument
!
! Revision 1.6  2005/11/17 14:37:24  haselbac
! Added more RVAL variables
!
! Revision 1.5  2005/09/13 21:36:45  haselbac
! Adapted to new init option
!
! Revision 1.4  2005/04/20 14:38:36  haselbac
! Added more int and real vals
!
! Revision 1.3  2005/03/29 22:28:31  haselbac
! Added setting of initFlowFlag for combo option
!
! Revision 1.2  2005/03/22 03:32:39  haselbac
! Added initialization of integer and real helper variables
!
! Revision 1.1  2004/12/01 16:50:26  haselbac
! Initial revision after changing case
!
! Revision 1.13  2004/11/14 19:35:42  haselbac
! Added initialization for incompressible fluid model
!
! Revision 1.12  2004/07/28 16:41:32  haselbac
! Bug fix: Initial values not assigned to all regions
!
! Revision 1.11  2004/04/08 03:15:24  wasistho
! nDummyCells in Rocflo read from INITFLOW section
!
! Revision 1.10  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.7  2003/09/15 00:36:19  haselbac
! Added hard-code option as input
!
! Revision 1.6  2003/05/16 02:27:43  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.5  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.4  2003/03/25 19:15:17  haselbac
! Fixed bug in RCSIdentString
!
! Revision 1.3  2003/03/15 16:27:39  haselbac
! Added KIND qualifyer
!
! Revision 1.2  2003/02/13 22:30:54  jferry
! removed RFLU_ prefix in RegisterFunction
!
! Revision 1.1  2003/01/28 16:12:56  haselbac
! Moved here from rfluprep
!
! Revision 1.4  2002/10/27 19:23:46  haselbac
! Removed tabs
!
! Revision 1.3  2002/10/07 14:11:29  haselbac
! Removed tabs
!
! Revision 1.2  2002/09/09 16:40:13  haselbac
! global and mixtInput now under regions
!
! Revision 1.1  2002/04/11 19:13:22  haselbac
! Initial revision
!
! ******************************************************************************

