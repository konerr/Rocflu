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
! Purpose: read in user input related to discrete particle module.
!
! Description: none.
!
! Input: user input file.
!
! Output: regions = total number of particles, drag model,injection model.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_ReadDisPartSection.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_ReadDisPartSection( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMixture, ONLY    : t_mixt_input
  USE ModPartLag, ONLY    : t_plag_input
  USE ModError
  USE ModParameters

  USE PLAG_ModParameters
  USE PLAG_ModInterfaces, ONLY : PLAG_ReadPdfFromFile

  USE ModInterfaces, ONLY: ReadSection

  IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================
  
  TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
!   Locals
! ==============================================================================

  INTEGER, PARAMETER :: NVALS_MAX=15

  CHARACTER(CHRLEN) :: RCSIdentString
  CHARACTER(20)     :: keys(NVALS_MAX)
  INTEGER           :: brbeg,brend,nVals
  LOGICAL           :: defined(NVALS_MAX)
  REAL(RFREAL)      :: vals(NVALS_MAX)

  TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ReadDisPartSection.F90,v $ $Revision: 1.1.1.1 $'

  global => regions(1)%global

  CALL RegisterFunction( global, 'PLAG_ReadDisPartSection',__FILE__ )

! ******************************************************************************
! Set variables and specify keywords 
! ******************************************************************************

  nVals   = NVALS_MAX

  defined(:) = .FALSE.

  keys(1)  = 'USED'
  keys(2)  = 'NPCLSMAX'  
  keys(3)  = 'INTRPLMIXTMODEL'
  keys(4)  = 'BREAKUPMODEL'
  keys(5)  = 'BREAKUPFAC'
  keys(6)  = 'BREAKUPWEBSWI'
  keys(7)  = 'FINDPCLMETHOD'
  keys(8)  = 'INITDIMENS'
  keys(9)  = 'NUNSTEADYDATA'
  keys(10) = 'COLLISION'
  keys(11) = 'COLLISIONPS'
  keys(12) = 'USESTABLEDT'
  keys(13) = 'LIMITFORCE'
  keys(14) = 'CFL'
  keys(15) = 'CRWFLAG'

! ******************************************************************************
! Read keywords 
! ******************************************************************************

  CALL ReadSection(global,IF_INPUT,nVals,keys(1:nVals),vals(1:nVals), &
                   defined(1:nVals))

  brbeg = LBOUND(regions,1)
  brend = UBOUND(regions,1)

! ******************************************************************************
! Copy values
! ******************************************************************************

  IF (defined(1)) THEN
     IF (vals(1) < 0.1) THEN
        regions(brbeg:brend)%plagInput%readStatus = 0
     ELSE
        regions(brbeg:brend)%plagInput%readStatus = 1
     ENDIF ! vals
  ENDIF ! defined

  IF (defined(2)) &
       regions(brbeg:brend)%plagInput%nPclsMax = MAX(1,INT(ABS(vals(2))+0.5_RFREAL))

  IF (defined(3)) THEN
    regions(brbeg:brend)%plagInput%intrplMixtModel    = NINT(ABS(vals(3)))
    IF ( vals(3) > -0.1 .AND. vals(3) < 0.1 ) &
      regions(brbeg:brend)%plagInput%intrplMixtModel  = ZEROTH_ORDER
    IF ( vals(3) > 0.9 .AND. vals(3) < 1.1 ) &
      regions(brbeg:brend)%plagInput%intrplMixtModel  = FIRST_ORDER
  ENDIF ! defined

  IF (defined(4)) THEN
    regions(brbeg:brend)%plagInput%breakupModel    = NINT(ABS(vals(4))) 
    IF ( vals(4) > -0.1 .AND. vals(4) < 0.1 ) &
      regions(brbeg:brend)%plagInput%breakupModel  = PLAG_BREAKUP_NOMODEL 
    IF ( vals(4) > 0.9 .AND. vals(4) < 1.1 ) &
      regions(brbeg:brend)%plagInput%breakupModel  = PLAG_BREAKUP_MODEL1
  ENDIF ! defined

  IF (defined(5)) &
    regions(brbeg:brend)%plagInput%breakupFac = ABS(vals(5))
              
  IF (defined(6)) THEN
    regions(brbeg:brend)%plagInput%breakupWebSwi   = NINT(ABS(vals(6))) 
    IF ( vals(6) > -0.1 .AND. vals(6) < 0.1 ) &
      regions(brbeg:brend)%plagInput%breakupWebSwi = PLAG_BREAKUP_NOWEBSWI 
    IF ( vals(6) > 0.9 .AND. vals(6) < 1.1 ) &
      regions(brbeg:brend)%plagInput%breakupWebSwi = PLAG_BREAKUP_WEBSWI1
  ENDIF ! defined

  IF ( defined(7) .EQV. .TRUE. ) THEN 
    regions(brbeg:brend)%plagInput%findPclMethod = vals(7)  
    ! Subbu - Check if the partmode=2 if HardCoded Pcl Tracking used
    !IF ( regions(brbeg:brend)%plagInput%findPclMethod == FIND_PCL_METHOD_HARDCODE ) THEN
    IF ( vals(7) == FIND_PCL_METHOD_HARDCODE ) THEN
     IF (global%prepPartMode .NE. PARTITION_MODE_IMPOSED) THEN
       CALL ErrorStop(global,ERR_PLAG_FINDPCL_PARTMODE_INCONSISTENT,__LINE__)
     END IF
    END IF
    ! Subbu - End Pcl Tracking Check 
  END IF ! defined

  IF ( defined(8) .EQV. .TRUE. ) THEN 
    regions(brbeg:brend)%plagInput%initDimens = vals(8)  
  END IF ! defined

  IF ( defined(9) .EQV. .TRUE. ) THEN 
    regions(brbeg:brend)%plagInput%nUnsteadyData = vals(9)  
  END IF ! defined

  IF ( defined(10) .EQV. .TRUE. ) THEN 
    IF ( NINT(ABS(vals(10))) > 0.9 .AND. NINT(ABS(vals(10))) < 1.1 ) THEN
      regions(brbeg:brend)%plagInput%flagCollision = .TRUE.
    END IF ! NINT(ABS(vals(10)))
  END IF ! defined

  IF ( defined(11) .EQV. .TRUE. ) THEN 
    regions(brbeg:brend)%plagInput%collisionPs = vals(11)  
  END IF ! defined

  IF ( defined(12) .EQV. .TRUE. ) THEN 
    IF ( NINT(ABS(vals(12))) > 0.9 .AND. NINT(ABS(vals(12))) < 1.1 ) THEN
    regions(brbeg:brend)%plagInput%flagStability = .TRUE.
    END IF ! NINT(ABS(vals(12)))
  END IF ! defined

  IF ( defined(13) .EQV. .TRUE. ) THEN 
    IF ( NINT(ABS(vals(13))) > 0.9 .AND. NINT(ABS(vals(13))) < 1.1 ) THEN
    regions(brbeg:brend)%plagInput%limitForce = .TRUE.
    END IF ! NINT(ABS(vals(13)))
  END IF ! defined

  IF ( defined(14) .EQV. .TRUE. ) THEN 
    regions(brbeg:brend)%plagInput%cfl = ABS(vals(14))
  END IF ! defined

  IF ( defined(15) .EQV. .TRUE. ) THEN 
    IF ( NINT(ABS(vals(15))) > 0.9 .AND. NINT(ABS(vals(15))) < 1.1 ) THEN
      regions(brbeg:brend)%plagInput%flagCRW = .TRUE.
    END IF ! NINT(ABS(vals(10)))
  END IF ! defined

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_ReadDisPartSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ReadDisPartSection.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.7  2008/12/06 08:43:52  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2008/01/12 17:41:08  haselbac
! Bug fix: Added reading of initDimens again
!
! Revision 1.4  2007/05/16 22:43:04  fnajjar
! Modified routines to be aligned with new bc datastructure and cosmetic changes
!
! Revision 1.3  2007/04/16 23:22:44  fnajjar
! Deleted all RFLO-pertinent calls and statements
!
! Revision 1.2  2007/04/15 02:39:59  haselbac
! Added reading of initDimens
!
! Revision 1.1  2007/04/09 18:50:27  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:34  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.5  2007/03/06 23:15:32  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.4  2005/05/17 18:35:41  fnajjar
! Further IO improvement
!
! Revision 1.3  2005/05/17 17:07:37  fnajjar
! Improved loading of values for proper error trapping
!
! Revision 1.2  2005/04/25 18:39:10  luca1
! Imposed PDF from file option for random particle ejection
!
! Revision 1.1  2004/12/01 20:58:05  fnajjar
! Initial revision after changing case
!
! Revision 1.13  2004/10/08 22:13:33  haselbac
! Added reading of findPclMethod
!
! Revision 1.12  2004/10/08 20:03:14  fnajjar
! Bug fix for nVals size and defined NVALS_MAX
!
! Revision 1.11  2004/07/29 16:58:07  fnajjar
! Inconsistent numbers in array vals leading to io clobbering
!
! Revision 1.10  2004/06/17 15:19:48  fnajjar
! Added infrastructure for ejection model and revamped DISPART input section
!
! Revision 1.9  2004/06/16 23:07:17  fnajjar
! Renamed variabled for CRE kernel
!
! Revision 1.8  2004/03/10 23:09:50  fnajjar
! Added maximum buffer size for corner-edge cells
!
! Revision 1.7  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.6  2004/02/26 21:02:19  haselbac
! Added RFLU support
!
! Revision 1.5  2003/09/17 21:05:13  fnajjar
! Added infrastructure for skewed Log distribution in injection model
!
! Revision 1.4  2003/09/13 20:14:22  fnajjar
! Added infrastructure for Breakup model
!
! Revision 1.3  2003/09/10 23:35:50  fnajjar
! Removed flags that are subsumed with Rocinteract
!
! Revision 1.2  2003/01/24 19:41:55  f-najjar
! Bug fix for nPclsBuffTot to read vals(14)
!
! Revision 1.1  2002/10/25 14:19:16  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************

