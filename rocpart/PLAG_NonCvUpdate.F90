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
! Purpose: update step for non conserved variables.
!
! Description: none.
!
! Input: region = current region
!        iReg   = current region number
!
! Output: region%plag = plag variables
!
! Notes: This corresponds to Part IX Step 26 in RocfluidMP framework.
!
!******************************************************************************
!
! $Id: PLAG_NonCvUpdate.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_NonCvUpdate( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMixture, ONLY    : t_mixt_input
  USE ModPartLag, ONLY    : t_plag

  USE PLAG_ModInterfaces, ONLY: PLAG_CalcDerivedVariables,   &
                                PLAG_IntrpMixtProperties
  USE ModInterfacesLagrangian, ONLY: PLAG_RFLU_ComputeVolFrac
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters         
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: region
  
! ... loop variables

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iPcls

  INTEGER      :: iFile
  REAL(RFREAL) :: diamL, massL, tauLR

  TYPE(t_plag),   POINTER :: pPlag  
  TYPE(t_global), POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_NonCvUpdate.F90,v $ $Revision: 1.1.1.1 $'
  
  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_NonCvUpdate',__FILE__ )
  
! get dimensions and set pointer ----------------------------------------------

  pPlag => region%plag

! Calculate derived variables and interpolate mixture properties --------------
!  Active for non-null number of particles in region --------------------------

  IF ( pPlag%nPcls > 0 ) THEN
            
! Get derived variables -------------------------------------------------------

!    WRITE(STDOUT,'(A)') '    Entering PLAG_CalcDerivedVariables'
    CALL PLAG_CalcDerivedVariables( region )
  
! Invoke interpolation for mixture properties ---------------------------------

!    WRITE(STDOUT,'(A)') '    Entering PLAG_IntrpMixtProperties'
    CALL PLAG_IntrpMixtProperties( region )

! Compute eulerian field ------------------------------------------------------
    IF ( global%moduleType == MODULE_TYPE_SOLVER ) THEN
      CALL PLAG_RFLU_ComputeVolFrac(region)
    END IF ! global%moduleType
! END TEMPORARY

  END IF ! nPcls

! TEMPORARY
! 
!    iFile = 700
!    DO iPcls = 1, 4 
!      massL = SUM( pPlag%cv(pPlag%cvPlagMass(:),iPcls) )
!      diamL = pPlag%dv(DV_PLAG_DIAM,iPcls)
!
!      IF ( diamL > 1.0E-14_RFREAL ) THEN
!        tauLR = 3.0_RFREAL*global%pi*pPlag%tv(TV_PLAG_MUELMIXT,iPcls)*diamL/massL
!      ELSE
!        tauLR = 0.0_RFREAL
!      ENDIF
! 
!      iFile = iFile+1
!      WRITE(iFile,'(1PE12.5,2X,I5,20(2X,1PE28.19))') global%currentTime+ global%dtMin, &
!            pPlag%aiv(AIV_PLAG_ICELLS,iPcls), &
!            pPlag%cv(CV_PLAG_XPOS:CV_PLAG_ZPOS,iPcls), &
!            pPlag%cv(pPlag%cvPlagMass(:),iPcls),&
!            massL,diamL,1.0_RFREAL/tauLR,pPlag%tv(TV_PLAG_MUELMIXT,iPcls),&
!            pPlag%cv(CV_PLAG_XMOM:CV_PLAG_ENER,iPcls), &
!            pPlag%dv(DV_PLAG_UVEL:DV_PLAG_WVEL,iPcls), &
!            pPlag%dv(DV_PLAG_UVELMIXT:DV_PLAG_WVELMIXT,iPcls)
!    END DO ! iPcls 
! 
! END TEMPORARY

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_NonCvUpdate

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_NonCvUpdate.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.4  2008/12/06 08:43:51  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/04/16 23:22:44  fnajjar
! Deleted all RFLO-pertinent calls and statements
!
! Revision 1.1  2007/04/09 18:50:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 20:57:53  fnajjar
! Initial revision after changing case
!
! Revision 1.9  2004/11/14 19:47:43  haselbac
! Changed interface
!
! Revision 1.8  2004/02/26 21:02:16  haselbac
! Removed iReg and iStage arguments, commented out writing to files
!
! Revision 1.7  2004/02/06 21:18:20  fnajjar
! Initial Integration of Rocpart with Rocflu
!
! Revision 1.6  2003/11/21 22:43:54  fnajjar
! Commented out PLAG output files
!
! Revision 1.5  2003/04/18 16:13:14  fnajjar
! Added iReg=10 and error trap for tauLR in I/O
!
! Revision 1.4  2003/04/14 20:16:51  fnajjar
! Included iReg=1 only to write debug file
!
! Revision 1.3  2003/04/14 19:46:52  fnajjar
! Added particle diameter and response time to file output
!
! Revision 1.2  2003/03/25 22:54:28  fnajjar
! Added dtMin in WRITE stamp to align with proper time stamp
!
! Revision 1.1  2003/02/04 19:10:11  f-najjar
! Initial Import
!
!******************************************************************************

