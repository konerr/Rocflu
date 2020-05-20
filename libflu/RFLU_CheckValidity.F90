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
! Purpose: Check validity of variables
!
! Description: None.
!
! Input: 
!   pRegion     Region data
!
! Output: None.
!
! Notes: 
!   1. Compute and check pressure here because, it being computed later 
!      through a call to MixtureProperties, an invalid value would not be
!      detected.
!
! ******************************************************************************
!
! $Id: RFLU_CheckValidity.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_CheckValidity(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModParameters
  USE ModMPI
  USE ModTools, ONLY: IsNan
  
  USE ModInterfaces, ONLY: MixtPerf_G_CpR, &
                           MixtPerf_P_DEoGVm2, &
                           MixtPerf_R_M, &
                           MixtPerf_T_DPR, &
                           RFLU_PrintLocInfo
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER, PARAMETER :: MAX_INVALID_LOCS = 10
  INTEGER :: icg,indCp,indMol,nLocs
  INTEGER :: loc(MAX_INVALID_LOCS,MIN_VAL:MAX_VAL)
  REAL(RFREAL) :: Eo,gamma,p,rgas,rho,rrho,t,u,v,Vm2,w
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,gv,dv
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_CheckValidity.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_CheckValidity',__FILE__)

  nLocs = 0

#ifdef ROCPROF
  CALL FPROFILER_BEGINS("RFLU::CheckValidity")
#endif

! ******************************************************************************
! Loop over cells and check for positivity
! ******************************************************************************

  cv => pRegion%mixt%cv
  gv => pRegion%mixt%gv
  dv => pRegion%mixt%dv

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol 

  DO icg = 1,pRegion%grid%nCells
    rho  = cv(CV_MIXT_DENS,icg)
    rrho = 1.0_RFREAL/rho
    u    = rrho*cv(CV_MIXT_XMOM,icg)
    v    = rrho*cv(CV_MIXT_YMOM,icg)
    w    = rrho*cv(CV_MIXT_ZMOM,icg)        
    Eo   = rrho*cv(CV_MIXT_ENER,icg)
    ! Subbu - Correction

    ! Pres & Temp already computed-Use these DV variables 
    ! Should not use perfect gas EoS to evalute p & T again
    ! when using real gas (JWL) EoS

    !rgas = MixtPerf_R_M(gv(GV_MIXT_MOL,icg*indMol))
    !gamma= MixtPerf_G_CpR(gv(GV_MIXT_CP,icg*indCp),rgas)
    !Vm2  = u*u + v*v + w*w
    !p    = MixtPerf_P_DEoGVm2(rho,Eo,gamma,Vm2)
    !t    = MixtPerf_T_DPR(rho,p,rgas)

    p    = dv(DV_MIXT_PRES,icg)
    t    = dv(DV_MIXT_TEMP,icg)

    ! Subbu - End Correction

    IF ( (IsNan(rho) .EQV. .TRUE.) .OR. & 
         (IsNan(u)   .EQV. .TRUE.) .OR. &
         (IsNan(v)   .EQV. .TRUE.) .OR. & 
         (IsNan(w)   .EQV. .TRUE.) .OR. & 
         (IsNan(p)   .EQV. .TRUE.) .OR. &
         (IsNan(t)   .EQV. .TRUE.) ) THEN
      nLocs = nLocs + 1   

      IF ( nLocs == 1 ) THEN 
        WRITE(STDOUT,'(A,1X,A,1X,I9)') SOLVER_NAME, & 
              'Invalid variables detected!'
              
        IF ( global%flowType == FLOW_UNSTEADY ) THEN 
          WRITE(STDOUT,'(A,3X,A,1X,1PE12.5)') SOLVER_NAME,'Current time:', &
                                              global%currentTime              
        ELSE 
          WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME, &
                                         'Current iteration number:', &
                                         global%currentIter           
        END IF ! global%flowType                 
                                            
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal 
        WRITE(STDOUT,'(A,6X,A,6(1X,A))') SOLVER_NAME,'#', &
                                         '   Density   ', &
                                         '  x-velocity ', &
                                         '  y-velocity ', &
                                         '  z-velocity ', &
                                         '   Pressure  ', &
                                         ' Temperature '       
      END IF ! nLocs

      IF ( nLocs <= MAX_INVALID_LOCS ) THEN 
        WRITE(STDOUT,'(A,4X,I3,6(1X,E13.6))') SOLVER_NAME,nLocs, & 
                                              rho,u,v,w,p,t
        loc(nLocs,MIN_VAL:MAX_VAL) = icg                                   
      END IF ! nLocs
    END IF ! cv    
  END DO ! icg

! ******************************************************************************
! Write out message and call error handling routine
! ******************************************************************************

  IF ( nLocs > 0 ) THEN 
    IF ( nLocs > MAX_INVALID_LOCS ) THEN 
       WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A,1X,I9,1X,A)') SOLVER_NAME, & 
             'Only wrote the first',MAX_INVALID_LOCS,'of',nLocs, & 
             'cells with invalid variables.'    
      CALL RFLU_PrintLocInfo(pRegion,loc,MAX_INVALID_LOCS, & 
                             LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
    ELSE 
      CALL RFLU_PrintLocInfo(pRegion,loc(1:nLocs,MIN_VAL:MAX_VAL),nLocs, & 
                             LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
    END IF ! nLocs
    
    CALL ErrorStop(global,ERR_INVALID_VALUE,__LINE__)   
  END IF ! nLocs

! ******************************************************************************
! End
! ******************************************************************************

#ifdef ROCPROF
  CALL FPROFILER_ENDS("RFLU::CheckValidity")
#endif

  CALL DeregisterFunction( global )

END SUBROUTINE RFLU_CheckValidity

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CheckValidity.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:34  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:52  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:48  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.5  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.4  2006/03/26 20:21:30  haselbac
! Removed dv declaration and definition
!
! Revision 1.3  2005/07/07 22:44:00  haselbac
! Added profiling calls, cosmetics
!
! Revision 1.2  2004/01/11 02:06:38  jiao
! Eliminated some redundant trailing spaces that made some lines too long.
! This changed was needed to compile with NAG F90 compiler.
!
! Revision 1.1  2003/12/04 03:23:33  haselbac
! Initial revision
!
! ******************************************************************************

