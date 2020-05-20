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
! Purpose: Write out user input for species.
!
! Description: None.
!
! Input:
!   region                Data associated with region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: SPEC_PrintUserInput.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_PrintUserInput(region)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModSpecies, ONLY: t_spec_type
  USE ModError
  USE ModParameters

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), INTENT(IN) :: region

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iSpec
  TYPE(t_global),     POINTER :: global
  TYPE(t_spec_type), POINTER :: pSpecType

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_PrintUserInput.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction(global,'SPEC_PrintUserInput',__FILE__)

! ******************************************************************************
! Print user input
! ******************************************************************************

  WRITE(STDOUT,'(A,7X,A,1X,I2)') SOLVER_NAME,'Number of species:', &
                                 region%specInput%nSpecies
  WRITE(STDOUT,'(A,7X,A,1X,I2)') SOLVER_NAME,'Number of EE species:', &
                                 region%specInput%nSpeciesEE                                 
  WRITE(STDOUT,'(A,7X,A,1X,L1)') SOLVER_NAME,'Source flag:', &
                                 region%specInput%sourceFlag
  WRITE(STDOUT,'(A,7X,A,1X,I1)') SOLVER_NAME,'Eq Eul velocity index:', &
                                 region%mixtInput%indSd

  DO iSpec = 1,region%specInput%nSpecies
    pSpecType => region%specInput%specType(iSpec)

    WRITE(STDOUT,'(A,7X,A,1X,I2)') SOLVER_NAME,'Species number:',iSpec
    WRITE(STDOUT,'(A,9X,A,1X,A)')  SOLVER_NAME,'Material:', &
                                   TRIM(pSpecType%pMaterial%name)
    WRITE(STDOUT,'(A,9X,A,1X,L1)') SOLVER_NAME,'Frozen flag:', &
                                   pSpecType%frozenFlag
    WRITE(STDOUT,'(A,9X,A,1X,I1)') SOLVER_NAME,'Source type:', &
                                   pSpecType%sourceType
    WRITE(STDOUT,'(A,9X,A,1X,ES13.6)') SOLVER_NAME,'Initial value:', &
                                       pSpecType%initVal
    WRITE(STDOUT,'(A,9X,A,1X,ES13.6)') SOLVER_NAME,'Schmidt number:', &
                                       pSpecType%schmidtNumber
                                       
    IF ( pSpecType%discreteFlag .EQV. .TRUE. ) THEN
      WRITE(STDOUT,'(A,9X,A)') SOLVER_NAME,'Discrete species:'     
      
      SELECT CASE ( pSpecType%velocityMethod )
        CASE ( SPEC_METHV_FLUIDVEL )
          WRITE(STDOUT,'(A,11X,A,1X,A)') SOLVER_NAME,'Velocity Method:', &
                                         'Fluid Velocity'
        CASE ( SPEC_METHV_EQEUL )
          WRITE(STDOUT,'(A,11X,A,1X,A)') SOLVER_NAME,'Velocity Method:', &
                                         'Equilibrium Eulerian Velocity'
          WRITE(STDOUT,'(A,11X,A,2(1X,I2))') SOLVER_NAME, & 
                                             'EE Species numbers:', &
                                             pSpecType%iSpec2iSpecEEv, & 
                                             pSpecType%iSpecEEv2iSpec                                        
        CASE DEFAULT
          CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )
      END SELECT ! velocityMethod

      WRITE(STDOUT,'(A,11X,A,1X,ES13.6)') SOLVER_NAME,'Diameter:', &
                                          pSpecType%diameter
      WRITE(STDOUT,'(A,11X,A,1X,ES13.6)') SOLVER_NAME,'Puff factor:', &
                                          pSpecType%puffFactor
      WRITE(STDOUT,'(A,11X,A,1X,ES13.6)') SOLVER_NAME,'Effective Density:', &
                                          pSpecType%effectiveDensity
      WRITE(STDOUT,'(A,11X,A,1X,ES13.6)') SOLVER_NAME,'Effective Volume:', &
                                          pSpecType%effectiveVolume
      WRITE(STDOUT,'(A,11X,A,1X,ES13.6)') SOLVER_NAME,'Material Volume:', &
                                          pSpecType%materialVolume
      WRITE(STDOUT,'(A,11X,A,1X,ES13.6)') SOLVER_NAME,'Coefficient for tau:',&
                                          pSpecType%tauCoefficient
      WRITE(STDOUT,'(A,11X,A,1X,L1)')     SOLVER_NAME,'Settling velocity:',&
                                          pSpecType%settlingFlag                                      
    ELSE
      WRITE(STDOUT,'(A,9X,A)') SOLVER_NAME,'Continuous species.'
    END IF ! discreteFlag
  END DO ! iSpec

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE SPEC_PrintUserInput

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_PrintUserInput.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:53  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:05  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:51:23  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.8  2005/11/27 01:55:11  haselbac
! Added printing of EEv variables
!
! Revision 1.7  2005/11/10 02:36:06  haselbac
! Added printing of settlingFlag, clean-up
!
! Revision 1.6  2004/07/30 22:47:37  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.5  2004/07/23 22:43:17  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.4  2004/04/01 21:31:58  haselbac
! Added printing of sourceFlag, cosmetics
!
! Revision 1.3  2004/02/02 22:41:50  haselbac
! Added printing of material name
!
! Revision 1.2  2004/01/29 22:59:36  haselbac
! Added printing of Schmidt number
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
! ******************************************************************************

