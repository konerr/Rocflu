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
! Purpose: Compute particle volume fraction for Eulerian and Lagrangian fields.
!
! Description: none.
!
! Input: pRegion = current region.
!
! Output: plag%vFracE and plag%vFracL
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_RFLU_ComputeVolFrac_Orig_17Feb2014.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_ComputeVolFrac(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModMixture, ONLY    : t_mixt
  USE ModError  
  USE ModParameters
  USE PLAG_ModParameters

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  LOGICAL :: flagCP,useInterpolation
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,fnDir,fnDirEnd,icg,icg2,icgMax,iCont,iPcl,isl,nCont, &
             nCount,nPcls 
  INTEGER, POINTER, DIMENSION(:) :: cvMass    
  REAL(RFREAL) :: alpha,densMixt,dx,dy,dz,maxVolFrac,spLoadL,sumTerms,term, &
                  volFrac,volL,volMixt
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_ComputeVolFrac_Orig_17Feb2014.F90,v $'

  global => pRegion%global
  
  CALL RegisterFunction( global, 'PLAG_RFLU_ComputeVolFrac',__FILE__ )

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

  nCont = pRegion%plagInput%nCont

  cvMass => pPlag%cvPlagMass

  alpha = 1.0_RFREAL

  nCount = 0
  flagCP = .FALSE. ! particle exceeding close-packing flag

  useInterpolation = .FALSE. ! TEMPORARY: Manoj: not using interpolation

! ******************************************************************************
! Compute volume fraction if there are particles in the region
! ******************************************************************************

  IF ( pPlag%nPcls > 0 ) THEN

! ==============================================================================
!   Reinitialize Eulerian field
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot
      pPlag%vFracE(1,icg) = 0.0_RFREAL
    END DO ! icg

! ==============================================================================
!   Compute volume of each particle
! ==============================================================================

    DO iPcl = 1,pPlag%nPcls
      volL = 0.0_RFREAL 
      DO iCont = 1, nCont
        volL = volL + pPlag%cv(cvMass(iCont),iPcl)/pRegion%plagInput%dens(iCont)
      END DO ! iCont   
        
      spLoadL = pPlag%arv(ARV_PLAG_SPLOAD,iPcl)

      pPlag%vFracL(1,iPcl) = volL*spLoadL
    END DO ! iPcl

! ==============================================================================
!   Compute Eulerian field for volume of particles
! ==============================================================================

    IF ( (useInterpolation) .AND. (pRegion%mixtInput%spaceOrder > 1) .AND. &
         (global%moduleType == MODULE_TYPE_SOLVER) ) THEN
      SELECT CASE ( pRegion%mixtInput%stencilDimensCells )
      CASE ( 1 ) 
        SELECT CASE ( pRegion%mixtInput%dimens ) 
          CASE ( 1 ) 
            fnDirEnd = XCOORD
          CASE ( 2 )
            fnDirEnd = YCOORD
          CASE ( 3 )       
            fnDirEnd = ZCOORD        
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pRegion%mixtInput%dimens

        DO iPcl = 1,pPlag%nPcls
          icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcl)
      
          dx = pGrid%cofg(XCOORD,icg) - pPlag%cv(CV_PLAG_XPOS,iPcl)
          dy = pGrid%cofg(YCOORD,icg) - pPlag%cv(CV_PLAG_YPOS,iPcl)
          dz = pGrid%cofg(ZCOORD,icg) - pPlag%cv(CV_PLAG_ZPOS,iPcl)

          IF ( (dx*dx + dy*dy + dz*dz) == 0.0_RFREAL ) THEN 
            pPlag%vFracE(1,icg) = pPlag%vFracE(1,icg) + pPlag%vFracL(1,iPcl)
          ELSE
            term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz) ! 1/distance**2 

            sumTerms =  term**alpha

            DO fnDir = XCOORD,fnDirEnd
              DO isl = 1,pGrid%c2cs1D(fnDir,icg)%nCellMembs
                icg2 = pGrid%c2cs1D(fnDir,icg)%cellMembs(isl)
        
                dx = pGrid%cofg(XCOORD,icg2) - pPlag%cv(CV_PLAG_XPOS,iPcl)
                dy = pGrid%cofg(YCOORD,icg2) - pPlag%cv(CV_PLAG_YPOS,iPcl)
                dz = pGrid%cofg(ZCOORD,icg2) - pPlag%cv(CV_PLAG_ZPOS,iPcl)

                term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz) ! 1/distance**2 

                sumTerms = sumTerms + term**alpha
              END DO ! isl
            END DO ! fnDir

            dx = pGrid%cofg(XCOORD,icg) - pPlag%cv(CV_PLAG_XPOS,iPcl)
            dy = pGrid%cofg(YCOORD,icg) - pPlag%cv(CV_PLAG_YPOS,iPcl)
            dz = pGrid%cofg(ZCOORD,icg) - pPlag%cv(CV_PLAG_ZPOS,iPcl)

            term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz) ! 1/distance**2 

            pPlag%vFracE(1,icg) = pPlag%vFracE(1,icg) &
                                + pPlag%vFracL(1,iPcl)*(term**alpha)/sumTerms

            DO fnDir = XCOORD,fnDirEnd
              DO isl = 1,pGrid%c2cs1D(fnDir,icg)%nCellMembs
                icg2 = pGrid%c2cs1D(fnDir,icg)%cellMembs(isl)
        
                dx = pGrid%cofg(XCOORD,icg2) - pPlag%cv(CV_PLAG_XPOS,iPcl)
                dy = pGrid%cofg(YCOORD,icg2) - pPlag%cv(CV_PLAG_YPOS,iPcl)
                dz = pGrid%cofg(ZCOORD,icg2) - pPlag%cv(CV_PLAG_ZPOS,iPcl)

                term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz) ! 1/distance**2 

                pPlag%vFracE(1,icg2) = pPlag%vFracE(1,icg2) &
                                     + pPlag%vFracL(1,iPcl)*(term**alpha)/sumTerms
              END DO ! isl
            END DO ! fnDir
          END IF ! (dx*dx + dy*dy + dz*dz)
        END DO ! iPcl
      CASE ( 2,3 ) 
        DO iPcl = 1,pPlag%nPcls
          icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcl)
       
          dx = pGrid%cofg(XCOORD,icg) - pPlag%cv(CV_PLAG_XPOS,iPcl)
          dy = pGrid%cofg(YCOORD,icg) - pPlag%cv(CV_PLAG_YPOS,iPcl)
          dz = pGrid%cofg(ZCOORD,icg) - pPlag%cv(CV_PLAG_ZPOS,iPcl)

          IF ( (dx*dx + dy*dy + dz*dz) == 0.0_RFREAL ) THEN 
            pPlag%vFracE(1,icg) = pPlag%vFracE(1,icg) + pPlag%vFracL(1,iPcl)
          ELSE
            term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz) ! 1/distance**2 

            sumTerms =  term**alpha

            DO isl = 1,pGrid%c2cs(icg)%nCellMembs 
              icg2 = pGrid%c2cs(icg)%cellMembs(isl)
        
              dx = pGrid%cofg(XCOORD,icg2) - pPlag%cv(CV_PLAG_XPOS,iPcl)
              dy = pGrid%cofg(YCOORD,icg2) - pPlag%cv(CV_PLAG_YPOS,iPcl)
              dz = pGrid%cofg(ZCOORD,icg2) - pPlag%cv(CV_PLAG_ZPOS,iPcl)

              term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz) ! 1/distance**2 

              sumTerms = sumTerms + term**alpha
            END DO ! isl

            dx = pGrid%cofg(XCOORD,icg) - pPlag%cv(CV_PLAG_XPOS,iPcl)
            dy = pGrid%cofg(YCOORD,icg) - pPlag%cv(CV_PLAG_YPOS,iPcl)
            dz = pGrid%cofg(ZCOORD,icg) - pPlag%cv(CV_PLAG_ZPOS,iPcl)

            term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz) ! 1/distance**2 

            pPlag%vFracE(1,icg) = pPlag%vFracE(1,icg) &
                                + pPlag%vFracL(1,iPcl)*(term**alpha)/sumTerms

            DO isl = 1,pGrid%c2cs(icg)%nCellMembs 
              icg2 = pGrid%c2cs(icg)%cellMembs(isl)

              dx = pGrid%cofg(XCOORD,icg2) - pPlag%cv(CV_PLAG_XPOS,iPcl)
              dy = pGrid%cofg(YCOORD,icg2) - pPlag%cv(CV_PLAG_YPOS,iPcl)
              dz = pGrid%cofg(ZCOORD,icg2) - pPlag%cv(CV_PLAG_ZPOS,iPcl)

              term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz) ! 1/distance**2 

              pPlag%vFracE(1,icg2) = pPlag%vFracE(1,icg2) &
                                   + pPlag%vFracL(1,iPcl)*(term**alpha)/sumTerms
            END DO ! isl
          END IF ! (dx*dx + dy*dy + dz*dz)
        END DO ! iPcl
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pRegion%mixtInput%stencilDimensCells
    ELSE
      DO iPcl = 1,pPlag%nPcls
        icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcl)

        pPlag%vFracE(1,icg) = pPlag%vFracE(1,icg) + pPlag%vFracL(1,iPcl)
      END DO ! iPcl
    END IF ! pRegion%mixtInput%spaceOrder

! ==============================================================================
!   Finalize volume fractions for Eulerian field
! ==============================================================================
    !Subbu - Debug
    DO iPcl = 1,pPlag%nPcls
      WRITE(*,*) 'iPcl,vFracL=',iPcl,pPlag%vFracL(1,iPcl)
    END DO

    maxVolFrac = 0.0_RFREAL
    icgMax = 0

    DO icg = 1,pGrid%nCellsTot
      pPlag%vFracE(1,icg) = pPlag%vFracE(1,icg)/pGrid%vol(icg)

      !Subbu - Debug
      WRITE (*,*) 'icg,vFracE=',icg,pPlag%vFracE(1,icg),pGrid%vol(icg)

      IF ( pPlag%vFracE(1,icg) > maxVolFrac ) THEN
        maxVolFrac = pPlag%vFracE(1,icg)
        icgMax = icg
      END IF ! pPlag%vFracE

      IF ( pPlag%vFracE(1,icg) > 0.6 ) THEN ! HARDCODE: Manoj
        nCount = nCount + 1
!        WRITE(*,'(I7,A,I7,1X,A,E12.6)') nCount," icg= ",icg, &
!                                     "  Particle vol frac=",pPlag%vFracE(1,icg)
        flagCP = .TRUE.
      END IF ! pPlag%vFracE
    END DO ! icg  

    IF ( flagCP .EQV. .TRUE. ) THEN 
      WRITE(*,'(A,I6,2(A,1X),E12.6,1X,A,1X,I6)') "Total ",nCount, &
                               " Cells with particle volFrac > 0.6.", &
                               " Maximum = ",maxVolFrac," in cell = ",icgMax
! TEMPORARY: Manoj: 2012-05-23: Letting code run beyong 60% volume fraction      
!      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
! END TEMPORARY
    END IF ! flagCP

! ==============================================================================
!   Reinitialize Lagrangian field
! ==============================================================================

    DO iPcl = 1,pPlag%nPcls
      pPlag%vFracL(1,iPcl) = 0.0_RFREAL
    END DO ! iPcl

! ==============================================================================
!   Compute Eulerian contribution to Lagrangian field
! ==============================================================================

    IF ( (useInterpolation) .AND. (pRegion%mixtInput%spaceOrder > 1) .AND. &
         (global%moduleType == MODULE_TYPE_SOLVER) ) THEN
      SELECT CASE ( pRegion%mixtInput%stencilDimensCells )
      CASE ( 1 ) 
        SELECT CASE ( pRegion%mixtInput%dimens ) 
          CASE ( 1 ) 
            fnDirEnd = XCOORD
          CASE ( 2 )
            fnDirEnd = YCOORD
          CASE ( 3 )       
            fnDirEnd = ZCOORD        
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pRegion%mixtInput%dimens

        DO iPcl = 1,pPlag%nPcls
          icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcl)
      
          dx = pGrid%cofg(XCOORD,icg) - pPlag%cv(CV_PLAG_XPOS,iPcl)
          dy = pGrid%cofg(YCOORD,icg) - pPlag%cv(CV_PLAG_YPOS,iPcl)
          dz = pGrid%cofg(ZCOORD,icg) - pPlag%cv(CV_PLAG_ZPOS,iPcl)

          IF ( (dx*dx + dy*dy + dz*dz) == 0.0_RFREAL ) THEN 
            pPlag%vFracL(1,iPcl) = pPlag%vFracE(1,icg)
          ELSE
            term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz) ! 1/distance**2 

            sumTerms =  term**alpha

            DO fnDir = XCOORD,fnDirEnd
              DO isl = 1,pGrid%c2cs1D(fnDir,icg)%nCellMembs
                icg2 = pGrid%c2cs1D(fnDir,icg)%cellMembs(isl)
        
                dx = pGrid%cofg(XCOORD,icg2) - pPlag%cv(CV_PLAG_XPOS,iPcl)
                dy = pGrid%cofg(YCOORD,icg2) - pPlag%cv(CV_PLAG_YPOS,iPcl)
                dz = pGrid%cofg(ZCOORD,icg2) - pPlag%cv(CV_PLAG_ZPOS,iPcl)

                term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz) ! 1/distance**2 

                sumTerms = sumTerms + term**alpha
              END DO ! isl
            END DO ! fnDir

            dx = pGrid%cofg(XCOORD,icg) - pPlag%cv(CV_PLAG_XPOS,iPcl)
            dy = pGrid%cofg(YCOORD,icg) - pPlag%cv(CV_PLAG_YPOS,iPcl)
            dz = pGrid%cofg(ZCOORD,icg) - pPlag%cv(CV_PLAG_ZPOS,iPcl)

            term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz) ! 1/distance**2 

            pPlag%vFracL(1,iPcl) = pPlag%vFracL(1,iPcl) &
                                 + pPlag%vFracE(1,icg)*(term**alpha)/sumTerms

            DO fnDir = XCOORD,fnDirEnd
              DO isl = 1,pGrid%c2cs1D(fnDir,icg)%nCellMembs
                icg2 = pGrid%c2cs1D(fnDir,icg)%cellMembs(isl)
        
                dx = pGrid%cofg(XCOORD,icg2) - pPlag%cv(CV_PLAG_XPOS,iPcl)
                dy = pGrid%cofg(YCOORD,icg2) - pPlag%cv(CV_PLAG_YPOS,iPcl)
                dz = pGrid%cofg(ZCOORD,icg2) - pPlag%cv(CV_PLAG_ZPOS,iPcl)

                term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz) ! 1/distance**2 

                pPlag%vFracL(1,iPcl) = pPlag%vFracL(1,iPcl) &
                                     + pPlag%vFracE(1,icg2)*(term**alpha)/sumTerms
              END DO ! isl
            END DO ! fnDir
          END IF ! (dx*dx + dy*dy + dz*dz)
        END DO ! iPcl
      CASE ( 2,3 ) 
        DO iPcl = 1,pPlag%nPcls
          icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcl)
       
          dx = pGrid%cofg(XCOORD,icg) - pPlag%cv(CV_PLAG_XPOS,iPcl)
          dy = pGrid%cofg(YCOORD,icg) - pPlag%cv(CV_PLAG_YPOS,iPcl)
          dz = pGrid%cofg(ZCOORD,icg) - pPlag%cv(CV_PLAG_ZPOS,iPcl)

          IF ( (dx*dx + dy*dy + dz*dz) == 0.0_RFREAL ) THEN 
            pPlag%vFracL(1,iPcl) = pPlag%vFracE(1,icg)
          ELSE
            term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz) ! 1/distance**2 

            sumTerms =  term**alpha

            DO isl = 1,pGrid%c2cs(icg)%nCellMembs 
              icg2 = pGrid%c2cs(icg)%cellMembs(isl)
        
              dx = pGrid%cofg(XCOORD,icg2) - pPlag%cv(CV_PLAG_XPOS,iPcl)
              dy = pGrid%cofg(YCOORD,icg2) - pPlag%cv(CV_PLAG_YPOS,iPcl)
              dz = pGrid%cofg(ZCOORD,icg2) - pPlag%cv(CV_PLAG_ZPOS,iPcl)

              term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz) ! 1/distance**2 

              sumTerms = sumTerms + term**alpha
            END DO ! isl

            dx = pGrid%cofg(XCOORD,icg) - pPlag%cv(CV_PLAG_XPOS,iPcl)
            dy = pGrid%cofg(YCOORD,icg) - pPlag%cv(CV_PLAG_YPOS,iPcl)
            dz = pGrid%cofg(ZCOORD,icg) - pPlag%cv(CV_PLAG_ZPOS,iPcl)

            term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz) ! 1/distance**2 

            pPlag%vFracL(1,iPcl) = pPlag%vFracL(1,iPcl) &
                                 + pPlag%vFracE(1,icg)*(term**alpha)/sumTerms

            DO isl = 1,pGrid%c2cs(icg)%nCellMembs 
              icg2 = pGrid%c2cs(icg)%cellMembs(isl)

              dx = pGrid%cofg(XCOORD,icg2) - pPlag%cv(CV_PLAG_XPOS,iPcl)
              dy = pGrid%cofg(YCOORD,icg2) - pPlag%cv(CV_PLAG_YPOS,iPcl)
              dz = pGrid%cofg(ZCOORD,icg2) - pPlag%cv(CV_PLAG_ZPOS,iPcl)

              term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz) ! 1/distance**2 

              pPlag%vFracL(1,iPcl) = pPlag%vFracL(1,iPcl) &
                                   + pPlag%vFracE(1,icg2)*(term**alpha)/sumTerms
            END DO ! isl
          END IF ! (dx*dx + dy*dy + dz*dz)
        END DO ! iPcl
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pRegion%mixtInput%stencilDimensCells
    ELSE
      DO iPcl = 1,pPlag%nPcls
        icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcl)

        pPlag%vFracL(1,iPcl) = pPlag%vFracE(1,icg)
      END DO ! iPcl
    END IF ! pRegion%mixtInput%spaceOrder

! TEMPORARY: Manoj    
!WRITE(*,*) "Stopping here ..."
!STOP
!WRITE(*,*) " After ComputeVolFrac ========================================"
!    DO iPcl = 1,pPlag%nPcls
!    DO iPcl = 1,1
!      icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcl)
!
!      WRITE(*,'(2(I5,1X),4(E12.6,1X))') iPcl,icg,pPlag%vFracE(1,icg),pPlag%vFracL(1,iPcl),pPlag%gradVFracE(XCOORD,1,icg),pPlag%gradVFracL(XCOORD,1,iPcl)
!    END DO ! iPcl
! END TEMPORARY: Manoj    
  END IF ! pPlag%nPcls

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLU_ComputeVolFrac

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_ComputeVolFrac_Orig_17Feb2014.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
!******************************************************************************

