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
! Purpose: define data types related to the whole gas mixture.
!
! Description: none
!
! Notes: none
!
! ******************************************************************************
!
! $Id: ModMixture.F90,v 1.7 2016/03/22 19:32:45 fred Exp $
!
! Copyright: (c) 2001-2005 by the University of Illinois
!
! ******************************************************************************

MODULE ModMixture

  USE ModDataTypes
  IMPLICIT NONE

! mixture input ---------------------------------------------------------------

  TYPE t_mixt_input

    INTEGER :: flowModel,fluidModel
    LOGICAL :: axiFlag,moveGrid,externalBc,computeTv,frozenFlag

    INTEGER :: dimens
    INTEGER :: nCv,nCvOld2,nDv,nGv,nGvAct,nTv
    INTEGER :: nGrad
    INTEGER :: indCp,indMfMixt,indMol,indSd
    REAL(RFREAL) :: prLam,prTurb,scnLam,scnTurb

! - turbulence modeling

    INTEGER :: turbModel

! - radiation

    LOGICAL :: radiUsed

! - numerics

    INTEGER      :: spaceDiscr, spaceOrder, pSwitchType
    INTEGER      :: timeScheme, ldiss(5)
    REAL(RFREAL) :: cfl, smoocf, vis2, vis4, pSwitchOmega, limFac, epsEntr
    REAL(RFREAL) :: ark(5), grk(5), trk(5), betrk(5)
    INTEGER :: spaceOrderBFaces
    INTEGER :: cReconstCells,cReconstFaces,reconst,stencilDimensCells, &
               stencilDimensFaces,stencilDimensBFaces
    REAL(RFREAL) :: cReconstCellsWeight,cReconstFacesWeight,dissFact,tolerICT

! - flow initialization 

    INTEGER :: prepIniCase
    INTEGER :: prepIntVal1,prepIntVal2,prepIntVal3,prepIntVal4,prepIntVal5, &
               prepIntVal6,prepIntVal7,prepIntVal8,prepIntVal9,prepIntVal10, &
               prepIntVal11,prepIntVal12,prepIntVal13,prepIntVal14,prepIntVal15, &
               prepIntVal16,prepIntVal17,prepIntVal18,prepIntVal19 !Fred - Coll model
    REAL(RFREAL) :: prepRealVal1,prepRealVal2,prepRealVal3,prepRealVal4, & 
                    prepRealVal5,prepRealVal6,prepRealVal7,prepRealVal8, & 
                    prepRealVal9,prepRealVal10,prepRealVal11,prepRealVal12, &
                    prepRealVal13,prepRealVal14,prepRealVal15,prepRealVal16, &
                    prepRealVal17,prepRealVal18,prepRealVal19,prepRealVal20, &
                    prepRealVal21,prepRealVal22,prepRealVal23, & ! Rahul-shktb case
                    prepRealVal24  !Fred - Adding 24 to allow for Variable JWL density   
    REAL(RFREAL) :: iniDens,iniPress,iniTemp,iniVelX,iniVelY,iniVelZ

! - grid motion

    INTEGER :: moveGridNIter,moveGridType
    REAL(RFREAL) :: moveGridSFact

! - gas model

    INTEGER :: gasModel

! - viscosity model

    INTEGER      :: viscModel
    REAL(RFREAL) :: refReNum,refVisc,refViscTemp,suthCoef

  END TYPE t_mixt_input

! mixture data ----------------------------------------------------------------

  TYPE t_mixt
    REAL(RFREAL), POINTER :: cv(:,:), cvOld(:,:)
    REAL(RFREAL), POINTER :: dv(:,:), tv(:,:), gv(:,:)
    REAL(RFREAL), DIMENSION(:,:), POINTER :: lim
    REAL(RFREAL), POINTER :: rhs(:,:), rhsOld(:,:), rhsSum(:,:), diss(:,:), & 
                             fterm(:,:)
    LOGICAL, POINTER :: set_p_JWL(:), set_e_JWL(:)
    REAL(RFREAL), POINTER :: JWL_p_prev(:,:)
    REAL(RFREAL), POINTER :: JWL_e_prev(:,:)
 
!Fred - Adding in flags/arrays to improve initial guesses for iterative JWL
!method 

#ifdef STATS
    REAL(RFREAL), POINTER :: tav(:,:), tavVert(:,:)
#endif
    INTEGER :: cvState
    INTEGER, DIMENSION(:), POINTER :: cvInfo
    REAL(RFREAL), DIMENSION(:), POINTER :: sigma
    REAL(RFREAL), DIMENSION(:), POINTER :: delP,mfMixt,vfMixt,vfMixtOld
    REAL(RFREAL), DIMENSION(:,:), POINTER :: cvVert,dvVert,gvVert
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: gradCell,gradCellOld, &
                                               gradCellOld2,gradFace
    REAL(RFREAL), POINTER :: sd(:,:)
    REAL(RFREAL), POINTER :: cvOld1(:,:),cvOld2(:,:)
    REAL(RFREAL), POINTER :: cvRef(:,:)

    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: gradCellE

    REAL(RFREAL), DIMENSION(:,:), POINTER :: levelSet
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: gradCellLevelSet
  END TYPE t_mixt

END MODULE ModMixture

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModMixture.F90,v $
! Revision 1.7  2016/03/22 19:32:45  fred
! Adding Glasser's model as a choice for collision modelling
!
! Revision 1.6  2016/02/06 17:23:01  fred
! Cleaning up order of RVAL variable declarations
!
! Revision 1.5  2016/02/05 20:25:27  fred
! Removing need to have SPEC flag on when compiling
!
! Revision 1.4  2016/02/04 22:17:52  fred
! Cleanup of unneeded variables
!
! Revision 1.3  2016/02/04 22:11:57  fred
! Adding JWL EOS iterative capabilities for cylindrical detonation problem
!
! Revision 1.2  2016/02/03 20:36:48  rahul
! prepRealVal21, prepRealVal22, prepRealVal23 have been added to read dx,dy
! and dz of shktb case from the .inp file. This is a temporary fix.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.6  2009/07/08 19:11:41  mparmar
! Added sigma and cvRef for absorbing layer
!
! Revision 1.5  2008/12/06 08:43:37  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:16:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2008/03/27 12:11:40  haselbac
! Added axiFlag
!
! Revision 1.2  2007/11/28 23:05:13  mparmar
! Added variables for SOLV_IMPLICIT_HM
!
! Revision 1.1  2007/04/09 18:49:11  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:17  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.63  2007/04/05 00:57:30  haselbac
! Added additional RVALxy params for 2p shocktube problems
!
! Revision 1.62  2007/03/19 21:40:47  haselbac
! Removed variables related to plotting, now in ModPlotting
!
! Revision 1.61  2006/08/19 15:38:54  mparmar
! Added spaceOrderBFaces and removed bGradFace
!
! Revision 1.60  2006/04/07 14:45:25  haselbac
! Added new stencilDimens params
!
! Revision 1.59  2006/03/26 20:21:57  haselbac
! Added initTemp
!
! Revision 1.58  2006/01/12 09:42:37  wasistho
! added tavVert
!
! Revision 1.57  2006/01/06 22:09:15  haselbac
! Added entry for stencilDimens
!
! Revision 1.56  2005/12/25 15:29:34  haselbac
! Added entries for constrained reconstruction
!
! Revision 1.55  2005/12/24 21:27:19  haselbac
! Added tolerICT
!
! Revision 1.54  2005/11/17 14:38:21  haselbac
! Added prepRealVal9 and prepRealVal10
!
! Revision 1.53  2005/11/10 22:23:05  fnajjar
! ACH: Added frozenFlag
!
! Revision 1.52  2005/11/10 02:21:18  haselbac
! Added nGvAct
!
! Revision 1.51  2005/10/31 21:08:53  haselbac
! Enabled gasModel for Rocflo
!
! Revision 1.50  2005/10/31 19:27:07  haselbac
! Added gasModel
!
! Revision 1.49  2005/10/27 18:58:07  haselbac
! Added parameter for constraints
!
! Revision 1.48  2005/09/22 17:08:31  hdewey2
! Added cvOld1 and cvOld2 to type t_mixt for transient implicit solver.
!
! Revision 1.47  2005/09/09 03:18:58  wasistho
! added prepIniCase
!
! Revision 1.46  2005/08/03 17:53:09  wasistho
! added domain-splitting feature to initial condition
!
! Revision 1.45  2005/07/25 12:22:00  haselbac
! Added pvInfo
!
! Revision 1.44  2005/07/11 19:27:16  mparmar
! Added reconst and lim
!
! Revision 1.43  2005/05/01 14:20:09  haselbac
! Added nPv, pv, pvVert
!
! Revision 1.42  2005/04/20 14:40:52  haselbac
! Added int and real prep vals, removed uniform flow checking vars
!
! Revision 1.41  2005/03/31 16:57:18  haselbac
! Removed nSd
!
! Revision 1.40  2005/03/22 03:33:46  haselbac
! Added prep init helper variables, cosmetics
!
! Revision 1.39  2005/03/09 14:54:40  haselbac
! Added dimensionality variable, cosmetics
!
! Revision 1.38  2005/01/13 21:38:00  haselbac
! Added rhsOld for incompressible solver
!
! Revision 1.37  2004/12/19 15:45:22  haselbac
! Added vfMixt for incompressible solver
!
! Revision 1.36  2004/11/02 02:29:19  haselbac
! Added fluidModel and nCv
!
! Revision 1.35  2004/09/02 02:33:37  wasistho
! added face-edge averaging input-option parameter in Rocflo
!
! Revision 1.34  2004/07/30 22:47:35  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.33  2004/07/28 15:41:00  jferry
! created global variable for spec use
!
! Revision 1.32  2004/07/08 02:18:10  haselbac
! Added dissFact
!
! Revision 1.31  2004/03/05 22:09:01  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.30  2004/03/03 23:55:39  jferry
! Allowed particles to be run with Euler case
!
! Revision 1.29  2004/03/02 21:49:21  jferry
! Added inrtUsed flag to mixture data structure
!
! Revision 1.28  2004/01/29 22:57:26  haselbac
! Added indMfMixt, mfMixt, and cvInfo, removed vf
!
! Revision 1.27  2003/11/25 21:03:13  haselbac
! Added specUsed and vf variables
!
! Revision 1.26  2003/07/03 21:48:45  jblazek
! Implemented dual-time stepping.
!
! Revision 1.25  2003/05/29 17:28:43  jblazek
! Implemented Roe scheme.
!
! Revision 1.24  2003/04/10 23:24:04  fnajjar
! Added infrastructure for viscosity models in t_mixt_input
!
! Revision 1.23  2003/03/31 16:14:11  haselbac
! Added moveGridType
!
! Revision 1.22  2003/03/15 17:48:58  haselbac
! Renamed cvVrtx, added dvVert and gvVert
!
! Revision 1.21  2003/01/28 16:44:39  haselbac
! Added grid motion input variables
!
! Revision 1.20  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.19  2002/09/09 14:57:31  haselbac
! mixtInput now under regions, added boundary face gradients
!
! Revision 1.18  2002/08/15 19:48:05  jblazek
! Implemented grid deformation capability.
!
! Revision 1.17  2002/08/01 01:30:20  wasistho
! Added gradFace for RFLU
!
! Revision 1.16  2002/07/25 15:15:27  haselbac
! Changed grad to gradCell
!
! Revision 1.15  2002/07/25 00:39:01  jblazek
! Option for TVD type pressure switch.
!
! Revision 1.14  2002/06/14 21:34:32  wasistho
! Added time avg statistics
!
! Revision 1.13  2002/05/04 16:59:10  haselbac
! Added variables for uniform flow preservation check
!
! Revision 1.12  2002/04/11 18:53:15  haselbac
! Added array for vertex-based solution vector
!
! Revision 1.11  2002/03/26 19:19:09  haselbac
! Added ROCFLU variables
!
! Revision 1.10  2002/02/27 18:38:20  jblazek
! Changed extrapol. to dummy cells at injection boundaries and slip walls.
!
! Revision 1.9  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.8  2002/02/06 00:15:39  jblazek
! Improved injection BC. Added pointers to gradients.
!
! Revision 1.7  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.6  2002/01/23 03:51:24  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.5  2002/01/16 22:03:35  jblazek
! Added time-stepping routines.
!
! Revision 1.4  2002/01/10 18:21:29  jblazek
! Added iteration number and initial residual to solution file.
!
! Revision 1.3  2002/01/02 16:20:19  jblazek
! Added flow initialization and dummy cell geometry.
!
! Revision 1.2  2001/12/11 21:59:29  jblazek
! memory allocation added.
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
! ******************************************************************************

