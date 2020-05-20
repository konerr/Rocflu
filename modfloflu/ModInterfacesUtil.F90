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
! Purpose: set explicit interfaces to utility subroutines and functions.
!
! Description: none
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ModInterfacesUtil.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModInterfacesUtil

  IMPLICIT NONE

  INTERFACE

  SUBROUTINE AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvar )
    USE ModDataTypes
    INTEGER               :: ijkD, ijkN1, ijkN2, iFBeg, iFEnd
    REAL(RFREAL), POINTER :: fvar(:,:)
  END SUBROUTINE AverageVecVar

  SUBROUTINE CentroidHexa( xyzNodes,cofgX,cofgY,cofgZ )
    USE ModDataTypes
    REAL(RFREAL) :: xyzNodes(3,8)
    REAL(RFREAL) :: cofgX, cofgY, cofgZ
  END SUBROUTINE CentroidHexa

  SUBROUTINE DescaleGridSpeeds(region)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region    
  END SUBROUTINE DescaleGridSpeeds

  SUBROUTINE FaceCentroidQuad( xyzNodes,fCenX,fCenY,fCenZ )
    USE ModDataTypes
    REAL(RFREAL) :: xyzNodes(3,4)
    REAL(RFREAL) :: fCenX, fCenY, fCenZ
  END SUBROUTINE FaceCentroidQuad

  SUBROUTINE FaceCentroidTria( xyzNodes,fCenX,fCenY,fCenZ )
    USE ModDataTypes
    REAL(RFREAL) :: xyzNodes(3,3)
    REAL(RFREAL) :: fCenX, fCenY, fCenZ
  END SUBROUTINE FaceCentroidTria

  SUBROUTINE FaceVectorQuad( xyzNodes,fVecX,fVecY,fVecZ )
    USE ModDataTypes
    REAL(RFREAL) :: xyzNodes(3,4)
    REAL(RFREAL) :: fVecX, fVecY, fVecZ
  END SUBROUTINE FaceVectorQuad

  SUBROUTINE FaceVectorTria( xyzNodes,fVecX,fVecY,fVecZ )
    USE ModDataTypes
    REAL(RFREAL) :: xyzNodes(3,3)
    REAL(RFREAL) :: fVecX, fVecY, fVecZ
  END SUBROUTINE FaceVectorTria

  SUBROUTINE ScaleGridSpeeds(region)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region    
  END SUBROUTINE ScaleGridSpeeds

  SUBROUTINE ScaleRotateVector( global,vect )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    REAL(RFREAL), POINTER   :: vect(:,:)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ScaleRotateVector

  SUBROUTINE SplitQuadFace( global,xyz1,xyz2,xyz3,xyz4,splitFlag )
    USE ModDataTypes
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
    INTEGER, INTENT(OUT)    :: splitFlag
    REAL(RFREAL), DIMENSION(3), INTENT(IN) :: xyz1,xyz2,xyz3,xyz4
  END SUBROUTINE SplitQuadFace

  SUBROUTINE VolumeHexa( xyzNodes,faceVecs,volume )
    USE ModDataTypes
    REAL(RFREAL) :: xyzNodes(3,8), faceVecs(3,6)
    REAL(RFREAL) :: volume
  END SUBROUTINE VolumeHexa

  END INTERFACE

END MODULE ModInterfacesUtil

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModInterfacesUtil.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:37  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:10  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:17  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2003/08/11 21:50:00  jblazek
! Splitted ModInterfaces into 4 sections.
!
!******************************************************************************

