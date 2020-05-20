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
! Purpose: set explicit interfaces to subroutines and functions
!          related to I/O.
!
! Description: none
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ModInterfacesIO.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModInterfacesIO

  IMPLICIT NONE

  INTERFACE

  INTEGER FUNCTION BuildPatchIdentifier(iRegion,iPatch)
    INTEGER, INTENT(IN) :: iPatch,iRegion
  END FUNCTION BuildPatchIdentifier

  SUBROUTINE MakeNumberedKeys(keys,indBegin,string,numBegin,numEnd,numSkip)
    CHARACTER(*), INTENT(inout) :: keys(:)
    CHARACTER(*), INTENT(in)    :: string
    INTEGER,      INTENT(in)    :: indBegin,numBegin,numEnd,numSkip
  END SUBROUTINE MakeNumberedKeys

  SUBROUTINE ReadBothRegionSection( global,fileID,nvals,nStrVals,keys, &
                                    strKeys,vals,strVals,brbeg,brend,  &
                                    defined,strDefined )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    INTEGER      :: fileID, nvals, nStrVals, brbeg, brend
    CHARACTER(*) :: keys(nvals), strKeys(nStrVals)
    LOGICAL      :: defined(nvals), strDefined(nStrVals)
    REAL(RFREAL) :: vals(nvals)
    CHARACTER(*) :: strVals(nStrVals)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadBothRegionSection

  SUBROUTINE ReadBothSection( global,fileID,nvals,nStrVals,keys,strKeys, &
                              vals,strVals,defined,strDefined )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    INTEGER      :: fileID, nvals, nStrVals
    CHARACTER(*) :: keys(nvals), strKeys(nStrVals)
    LOGICAL      :: defined(nvals), strDefined(nStrVals)
    REAL(RFREAL) :: vals(nvals)
    CHARACTER(*) :: strVals(nStrVals)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadBothSection

  SUBROUTINE ReadABCSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadABCSection

  SUBROUTINE ReadGFMSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadGFMSection

  SUBROUTINE ReadAccelerationSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadAccelerationSection

  SUBROUTINE ReadFlowSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadFlowSection

  SUBROUTINE ReadForcesSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadForcesSection

  SUBROUTINE ReadFormatsSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadFormatsSection

  SUBROUTINE ReadGridMotionSection(regions)
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadGridMotionSection

  SUBROUTINE ReadInitFlowSection(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), DIMENSION(:), POINTER :: regions
  END SUBROUTINE ReadInitFlowSection

  SUBROUTINE ReadInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadInputFile

  SUBROUTINE ReadListSection( global,fileID,key,nCols,nRows,vals,defined )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    INTEGER      :: fileID, nCols, nRows
    CHARACTER(*) :: key
    LOGICAL      :: defined
    REAL(RFREAL), POINTER :: vals(:,:)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadListSection
  
  SUBROUTINE ReadMiscSection(global)
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadMiscSection

  SUBROUTINE ReadMixtureSection(regions)
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadMixtureSection

  SUBROUTINE ReadMultigridSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadMultigridSection

  SUBROUTINE ReadMvFrameSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadMvFrameSection

  SUBROUTINE ReadNumericsSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadNumericsSection

  SUBROUTINE ReadPatchSection( global,fileID,nvals,keys,vals, &
                               prbeg,prend,distrib,fname,bcName,defined )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    INTEGER      :: fileID, nvals, prbeg, prend, distrib
    CHARACTER(*) :: keys(nvals), fname
    CHARACTER(*) :: bcName
    LOGICAL      :: defined(nvals)
    REAL(RFREAL) :: vals(nvals)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadPatchSection

  SUBROUTINE ReadPostSection(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadPostSection

  SUBROUTINE ReadPrefixedListSection( global,fileID,key,nCols,nRows, &
                                      vals,strVals,defined )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    INTEGER      :: fileID, nCols, nRows
    CHARACTER(*) :: key
    LOGICAL      :: defined
    REAL(RFREAL), POINTER :: vals(:,:)
    CHARACTER(*), POINTER :: strVals(:)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadPrefixedListSection

  SUBROUTINE ReadPrepSection(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadPrepSection

  SUBROUTINE ReadProbeSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadProbeSection

  SUBROUTINE ReadRandomSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadRandomSection

  SUBROUTINE ReadReferenceSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadReferenceSection

  SUBROUTINE ReadTimeZoomingSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadTimeZoomingSection

  SUBROUTINE ReadRocketSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadRocketSection

  SUBROUTINE ReadRegionSection( global,fileID,nvals,keys,vals, &
                                brbeg,brend,defined )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    INTEGER      :: fileID, nvals, brbeg, brend
    CHARACTER(*) :: keys(nvals)
    LOGICAL      :: defined(nvals)
    REAL(RFREAL) :: vals(nvals)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadRegionSection

  SUBROUTINE ReadSection( global,fileID,nvals,keys,vals,defined )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    INTEGER      :: fileID, nvals
    CHARACTER(*) :: keys(nvals)
    LOGICAL      :: defined(nvals)
    REAL(RFREAL) :: vals(nvals)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadSection

  SUBROUTINE ReadStringSection( global,fileID,nvals,keys,vals,defined )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    INTEGER        :: fileID, nvals
    CHARACTER(*)   :: keys(nvals)
    LOGICAL        :: defined(nvals)
    CHARACTER(*)   :: vals(nvals)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadStringSection

  SUBROUTINE ReadThrustSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadThrustSection

  SUBROUTINE ReadTimestepSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadTimestepSection

  SUBROUTINE ReadTransformSection( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE ReadTransformSection

  SUBROUTINE ReadViscositySection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE ReadViscositySection

  SUBROUTINE WriteConvergence( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE WriteConvergence

  SUBROUTINE WriteTotalMass(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE WriteTotalMass

  SUBROUTINE WriteProbe( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE WriteProbe

  SUBROUTINE WriteThrust( global )
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE WriteThrust

  END INTERFACE

END MODULE ModInterfacesIO

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModInterfacesIO.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.7  2009/08/28 18:29:46  mtcampbe
! RocfluMP integration with Rocstar and some makefile tweaks.  To build
! Rocstar with new Rocflu:
! make ROCFLU=RocfluMP
! To build Rocstar with the new RocfluND:
! make ROCFLU=RocfluMP HYPRE=/the/hypre/install/path
!
! Revision 1.6  2009/07/08 19:11:40  mparmar
! Added ReadABCSection
!
! Revision 1.5  2008/12/06 08:43:37  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:16:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/11/28 23:05:12  mparmar
! Added declaration of ReadMiscSection function
!
! Revision 1.2  2007/06/18 17:54:53  mparmar
! Added declaration of ReadMvFrameSection
!
! Revision 1.1  2007/04/09 18:49:10  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:17  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.6  2005/10/31 19:26:48  haselbac
! Added interface for ReadMixtureSection
!
! Revision 1.5  2004/06/16 20:00:54  haselbac
! Removed buildFileNameXXX routines
!
! Revision 1.4  2004/04/08 03:17:07  wasistho
! nDummyCells in Rocflo read from INITFLOW section
!
! Revision 1.3  2003/11/21 22:35:51  fnajjar
! Update Random Number Generator
!
! Revision 1.2  2003/08/28 20:05:39  jblazek
! Added acceleration terms.
!
! Revision 1.1  2003/08/11 21:50:00  jblazek
! Splitted ModInterfaces into 4 sections.
!
!******************************************************************************

