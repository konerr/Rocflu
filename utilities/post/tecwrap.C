/* *******************************************************************
 * Illinois Open Source License                                      *
 *                                                                   *
 * University of Illinois/NCSA                                       * 
 * Open Source License                                               *
 *                                                                   *
 * Copyright@2008, University of Illinois.  All rights reserved.     *
 *                                                                   *
 *  Developed by:                                                    *
 *                                                                   *
 *     Center for Simulation of Advanced Rockets                     *
 *                                                                   *
 *     University of Illinois                                        *
 *                                                                   *
 *     www.csar.uiuc.edu                                             *
 *                                                                   *
 * Permission is hereby granted, free of charge, to any person       *
 * obtaining a copy of this software and associated documentation    *
 * files (the "Software"), to deal with the Software without         *
 * restriction, including without limitation the rights to use,      *
 * copy, modify, merge, publish, distribute, sublicense, and/or      *
 * sell copies of the Software, and to permit persons to whom the    *
 * Software is furnished to do so, subject to the following          *
 * conditions:                                                       *
 *                                                                   *
 *                                                                   *
 * @ Redistributions of source code must retain the above copyright  * 
 *   notice, this list of conditions and the following disclaimers.  *
 *                                                                   * 
 * @ Redistributions in binary form must reproduce the above         *
 *   copyright notice, this list of conditions and the following     *
 *   disclaimers in the documentation and/or other materials         *
 *   provided with the distribution.                                 *
 *                                                                   *
 * @ Neither the names of the Center for Simulation of Advanced      *
 *   Rockets, the University of Illinois, nor the names of its       *
 *   contributors may be used to endorse or promote products derived * 
 *   from this Software without specific prior written permission.   *
 *                                                                   *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
 * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
 * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
 *********************************************************************
 * Please acknowledge The University of Illinois Center for          *
 * Simulation of Advanced Rockets in works and publications          *
 * resulting from this software or its derivatives.                  *
 *********************************************************************/
// $Id: tecwrap.C,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
//
// Purpose: Wrappers for TECPLOT functions so can use TECPLOT even when have 
//   conflict with lower/uppercase and trailing underscores.
//
// Description: None
//
// Notes: None
//
// Copyright: (c) 2005 by the University of Illinois

#include <string>

using namespace std;

#include "TECIO.h"

#undef TECINI100
#undef TECZNE100
#undef TECDAT100
#undef TECNOD100
#undef TECEND100
#undef TECFIL100


extern "C" {
  int TECINI100(char *Title,
                char *Variables,
                char *FName,
                char *ScratchDir,
                int  *Debug,
                int  *VIsDouble, 
                int TitleLen,
                int VariablesLen, 
                int FNameLen, 
                int ScratchDirLen);
  int TECINI100_(char *Title,
                 char *Variables,
                 char *FName,
                 char *ScratchDir,
                 int  *Debug,
                 int  *VIsDouble, 
                 int TitleLen,
                 int VariablesLen, 
                 int FNameLen, 
                 int ScratchDirLen);
  int TECZNE100(char *ZoneTitle,
                int *ZoneType,
                int *IMxOrNumPts,
                int *JMxOrNumElements,
                int *KMx,
                int *ICellMx,
                int *JCellMx,
                int *KCellMx,
                int *IsBlock,
                int *NumFaceConnections,
                int *FaceNeighborMode,
                int *ValueLocation,
                int *ShareVarFromZone,
                int *ShareConnectivityFromZone, 
                int ZoneTitleLen);
  int TECZNE100_(char *ZoneTitle,
                 int *ZoneType,
                 int *IMxOrNumPts,
                 int *JMxOrNumElements,
                 int *KMx,
                 int *ICellMx,
                 int *JCellMx,
                 int *KCellMx,
                 int *IsBlock,
                 int *NumFaceConnections,
                 int *FaceNeighborMode,
                 int *ValueLocation,
                 int *ShareVarFromZone,
                 int *ShareConnectivityFromZone, 
                 int ZoneTitleLen);
  int TECDAT100(int  *N,
                void *FieldData,
                int  *IsDouble);
  int TECDAT100_(int  *N,
                 void *FieldData,
                 int  *IsDouble);
  int TECNOD100(int *NData);
  int TECNOD100_(int *NData);
  int TECEND100(void);
  int TECEND100_(void);
  int TECFIL100(int *F);
  int TECFIL100_(int *F);
};


// 
// TECINI100
//

int TECINI100(char *Title,
              char *Variables,
              char *FName,
              char *ScratchDir,
              int  *Debug,
              int  *VIsDouble, 
              int TitleLen,
              int VariablesLen, 
              int FNameLen, 
              int ScratchDirLen) {
              
  string TitleS(Title,TitleLen);
  string VariablesS(Variables,VariablesLen);
  string FNameS(FName,FNameLen);
  string ScratchDirS(ScratchDir,ScratchDirLen);                 
              
  return (tecini100((char *)TitleS.c_str(),
                    (char *)VariablesS.c_str(),
                    (char *)FNameS.c_str(),
                    (char *)ScratchDirS.c_str(),
                    Debug,
                    VIsDouble));                        
}

int TECINI100_(char *Title,
               char *Variables,
               char *FName,
               char *ScratchDir,
               int  *Debug,
               int  *VIsDouble, 
               int TitleLen,
               int VariablesLen, 
               int FNameLen, 
               int ScratchDirLen) {

  string TitleS(Title,TitleLen);
  string VariablesS(Variables,VariablesLen);
  string FNameS(FName,FNameLen);
  string ScratchDirS(ScratchDir,ScratchDirLen);                 
              
  return (tecini100((char *)TitleS.c_str(),
                    (char *)VariablesS.c_str(),
                    (char *)FNameS.c_str(),
                    (char *)ScratchDirS.c_str(),
                    Debug,
                    VIsDouble));
}

//
// TECZNE100
//

int TECZNE100(char *ZoneTitle,
              int *ZoneType,
              int *IMxOrNumPts,
              int *JMxOrNumElements,
              int *KMx,
              int *ICellMx,
              int *JCellMx,
              int *KCellMx,
              int *IsBlock,
              int *NumFaceConnections,
              int *FaceNeighborMode,
              int *ValueLocation,
              int *ShareVarFromZone,
              int *ShareConnectivityFromZone, 
              int ZoneTitleLen) {

  string ZoneTitleS(ZoneTitle,ZoneTitleLen);
  
  return (teczne100((char *)ZoneTitle,
                    ZoneType,
                    IMxOrNumPts,
                    JMxOrNumElements,
                    KMx,
                    ICellMx,
                    JCellMx,
                    KCellMx,
                    IsBlock,
                    NumFaceConnections,
                    FaceNeighborMode,
                    ValueLocation,
                    ShareVarFromZone,
                    ShareConnectivityFromZone));

}

int TECZNE100_(char *ZoneTitle,
               int *ZoneType,
               int *IMxOrNumPts,
               int *JMxOrNumElements,
               int *KMx,
               int *ICellMx,
               int *JCellMx,
               int *KCellMx,
               int *IsBlock,
               int *NumFaceConnections,
               int *FaceNeighborMode,
               int *ValueLocation,
               int *ShareVarFromZone,
               int *ShareConnectivityFromZone, 
               int ZoneTitleLen) {

  string ZoneTitleS(ZoneTitle,ZoneTitleLen);
  
  return (teczne100((char *)ZoneTitle,
                    ZoneType,
                    IMxOrNumPts,
                    JMxOrNumElements,
                    KMx,
                    ICellMx,
                    JCellMx,
                    KCellMx,
                    IsBlock,
                    NumFaceConnections,
                    FaceNeighborMode,
                    ValueLocation,
                    ShareVarFromZone,
                    ShareConnectivityFromZone));
}

//
// TECDAT100
//

int TECDAT100(int  *N,
              void *FieldData,
              int  *IsDouble) {
  
  return (tecdat100(N,
                    FieldData,
                    IsDouble));
}

int TECDAT100_(int  *N,
               void *FieldData,
               int  *IsDouble) {

  return (tecdat100(N,
                    FieldData,
                    IsDouble));
}

//
// TECNOD100
//

int TECNOD100(int *NData) {

  return(tecnod100(NData));

}

int TECNOD100_(int *NData) {

  return(tecnod100(NData));

}

//
// TECEND100
//

int TECEND100(void) {

  return(tecend100());

}

int TECEND100_(void) {

  return(tecend100());

}


//
// TECFIL100
//

int TECFIL100(int *F) {

  return(tecfil100(F));

}

int TECFIL100_(int *F) {

  return(tecfil100(F));

}



// RCS Revision history:
//
// $Log: tecwrap.C,v $
// Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
// merged rocflu micro and macro
//
// Revision 1.1.1.1  2014/07/15 14:31:37  brollin
// New Stable version
//
// Revision 1.3  2008/12/06 08:43:59  mtcampbe
// Updated license.
//
// Revision 1.2  2008/11/19 22:17:13  mtcampbe
// Added Illinois Open Source License/Copyright
//
// Revision 1.1  2007/04/09 18:58:09  haselbac
// Initial revision after split from RocfloMP
//
// Revision 1.1  2005/05/03 20:35:32  haselbac
// Initial revision
//

