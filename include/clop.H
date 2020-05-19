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

// Utility Functions
vector<string> Vectize(const char **);
vector<string> Vectize(const char **,int);
const string stripdir(const string &);

class
clop
{
protected:
  string       _longop;
  char         _shortop;
public:
  clop()
    : _shortop((char)0)
  {};
  clop(const clop &op)
    : _longop(op._longop),_shortop(op._shortop)
  {};
  clop(const string &lo,const char &so)
    : _longop(lo),_shortop(so)
  {};
  bool operator==(const clop &op) const
  {
    return(_longop  == op._longop  && 
           _shortop == op._shortop);
  };
  bool operator==(const string &s) const
  {
    return(s == _longop);
  };
  clop &operator=(const clop &op)
  {
    _longop = op._longop;
    _shortop = op._shortop;
    return(*this);
  };
  string &longop()
  {
    return(_longop);
  };
  const string &longop() const
  {
    return(_longop);
  };
  char &shortop()
  {
    return(_shortop);
  };
  const char &shortop() const
  {
    return(_shortop);
  };
};

bool
GetOp(const string &ops,const vector<string> &args);

bool
GetOp(const string &ops,string &rv,const vector<string> &args);

bool
CheckOps(const vector<string> &args,string &badop);

void
AddOp(const string &lo,const char &so);

// ADD UNDERSCORE IF NEEDED
#ifdef FORTRAN_UNDERSCORE
#  define F_FUNC(x)  x ## _
#else
#  define F_FUNC(x)  x
#endif

#ifdef FORTRAN_LOWERCASE 

// FUNCTIONS FOR LOWERCASE FORTRAN SYMBOLS
#define RFLUCONV            F_FUNC(rfluconv)
#define RFLUCLONE           F_FUNC(rfluclone)
#define RFLUEXTR            F_FUNC(rfluextr)
#define RFLUINIT            F_FUNC(rfluinit)
#define RFLUMAP             F_FUNC(rflumap)
#define RFLUMP              F_FUNC(rflump)
#define RFLUPART            F_FUNC(rflupart)
#define RFLUPICK            F_FUNC(rflupick)
#define RFLUPOST            F_FUNC(rflupost)

#else 

// FUNCTIONS FOR UPPERCASE FORTRAN SYMBOLS
#define RFLUCONV            F_FUNC(RFLUCONV)
#define RFLUCLONE           F_FUNC(RFLUCLONE)
#define RFLUEXTR            F_FUNC(RFLUEXTR)
#define RFLUINIT            F_FUNC(RFLUINIT)
#define RFLUMAP             F_FUNC(RFLUMAP)
#define RFLUMP              F_FUNC(RFLUMP)
#define RFLUPART            F_FUNC(RFLUPART)
#define RFLUPICK            F_FUNC(RFLUPICK)
#define RFLUPOST            F_FUNC(RFLUPOST)

#endif


