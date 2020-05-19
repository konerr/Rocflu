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
#include <iostream>
#ifndef list
#include <list>
#endif
#ifndef vector
#include <vector>
#endif
#ifndef string
#include <string>
#endif

using namespace std;

#include "clop.H"

list<clop> cloption_list;

const string
stripdir(const string &pname)
{
  string retval;
  string::size_type x = pname.find("/");
  if(x == string::npos)
    return(pname);
  return(pname.substr(pname.find_last_of("/")+1));
}

// Adds a command line option to the list of possibilities
void
AddOp(const string &lo,const char &so)
{
  cloption_list.push_back(clop(lo,so));
}

// Checks to see if the command line options are
// correct - if there are any bad options, it reports
// the first one encountered
bool
CheckOps(const vector<string> &args,string &badop)
{
  vector<string>::const_iterator ai = args.begin();
  ai++;
  while(ai != args.end()){
    string option_string;
    bool found = false;
    bool longo = false;
    if((*ai)[0] == '-'){
      if((*ai)[1] == '-'){
	longo = true;
	option_string = (*ai).substr(2);
      }
      list<clop>::const_iterator oi = cloption_list.begin();
      while(oi != cloption_list.end() && !found){
	if(longo){
	  if(oi->longop() == option_string)
	    found = true;
	}
	else
	  if((*ai)[1] == oi->shortop())
	    found =  true;
	oi++;
      }
      if(!found){
	badop = *ai;
	return(false);
      }
    }
    ai++;
  } 
  return(true);
}
    
// Checks to see if a boolean option is set, returns true if so
bool
GetOp(const string &ops,const vector<string> &args)
{
  list<clop>::const_iterator oi = cloption_list.begin();
  bool found = false;
  while(oi != cloption_list.end() && !found){
    if(*oi == ops)
      found = true;
    else
      oi++;
  }
  if(!found)
    return(false);
  vector<string>::const_iterator ai = args.begin();
  while(ai != args.end()){
    if((*ai)[0] == '-'){
      if((*ai)[1] != '-'){
	if((*ai)[1] == oi->shortop()){
	  return true;
	}
      }
      else if((*ai)[1] == '-'){
	string op((*ai).substr(2));
	if(op == oi->longop()){
	  return true;
	}
      }
    }
    ai++;
  }
  return false;
}

// Checks to see if an option is set, and sets the option argument
// if it exists
bool
GetOp(const string &ops,string &rv,const vector<string> &args)
{
  rv.erase();
  list<clop>::const_iterator oi = cloption_list.begin();
  bool found = false;
  while(oi != cloption_list.end() && !found){
    if(*oi == ops)
      found = true;
    else
      oi++;
  }
  if(!found){
    return(false);
  }
  vector<string>::const_iterator ai = args.begin();
  while(ai != args.end()){
    if((*ai)[0] == '-'){
      if((*ai)[1] != '-'){
	if((*ai)[1] == oi->shortop()){
	  ai++;
	  if(ai == args.end())
	    return true;
	  else {
	    if((*ai)[0] != '-')
	      rv = *ai;
	    return true;
	  }
	}
      }
      else if((*ai)[1] == '-'){
	string op((*ai).substr(2));
	if(op == oi->longop()){
	  ai++;
	  if(ai == args.end())
	    return true;
	  else {
	    if((*ai)[0] != '-')
	      rv = *ai;
	    return true;
	  }
	}
      }
    }
    ai++;
  }
  return false;
}

// Utilities
vector<string>
Vectize(const char **in)
{
  vector<string> retVal;
  int i = 0;
  while(in[i] != NULL)
    retVal.push_back(in[i++]);
  return retVal;
}

vector<string>
Vectize(const char **in,int n)
{
  vector<string> retVal;
  if(n <= 0) return retVal;
  int i = 0;
  while((in[i] != NULL) && i < n)
    retVal.push_back(in[i++]);
  return retVal;
}

