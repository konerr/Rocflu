################################################################################
#
# $Id: Makefile,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
#
# Purpose: Makefile for Rocflu-MP pre- and postprocessing utilities.
#
# Description: none.
#
# Notes: none.
#
# Copyright: (c) 2001-2007 by the University of Illinois
#
################################################################################

all:
	$(MAKE) -C conv UTIL=yes
	$(MAKE) -C clone UTIL=yes
	$(MAKE) -C extr UTIL=yes        
	$(MAKE) -C init UTIL=yes
	$(MAKE) -C map UTIL=yes
	$(MAKE) -C part UTIL=yes
	$(MAKE) -C pick UTIL=yes
	$(MAKE) -C post UTIL=yes

clean install:
	$(MAKE) -C conv $@
	$(MAKE) -C clone $@
	$(MAKE) -C extr $@        
	$(MAKE) -C init $@
	$(MAKE) -C map $@
	$(MAKE) -C part $@
	$(MAKE) -C pick $@
	$(MAKE) -C post $@

################################################################################
#
# RCS Revision history:
#
# $Log: Makefile,v $
# Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
# merged rocflu micro and macro
#
# Revision 1.1.1.1  2014/07/15 14:31:37  brollin
# New Stable version
#
# Revision 1.3  2007/11/27 13:11:05  haselbac
# Added entry for rfluextr
#
# Revision 1.2  2007/08/07 17:15:50  haselbac
# Added commands to compile rfluclone
#
# Revision 1.1  2007/04/09 18:45:47  haselbac
# Initial revision after split from RocfloMP
#
# Revision 1.1  2007/04/09 17:56:41  haselbac
# Initial revision after split from RocfloMP
#
# Revision 1.14  2006/01/10 21:53:07  wasistho
# make rflopost as default
#
# Revision 1.13  2005/04/15 15:05:48  haselbac
# Removed Charm/FEM stuff, cosmetics
#                                                             
# Revision 1.12  2004/10/19 19:22:24  haselbac                                   
# Changed so rflumap gets built and cleaned out if CHARM == 1                    
#
# Revision 1.11  2004/08/17 00:59:05  wasistho                                   
# prepared for utilities/rocflo/toflu                                            
#
# Revision 1.10  2004/03/19 21:16:56  haselbac                                   
# Added rinit for RFLU                                                           
#
# Revision 1.9  2003/09/25 15:40:22  jferry                                      
# Implented Rocsmoke post-processing                                             
#
# Revision 1.8  2003/06/04 21:53:15  haselbac                                    
# Removed cells and added pick                                                   
#
# Revision 1.7  2003/05/06 16:20:36  fnajjar                                     
# Included calls for Rocpart postprocessing tool                                 
#
# Revision 1.6  2003/04/11 15:15:18  jiao                                        
# Changed RFLU_CHARM back to CHARM.                                              
#
# Revision 1.5  2003/04/01 16:35:08  haselbac                                    
# Added rflucells, restriction for serial utils                                  
#
# Revision 1.4  2003/03/25 19:09:13  haselbac                                    
# Removed POST for Rocflu                                                        
#
# Revision 1.3  2003/03/21 15:36:07  haselbac                                    
# Added install target                                                           
#
# Revision 1.2  2003/03/20 23:45:48  jiao                                        
# ACH: Added rocflu/conv.                                                        
#
# Revision 1.1.1.1  2003/03/20 22:48:04  haselbac                                
# Initial revision                                                               
#
# Revision 1.14  2002/12/20 19:38:02  jblazek                                    
# Added tool to count blocks in a box.                                           
#
# Revision 1.13  2002/10/19 00:40:31  jblazek                                    
# Added utility (rflosurf) to write out surface grids for GenX.                  
#
# Revision 1.12  2002/10/16 21:20:30  haselbac                                   
# Added rflump                                                                   
#
# Revision 1.11  2002/09/17 13:43:00  jferry                                     
# Added Time-dependent boundary conditions                                       
#
# Revision 1.10  2002/09/10 00:00:52  f-najjar                                   
# Reactivate RocVAV compliation flag                                             
#
# Revision 1.9  2002/09/05 17:40:22  jblazek                                     
# Variable global moved into regions().                                          
#
# Revision 1.8  2002/07/12 21:50:07  jblazek                                     
# Added tool to split single grid into multiple regions.                         
#
# Revision 1.7  2002/06/22 01:13:38  jblazek                                     
# Modified interfaces to BC routines.                                            
#
# Revision 1.6  2002/06/10 22:30:13  haselbac                                    
# Fixed spacing problem, moved vav inside ifdef RFLO                             
#
# Revision 1.5  2002/06/10 18:03:55  jblazek                                     
# Added V&V to utilities/Makefile. Corrected path to libraries in vav/Makefile.  
#
# Revision 1.4  2002/01/12 00:02:49  jblazek                                     
# Added postprocessor.                                                           
#
# Revision 1.3  2001/12/22 00:09:39  jblazek                                     
# Added routines to store grid and solution.                                     
#
# Revision 1.2  2001/12/21 23:56:51  jblazek                                     
# Added utility to convert 2D grids to 3D.                                       
#
# Revision 1.1.1.1  2001/12/03 21:44:05  jblazek                                 
# Import of RocfluidMP                                                           
#
################################################################################
