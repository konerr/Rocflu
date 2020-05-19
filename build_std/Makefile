################################################################################
#
# $Id: Makefile,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
#
# Purpose: Makefile for RocfluidMP standalone codes
#
# Copyright: (c) 2003 by the University of Illinois
#
################################################################################

LIBINT=../$(BUILDLIB_DIR)/libflu.a
EXEC=../rflump
FLUIDHOME=../
MODDIR=../$(BUILDLIB_DIR)
SRCF77=
SRCF90=
SYSTEM=$(shell uname -s)

FPPFLAGS+= 

### Include common and dependent makefiles #####################################
include ../Makefile.common
include ../Makefile.$(SYSTEM)
LIBEXT= $(LIBSYS) $(LIBMPI) $(LIBLAP) $(LIBPET) $(LIBX11)

### vpath definitions ##########################################################
vpath Mod%.F90 		../modfloflu
vpath RFLU_Mod%.F90	../modflu
vpath RFLU_%.F90	../libflu
vpath %.F90		../libfloflu:../rocflu:../standalone
vpath %.f		../libflu
vpath %.C		../standalone

vpath %.o .:../$(BUILDLIB_DIR)

### Define SRCF90 and SRCF77 ###################################################
include ../standalone/$(FILELIST)

include ../Makefile.dep

### Rules ######################################################################
OBJF77=$(SRCF77:.f=.o)
OBJF90=$(SRCF90:.F90=.o)
OBJCXX=$(SRCCXX:.C=.o)

### Targets ####################################################################
.PHONY:	all clean

all:    $(EXEC)

$(EXEC): $(OBJF77) $(OBJF90) $(LIBINT) $(OBJCXX)
	$(LINKER) $(OBJF77) $(OBJF90) $(OBJCXX) -o $(EXEC) $(LDFLAGS) $(LIBINT) $(LIBEXT) $(CXXLINKFLAGS)
clean:  
	$(RM) $(EXEC)
	$(RM) $(MODULES)
	$(RM) $(RIFFRAFF)

install:
	-$(CP) $(EXEC) $(BINDIR)/.

################################################################################
#
# RCS Revision history:
#
# $Log: Makefile,v $
# Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
# merged rocflu micro and macro
#
# Revision 1.1.1.1  2014/07/15 14:31:38  brollin
# New Stable version
#
# Revision 1.3  2008/04/06 13:31:14  haselbac
# Removed -I../modflo; fix bad prior checkin
#
# Revision 1.1  2007/04/09 18:45:36  haselbac
# Initial revision after split from RocfloMP
#
# Revision 1.1  2007/04/09 17:56:27  haselbac
# Initial revision after split from RocfloMP
#
# Revision 1.12  2005/07/05 16:15:18  haselbac
# Removed CHARM ifdef, cosmetics
#			        
# Revision 1.11  2005/05/03 03:03:21  haselbac  
# Changed to linking with C++		        
#
# Revision 1.10  2004/12/19 15:42:38  haselbac  
# Expanded LIBEXT list  		        
#
# Revision 1.9  2003/12/04 03:19:18  haselbac   
# Modified clean target 		        
#
# Revision 1.8  2003/11/20 16:40:32  mdbrandy   
# Backing out RocfluidMP changes from 11-17-03  
#
# Revision 1.5  2003/05/15 02:56:59  jblazek    
# Inlined index function.		        
#
# Revision 1.4  2003/04/11 15:15:18  jiao       
# Changed RFLU_CHARM back to CHARM.	        
#
# Revision 1.3  2003/03/21 23:05:34  haselbac   
# Added LIBINT to deps, made CP robust          
#
# Revision 1.2  2003/03/21 17:09:09  haselbac   
# Added missing header and footer	        
#
################################################################################
