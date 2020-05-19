################################################################################
#
# $Id: Makefile,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
#
# Purpose: Makefile for ROCFLU-MP's extractor.
#
# Description: None
#
# Notes: None
#
# Copyright: (c) 2007 by the University of Illinois
#
################################################################################

FLUIDHOME=../..
SRCDIR = $(FLUIDHOME)/utilities/extr

EXEC= $(FLUIDHOME)/rfluextr

LIBINT=$(FLUIDHOME)/$(BUILDLIB_DIR)/libflu.a

MODDIR=$(FLUIDHOME)/$(BUILDLIB_DIR)
SRCF77=
SRCF90=
SYSTEM=$(shell uname -s)

### Include common and dependent makefiles #####################################
include $(FLUIDHOME)/Makefile.common
include $(FLUIDHOME)/Makefile.$(SYSTEM)
LIBEXT= $(LIBSYS) 

### Define SRCF90 and SRCF77 ###################################################
include $(SRCDIR)/$(FILELIST)

vpath %.F90 $(SRCDIR)
vpath Mod%.o .:$(FLUIDHOME)/$(BUILDLIB_DIR)
vpath RFLU_Mod%.o .:$(FLUIDHOME)/$(BUILDLIB_DIR)
vpath PLAG_Mod%.o .:$(FLUIDHOME)/$(BUILDLIB_DIR)
vpath %.C $(SRCDIR)

include $(FLUIDHOME)/Makefile.dep

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
# Revision 1.1.1.1  2014/07/15 14:31:37  brollin
# New Stable version
#
# Revision 1.1  2007/11/27 13:13:10  haselbac
# Initial revision
#
################################################################################
