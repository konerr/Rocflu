################################################################################
#
# $Id: Makefile,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
#
# Purpose: Makefile for ROCFLU-MP's converter.
#
# Description: none
#
# Notes: none
#
# Copyright: (c) 2003 by the University of Illinois
#
################################################################################

FLUIDHOME=../..
SRCDIR = $(FLUIDHOME)/utilities/conv

EXEC= $(FLUIDHOME)/rfluconv

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
#   $Log: Makefile,v $
#   Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
#   merged rocflu micro and macro
#
#   Revision 1.1.1.1  2014/07/15 14:31:37  brollin
#   New Stable version
#
#   Revision 1.1  2007/04/09 18:46:33  haselbac
#   Initial revision after split from RocfloMP
#
#   Revision 1.1  2007/04/09 17:57:37  haselbac
#   Initial revision after split from RocfloMP
#
#   Revision 1.5  2005/05/03 03:04:03  haselbac
#   Changed to linking with C++
#
#   Revision 1.4  2003/12/04 03:19:27  haselbac
#   Removed LIBTEC and modified clean target
#
#   Revision 1.3  2003/03/31 19:51:22  haselbac
#   Removed Charm version of executable
#
#   Revision 1.2  2003/03/21 23:07:39  haselbac
#   Added LIBINT to deps, made CP robust, added _charm to EXEC
#
#   Revision 1.1.1.1  2003/03/20 22:51:04  haselbac
#   Initial revision
#
################################################################################
