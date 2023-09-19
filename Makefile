# Target variables
MODE ?= Debug
DISTRO := $(shell lsb_release -r 2>/dev/null | grep Release | awk '{ print $$2 }')

# Build cores.
CORES = $(shell nproc)

# The Directories, Source, Includes, Objects, Binary 
ROOT   	        := .
SRC             := $(ROOT)/src
3RD_PARTY_DIR   := $(SRC)/3rd-party
SBG_LIB         := $(3RD_PARTY_DIR)/sbg
SBG_DEV         := sb-graph-dev

all: sbg-partitioner 

sbg-partitioner-doc: sbg-partitioner
	@echo BUILDING PROJECT DOCUMENTATION
	@cd src && $(MAKE) doc 
	@echo Done

lib-gtest: | create-folders 
ifeq ("$(wildcard $(3RD_PARTY_DIR)/gtest/usr/lib)","")
	mkdir -p $(3RD_PARTY_DIR)/gtest/usr
	cp -r $(3RD_PARTY_DIR)/gtest/ubuntu-$(DISTRO)/* $(3RD_PARTY_DIR)/gtest/usr/
endif

lib-sbg:  
ifeq ("$(wildcard $(SBG_DEV)/lib)","")
	cd $(SBG_LIB); tar xvzf $(SBG_DEV).tar.gz
	cd $(SBG_LIB)/$(SBG_DEV); autoconf 
	cd $(SBG_LIB)/$(SBG_DEV); ./configure
	cd $(SBG_LIB)/$(SBG_DEV); make
	cd $(SBG_LIB)/$(SBG_DEV); make install
endif

test:
	@echo COMPILE AND RUN TESTS
	@cd src && $(MAKE) test
	@echo Done

create-folders:
	@mkdir -p $(ROOT)/bin
	@mkdir -p $(ROOT)/lib

install:	
	@echo "Installing GTEST libraries."
	rm $(ROOT)/lib/*
	ln -rsvf $(ROOT)/src/3rd-party/gtest/usr/lib/*.a $(ROOT)/lib/
	ln -rsvf $(SBG_DEV)/lib/*.a $(ROOT)/lib/
	
.PHONY: clean all 

clean: 
	@cd src && $(MAKE) clean 

help:
	@echo "make MODE=<Debug|Release> "
	@echo "Default values:"
	@echo ""
	@echo "MODE=Debug"
