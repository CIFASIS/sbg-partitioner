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

sbg-partitioner:
	@echo BUILDING PROJECT SBG PARTITIONER
	@cd src && $(MAKE) MODE=$(MODE)
	@echo Done

test:
	@echo COMPILE AND RUN TESTS
	@cd src && $(MAKE) test
	@echo Done

.PHONY: clean all 

clean: 
	@cd src && $(MAKE) clean 

help:
	@echo "make MODE=<Debug|Release> "
	@echo "Default values:"
	@echo ""
	@echo "MODE=Debug"
