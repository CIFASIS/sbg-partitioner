# Target variables
MODE ?= Debug
sbg_branch ?= sb-graph-dev
build_sbg ?= True
repo_checkout=ssh

# Build cores.
CORES = $(shell nproc)

# The Directories, Source, Includes, Objects, Binary 
ROOT   	        := .
SRC             := $(ROOT)/src
3RD_PARTY_DIR   := $(SRC)/3rd-party
SBG_LIB         := $(3RD_PARTY_DIR)/sbg
SBG_DEV         := sb-graph-dev

all: sbg-partitioner 

update-sbg-lib:
	@echo UPDATING SBG LIBRARY.
	@cd src && $(MAKE) update-sbg MODE=$(MODE) sbg_branch=$(sbg_branch) repo_checkout=$(repo_checkout)
	@echo Done

sbg-partitioner:
	@echo BUILDING PROJECT SBG PARTITIONER
	@cd src && $(MAKE) MODE=$(MODE) build_sbg=$(build_sbg)
	@echo Done

test:
	@echo COMPILE AND RUN TESTS
	@cd src && $(MAKE) test
	@echo Done

.PHONY: clean all 

clean: 
	@cd src && $(MAKE) clean 

help:
	@echo "make MODE=<Debug|Release> sbg_branch=<BRANCH_NAME> build_sbg=<True|False> repo_checkout=<ssh|https>"
	@echo "Default values:"
	@echo ""
	@echo "MODE=Debug"
	@echo "sbg_branch=sb-graph-dev"
	@echo "build_sbg=True"
	@echo "repo_checkout=ssh"
