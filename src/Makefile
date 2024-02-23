# Target variables
MODE ?= Debug
DISTRO := $(shell lsb_release -r 2>/dev/null | grep Release | awk '{ print $$2 }')
sbg_branch ?= sb-graph-dev
build_sbg ?= True
repo_checkout ?= ssh

# Compiler and Linker
CXX     := g++

# The Directories, Source, Includes, Objects, Binary 
ROOT   	        := ..
SCR_DIR         := .
BIN_DIR         := $(ROOT)/bin
3RD_PARTY_DIR   := $(SCR_DIR)/3rd-party
BUILD_DIR    	:= $(SCR_DIR)/obj/release
ifeq ($(MODE),Debug)
BUILD_DIR    	:= $(SCR_DIR)/obj/debug
endif
SBG_LIB_PATH    := $(3RD_PARTY_DIR)/sbg
SBG_DEV         := sb-graph-dev
SBG_LIB         := $(SBG_DEV).tar.gz
BOOST_LIB_PATH  := $(3RD_PARTY_DIR)/boost
BOOST_LIB       := boost-1.81.0.tar.xz
RM = rm -rf

# The Target Binary Program
TARGET      := $(BIN_DIR)/sbg-partitioner

# Source files
SOURCES := build_sb_graph.cpp \
		   dfs_on_sbg.cpp \
		   main.cpp \
		   partition_graph.cpp \
		   partition_strategy.cpp

# Flags, Libraries and Includes
INCLUDES := -I. -I$(SBG_LIB_PATH)/$(SBG_DEV)/usr/include -I$(BOOST_LIB_PATH)/include -I$(3RD_PARTY_DIR)/rapidjson/include
CXXFLAGS := -std=c++17 -Wall -Werror -Wno-reorder -O3 
ifeq ($(MODE),Debug)
CXXFLAGS 	+= -ggdb  
endif
LIB_SBG_PARTITIONER = lib/libsbg-partitioner.a
LIBS     := -L./lib -L$(SBG_LIB_PATH)/$(SBG_DEV)/usr/lib -lsbgraph

default: sbg-partitioner

# Sources

$(BUILD_DIR)/%.o : %.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) -MM -MT $@ -MF $(patsubst %.o,%.d,$@) $<
	$(CXX) $(INCLUDES) -c $< -o $@ $(CXXFLAGS)

lib-gtest: | create-folders 
ifeq ("$(wildcard $(3RD_PARTY_DIR)/gtest/usr/lib)","") 
	mkdir -p $(3RD_PARTY_DIR)/gtest/usr
	cp -r $(3RD_PARTY_DIR)/gtest/ubuntu-$(DISTRO)/* $(3RD_PARTY_DIR)/gtest/usr/
endif

lib-boost:
ifeq ("$(wildcard $(BOOST_LIB_PATH)/include)","")
	cd $(BOOST_LIB_PATH); tar -xvf $(BOOST_LIB)
endif

update-sbg:
	cd $(SBG_LIB_PATH); ./update.py --branch_name $(sbg_branch) --repo_checkout $(repo_checkout)
	cd $(SBG_LIB_PATH); tar -zxvf $(SBG_LIB)

lib-sbg: lib-boost  
ifeq ("$(wildcard $(SBG_LIB_PATH)/$(SBG_DEV))","")
	make update-sbg
endif
ifeq ($(build_sbg), True)
	cd $(SBG_LIB_PATH)/$(SBG_DEV); autoconf
	cd $(SBG_LIB_PATH)/$(SBG_DEV); ./configure
	cd $(SBG_LIB_PATH)/$(SBG_DEV); make boost_libdir=../../boost
	cd $(SBG_LIB_PATH)/$(SBG_DEV); mkdir -p usr
	cd $(SBG_LIB_PATH)/$(SBG_DEV); make install prefix=./usr
endif

sbg-partitioner: lib-boost lib-sbg | create-folders
	$(CXX) $(INCLUDES) $(CXXFLAGS) $(SOURCES) -o $(TARGET) $(LIBS)

create-folders::
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(BUILD_DIR)

test: lib-gtest
	@cd test && $(MAKE)

clean:
	$(RM) $(BUILD_DIR)
	@cd test && $(MAKE) clean

help:
	@echo "make MODE=<Debug|Release> sbg_branch=<BRANCH_NAME> build_sbg=<True|False> repo_checkout=<ssh|https>"
	@echo "Default values:"
	@echo ""
	@echo "MODE=Debug"
	@echo "sbg_branch=sb-graph-dev"
	@echo "build_sbg=True"
	@echo "repo_checkout=ssh"
