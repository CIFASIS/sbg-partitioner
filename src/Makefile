# Target variables
MODE ?= Debug

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
SBG_LIB         := $(SBG_DEV).tar.xz
BOOST_LIB_PATH  := $(3RD_PARTY_DIR)/boost
BOOST_LIB       := boost-1.81.0.tar.xz
RM = rm -rf

# The Target Binary Program
TARGET      := $(BIN_DIR)/sbg-partitioner

# Flags, Libraries and Includes
INCLUDES := -I. -I$(SBG_LIB)/$(SBG_DEV)/usr/include -I$(BOOST_LIB_PATH)/include 
CXXFLAGS := -std=c++17 -Wall -Werror -Wno-reorder -O3 
ifeq ($(MODE),Debug)
CXXFLAGS 	+= -ggdb  
endif
LIB_SBG_PARTITIONER = lib/libsbg-partitioner.a
LIBS     := -L./lib -L$(SBG_LIB)/$(SBG_DEV)/usr/lib -lsbgraph

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

# Temp use precompiled lib for sbg.
lib-sbg:  
ifeq ("$(wildcard $(SBG_LIB_PATH)/$(SBG_DEV)/usr/lib)","")
	cd $(SBG_LIB_PATH); tar -xvf $(SBG_LIB)
endif

sbg-partitioner: lib-boost lib-sbg | create-folders
	$(CXX) $(INCLUDE) $(CXXFLAGS) main.cpp -o $(TARGET) $(LIBS) 

create-folders::
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(BUILD_DIR)

test: lib-gtest
	@cd test && $(MAKE)

clean:
	$(RM) $(BUILD_DIR)
	@cd test && $(MAKE) clean

help:
	@echo "make MODE=<Debug|Release> "
	@echo "Default values:"
	@echo ""
	@echo "MODE=Debug"
