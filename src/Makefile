# Target variables
MODE ?= Release
DISTRO := $(shell lsb_release -r 2>/dev/null | grep Release | awk '{ print $$2 }')
sbg_branch ?= sbg-partioner-related-branch
build_sbg ?= True
repo_checkout ?= ssh

# Compiler and Linker
CXX     := g++

# The Directories, Source, Includes, Objects, Binary 
ROOT   	        := ..
SCR_DIR         := .
BIN_DIR         := $(ROOT)/bin
LIB_DIR         := $(ROOT)/lib
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
TARGET_METRICS      := $(BIN_DIR)/sbg-partitioner-metrics

# Source files
SOURCES := build_sb_graph.cpp \
		   dfs_on_sbg.cpp \
		   partition_graph.cpp \
		   kernighan_lin_partitioner.cpp \
		   partition_metrics_api.cpp \
		   partition_strategy.cpp \
		   weighted_sb_graph.cpp
OSOURCES := $(SOURCES:.cpp=.o)
MAIN_SRC := main.cpp
MAIN_OBJ := $(MAIN_SRC:.cpp=.o)
METRICS_SRC := partition_metrics.cpp
METRICS_OBJ := $(METRICS_SRC:.cpp=.o)

# Flags, Libraries and Includes
INCLUDES := -I. -I$(SBG_LIB_PATH)/$(SBG_DEV)/usr/include -I$(BOOST_LIB_PATH)/include -I$(3RD_PARTY_DIR)/rapidjson/include
CXXFLAGS := -std=c++17 -Wall -Werror -Wno-reorder -O3 -pthread
ifeq ($(MODE),Debug)
CXXFLAGS 	+= -ggdb  
CXXFLAGS    += -D PARTITION_SANITY_CHECK
CXXFLAGS    += -D SBG_PARTITIONER_LOGGING
endif
LIB_SBG_PARTITIONER = $(LIB_DIR)/libsbg-partitioner.a
LIBS     := -L./lib -L$(SBG_LIB_PATH)/$(SBG_DEV)/usr/lib -lsbgraph -lstdc++fs

default: sbg-partitioner

# Sources

$(SCR_DIR)/%.o : %.cpp
	# $(CXX) $(INCLUDES) $(CXXFLAGS) -MM -MT $@ -MF $(patsubst %.o,%.d,$@) $<
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

sbg-partitioner-lib: build-sbg-partitioner-lib lib-boost lib-sbg | create-folders


build-sbg-partitioner-lib: $(LIB_SBG_PARTITIONER)

$(LIB_SBG_PARTITIONER): $(OSOURCES) | create-folders
	ar rcs $(LIB_SBG_PARTITIONER) $(SCR_DIR)/*.o $(SBG_LIB_PATH)/$(SBG_DEV)/usr/lib/libsbgraph.a

sbg-partitioner-main: $(MAIN_OBJ)
	$(CXX) $(INCLUDES) -c $(MAIN_SRC) -o $(MAIN_OBJ) $(CXXFLAGS)


sbg-partitioner: lib-sbg lib-boost sbg-partitioner-main sbg-partitioner-lib | create-folders
	$(CXX) $(MAIN_OBJ) -L$(LIB_DIR) -lsbg-partitioner -o $(TARGET) $(CXXFLAGS) $(LIBS)

sbg-partitioner-metrics-main: $(METRICS_OBJ)
	$(CXX) $(INCLUDES) -c $(METRICS_SRC) -o $(METRICS_OBJ) $(CXXFLAGS)

sbg-partitioner-metrics: sbg-partitioner-metrics-main sbg-partitioner-lib lib-boost lib-sbg | create-folders
	$(CXX) $(METRICS_OBJ) -L$(LIB_DIR) -lsbg-partitioner -o $(TARGET_METRICS) $(CXXFLAGS) $(LIBS)


create-folders::
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(LIB_DIR)

test: lib-gtest
	@cd test && $(MAKE)

clean:
	$(RM) $(BUILD_DIR)
	$(RM) $(OSOURCES)
	$(RM) $(MAIN_OBJ)
	$(RM) $(METRICS_OBJ)
	@cd test && $(MAKE) clean

help:
	@echo "make MODE=<Debug|Release> sbg_branch=<BRANCH_NAME> build_sbg=<True|False> repo_checkout=<ssh|https>"
	@echo "Default values:"
	@echo ""
	@echo "MODE=Debug"
	@echo "sbg_branch=sbg-partioner-related-branch"
	@echo "build_sbg=True"
	@echo "repo_checkout=ssh"
