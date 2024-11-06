# Directories
SBG_PARTITIONER := ../..
SRC_DIR         := .
3RD_PARTY_DIR   := $(SRC_DIR)/3rd-party
USR_DIR         := $(SRC_DIR)/usr
BIN_DIR         := $(USR_DIR)/bin
LIB_DIR         := $(USR_DIR)/lib
INC_DIR         := $(USR_DIR)/include
BUILD_DIR       := $(USR_DIR)/obj

# Libraries
SBG_PARTITIONER_LIB := $(LIB_DIR)/

# Source includes
SBG_PARTITIONER_INC := $(SBG_PARTITIONER)/src
SB_GRAPH_INC        := $(SBG_PARTITIONER)/src/3rd-party/sbg/sb-graph-dev/usr/include
BOOST_INC           := $(SBG_PARTITIONER)/src/3rd-party/boost/include

# Link libraries
LINK_LIBS := -lmetis -lscotch -lscotcherr -lpatoh -lkahip -lsbg-partitioner -lsbgraph

# Default targets
default: graph_partitioner

# Includes
include 3rd-party/Makefile.include

create-folders:
	@mkdir -p $(BIN_DIR) $(BUILD_DIR) $(LIB_DIR) $(INC_DIR)

# Compiler
CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall -I$(INC_DIR) -I$(SBG_PARTITIONER_INC) -I$(SB_GRAPH_INC) -I$(BOOST_INC)

# Source files
SRC = graph_partitioner.cpp

# Object files
OBJ=$(addprefix $(BUILD_DIR)/, $(SRC:.cpp=.o))

sbg_partitioner_lib:
	cd $(SBG_PARTITIONER) && make
	cd $(SBG_PARTITIONER) && cp ./lib/* src/external_tools/$(LIB_DIR)/
	cd $(SBG_PARTITIONER) && cp ./src/3rd-party/sbg/sb-graph-dev/usr/lib/* src/external_tools/$(LIB_DIR)/

graph_partitioner: kahip patoh sbg_partitioner_lib create-folders $(OBJ) 
	$(CXX) main.cpp -o $(BIN_DIR)/$@ $(OBJ) -L$(LIB_DIR) -L$(LIB_KAHIP_DIR) -L$(LIB_PATOH_DIR) $(LINK_LIBS) 

$(BUILD_DIR)/%.o: %.cpp create-folders
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean::
	rm -rf $(TARGET) $(USR_DIR)

# Phony targets
.PHONY: all clean