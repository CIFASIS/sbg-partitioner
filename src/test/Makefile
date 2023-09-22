# The Directories, Source, Includes, Objects, Binary 
ROOT_DIR       := ../
SRC_DIR        := .
INT_DIR	   	   := $(SRC_DIR)/integration
SYS_DIR   	   := $(SRC_DIR)/system
DATA_DIR   	   := $(SYS_DIR)/test_data
USR_DIR        := $(SRC_DIR)/usr
BUILD_DIR  	   := $(USR_DIR)/obj
BIN_DIR        := $(USR_DIR)/bin

# Flags, Libraries and Includes
GOOGLE_TEST_INSTALL = $(ROOT_DIR)/3rd-party/gtest
GOOGLE_TEST_LIB = gtest
GOOGLE_MOCK_LIB = gmock
GOOGLE_TEST_INCLUDE = $(GOOGLE_TEST_INSTALL)/usr/include
G++ = g++
G++_FLAGS = -c -Wall -I $(GOOGLE_TEST_INCLUDE) -std=c++17
LD_FLAGS = -L $(GOOGLE_TEST_INSTALL)/usr/lib -l $(GOOGLE_TEST_LIB) -l $(GOOGLE_MOCK_LIB) -l pthread 
RM = rm -rf

# The Target Binary Program
INT_TEST    := $(BIN_DIR)/int-test-sbg-partitioner
SYS_TEST	:= $(BIN_DIR)/sys-test-sbg-partitioner
RUN_TESTS   := int-test-sbg-partitioner

# Source files.
MAIN_SRC = $(SRC_DIR)/main.cpp

INT_SRC  =	$(SRC_DIR)/dummy_test.cpp

SYS_SRC = $(SRC_DIR)/sbg_part_test.cpp

# Objects
INT_OBJ=$(addprefix $(BUILD_DIR)/int_, $(notdir $(INT_SRC:.cpp=.o)))
SYS_OBJ=$(addprefix $(BUILD_DIR)/sys_, $(notdir $(SYS_SRC:.cpp=.o)))
MAIN_OBJ=$(addprefix $(BUILD_DIR)/, $(notdir $(MAIN_SRC:.cpp=.o)))

$(BUILD_DIR)/int_%.o : $(INT_DIR)/%.cpp
	$(G++) $(G++_FLAGS) $< -o $@

$(BUILD_DIR)/sys_%.o : $(SYS_DIR)/%.cpp
	$(G++) $(G++_FLAGS) $< -o $@

$(BUILD_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(G++) $(G++_FLAGS) $< -o $@ 

all: $(RUN_TESTS)

$(INT_TEST): $(MAIN_OBJ) $(INT_OBJ) $(BUILD_DIR)
		g++ -o $(INT_TEST) $(MAIN_OBJ) $(INT_OBJ) $(LD_FLAGS) $(LIB)

$(SYS_TEST): $(MAIN_OBJ) $(SYS_OBJ) $(BUILD_DIR)
		g++ -o $(SYS_TEST) $(MAIN_OBJ) $(SYS_OBJ) $(LD_FLAGS) $(LIB)

$(RUN_TESTS): $(INT_TEST) $(SYS_TEST) $(BUILD_DIR)
		@echo Setup test data dir.
		@rm -rf $(DATA_DIR)/*
		@echo Done
		@echo Running integration tests.
		$(INT_TEST)
		@echo Done.
		@echo Running system tests.
		$(SYS_TEST)
		@echo Done.

$(INT_OBJ): | $(BUILD_DIR)

$(SYS_OBJ): | $(BUILD_DIR)

$(MAIN_OBJ): | $(BUILD_DIR)

$(BUILD_DIR):
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(DATA_DIR)

clean:
	$(RM) $(USR_DIR) 
                    
.PHONY: all clean