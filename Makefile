###############################################################################
#                                                                             #
#        Generic and Simple Makefile for this C++ file structure              #
#                                                                             #
###############################################################################

#------------ Compilers ------------#
XC := icpc 
FC := ifort
CC := icc

#------------ File Tree ------------#
SRC_DIR := src
OBJ_DIR := build
TARGET := bin/main.x

#------------ Variables ------------#
SRCEXT := cpp
SRC := $(shell find $(SRC_DIR) -type f -name *.$(SRCEXT))
OBJ := $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(SRC:.$(SRCEXT)=.o))

#---------- Compile Flags ----------#
XFLAGS := -std=c++11 -fast
FFLAGS := -std=15 -march=native -O2

#---------- Linker Flags -----------#
XINC := -I include
FINC := -I ${MKLROOT}/include

###############################################################################
#                                                                             #
#                            End of Definitions                               #
#                                                                             #
###############################################################################

$(TARGET): $(OBJ)
	@echo " Linking...";
	@echo " $(XC) $^ -o $(TARGET) "; $(XC) $^ -o $(TARGET) 

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.$(SRCEXT)
	@mkdir -p $(OBJ_DIR)
	@echo " $(XC) $(XFLAGS) $(XINC) -c -o $@ $<"; $(XC) $(XFLAGS) $(XINC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 	
	@echo " $(RM) -r $(OBJ_DIR) $(TARGET)"; $(RM) -r $(OBJ_DIR) $(TARGET)

test:
	$(XC) $(XFLAGS) test/test.cpp $(XINC) -o bin/test

.PHONY: clean
