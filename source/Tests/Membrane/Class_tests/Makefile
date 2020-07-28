TARGET=Class_test
CXXFLAGS=-std=c++14 -O3
CXX=icpc

OpenMM_INSTALL_DIR=/usr/local/openmm
BINDIR=.
SRCDIR=../../../../source
INCDIR=../../../../include
OBJDIR=../../../../objects


INCDIRS=-I$(INCDIR) -I$(OpenMM_INSTALL_DIR)/include -I/usr/local/include/
#INCDIRS=-I$(INCDIR) -I/usr/local/include/
LIB_DIR=-L$(OpenMM_INSTALL_DIR)/lib -L/usr/local/lib
#LIB_DIR= -L/usr/local/lib
LIBS= -lOpenMM -lgtest -lgtest_main
#LIBS= -lgtest -lgtest_main

#SRCFILES=$(wildcard $(SRCDIR)/Tests/*.cpp) $(wildcard $(SRCDIR)/Membrane/*.cpp) $(wildcard $(SRCDIR)/Chromatin/*.cpp) $(wildcard $(SRCDIR)/Actin/*.cpp) $(wildcard $(SRCDIR)/ECM/*.cpp) $(wildcard $(SRCDIR)/Membrane_Actin/*.cpp) $(wildcard $(SRCDIR)/Membrane_Chromatin/*.cpp) $(wildcard $(SRCDIR)/Membrane_Particle/*.cpp) $(wildcard $(SRCDIR)/OpenMM/*.cpp) $(wildcard $(SRCDIR)/Genfuncs/*.cpp)
SRCFILES=$(wildcard $(SRCDIR)/Tests/Membrane/Class_tests/*.cpp) $(wildcard $(SRCDIR)/Membrane/*.cpp) $(SRCDIR)/Genfuncs/General_functions.cpp $(SRCDIR)/Genfuncs/resamplekin.cpp $(SRCDIR)/Genfuncs/map_constructors.cpp $(SRCDIR)/Genfuncs/interaction_map_reader.cpp 
OBJFILES=$(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCFILES))
DEPFILES=$(wildcard $(INCDIR)/*.hpp) $(wildcard $(INCDIR)/*.h)


INC=-I$(DEPFILES)

all: $(BINDIR)/$(TARGET)
	@echo Finished!

$(BINDIR)/$(TARGET): $(OBJFILES)
	@$(CXX) $(CXXFLAGS) $(LIB_DIR) $(LIBS) $? -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCDIRS) -c $< -o $@

SUBDIR_ROOTS := ../../../../objects
DIRS := . $(shell find $(SUBDIR_ROOTS) -type d)
GARBAGE_PATTERNS := *.o
GARBAGE := $(foreach DIR,$(DIRS),$(addprefix $(DIR)/,$(GARBAGE_PATTERNS)))

clean:
	rm -rf $(GARBAGE)
