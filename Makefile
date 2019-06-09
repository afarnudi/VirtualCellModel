TARGET=VCM
CXXFLAGS=-std=c++11 -O2
CXX=g++

OpenMM_INSTALL_DIR=/usr/local/openmm
BINDIR=bin
SRCDIR=source
INCDIR=include
OBJDIR=objects


INCDIRS=-I$(INCDIR) -I$(OpenMM_INSTALL_DIR)/include
LIB_DIR=$(OpenMM_INSTALL_DIR)/lib
LIBS= -lOpenMM

SRCFILES=$(wildcard $(SRCDIR)/*.cpp) $(wildcard $(SRCDIR)/Membrane/*.cpp) $(wildcard $(SRCDIR)/Chromatin/*.cpp) $(wildcard $(SRCDIR)/Actin/*.cpp) $(wildcard $(SRCDIR)/ECM/*.cpp) $(wildcard $(SRCDIR)/Membrane_Actin/*.cpp) $(wildcard $(SRCDIR)/Membrane_Chromatin/*.cpp) $(wildcard $(SRCDIR)/Membrane_Particle/*.cpp)
OBJFILES=$(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCFILES)) 
DEPFILES=$(wildcard $(INCDIR)/*.hpp) $(wildcard $(INCDIR)/*.h)


INC=-I$(DEPFILES)

all: $(BINDIR)/$(TARGET)
	Mac: export DYLD_LIBRARY_PATH=/usr/local/openmm/lib
	Lin: export LD_LIBRARY_PATH=/usr/local/openmm/lib 

$(BINDIR)/$(TARGET): $(OBJFILES)
	$(CXX) $(CXXFLAGS) -L$(LIB_DIR) $(LIBS) $? -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp 
	$(CXX) $(CXXFLAGS) $(INCDIRS) -c $< -o $@

SUBDIR_ROOTS := objects 
DIRS := . $(shell find $(SUBDIR_ROOTS) -type d)
GARBAGE_PATTERNS := *.o
GARBAGE := $(foreach DIR,$(DIRS),$(addprefix $(DIR)/,$(GARBAGE_PATTERNS)))

clean: 
	rm -rf $(GARBAGE)
