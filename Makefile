TARGET=VCM
CXXFLAGS=-std=c++14 -O2
CXX=g++

# CXXFLAGS=-std=c++14 -O3
# CXX=icpc

#OpenMM_INSTALL_DIR=/scratch/alifarnudi/local/openmm
OpenMM_INSTALL_DIR=/usr/local/openmm
BINDIR=bin
SRCDIR=source
INCDIR=include
OBJDIR=objects


INCDIRS=-I$(INCDIR) -I$(OpenMM_INSTALL_DIR)/include
LIB_DIR=-L$(OpenMM_INSTALL_DIR)/lib  
LIBS= -lOpenMM 

SRCFILES=$(wildcard $(SRCDIR)/main.cpp) $(wildcard $(SRCDIR)/Membrane/*.cpp) $(wildcard $(SRCDIR)/Chromatin/*.cpp) $(wildcard $(SRCDIR)/Actin/*.cpp) $(wildcard $(SRCDIR)/ECM/*.cpp) $(wildcard $(SRCDIR)/Membrane_Actin/*.cpp) $(wildcard $(SRCDIR)/OpenMM/*.cpp) $(wildcard $(SRCDIR)/Genfuncs/*.cpp) 
OBJFILES=$(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCFILES)) 
DEPFILES=$(wildcard $(INCDIR)/*.hpp) $(wildcard $(INCDIR)/*.h)


INC=-I$(DEPFILES)

all: $(BINDIR)/$(TARGET)
	@mkdir $(BINDIR)/Results; true
	@mkdir $(BINDIR)/Results/Relaxation; true
	@mkdir $(BINDIR)/Results/Reports; true
	@mkdir $(BINDIR)/Results/Resumes; true
	@mkdir $(BINDIR)/Results/Resumes/OpenMM; true
	@echo Finished!
	@echo 
	@echo Don\'t forget to export OpenMM\'s Dynamic Library before running the executable. 
	@echo Default paths are:
	@echo     Mac: export DYLD_LIBRARY_PATH=/usr/local/openmm/lib
	@echo     Lin: export LD_LIBRARY_PATH=/usr/local/openmm/lib 
	@echo

$(BINDIR)/$(TARGET): $(OBJFILES)
	@$(CXX) $(CXXFLAGS) $(LIB_DIR) $(LIBS) $? -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)	
	$(CXX) $(CXXFLAGS) $(INCDIRS) -c $< -o $@

SUBDIR_ROOTS := objects 
DIRS := . $(shell find $(SUBDIR_ROOTS) -type d)
GARBAGE_PATTERNS := *.o
GARBAGE := $(foreach DIR,$(DIRS),$(addprefix $(DIR)/,$(GARBAGE_PATTERNS)))

clean: 
	rm -rf $(GARBAGE)
