TARGET=Analysis_test
CXXFLAGS=-std=c++14 -O3
CXX=g++   #icpc

#OpenMM_INSTALL_DIR=/scratch/alifarnudi/local/openmm
OpenMM_INSTALL_DIR=/usr/local/openmm

OSNAME := $(shell uname -s)

Boost_LIB_Mac_Dir = /usr/local/lib
Boost_LIB_Lin_Dir = /usr/lib/x86_64-linux-gnu

MacLIBS= -lOpenMM -lboost_filesystem -lgtest -lgtest_main
LinLIBS= -lOpenMM -lboost_filesystem -lboost_system -lgtest -lgtest_main

ifeq ($(OSNAME),Darwin)
  Boost_LIB_Dir = $(Boost_LIB_Mac_Dir)
  LIBS = $(MacLIBS)
  CXX = g++
else
  Boost_LIB_Dir = $(Boost_LIB_Mac_Dir)
  LIBS = $(LinLIBS)
  CXX = clang++
endif


BINDIR=.
SRCDIR=../../../../source
INCDIR=../../../../include
OBJDIR=../../../../objects


INCDIRS=-I$(INCDIR) -I$(OpenMM_INSTALL_DIR)/include 
LIB_DIR=-L$(OpenMM_INSTALL_DIR)/lib -L$(Boost_LIB_Dir)

SRCFILES=$(wildcard $(SRCDIR)/Tests/Membrane/Analysis_tests/Analysis_test.cpp)  $(wildcard $(SRCDIR)/Membrane/*.cpp) $(wildcard $(SRCDIR)/Chromatin/*.cpp) $(wildcard $(SRCDIR)/Actin/*.cpp) $(wildcard $(SRCDIR)/ECM/*.cpp) $(wildcard $(SRCDIR)/Membrane_Actin/*.cpp) $(wildcard $(SRCDIR)/OpenMM/*.cpp) $(wildcard $(SRCDIR)/Genfuncs/*.cpp)  
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
	rm $(TARGET)
