TARGET=VCM
CXXFLAGS=-std=c++11 -O2
CXX=g++

BINDIR=bin
SRCDIR=source
INCDIR=include
OBJDIR=objects


INCDIRS=-I$(INCDIR) 
INCOPENMM=-I/usr/local/openmm/include

SRCFILES=$(wildcard $(SRCDIR)/*.cpp) $(wildcard $(SRCDIR)/Membrane/*.cpp) $(wildcard $(SRCDIR)/Chromatin/*.cpp) $(wildcard $(SRCDIR)/Actin/*.cpp) $(wildcard $(SRCDIR)/ECM/*.cpp) $(wildcard $(SRCDIR)/Membrane\-Actin/*.cpp) $(SRCDIR)/Membrane\-Chromatin/*.cpp)
OBJFILES=$(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCFILES)) 
DEPFILES=$(wildcard $(INCDIR)/*.hpp) $(wildcard $(INCDIR)/*.h)


INC=-I$(DEPFILES)

$(BINDIR)/$(TARGET): $(OBJFILES)
	$(CXX) $(CXXFLAGS) $? -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp 
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(INCOPENMM) -c $< -o $@

.PHONY: clean
#clean:
#	rm -f $(OBJDIR)/$(OBJFILES) $(BINDIR)/$(TARGET)
