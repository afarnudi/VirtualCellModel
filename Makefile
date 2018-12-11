TARGET=Memb
CXXFLAGS=-std=c++11 -O2

BINDIR=bin
SRCDIR=source
INCDIR=include
OBJDIR=objects


INCDIRS=-I$(INCDIR)

SRCFILES=$(wildcard $(SRCDIR)/*.cpp) $(wildcard $(SRCDIR)/Membrane/*.cpp) $(wildcard $(SRCDIR)/Chromatin/*.cpp)
OBJFILES=$(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCFILES)) 
DEPFILES=$(wildcard $(INCDIR)/*.hpp) $(wildcard $(INCDIR)/*.h)


INC=-I$(DEPFILES)

$(BINDIR)/$(TARGET): $(OBJFILES)
	$(CXX) $(CXXFLAGS) $? -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp 
	$(CXX) $(CXXFLAGS) $(INCDIRS) -c $< -o $@

.PHONY: clean
#clean:
#	rm -f $(OBJDIR)/$(OBJFILES) $(BINDIR)/$(TARGET)