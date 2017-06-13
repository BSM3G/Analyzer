ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

CXX = g++
CXXFLAGS += -Wall $(ROOTCFLAGS) -I./ -O3

LD = g++
LDFLAGS += -Wall $(ROOTLIBS) -lGenVector -O3

ifdef DEBUG
CXXFLAGS = -O0 -g -pg -Wall $(ROOTCFLAGS) -I./
LDFLAGS = -O0 -g -Wall $(ROOTLIBS) -lGenVector
endif

LIBS=

SRCDIR = src
BTAGDIR = $(SRCDIR)/btagging
OBJDIR = obj
EXE = Analyzer

#------------------------------------------------------------------------------
SOURCES = $(wildcard src/*.cc)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)

BTAGSRC = $(wildcard $(BTAGDIR)/*.cpp)
BTAGOBJ = $(BTAGSRC:$(BTAGDIR)/%.cpp=$(OBJDIR)/%.o)

#------------------------------------------------------------------------------

all: $(EXE)


$(EXE): $(OBJECTS) $(BTAGOBJ)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)


obj/main.o: src/main.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cc $(SRCDIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

%: $(OBJDIR)/%.o
	$(LD) -o $@ $(LDFLAGS) $<  $(LIBS)

$(OBJDIR)/%.o: $(BTAGDIR)/%.cpp $(BTAGDIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	rm $(OBJDIR)/*

