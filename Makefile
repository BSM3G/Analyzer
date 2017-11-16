ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)


# Paths for CMSSW libraries:
ifndef CMSSW_RELEASE_BASE
# If you're using a patched CMSSW release, some of the libs are still in the base release, so you also have to look there.
CMSSW_RELEASE_BASE_NOPATCH := $(shell echo $(CMSSW_RELEASE_BASE) | sed -e 's/-patch//' -e 's/_patch.//')
CMSSW_BOOST_BASE := $(shell cat $(CMSSW_RELEASE_BASE)/config/toolbox/$(SCRAM_ARCH)/tools/selected/boost.xml | grep 'name="BOOST_BASE"' | sed -e 's/.*default="//' | sed -e 's/"\/>//')

CMSSW_LIB_PATHS := -L$(CMSSW_BASE)/lib/$(SCRAM_ARCH)
CMSSW_LIB_PATHS += -L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH)
CMSSW_LIB_PATHS += -L$(CMSSW_RELEASE_BASE_NOPATCH)/lib/$(SCRAM_ARCH)
CMSSW_LIB_PATHS += -L$(CMSSW_BOOST_BASE)/lib

# For the headers there are symlinks.
CMSSW_INC_PATHS := -isystem$(CMSSW_BASE)/src
CMSSW_INC_PATHS += -isystem$(CMSSW_RELEASE_BASE)/src
CMSSW_INC_PATHS += -isystem$(CMSSW_BOOST_BASE)/include
endif


CXX = g++
CXXFLAGS += -Wall $(ROOTCFLAGS) -I./
CXXSPEED = -O3

LD = g++
LDFLAGS += -Wall $(ROOTLIBS) -lGenVector
LDSPEED = -O3

# Gather all additional flags
ifndef CMSSW_RELEASE_BASE
EXTRA_CFLAGS  := $(CMSSW_INC_PATHS)
EXTRA_LDFLAGS := $(CMSSW_LIB_PATHS) $(CMSSW_LIBS)
endif


ifdef FAST
CXXSPEED= -Ofast
LDSPEED= -Ofast
endif


CXXFLAGS+=$(CXXSPEED)
LDFLAGS+=$(LDSPEED)


##This will make it very slow
ifdef DEBUG
CXXFLAGS = -Og -g -pg -Wall $(ROOTCFLAGS) -I./
LDFLAGS = -Og -g -Wall $(ROOTLIBS) -lGenVector
endif

CXXFLAGS+=$(EXTRA_CFLAGS) -Wno-deprecated
LDFLAGS+=$(EXTRA_LDFLAGS)
LIBS=

SRCDIR = src
BTAGDIR = $(SRCDIR)/btagging
MT2DIR = $(SRCDIR)/mt2
OBJDIR = obj
EXE = Analyzer

#------------------------------------------------------------------------------
SOURCES = $(wildcard src/*.cc)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)

BTAGSRC = $(wildcard $(BTAGDIR)/*.cpp)
BTAGOBJ = $(BTAGSRC:$(BTAGDIR)/%.cpp=$(OBJDIR)/%.o)

MT2SRC = $(wildcard $(MT2DIR)/*.cc)
MT2OBJ = $(MT2SRC:$(MT2DIR)/%.cc=$(OBJDIR)/%.o)

#------------------------------------------------------------------------------

all: $(EXE)


$(EXE): $(OBJECTS) $(BTAGOBJ) $(MT2OBJ)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)


obj/main.o: src/main.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cc $(SRCDIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

%: $(OBJDIR)/%.o
	$(LD) -o $@ $(LDFLAGS) $<  $(LIBS)

$(OBJDIR)/%.o: $(BTAGDIR)/%.cpp $(BTAGDIR)/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
$(OBJDIR)/%.o: $(MT2DIR)/%.cc $(MT2DIR)/%.hh
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	rm $(OBJDIR)/*

