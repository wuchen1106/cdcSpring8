CXX = g++
BINARIESDIR = BinaryFiles
BINDIR = $(BINARIESDIR)/bin
LIBDIR = $(BINARIESDIR)/lib
APPDIR = app
SRCDIR = src
SUFFIX = .cxx
HEADER = .h
EXCUTE = 

CXXFLAGS += -g -W -Wall -O2
CXXFLAGS += $(shell root-config --cflags)
LIBS     += $(shell $(ROOTSYS)/bin/root-config --glibs) -lMinuit
LIBS     += -pthread -lm -ldl -rdynamic -lGeom -lEG
INCS     += -I$(shell $(ROOTSYS)/bin/root-config --incdir)
INCS     += -Isrc -Iapp

SRCA := $(wildcard $(APPDIR)/*$(SUFFIX))
SRCS := $(wildcard $(SRCDIR)/*$(SUFFIX))
HEADERS := $(wildcard $(SRCDIR)/*$(HEADER))

TGTS = $(addprefix $(BINDIR)/, $(notdir $(basename $(SRCA))))
OBJS = $(addprefix $(LIBDIR)/, $(notdir $(SRCS:$(SUFFIX)=.o)))

.PHONY: all
all: $(TGTS) $(OBJS)

$(BINDIR)/%: $(OBJS) $(APPDIR)/%$(SUFFIX) $(HEADERS)
	mkdir -p $(BINDIR); \
	$(CXX) $(CXXFLAGS) $(LIBS) $(INCS) $(APPDIR)/$(notdir $@)$(SUFFIX) -o $@${EXCUTE} $(filter-out $(APPDIR)/%$(SUFFIX), $^)

$(LIBDIR)/%.o: $(SRCDIR)/%$(SUFFIX) $(HEADERS)
	mkdir -p $(LIBDIR); \
	$(CXX) $(INCS) -c $(CXXFLAGS) $< -o $@

.PHONY: clean
clean:	
	rm -rf $(BINARIESDIR)
