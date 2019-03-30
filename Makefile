CXX = g++
BINARIESDIR = BinaryFiles
BINDIR = $(BINARIESDIR)/bin
LIBDIR = $(BINARIESDIR)/lib
APPDIR = app
SRCDIR = src
SUFFIX = .cxx
HEADER = .hxx
EXCUTE = 

LDFLAGS  += -shared
CXXFLAGS += -g -W -Wall -O2 -fPIC
CXXFLAGS += $(shell root-config --cflags)
LIBS     += $(shell $(ROOTSYS)/bin/root-config --glibs) -lMinuit -L $(LIBDIR) -lTarget
LIBS     += -pthread -lm -ldl -rdynamic -lGeom -lEG
INCS     += -I$(shell $(ROOTSYS)/bin/root-config --incdir)
INCS     += -Isrc -Iapp

SRCA := $(wildcard $(APPDIR)/*$(SUFFIX))
SRCS := $(wildcard $(SRCDIR)/*$(SUFFIX))
HEADERS := $(wildcard $(SRCDIR)/*$(HEADER))

TGTS = $(addprefix $(BINDIR)/, $(notdir $(basename $(SRCA))))
SHLS = $(LIBDIR)/libTarget.so
OBJS = $(addprefix $(LIBDIR)/, $(notdir $(SRCS:$(SUFFIX)=.o)))

.PHONY: all
all: $(TGTS) $(SHLS) $(OBJS)

$(BINDIR)/%: $(APPDIR)/%$(SUFFIX) $(HEADERS) $(SHLS)
	mkdir -p $(BINDIR); \
	$(CXX) $(CXXFLAGS) $(LIBS) $(INCS) $(APPDIR)/$(notdir $@)$(SUFFIX) -o $@${EXCUTE} $(filter-out $(APPDIR)/%$(SUFFIX), $^)

$(SHLS): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

$(LIBDIR)/%.o: $(SRCDIR)/%$(SUFFIX) $(HEADERS)
	mkdir -p $(LIBDIR); \
	$(CXX) $(INCS) -c $(CXXFLAGS) $< -o $@

.PHONY: clean
clean:	
	rm -rf $(BINARIESDIR)
