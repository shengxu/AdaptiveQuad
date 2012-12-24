TARGETS := adapsimpson

CC := g++
CXX := g++
# You can add -Werr to GCC to force all warnings to turn into errors
CFLAGS := -g -Wall
CXXFLAGS := -g -Wall -std=c++0x
LDFLAGS := -lpthread

HEADERS := \
	parameters.h \
	interp.h \
	adapsimpsint.h \
	xsdata.h \
	util.h \

# Blank line ends list.

# If you add a new file called "filename.c", you should
# add "filename.o \" to this list.
OBJS := \
	xsdata.o \
	adapsimpsint.o \
	main.o \


# Blank line ends list.

OLDMODE := $(shell cat .buildmode 2> /dev/null)
ifeq ($(DEBUG),1)
CFLAGS := -DDEBUG -O0 $(CFLAGS)
CXXFLAGS := -DDEBUG -O0 $(CXXFLAGS)
ifneq ($(OLDMODE),debug)
$(shell echo debug > .buildmode)
endif
else
CFLAGS := -DNDEBUG -O3 $(CFLAGS)
CXXFLAGS := -DNDEBUG -O3 $(CXXFLAGS)
ifneq ($(OLDMODE),nodebug)
$(shell echo nodebug > .buildmode)
endif
endif

# make all targets specified
all: $(TARGETS)

.PHONY: pintool
pintool:
	$(MAKE) -C pintool

adapsimpson: $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)
	
test: xsdata.o test.o
	$(CXX) $^ -o $@ $(LDFLAGS)

# compile objects

# pattern rule for building objects
%.o: %.cxx %.h $(HEADERS) .buildmode Makefile
	$(CXX) $(CXXFLAGS) -c $< -o $@
%.o: %.c %.h $(HEADERS) .buildmode Makefile
	$(CC) $(CFLAGS) -c $< -o $@

# remove targets and .o files as well as output generated by CQ
clean:
	$(RM) $(TARGETS) $(OBJS) test.o test adapsimpson *.std* .buildmode
