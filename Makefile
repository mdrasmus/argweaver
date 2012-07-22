#
# ArgHmm
# Matt Rasmussen
# Copyright 2012
#
# Makefile
#

# install prefix paths
prefix = /usr


# C++ compiler options
CXX = g++

CFLAGS := $(CFLAGS) \
    -Wall -fPIC \
    -Isrc


#=============================================================================
# optional CFLAGS

# profiling
ifdef PROFILE
	CFLAGS := $(CFLAGS) -pg
endif

# debugging
ifdef DEBUG	
	CFLAGS := $(CFLAGS) -g
else
	CFLAGS := $(CFLAGS) -O3 -funroll-loops
endif


#=============================================================================
# ArgHmm files

# package
PKG_VERSION:=$(shell python -c 'import arghmm; print arghmm.PROGRAM_VERSION_TEXT')
PKG_NAME=arghmm
PKG=dist/$(PKG_NAME)-$(PKG_VERSION).tar.gz
PKG_DIR=dist/$(PKG_NAME)-$(PKG_VERSION)

# program files
SCRIPTS =  
PROGS = bin/argsample
BINARIES = $(PROGS) $(SCRIPTS)

ARGHMM_SRC = \
    src/emit.cpp \
    src/est_popsize.cpp \
    src/hmm.cpp \
    src/itree.cpp \
    src/local_tree.cpp \
    src/logging.cpp \
    src/newick.cpp \
    src/parsing.cpp \
    src/ptree.cpp \
    src/recomb.cpp \
    src/sample_arg.cpp \
    src/sample_thread.cpp \
    src/seq.cpp \
    src/states.cpp \
    src/sequences.cpp \
    src/thread.cpp \
    src/total_prob.cpp \
    src/trans.cpp \
    src/Tree.cpp \
    src/t2exp.cpp

ALL_SRC = \
    $(ARGHMM_SRC) \
    src/argsample.cpp


ARGHMM_OBJS = $(ARGHMM_SRC:.cpp=.o)
ALL_OBJS = $(ALL_SRC:.cpp=.o) 

LIBS =
# `gsl-config --libs`
#-lgsl -lgslcblas -lm


#=======================
# ArgHmm C-library files
LIBARGHMM = lib/libarghmm.a
LIBARGHMM_SHARED_NAME = libarghmm.so
LIBARGHMM_SHARED = lib/$(LIBARGHMM_SHARED_NAME)
LIBARGHMM_SHARED_INSTALL = $(prefix)/lib/$(LIBARGHMM_SHARED_NAME)
LIBARGHMM_OBJS = $(ARGHMM_OBJS)


#=============================================================================
# targets

# default targets
all: $(PROGS) $(LIBARGHMM) $(LIBARGHMM_SHARED)

bin/argsample: src/argsample.o $(LIBARGHMM)
	$(CXX) $(CFLAGS) -o bin/argsample src/argsample.o $(LIBARGHMM)


#-----------------------------
# ARGHMM C-library
lib: $(LIBARGHMM) $(LIBARGHMM_SHARED)

$(LIBARGHMM): $(LIBARGHMM_OBJS)
	mkdir -p lib
	$(AR) -r $(LIBARGHMM) $(LIBARGHMM_OBJS)

$(LIBARGHMM_SHARED): $(LIBARGHMM_OBJS) 
	mkdir -p lib
	$(CXX) -o $(LIBARGHMM_SHARED) -shared $(LIBARGHMM_OBJS) $(LIBS)


#-----------------------------
# packaging

pkg:
	python make-pkg.py $(PKG_DIR)

$(PKG):
	python make-pkg.py $(PKG_DIR)

#-----------------------------
# install

install: $(BINARIES) $(LIBARGHMM_SHARED_INSTALL)
	mkdir -p $(prefix)/bin
	cp $(BINARIES) $(prefix)/bin
	echo $(LIBARGHMM_SHARED_INSTALL)
	python setup.py install --prefix=$(prefix)

pylib: $(LIBARGHMM_SHARED_INSTALL)
	python setup.py install --prefix=$(prefix)


$(LIBARGHMM_SHARED_INSTALL): $(LIBARGHMM_SHARED)
	mkdir -p $(prefix)/lib
	cp $(LIBARGHMM_SHARED) $(LIBARGHMM_SHARED_INSTALL)

#=============================================================================
# basic rules

$(ALL_OBJS): %.o: %.cpp
	$(CXX) -c $(CFLAGS) -o $@ $<


clean:
	rm -f $(ARGHMM_OBJS) $(LIBARGHMM) $(LIBARGHMM_SHARED)

clean-obj:
	rm -f $(ARGHMM_OBJS)


#=============================================================================
# dependencies

dep:
	touch Makefile.dep
	makedepend -f Makefile.dep -Y src/*.cpp src/*.h

Makefile.dep:
	touch Makefile.dep
	makedepend -f Makefile.dep -Y src/*.cpp src/*.h

include Makefile.dep
# DO NOT DELETE
