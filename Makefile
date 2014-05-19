#
# ARGweaver
# Matt Rasmussen
# Copyright 2012-2013
#
# Makefile
#

# install prefix paths
prefix = /usr


# programs
CXX = g++
PYTHON = python

# C++ compiler options
CFLAGS := $(CFLAGS) \
    -Wall -fPIC \
    -Isrc

GTEST_URL = 'http://googletest.googlecode.com/files/gtest-1.7.0.zip'
GTEST_DIR = gtest-1.7.0

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
# ARGweaver files

# package
PKG_VERSION:=$(shell $(PYTHON) -c 'import argweaver; print argweaver.PROGRAM_VERSION_TEXT' 2>/dev/null)
PKG_NAME=argweaver
DIST=dist
PKG=$(DIST)/$(PKG_NAME)-$(PKG_VERSION).tar.gz
PKG_DIR=$(DIST)/$(PKG_NAME)-$(PKG_VERSION)

# program files
SCRIPTS = bin/*
PROGS = bin/arg-sample bin/arg-summarize bin/smc2bed
BINARIES = $(PROGS) $(SCRIPTS)

ARGWEAVER_SRC = \
    src/compress.cpp \
    src/emit.cpp \
    src/est_popsize.cpp \
    src/fs.cpp \
    src/hmm.cpp \
    src/IntervalIterator.cpp \
    src/itree.cpp \
    src/local_tree.cpp \
    src/logging.cpp \
    src/matrices.cpp \
    src/mem.cpp \
    src/model.cpp \
    src/newick.cpp \
    src/parsing.cpp \
    src/ptree.cpp \
    src/recomb.cpp \
    src/sample_arg.cpp \
    src/sample_thread.cpp \
    src/seq.cpp \
    src/states.cpp \
    src/sequences.cpp \
    src/tabix.cpp \
    src/thread.cpp \
    src/total_prob.cpp \
    src/track.cpp \
    src/trans.cpp \
    src/Tree.cpp \
    src/t2exp.cpp

#    src/prior.cpp
#    src/expm/matrix_exponential.cpp
#    src/expm/r8lib.cpp


ALL_SRC = \
    $(ARGWEAVER_SRC) \
    src/arg-sample.cpp \
    src/arg-summarize.cpp \
    src/smc2bed.cpp


ARGWEAVER_OBJS = $(ARGWEAVER_SRC:.cpp=.o)
ALL_OBJS = $(ALL_SRC:.cpp=.o)

LIBS =
# `gsl-config --libs`
#-lgsl -lgslcblas -lm


# ARGweaver C-library files
LIBARGWEAVER = lib/libargweaver.a
LIBARGWEAVER_SHARED_NAME = libargweaver.so
LIBARGWEAVER_SHARED = lib/$(LIBARGWEAVER_SHARED_NAME)
LIBARGWEAVER_SHARED_INSTALL = $(prefix)/lib/$(LIBARGWEAVER_SHARED_NAME)
LIBARGWEAVER_OBJS = $(ARGWEAVER_OBJS)


# ARGweaver C++ tests
CFLAGS_TEST = -I $(GTEST_DIR)/include
LIBS_TEST = -Llib -lgtest -lgtest_main -lpthread
GTEST_SRC = gtest-1.7.0
TEST_SRC = \
	src/tests/test.cpp \
	src/tests/test_local_tree.cpp \
	src/tests/test_prob.cpp

TEST_OBJS = $(TEST_SRC:.cpp=.o)


#=============================================================================
# targets

.PHONY: all pkg test ctest cq install clean cleanobj lib pylib gtest

# default targets
all: $(PROGS) $(LIBARGWEAVER) $(LIBARGWEAVER_SHARED)

bin/arg-sample: src/arg-sample.o $(LIBARGWEAVER)
	$(CXX) $(CFLAGS) -o bin/arg-sample src/arg-sample.o $(LIBARGWEAVER)

bin/smc2bed: src/smc2bed.o $(LIBARGWEAVER)
	$(CXX) $(CFLAGS) -o bin/smc2bed src/smc2bed.o $(LIBARGWEAVER)


bin/arg-summarize: src/arg-summarize.o $(LIBARGWEAVER)
	$(CXX) $(CFLAGS) -o bin/arg-summarize src/arg-summarize.o $(LIBARGWEAVER)


#-----------------------------
# ARGWEAVER C-library
lib: $(LIBARGWEAVER) $(LIBARGWEAVER_SHARED)

$(LIBARGWEAVER): $(LIBARGWEAVER_OBJS)
	mkdir -p lib
	$(AR) -r $(LIBARGWEAVER) $(LIBARGWEAVER_OBJS)

$(LIBARGWEAVER_SHARED): $(LIBARGWEAVER_OBJS)
	mkdir -p lib
	$(CXX) -o $(LIBARGWEAVER_SHARED) -shared $(LIBARGWEAVER_OBJS) $(LIBS)


#-----------------------------
# packaging

pkg: $(PKG)

$(PKG):
	mkdir -p $(DIST)
	git archive --format=tar --prefix=$(PKG_NAME)-$(PKG_VERSION)/ HEAD | \
	gzip > $(PKG)

#-----------------------------
# testing

test: $(LIBARGWEAVER_SHARED)
	nosetests -v test

cq:
	nosetests -v test/test_codequality.py

ctest: src/tests/test
	src/tests/test

src/tests/test: $(TEST_OBJS) $(LIBARGWEAVER)
	$(CXX) -o src/tests/test $(TEST_OBJS) $(LIBS_TEST) $(LIBARGWEAVER)

$(TEST_OBJS): %.o: %.cpp
	$(CXX) -c $(CFLAGS) $(CFLAGS_TEST) -o $@ $<

# Download and install gtest unit-testing framework.
gtest:
	wget $(GTEST_URL) -O gtest.zip
	rm -rf $(GTEST_SRC)
	unzip gtest
	cd $(GTEST_SRC) && ./configure && make
	cp $(GTEST_SRC)/lib/.libs/libgtest_main.a lib/
	cp $(GTEST_SRC)/lib/.libs/libgtest.a lib/

#-----------------------------
# install

install: $(BINARIES) $(LIBARGWEAVER_SHARED_INSTALL)
	mkdir -p $(prefix)/bin
	cp $(BINARIES) $(prefix)/bin
	$(PYTHON) setup.py install --prefix=$(prefix)

pylib: $(LIBARGWEAVER_SHARED_INSTALL)
	$(PYTHON) setup.py install --prefix=$(prefix)


$(LIBARGWEAVER_SHARED_INSTALL): $(LIBARGWEAVER_SHARED)
	mkdir -p $(prefix)/lib
	cp $(LIBARGWEAVER_SHARED) $(LIBARGWEAVER_SHARED_INSTALL)

#=============================================================================
# basic rules

$(ALL_OBJS): %.o: %.cpp
	$(CXX) -c $(CFLAGS) -o $@ $<

clean:
	rm -f $(ALL_OBJS) $(LIBARGWEAVER) $(LIBARGWEAVER_SHARED) $(TEST_OBJS)

clean-test:
	rm -f $(TEST_OBJS)

clean-obj:
	rm -f $(ALL_OBJS) $(TEST_OBJS)
