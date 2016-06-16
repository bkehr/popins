#-*-makefile-*-

-include popins.config

CXX=g++ -std=c++11
CC=$(CXX)

TOOLS=-DSAMTOOLS=\"$(SAMTOOLS)\" -DBWA=\"$(BWA)\" -DSICKLE=\"$(SICKLE)\" -DVELVETH=\"$(VELVETH)\" -DVELVETG=\"$(VELVETG)\"

GIT_DATE := $(shell git log --pretty=format:"%cd" | head -n 1)
GIT_VERSION := $(shell git describe --always)

WARN= -W -Wall
FLAGS=-I$(SEQAN_LIB) -DSEQAN_HAS_ZLIB=1 -DVERSION=\"$(GIT_VERSION)\" -DVERSION_DATE=\""$(GIT_DATE)"\"

# Flags for optimization OR for debugging
RELEASE_FLAGS=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0
DEBUG_FLAGS=-g -fno-inline -DDEBUG -DBOUNDS_CHECK -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1

# Replace RELEASE_FLAGS by DEBUG_FLAGS for debugging
CXXFLAGS = $(WARN) $(RELEASE_FLAGS) $(INCLUDE)
#CC = g++ -fno-merge-constants

LDLIBS=-lz -lrt -pthread

.cpp.o:popins.cpp
	$(CC) -c $(CXXFLAGS) $(FLAGS) $(TOOLS) $< -o $@

all:popins

popins:popins.o all.dep
	$(CC) -o $@ $(CXXFLAGS) popins.o $(LDLIBS)

all.dep:
	$(CC) -c -MM popins.cpp > all.dep

-include all.dep

depend:
	rm all.dep
	make all.dep

clean:
	rm -f all.dep *.o popins

default:
	all
