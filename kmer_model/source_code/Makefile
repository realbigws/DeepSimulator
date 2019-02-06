

# make the subdir first
CXX = g++
#FLAG = -O3 -static
FLAG = -O3
SHELL :=/bin/bash

all:
	$(CXX) kmer_poremodel.cpp kmer_index.cpp $(FLAG) -o kmer_poremodel
	mv kmer_poremodel ../

# clean all the temp files
clean:
