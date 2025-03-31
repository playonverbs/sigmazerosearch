all: bdt_prep bdt_apply

.PHONY: clean

CXXFLAGS := $(shell root-config --cflags) -O3 -ggdb -Wall -Wextra -pedantic
LDFLAGS := $(shell root-config --libs --ldflags) -lTMVA -lTMVAGui -lboost_system -lboost_program_options -lboost_filesystem

bdt_prep.o: bdt_prep.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

bdt_prep: bdt_prep.o
	$(CXX) $< $(LDFLAGS) -o $@

bdt_apply.o: bdt_apply.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

bdt_apply: bdt_apply.o
	$(CXX) $< $(LDFLAGS) -o $@

clean:
	rm bdt_prep.o bdt_prep
