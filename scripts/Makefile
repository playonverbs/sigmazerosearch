all: build build/bdt_prep build/bdt_apply build/numi_beamline build/numi_beamline_weights_apply

.PHONY: clean build

CXXFLAGS := $(shell root-config --cflags) -O3 -Wall -Wextra -pedantic
LDFLAGS := $(shell root-config --libs --ldflags) -lTMVA -lTMVAGui -lboost_system -lboost_program_options -lboost_filesystem

build/bdt_prep.o: bdt_prep.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

build/bdt_prep: build/bdt_prep.o build
	$(CXX) $< $(LDFLAGS) -o $@

build/bdt_apply.o: bdt_apply.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

build/bdt_apply: build/bdt_apply.o build
	$(CXX) $< $(LDFLAGS) -o $@

build/numi_beamline.o: numi_beamline.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

build/numi_beamline: build/numi_beamline.o build
	$(CXX) $< $(LDFLAGS) -o $@

build/numi_beamline_weights_apply.o: numi_beamline_weights_apply.cc
	$(CXX) -c $(CXXFLAGS) $< -o $@

build/numi_beamline_weights_apply: build/numi_beamline_weights_apply.o build
	$(CXX) $< $(LDFLAGS) -o $@

build:
	mkdir -p build

clean:
	rm build/bdt_prep.o build/bdt_prep build/bdt_apply.o build/bdt_apply build/numi_beamline_weights_apply build/numi_beamline_weights_apply.o
