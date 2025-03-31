#include <cstdlib>
#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TROOT.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Config.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
/* namespace fs = boost::filesystem; */

void bdt_prep(std::string in_file, std::string out_file) {
	std::unique_ptr<TFile> input{TFile::Open(in_file.c_str(), "READ")};

	if (!input || input->IsZombie()) {
		throw std::runtime_error("ERROR: could not open file");
	}

	std::unique_ptr<TFile> output{TFile::Open(out_file.c_str(), "RECREATE")};
	if (!output || output->IsZombie()) {
		throw std::runtime_error("ERROR: could not open file");
	}

	TTree* t_Signal = input->Get<TTree>("bdt/SignalTree");
	TTree* t_Background = input->Get<TTree>("bdt/BackgroundTree");

	double SignalWeight = 1.0;
	double BackgroundWeight = 1.0;

	TMVA::Factory* factory = new TMVA::Factory("sigma0_photon_bdt", output.get());
	TMVA::DataLoader* loader = new TMVA::DataLoader("dataset");

	TMVA::gConfig().GetVariablePlotting().fPlotFormat = TMVA::Config::VariablePlotting::kPDF;

	loader->AddSignalTree(t_Signal, SignalWeight);
	loader->AddBackgroundTree(t_Background, BackgroundWeight);

	loader->AddVariable("pfp_trk_shr_score", 'F');
	loader->AddVariable("trk_llrpid", 'F');
	/* loader->AddVariable("trk_three_plane_dedx", 'F'); */
	loader->AddVariable("shr_length", 'F');
	loader->AddVariable("shr_open_angle", 'F');
	loader->AddVariable("nusep", "Separation from nu vertex", "cm", 'F');

	// Remove entries with error values
	TCut preselection = "pfp_trk_shr_score < 0.5 && pfp_trk_shr_score >= 0 && trk_llrpid > -50 && shr_length > -5 && shr_open_angle > -5";
	loader->PrepareTrainingAndTestTree(preselection, "");

	factory->BookMethod(loader, TMVA::Types::kBDT, "BDT");
	factory->BookMethod(loader, TMVA::Types::kBDT, "BDTG",
						"V:BoostType=Grad");
	factory->BookMethod(loader, TMVA::Types::kLikelihood, "Likelihood");

	factory->TrainAllMethods();

	factory->TestAllMethods();

	factory->EvaluateAllMethods();

	output->Write();

	if (!gROOT->IsBatch())
		TMVA::TMVAGui(out_file.c_str());
}

int main(int argc, char** argv) {
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "output help message")
		("input,i", po::value<std::string>(), "input file name")
		("output,o", po::value<std::string>(), "output file name")
		(
		 "presel,p",
		 po::value<std::string>()->default_value(
			 "pfp_trk_shr_score > -0.5 && trk_llrpid > -50 && shr_length > -5 && shr_open_angle > -5"
			 ),
		"preselection cut"
		);

	po::positional_options_description p;
	p.add("input", 1);
	p.add("output", 1);

	po::variables_map vm;
	/* po::store(po::parse_command_line(argc, argv, desc), vm); */
	po::store(po::command_line_parser(argc, argv).
			options(desc).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << "bdt_prep [OPTIONS] <input> <output>\n";
		std::cout << desc <<"\n";
		return 1;
	}

	if (vm.count("input") && vm.count("output")) {
		bdt_prep(vm["input"].as<std::string>(), vm["output"].as<std::string>());
	} else {
		std::cout << "ERROR: provide an input and output file\n";

		std::cout << "bdt_prep [OPTIONS] <input> <output>\n";
		std::cout << desc << "\n";
		return 1;
	}

	return 0;
}
