#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TVector3.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;

static const std::vector<std::string> branch_names = {
	"reco_primary_vtx_x",
	"reco_primary_vtx_y",
	"reco_primary_vtx_z",
	"pfp_trk_shr_score",
	"trk_llrpid",
	"shr_start_x",
	"shr_start_y",
	"shr_start_z",
	"shr_length",
	"shr_open_angle",
	"pfp_photon_bdt_score"
};

inline double calculate_nusep(TVector3 nu_vertex, TVector3 shower_vertex) {
	return (shower_vertex - nu_vertex).Mag();
}


void bdt_apply(std::string in_file, std::string weights_file) {
	Double_t reco_primary_vtx_x,
			 reco_primary_vtx_y,
			 reco_primary_vtx_z;
	std::vector<double>* pfp_trk_shr_score = 0;
	std::vector<double>* trk_llrpid = 0;
	std::vector<double>* shr_length = 0;
	std::vector<double>* shr_open_angle = 0;
	std::vector<double>* shr_start_x = 0;
	std::vector<double>* shr_start_y = 0;
	std::vector<double>* shr_start_z = 0;

	std::vector<double>* pfp_photon_bdt_score = 0;

	std::unique_ptr<TFile> input(TFile::Open(in_file.c_str(), "UPDATE"));
	if (!input || input->IsZombie()) {
		throw std::runtime_error("ERROR: could not open input file");
	}

	/* std::unique_ptr<TFile> output(TFile::Open(out_file.c_str(), "UPDATE")); */
	/* if (!output || output->IsZombie()) { */
	/* 	throw std::runtime_error("ERROR: could not open output file"); */
	/* } */

	TTree* t_Output = input->Get<TTree>("ana/OutputTree");
	if (!t_Output || t_Output->IsZombie()) {
		throw std::runtime_error("ERROR: could not get ana/OutputTree from input file");
	}

	std::unique_ptr<TMVA::Reader> reader =
		std::make_unique<TMVA::Reader>("V:Color");

	// declare per-pfp variables for the MVA
	Float_t _pfp_trk_shr_score,
			_trk_llrpid,
			/* _trk_three_plane_dedx, */
			_shr_length,
			_shr_open_angle,
			_nusep;
	reader->AddVariable("pfp_trk_shr_score", &_pfp_trk_shr_score);
	reader->AddVariable("trk_llrpid", &_trk_llrpid);
	/* reader->AddVariable("trk_three_plane_dedx", &_trk_three_plane_dedx); */
	reader->AddVariable("shr_length", &_shr_length);
	reader->AddVariable("shr_open_angle", &_shr_open_angle);
	reader->AddVariable("nusep", &_nusep);
	reader->BookMVA("BDT", weights_file.c_str());

	t_Output->SetBranchStatus("*", 1);

	/// NOTE: use a raw pointer here as std::unique_ptr clashes with ROOTs memory mgmt.
	TBranch* bdt_branch = t_Output->Branch("pfp_photon_bdt_score", &pfp_photon_bdt_score);
	/* std::unique_ptr<TBranch> bdt_branch(t_Output->Branch("pfp_photon_bdt_score", &pfp_photon_bdt_score)); */

	t_Output->SetBranchAddress("reco_primary_vtx_x", &reco_primary_vtx_x);
	t_Output->SetBranchAddress("reco_primary_vtx_y", &reco_primary_vtx_y);
	t_Output->SetBranchAddress("reco_primary_vtx_z", &reco_primary_vtx_z);
	t_Output->SetBranchAddress("pfp_trk_shr_score",  &pfp_trk_shr_score);
	t_Output->SetBranchAddress("trk_llrpid",         &trk_llrpid);
	t_Output->SetBranchAddress("shr_length",         &shr_length);
	t_Output->SetBranchAddress("shr_open_angle",     &shr_open_angle);
	t_Output->SetBranchAddress("shr_start_x",        &shr_start_x);
	t_Output->SetBranchAddress("shr_start_y",        &shr_start_y);
	t_Output->SetBranchAddress("shr_start_z",        &shr_start_z);

	Long64_t n_entries = t_Output->GetEntries();

	for (Long64_t i_evt = 0; t_Output->LoadTree(i_evt) >= 0 ; ++i_evt) {
		t_Output->GetEntry(i_evt);

		pfp_photon_bdt_score->clear();
		
		TVector3 nu_vertex(reco_primary_vtx_x,
						   reco_primary_vtx_y,
						   reco_primary_vtx_z);

		// iterate through all PFPs
		for (size_t i_pfp = 0; i_pfp < pfp_trk_shr_score->size(); ++i_pfp) {
			_pfp_trk_shr_score    = pfp_trk_shr_score->at(i_pfp);
			_trk_llrpid           = trk_llrpid->at(i_pfp);
		/* 	/1* _trk_three_plane_dedx = trk_three_plane_dedx->at(i_pfp); *1/ */
			_shr_length           = shr_length->at(i_pfp);
			_shr_open_angle       = shr_open_angle->at(i_pfp);

			TVector3 shr_start(shr_start_x->at(i_pfp),
							   shr_start_y->at(i_pfp),
							   shr_start_z->at(i_pfp));

			_nusep = calculate_nusep(nu_vertex, shr_start);

			// NOTE: here we evaluate the bdt response for all PFPs, the cuts
			// fed into the BDT training should also be applied later.
			Double_t bdt_value = reader->EvaluateMVA("BDT");
			pfp_photon_bdt_score->emplace_back(bdt_value);
		}

		bdt_branch->Fill();
		if (i_evt % 100 == 0)
			std::printf("%lld / %lld (%.1f %%) events processed\n",
					i_evt, n_entries, 100.0 * (static_cast<float>(i_evt) / static_cast<float>(n_entries)));
	}

	// catch write errors for the TTree and the TFile.
	if (t_Output->Write() == 0) {
		std::runtime_error("ERROR: could not write updates to file");
	};

	if (input->Write() == 0) {
		std::runtime_error("ERROR: could not write updates to file");
	};
	input->Close();

	std::printf("BDT Response branch successfully added to %s.\n", in_file.c_str());
}

int main(int argc, char** argv) {
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "output help messaege")
		("input,i", po::value<std::string>(), "input file name")
		("weights,w",
		 po::value<std::string>()->default_value("./dataset/weights/sigma0_photon_bdt_BDT.weights.xml"),
		 "weights file path created by bdt_prep (via TMVA).");

	po::positional_options_description p;
	p.add("input", 1);

	po::variables_map vm;

	/* po::store(po::parse_command_line(argc, argv, desc), vm); */
	po::store(po::command_line_parser(argc, argv).
			  options(desc).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << "bdt_apply [OPTIONS] <input>\n";
		std::cout << desc << "\n";
		return 1;
	}

	if (vm.count("input") && vm.count("weights")) {
		bdt_apply(vm["input"].as<std::string>(), vm["weights"].as<std::string>());
	} else {
		std::cout << "ERROR: provide an input file\n";

		std::cout << "bdt_apply [OPTIONS] <input>\n";
		std::cout << desc << "\n";
		return 1;
	}

	return 0;
}
