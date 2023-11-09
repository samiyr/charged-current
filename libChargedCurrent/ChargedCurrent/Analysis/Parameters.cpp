#ifndef ANALYSIS_PARAMETERS_H
#define ANALYSIS_PARAMETERS_H

#include <string>
#include <any>
#include <filesystem>

#include "Common/Process.cpp"
#include "Common/Constants.cpp"
#include "Common/ScaleDependence.cpp"
#include "Common/PerturbativeQuantity.cpp"

#include "Decay/DecayParametrization.cpp"

#include "Integration/IntegrationParameters.cpp"

struct AnalysisParameters {
	std::string pdf_set = "EPPS21nlo_CT18Anlo_Fe56";
	bool explicit_isospin = false;
	double Z = 1.0;
	double A = 1.0;

	Process process = Process(Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino);
	Particle decay_lepton = Constants::Particles::Muon;

	double charm_mass = Constants::Charm::Mass;
	double minimum_lepton_momentum = 5.0;

	ScaleVariation scale_variation = ScaleVariation::None;
	bool pdf_error_sets = false;

	PerturbativeOrder order = PerturbativeOrder::NLO;
	bool use_nlp_nlo = false;

	IntegrationParameters integration;

	bool parallelize = true;
	unsigned int number_of_threads = Utility::get_default_thread_count();

	DecayParametrization decay_parametrization = DecayParametrization::fit1();
	std::vector<DecayParametrization> decay_parametrization_set = {};
	bool decay_variation = false;
	bool use_decay_grid = false;
	std::filesystem::path decay_grid_folder = std::filesystem::current_path() / "DecayGrids";

	bool freeze_factorization = false;
	bool freeze_fragmentation = true;

	double Q2_min = 1.0;

	double primary_muon_min_energy = 0.0;
	double hadronic_min_energy = 0.0;
};

#endif