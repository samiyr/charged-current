#ifndef CHARGED_CURRENT_DIS_H
#define CHARGED_CURRENT_DIS_H

#define RUN_TESTS false

#include <iostream>
#include <array>

#include <ChargedCurrent/DIS/DIS.cpp>

#include <ChargedCurrent/SIDIS/SIDIS.cpp>

#include <ChargedCurrent/PDF/Interfaces/LHAInterface.cpp>

#include <ChargedCurrent/Common/TRFKinematics.cpp>

#include <ChargedCurrent/Analysis/Analysis.cpp>

int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv) {
	LHAInterface<>::disable_verbosity();

	const unsigned int number_of_threads = 4;

	const std::vector<double> x_bins = {
		0.02, 0.02125, 0.0225, 0.02375, 0.025, 0.02625, 0.0275, 0.02875, 0.03, 0.03125, 0.0325, 0.03375, 0.035, 0.03625, 0.0375, 0.03875, 
		0.04, 0.04125, 0.0425, 0.04375, 0.045, 0.04625, 0.0475,
		0.05, 0.0625, 0.075, 0.0875, 0.1, 0.1125, 0.125, 0.1375, 0.15, 0.1625, 0.175, 0.1875, 0.2, 0.2125, 0.225, 0.2375, 0.25, 
		0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 
		0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9
	};
	const std::vector<double> x_bins_2 = {
		0.02, 0.0225, 0.025, 0.0275, 0.03, 0.0325, 0.035, 0.0375, 
		0.04, 0.0425, 0.045, 0.0475,
		0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 
		0.275, 0.3, 0.325, 0.35, 0.375, 0.4
	};

	const auto bar = []<
		is_scale_dependence RenormalizationScale, is_scale_dependence FactorizationScale, is_scale_dependence FragmentationScale
	>(Analysis<RenormalizationScale, FactorizationScale, FragmentationScale> analysis) {
		analysis.params.process.type = Process::Type::AntiNeutrinoToAntiLepton;
		return analysis;
	};

	Analysis base;
	base.params.charm_mass = 1.3;
	base.params.scale_variation = ScaleVariation::None;
	base.params.integration.cuba.maximum_relative_error = 1e-2;
	base.params.order = PerturbativeOrder::NLO;
	base.params.process.type = Process::Type::NeutrinoToLepton;
	base.params.number_of_threads = number_of_threads;

	Analysis nlo = base;

	Analysis nnlo = base;
	nnlo.params.order = PerturbativeOrder::NNLO;

	Analysis scale = base;
	scale.params.scale_variation = ScaleVariation::All;

	Analysis errors = base;
	errors.params.pdf_error_sets = true;
	errors.params.order = PerturbativeOrder::NNLO;

	// nlo.dis().charm_production(AnalysisSet::NuTeV_old, x_bins, "Data/DIS/CharmProduction/nutev_old_neutrino.csv");
	// nlo.dis().charm_production(AnalysisSet::NuTeV, x_bins, "Data/DIS/CharmProduction/nutev_new_neutrino.csv");
	// nlo.dis().charm_production(AnalysisSet::CCFR, x_bins, "Data/DIS/CharmProduction/ccfr_neutrino.csv");

	// bar(nlo).dis().charm_production(AnalysisSet::NuTeV_old, x_bins, "Data/DIS/CharmProduction/nutev_old_antineutrino.csv");
	// bar(nlo).dis().charm_production(AnalysisSet::NuTeV, x_bins, "Data/DIS/CharmProduction/nutev_new_antineutrino.csv");
	// bar(nlo).dis().charm_production(AnalysisSet::CCFR, x_bins, "Data/DIS/CharmProduction/ccfr_antineutrino.csv");

	// nnlo.sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/nutev_old_neutrino.csv");
	// nnlo.sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/nutev_new_neutrino.csv");
	// nnlo.sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ccfr_neutrino.csv");

	// bar(nnlo).sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/nutev_old_antineutrino.csv");
	// bar(nnlo).sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/nutev_new_antineutrino.csv");
	// bar(nnlo).sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ccfr_antineutrino.csv");

	// nnlo.sidis().muon_pair_production_only_D0(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/D0/nutev_old_neutrino.csv");
	// nnlo.sidis().muon_pair_production_only_D0(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/D0/nutev_new_neutrino.csv");
	// nnlo.sidis().muon_pair_production_only_D0(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/D0/ccfr_neutrino.csv");

	// bar(nnlo).sidis().muon_pair_production_only_D0(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/D0/nutev_old_antineutrino.csv");
	// bar(nnlo).sidis().muon_pair_production_only_D0(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/D0/nutev_new_antineutrino.csv");
	// bar(nnlo).sidis().muon_pair_production_only_D0(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/D0/ccfr_antineutrino.csv");

	// scale.sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ScaleDependence/nutev_old_neutrino.csv");
	// scale.sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ScaleDependence/nutev_new_neutrino.csv");
	// scale.sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ScaleDependence/ccfr_neutrino.csv");

	// bar(scale).sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ScaleDependence/nutev_old_antineutrino.csv");
	// bar(scale).sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ScaleDependence/nutev_new_antineutrino.csv");
	// bar(scale).sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ScaleDependence/ccfr_antineutrino.csv");

	errors.sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins_2, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ErrorSets/nutev_old_neutrino.csv");
	errors.sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins_2, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ErrorSets/nutev_new_neutrino.csv");
	errors.sidis().muon_pair_production(AnalysisSet::CCFR, x_bins_2, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ErrorSets/ccfr_neutrino.csv");

	bar(errors).sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins_2, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ErrorSets/nutev_old_antineutrino.csv");
	bar(errors).sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins_2, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ErrorSets/nutev_new_antineutrino.csv");
	bar(errors).sidis().muon_pair_production(AnalysisSet::CCFR, x_bins_2, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ErrorSets/ccfr_antineutrino.csv");

	nlo.sidis().muon_pair_production_quark_gluon_channels(x_bins, AnalysisConstants::NuTeV::New::Neutrino::y_bins, AnalysisConstants::NuTeV::New::Neutrino::E_bins, 
		{
			"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_quark_to_quark.csv",
			"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_quark_to_gluon.csv",
			"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_gluon_to_quark.csv",
			"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_gluon_to_gluon.csv"
		}
	);

	nlo.sidis().muon_pair_production_flavor_decomposition_quark_to_quark(
		x_bins, AnalysisConstants::NuTeV::New::Neutrino::y_bins, AnalysisConstants::NuTeV::New::Neutrino::E_bins,
		"Data/SIDIS/MuonPairProduction/CharmedHadrons/FlavorDecomposition/nutev_new"
	);
	nlo.sidis().muon_pair_production_flavor_decomposition_gluon_to_quark(
		x_bins, AnalysisConstants::NuTeV::New::Neutrino::y_bins, AnalysisConstants::NuTeV::New::Neutrino::E_bins,
		"Data/SIDIS/MuonPairProduction/CharmedHadrons/FlavorDecomposition/nutev_new"
	);

	nlo.sidis().muon_pair_production_fragmentation_decomposition(
		x_bins, AnalysisConstants::NuTeV::New::Neutrino::y_bins, AnalysisConstants::NuTeV::New::Neutrino::E_bins,
		"Data/SIDIS/MuonPairProduction/CharmedHadrons/FragmentationDecomposition/nutev_new"
	);

	Analysis proton_nlo;
	proton_nlo.params.pdf_set = "CT14nnlo_NF3";
	nlo.params.integration.cuba.maximum_relative_error = 1e-2;
	nlo.params.order = PerturbativeOrder::NLO;
	nlo.params.number_of_threads = number_of_threads;

	proton_nlo.dis().charm_production_mass_scaling_comparison(
		1.3, 1.5, AnalysisSet::NuTeV_old, x_bins, "Data/DIS/CharmProduction/MassScaling/", "nutev_old_13.csv", "nutev_old_15.csv"
	);
	proton_nlo.dis().charm_production_mass_scaling_comparison(
		1.3, 1.5, AnalysisSet::NuTeV, x_bins, "Data/DIS/CharmProduction/MassScaling/", "nutev_new_13.csv", "nutev_new_15.csv"
	);
	proton_nlo.dis().charm_production_mass_scaling_comparison(
		1.3, 1.5, AnalysisSet::CCFR, x_bins, "Data/DIS/CharmProduction/MassScaling/", "ccfr_13.csv", "ccfr_15.csv"
	);

	return 0;
}

#endif