#ifndef CHARGED_CURRENT_DIS_H
#define CHARGED_CURRENT_DIS_H

#define RUN_TESTS false

#include <iostream>
#include "DIS/DIS.cpp"
#include "SIDIS/SIDIS.cpp"
#include "PDF/Interfaces/LHAInterface.cpp"
#include "Common/TRFKinematics.cpp"
#include "Analysis/Analysis.cpp"
#include <array>

int enable_testing(int, char **);
#if RUN_TESTS
#include "Tests/Tests.cpp"
#include <gtest/gtest.h>
#endif

int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv) {
	LHAInterface<>::disable_verbosity();

	#if RUN_TESTS
	enable_testing(argc, argv);
	return 0;
	#endif

	const std::vector<double> x_bins = {0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3};

	Analysis analysis;

	analysis.dis().charm_production(AnalysisSet::NuTeV_old, x_bins, "Data/DIS/CharmProduction/nutev_old.csv");
	analysis.dis().charm_production(AnalysisSet::NuTeV, x_bins, "Data/DIS/CharmProduction/nutev_new.csv");
	analysis.dis().charm_production(AnalysisSet::CCFR, x_bins, "Data/DIS/CharmProduction/ccfr.csv");

	analysis.sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/nutev_old.csv");
	analysis.sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/nutev_new.csv");
	analysis.sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ccfr.csv");

	analysis.sidis().muon_pair_production_only_D0(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/D0/nutev_old.csv");
	analysis.sidis().muon_pair_production_only_D0(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/D0/nutev_new.csv");
	analysis.sidis().muon_pair_production_only_D0(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/D0/ccfr.csv");

	analysis.sidis().muon_pair_production_quark_gluon_channels(x_bins, AnalysisConstants::NuTeV::New::y_bins, AnalysisConstants::NuTeV::New::E_bins, 
		{
			"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_quark_to_quark.csv",
			"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_quark_to_gluon.csv",
			"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_gluon_to_quark.csv",
			"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_gluon_to_gluon.csv"
		}
	);

	analysis.sidis().muon_pair_production_flavor_decomposition_quark_to_quark(x_bins, AnalysisConstants::NuTeV::New::y_bins, AnalysisConstants::NuTeV::New::E_bins,
		"Data/SIDIS/MuonPairProduction/CharmedHadrons/FlavorDecomposition/nutev_new"
	);
	analysis.sidis().muon_pair_production_flavor_decomposition_gluon_to_quark(x_bins, AnalysisConstants::NuTeV::New::y_bins, AnalysisConstants::NuTeV::New::E_bins,
		"Data/SIDIS/MuonPairProduction/CharmedHadrons/FlavorDecomposition/nutev_new"
	);

	analysis.sidis().muon_pair_production_fragmentation_decomposition(x_bins, AnalysisConstants::NuTeV::New::y_bins, AnalysisConstants::NuTeV::New::E_bins,
		"Data/SIDIS/MuonPairProduction/CharmedHadrons/FragmentationDecomposition/nutev_new"
	);


	return 0;
}

int enable_testing([[maybe_unused]] int argc, [[maybe_unused]] char **argv) {
	#if RUN_TESTS
	::testing::InitGoogleTest(&argc, argv);
  	return RUN_ALL_TESTS();
	#else
	return 0;
	#endif
}

#endif