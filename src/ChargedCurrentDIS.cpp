#ifndef CHARGED_CURRENT_DIS_H
#define CHARGED_CURRENT_DIS_H

#include <iostream>
#include "DIS/DIS.cpp"
#include "SIDIS/SIDIS.cpp"
#include "PDF/Interfaces/LHAInterface.cpp"
#include "Utility/Tests.cpp"
#include "Common/TRFKinematics.cpp"
#include "Analysis/Analysis.cpp"
// #include "PythiaSIDIS.cpp"
#include <array>

int main() {
	static_assert(std::is_nothrow_move_constructible<LHAInterface>::value, "MyType should be noexcept MoveConstructible"); 
	LHAInterface::disable_verbosity();

	const std::vector<double> x_bins = {0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3};

	// Analysis::Inclusive::charm_production_nutev_old(x_bins, "Data/DIS/CharmProduction/nutev_old.csv");
	// Analysis::Inclusive::charm_production_nutev_new(x_bins, "Data/DIS/CharmProduction/nutev_new.csv");
	// Analysis::Inclusive::charm_production_ccfr(x_bins, "Data/DIS/CharmProduction/ccfr.csv");

	Analysis::SemiInclusive::muon_pair_production_nutev_old(x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/nutev_old.csv");
	// Analysis::SemiInclusive::muon_pair_production_nutev_new(x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/nutev_new.csv");
	// Analysis::SemiInclusive::muon_pair_production_ccfr(x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ccfr.csv");

	// Analysis::SemiInclusive::muon_pair_production_nutev_old_only_D0(x_bins, "Data/SIDIS/MuonPairProduction/D0/nutev_old.csv");
	// Analysis::SemiInclusive::muon_pair_production_nutev_new_only_D0(x_bins, "Data/SIDIS/MuonPairProduction/D0/nutev_new.csv");
	// Analysis::SemiInclusive::muon_pair_production_ccfr_only_D0(x_bins, "Data/SIDIS/MuonPairProduction/D0/ccfr.csv");

	// Analysis::SemiInclusive::ChannelDecomposition::muon_pair_production_quark_gluon_channels(x_bins, Analysis::NuTeV::New::y_bins, Analysis::NuTeV::New::E_bins,
	// 	{
	// 		"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_quark_to_quark.csv",
	// 		"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_quark_to_gluon.csv",
	// 		"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_gluon_to_quark.csv",
	// 		"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_gluon_to_gluon.csv"
	// 	}
	// );

	// Analysis::SemiInclusive::FlavorDecomposition::muon_pair_production_quark_to_quark(x_bins, Analysis::NuTeV::New::y_bins, Analysis::NuTeV::New::E_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/FlavorDecomposition/nutev_new");
	// Analysis::SemiInclusive::FlavorDecomposition::muon_pair_production_gluon_to_quark(x_bins, Analysis::NuTeV::New::y_bins, Analysis::NuTeV::New::E_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/FlavorDecomposition/nutev_new");

	// Analysis::SemiInclusive::FragmentationDecomposition::muon_pair_production(x_bins, Analysis::NuTeV::New::y_bins, Analysis::NuTeV::New::E_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/FragmentationDecomposition/nutev_new");

	// Tests::dis_cross_section_tests();
	// return 0;

	// const double minimum_lepton_momentum = 5.0;
	// const Particle target = Constants::Particles::Proton;
	// const auto decay_function = DecayFunctions::decay_function;

	// const DecayParametrization parametrization = DecayParametrization::fit1();

	// SIDIS sidis(
	// 	{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
	// 	LHAInterface("EPPS21nlo_CT18Anlo_Fe56"),
	// 	FragmentationConfiguration(
	// 		{
	// 			LHAInterface("kkks08_opal_d0___mas"), 
	// 			LHAInterface("kkks08_opal_d+___mas"), 
	// 			LHAInterface("bkk05_D3_d_s_nlo"), 
	// 			LHAInterface("bkk05_D3_lambda_c_nlo")
	// 		},
	// 		{
	// 			Decay(parametrization, Constants::Particles::D0, target, decay_function, minimum_lepton_momentum),
	// 			Decay(parametrization, Constants::Particles::Dp, target, decay_function, minimum_lepton_momentum),
	// 			Decay(parametrization, Constants::Particles::Ds, target, decay_function, minimum_lepton_momentum),
	// 			Decay(parametrization, Constants::Particles::LambdaC, target, decay_function, minimum_lepton_momentum)
	// 		}
	// 	),
	// 	500'000,
	// 	Process {Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino }
	// 	// ScaleDependence::multiplicative(2.0),
	// 	// ScaleDependence::multiplicative(2.0)
	// );

	// sidis.charm_mass = 1.3;

	// sidis.max_chi_squared_deviation = 0.5;
	// sidis.max_relative_error = 1e-2;
	// sidis.iter_max = 5;
	// sidis.combine_integrals = true;

	// const std::vector<double> y_bins_nutev = {0.334, 0.573, 0.790};
	// const std::vector<double> y_bins_ccfr = {0.32, 0.57, 0.795};
	
	// const std::vector<double> E_bins_nutev = {90.18, 174.37, 244.72};
	// const std::vector<double> E_bins_ccfr = {109.46, 209.89, 332.7};

	// sidis.lepton_pair_cross_section_xy(
	// 	{0.017970, 0.057648, 0.10725, 0.12700, 0.15684, 0.20644}, 
	// 	{0.776},
	// 	{143.74},
	// 	"garbage.csv"
	// );
	
	return 0;
}

#endif