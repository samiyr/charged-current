#ifndef CHARGED_CURRENT_DIS_H
#define CHARGED_CURRENT_DIS_H

#include <iostream>
#include "DIS.cpp"
#include "SIDIS.cpp"
#include "LHAInterface.cpp"
#include "Tests.cpp"
#include "TRFKinematics.cpp"
#include "Analysis.cpp"
// #include "PythiaSIDIS.cpp"
#include <array>

int main() {
	LHAInterface::disable_verbosity();

	const std::vector<double> x_bins = {0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35};

	return 0;

	Analysis::Inclusive::charm_production_nutev_old(x_bins, "Data/DIS/charm_production_nutev_old.csv");
	Analysis::Inclusive::charm_production_nutev_new(x_bins, "Data/DIS/charm_production_nutev_new.csv");
	Analysis::Inclusive::charm_production_ccfr(x_bins, "Data/DIS/charm_production_ccfr.csv");

	// Analysis::SemiInclusive::muon_pair_production_nutev_old(x_bins, "Data/SIDIS/CharmedHadrons/muon_pair_production_nutev_old.csv");
	// Analysis::SemiInclusive::muon_pair_production_nutev_new(x_bins, "Data/SIDIS/CharmedHadrons/muon_pair_production_nutev_new.csv");
	// Analysis::SemiInclusive::muon_pair_production_ccfr(x_bins, "Data/SIDIS/CharmedHadrons/muon_pair_production_ccfr.csv");

	// Analysis::SemiInclusive::muon_pair_production_nutev_old_only_D0(x_bins, "Data/SIDIS/D0/muon_pair_production_nutev_old.csv");
	// Analysis::SemiInclusive::muon_pair_production_nutev_new_only_D0(x_bins, "Data/SIDIS/D0/muon_pair_production_nutev_new.csv");
	// Analysis::SemiInclusive::muon_pair_production_ccfr_only_D0(x_bins, "Data/SIDIS/D0/muon_pair_production_ccfr.csv");

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