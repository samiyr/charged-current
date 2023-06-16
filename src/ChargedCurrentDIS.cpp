#ifndef CHARGED_CURRENT_DIS_H
#define CHARGED_CURRENT_DIS_H

#include <iostream>
#include "DIS.cpp"
#include "SIDIS.cpp"
#include "LHAInterface.cpp"
#include "Tests.cpp"
#include "TRFKinematics.cpp"
// #include "PythiaSIDIS.cpp"

int main() {
	LHAInterface::disable_verbosity();

	Tests::dis_cross_section_tests();
	return 0;

	const double N = 7.365;
	const double alpha = 1.4;
	const double beta = 2.276;
	const double gamma = 2.040;
	const double minimum_lepton_momentum = 5.0;
	const Particle target = Constants::Particles::Proton;
	const auto decay_function = DecayFunctions::decay_function;

	const DecayParametrization parametrization(N, alpha, beta, gamma);

	SIDIS sidis(
		{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
		LHAInterface("EPPS21nlo_CT18Anlo_Fe56"),
		FragmentationConfiguration(
			{
				LHAInterface("kkks08_opal_d0___mas"), 
				LHAInterface("kkks08_opal_d+___mas"), 
				LHAInterface("bkk05_D3_d_s_nlo"), 
				LHAInterface("bkk05_D3_lambda_c_nlo")
			},
			{
				Decay(parametrization, Constants::Particles::D0, target, decay_function, minimum_lepton_momentum),
				Decay(parametrization, Constants::Particles::Dp, target, decay_function, minimum_lepton_momentum),
				Decay(parametrization, Constants::Particles::Ds, target, decay_function, minimum_lepton_momentum),
				Decay(parametrization, Constants::Particles::LambdaC, target, decay_function, minimum_lepton_momentum)
			}
		),
		200'000,
		Process {Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino }
		// ScaleDependence::multiplicative(2.0),
		// ScaleDependence::multiplicative(2.0)
	);

	sidis.charm_mass = 1.3;

	sidis.max_chi_squared_deviation = 0.5;
	sidis.max_relative_error = 1e-2;
	sidis.iter_max = 5;
	sidis.combine_integrals = true;

	std::cout << "Remember alpha_s!" << std::endl;

	sidis.lepton_pair_cross_section_xy(
		{0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35}, 
		{0.334, 0.573, 0.790},
		{90.18, 174.37, 244.72},
		"neutrino_sidis_nutev.csv"
	);
	sidis.lepton_pair_cross_section_xy(
		{0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35}, 
		{0.32, 0.57, 0.795},
		{109.46, 209.89, 332.7},
		"neutrino_sidis_ccfr.csv"
	);
	// sidis.lepton_pair_cross_section_xy(
	// 	{0.02245, 0.03233, 0.0422, 0.10145, 0.15082, 0.20019}, 
	// 	{0.790},
	// 	{90.18},
	// 	"garbage.csv"
	// );
	// sidis.lepton_pair_cross_section(
	// 	{0.3, 0.35}, 
	// 	{0.32},
	// 	{109.46},
	// 	"garbage.csv"
	// );
	
	return 0;
}

#endif