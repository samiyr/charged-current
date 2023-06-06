#ifndef CHARGED_CURRENT_DIS_H
#define CHARGED_CURRENT_DIS_H

#include <iostream>
#include "DIS.cpp"
#include "SIDIS.cpp"
#include "LHAInterface.cpp"
#include "TrivialInterface.cpp"
#include "DeltaInterface.cpp"
#include "SIDISComparison.cpp"
#include "Tests.cpp"
#include "TRFKinematics.cpp"
// #include "PythiaSIDIS.cpp"

int main() {
	// LHAInterface pdf("CT18ANLO");
	// LHAInterface ff("kkks08_global_d0_mas");

	// const double x = 0.4;
	// const double z = 0.1;
	// const double Q2 = 10.0;

	// pdf.evaluate(x, Q2);
	// ff.evaluate(z, Q2);

	// const double value = PDFCommon::xq_zq_sum(pdf, ff, FlavorInfo({Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom}), false, Process {Process::Type::NeutrinoToLepton});
	// std::cout << 2 * value / z << std::endl;
	// const double value3 = PDFCommon::xq_zq_sum(pdf, ff, FlavorInfo({Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom}), true, Process {Process::Type::NeutrinoToLepton});
	// std::cout << 2 * value3 / z << std::endl;


	// LHAInterface pdf2("CT18ANLO");
	// LHAInterface ff2("kkks08_global_d0_mas");

	// using namespace Constants;
	// double value_2 = 0.0;
	// value_2 += V_ud * pdf2.xf_evaluate(Flavor::Down, x, Q2) * ff2.xf_evaluate(Flavor::Up, z, Q2);
	// value_2 += V_us * pdf2.xf_evaluate(Flavor::Strange, x, Q2) * ff2.xf_evaluate(Flavor::Up, z, Q2);
	// value_2 += V_ub * pdf2.xf_evaluate(Flavor::Bottom, x, Q2) * ff2.xf_evaluate(Flavor::Up, z, Q2);
	// value_2 += V_cd * pdf2.xf_evaluate(Flavor::Down, x, Q2) * ff2.xf_evaluate(Flavor::Charm, z, Q2);
	// value_2 += V_cs * pdf2.xf_evaluate(Flavor::Strange, x, Q2) * ff2.xf_evaluate(Flavor::Charm, z, Q2);
	// value_2 += V_cb * pdf2.xf_evaluate(Flavor::Bottom, x, Q2) * ff2.xf_evaluate(Flavor::Charm, z, Q2);
	// value_2 += V_ud * pdf2.xf_evaluate(Flavor::AntiUp, x, Q2) * ff2.xf_evaluate(Flavor::AntiDown, z, Q2);
	// value_2 += V_us * pdf2.xf_evaluate(Flavor::AntiUp, x, Q2) * ff2.xf_evaluate(Flavor::AntiStrange, z, Q2);
	// value_2 += V_ub * pdf2.xf_evaluate(Flavor::AntiUp, x, Q2) * ff2.xf_evaluate(Flavor::AntiBottom, z, Q2);
	// value_2 += V_cd * pdf2.xf_evaluate(Flavor::AntiCharm, x, Q2) * ff2.xf_evaluate(Flavor::AntiDown, z, Q2);
	// value_2 += V_cs * pdf2.xf_evaluate(Flavor::AntiCharm, x, Q2) * ff2.xf_evaluate(Flavor::AntiStrange, z, Q2);
	// value_2 += V_cb * pdf2.xf_evaluate(Flavor::AntiCharm, x, Q2) * ff2.xf_evaluate(Flavor::AntiBottom, z, Q2);

	// double value3_2 = V_ud * pdf2.xf_evaluate(Flavor::Down, x, Q2) * ff2.xf_evaluate(Flavor::Up, z, Q2);
	// value3_2 += V_us * pdf2.xf_evaluate(Flavor::Strange, x, Q2) * ff2.xf_evaluate(Flavor::Up, z, Q2);
	// value3_2 += V_ub * pdf2.xf_evaluate(Flavor::Bottom, x, Q2) * ff2.xf_evaluate(Flavor::Up, z, Q2);
	// value3_2 += V_cd * pdf2.xf_evaluate(Flavor::Down, x, Q2) * ff2.xf_evaluate(Flavor::Charm, z, Q2);
	// value3_2 += V_cs * pdf2.xf_evaluate(Flavor::Strange, x, Q2) * ff2.xf_evaluate(Flavor::Charm, z, Q2);
	// value3_2 += V_cb * pdf2.xf_evaluate(Flavor::Bottom, x, Q2) * ff2.xf_evaluate(Flavor::Charm, z, Q2);
	// value3_2 -= V_ud * pdf2.xf_evaluate(Flavor::AntiUp, x, Q2) * ff2.xf_evaluate(Flavor::AntiDown, z, Q2);
	// value3_2 -= V_us * pdf2.xf_evaluate(Flavor::AntiUp, x, Q2) * ff2.xf_evaluate(Flavor::AntiStrange, z, Q2);
	// value3_2 -= V_ub * pdf2.xf_evaluate(Flavor::AntiUp, x, Q2) * ff2.xf_evaluate(Flavor::AntiBottom, z, Q2);
	// value3_2 -= V_cd * pdf2.xf_evaluate(Flavor::AntiCharm, x, Q2) * ff2.xf_evaluate(Flavor::AntiDown, z, Q2);
	// value3_2 -= V_cs * pdf2.xf_evaluate(Flavor::AntiCharm, x, Q2) * ff2.xf_evaluate(Flavor::AntiStrange, z, Q2);
	// value3_2 -= V_cb * pdf2.xf_evaluate(Flavor::AntiCharm, x, Q2) * ff2.xf_evaluate(Flavor::AntiBottom, z, Q2);

	// std::cout << 2 * value_2 / z << std::endl;
	// std::cout << 2 * value3_2 / z << std::endl;

	// Tests::sidis_lo_cross_section_integration_test();
	// return 0;
	// Tests::run_tests();
	// return 0;
	// const std::vector<double> Q2_list = {2.25, 10, 20, 50, 100, 500, 1000, 10'000, 50'000, 100'000, 1'000'000};
	// const std::vector<std::string> ff_sets = {
	// 	"kkks08_alep_d+st_m00",
	// 	"kkks08_alep_d+st_mas",
	// 	"kkks08_bell_d+st_m00",
	// 	"kkks08_bell_d+st_mas",
	// 	"kkks08_belle_d+__m00",
	// 	"kkks08_belle_d+__mas",
	// 	"kkks08_belle_d0__m00",
	// 	"kkks08_cleo_d0___mas",
	// 	"kkks08_glob_d+st_m00",
	// 	"kkks08_glob_d+st_mas",
	// 	"kkks08_global_d+_m00",
	// 	"kkks08_global_d+_mas",
	// 	"kkks08_global_d0_m00",
	// 	"kkks08_global_d0_mas",
	// 	"kkks08_opal_d+___m00",
	// 	"kkks08_opal_d+___mas",
	// 	"kkks08_opal_d+st_m00",
	// 	"kkks08_opal_d+st_mas",
	// 	"kkks08_opal_d0___m00",
	// 	"kkks08_opal_d0___mas"
	// };

	// for (auto const &set : ff_sets) {
	// 	Tests::evaluate_pdf(LHAInterface(set), Q2_list, "kkks08/" + set + "_data.csv", true);
	// }

	// const std::vector<double> Q2_list = {2.25, 10, 20, 50, 100, 500, 1000, 10'000, 50'000, 100'000, 1'000'000};
	// const std::vector<std::string> ff_sets = {
	// 	"bkk05_D3_d_s_lo",
	// 	"bkk05_D3_d_s_nlo",
	// 	"bkk05_D3_d+_lo",
	// 	"bkk05_D3_d+_nlo",
	// 	"bkk05_D3_d0_lo",
	// 	"bkk05_D3_d0_nlo",
	// 	"bkk05_D3_lambda_c_lo",
	// 	"bkk05_D3_lambda_c_nlo",
	// 	"bkk05_D3_opal_d_st_lo",
	// 	"bkk05_D3_opal_d_st_nlo"
	// };

	// for (auto const &set : ff_sets) {
	// 	Tests::evaluate_pdf(LHAInterface(set), Q2_list, "bkk05/" + set + "_data.csv", true);
	// }
	// return 0;

	// Tests::evaluate_pdf(LHAInterface("JAM20-SIDIS_FF_pion_nlo"), Q2_list, "jam.csv");

	// return 0;

	// PythiaSIDIS::cross_section(200'000'000, "pythia_sidis.csv");

	// return 0;

	// Tests::sidis_lo_cross_section_integration_test();
	// Tests::lepton_pair_lo_cross_section_integration_tests();
	// return 0;

	const double N = 7.365;
	const double alpha = 1.4;
	const double beta = 2.276;
	const double gamma = 2.0;
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
		100'000,
		Process {Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino }
	);
	sidis.global_sqrt_s = 21.5465;
	sidis.max_chi_squared_deviation = 0.5;
	sidis.max_relative_error = 1e-2;
	sidis.iter_max = 5;

	sidis.lepton_pair_cross_section(
		{0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35}, 
		{0.334},
		{90.2},
		"lepton_pair_data.csv"
	);
	
	return 0;
}

#endif