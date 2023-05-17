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
// #include "PythiaSIDIS.cpp"

int main(int argc, char const *argv[]) {
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

	// Tests::decay_function_tests_2();
	// return 0;
	Tests::run_tests();
	return 0;
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

	// Tests::evaluate_pdf(LHAInterface("JAM20-SIDIS_FF_pion_nlo"), Q2_list, "jam.csv");

	// return 0;

	// PythiaSIDIS::cross_section(200'000'000, "pythia_sidis.csv");

	// return 0;

	SIDIS sidis(
		318,
		{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
		LHAInterface("CT18ANLO"),
		// LHAInterface("JAM20-SIDIS_FF_pion_nlo"),
		// LHAInterface("JAM20-SIDIS_FF_kaon_nlo"),
		// LHAInterface("JAM20-SIDIS_FF_hadron_nlo"),
		LHAInterface("kkks08_global_d0_mas"),
		// LHAInterface("kkks08_opal_d0___mas"),
		// LHAInterface("kkks08_cleo_d0___mas"),
		// LHAInterface("kkks08_belle_d0__m00"),
		20'000,
		Process {Process::Type::NeutrinoToLepton}
	);
	sidis.max_chi_squared_deviation = 0.2;
	sidis.iter_max = 10;

	// std::cout << sidis.F2(x, z, Q2).lo << std::endl;
	// std::cout << sidis.FL(x, z, Q2).lo << std::endl;
	// std::cout << sidis.xF3(x, z, Q2).lo << std::endl;
	// return 0;

	sidis.differential_cross_section({0.0001, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}, {0.0001, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}, {100}, "sidis_cross_sections.csv");
	// sidis.differential_cross_section({0.002, 0.1, 0.2}, {0.1, 0.2, 0.3}, {10, 100, 1000}, "sidis_cross_sections.csv");

	// std::cout << sidis.F2(0.2, 0.3, 10).nlo << std::endl;
	// std::cout << sidis.FL(0.2, 0.3, 10).nlo << std::endl;
	// std::cout << sidis.xF3(0.2, 0.3, 10).nlo << std::endl;
	
	// const PerturbativeResult f2 = sidis.compute_structure_function(StructureFunction::F2, 0.01, 1 - 1e-12, 10);
	// std::cout << f2.nlo << std::endl;

	// DIS dis(318,
	// 	{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
	// 	// TrivialInterface(),
	// 	LHAInterface("CT18ANLO"),
	// 	// FunctionalFormInterface([](const FlavorType flavor, const double x, const double Q2) {
	// 	// 	return x * (1 - x);
	// 	// }),
	// 	100'000,
	// 	Process {Process::Type::NeutrinoToLepton}
	// );

	// std::cout << dis.F2(0.2, 10).nlo << std::endl;

	// const double sidis_result = SIDISComparison::F2_z_integral(0.2, 10, 1e-3);
	// std::cout << sidis_result << std::endl;
	// // dis.differential_cross_section({0.002}, {200}, "differential_cross_sections.csv");
	// dis.differential_cross_section({0.002, 0.1, 0.2, 0.5}, {50, 500, 1'000, 10'000}, "differential_cross_sections.csv");

	// dis.y_max = 0.9;	
	// dis.cross_section({2e2, 5e2, 1e3, 2e3, 3e3, 5e3, 1e4, 2e4, 3e4}, "output2.csv");

	// dis.differential_cross_section({1e-4, 2e-4, 3e-4, 5e-4, 1e-3, 2e-3, 3e-3, 5e-3, 1e-2}, {10}, "dis_cross_section_small_x.csv");
	// dis.compute_all_structure_function({1e-4, 2e-4, 3e-4, 5e-4, 1e-3, 2e-3, 3e-3, 5e-3, 1e-2}, {100}, "dis_structure_functions_small_x.csv");
	// dis.compute_all_structure_function({0.002, 0.1, 0.2, 0.5}, {200, 500, 1'000, 10'000}, "structure_functions.csv");
	// dis.compute_all_structure_function({0.002}, {200}, "structure_functions.csv");

	return 0;
}

#endif