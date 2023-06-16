#ifndef TESTS_H
#define TESTS_H

#include "DIS.cpp"
#include "SIDIS.cpp"
#include <iostream>
#include "Utility.cpp"
#include "LHAInterface.cpp"
#include "DecayFunctions.cpp"
#include "FunctionalFormInterface.cpp"

namespace Tests {
	/// Checks that LHAInterface evaluates to correct values at various values of x using the CT18NLO PDF set.
	/// Checks are done against hard-coded values.
	bool pdf_evaluation_tests() {
		bool flag = true;
		LHAInterface pdf("CT18NLO");

		std::cout << "CT18NLO central PDF test with up at Q^2 = 4 GeV^2:" << std::endl;
		flag &= double_comparison(pdf.xf_evaluate(Flavor::Up, 1e-3, 4), 0.488643);
		flag &= double_comparison(pdf.xf_evaluate(Flavor::Up, 1e-2, 4), 0.444769);
		flag &= double_comparison(pdf.xf_evaluate(Flavor::Up, 1e-1, 4), 0.64825);
		flag &= double_comparison(pdf.xf_evaluate(Flavor::Up, 0.2, 4), 0.700196);
		flag &= double_comparison(pdf.xf_evaluate(Flavor::Up, 0.5, 4), 0.287891);
		flag &= double_comparison(pdf.xf_evaluate(Flavor::Up, 0.9, 4), 0.00186158);
		
		std::cout << "CT18NLO central PDF test with down at Q^2 = 4 GeV^2:" << std::endl;

		flag &= double_comparison(pdf.xf_evaluate(Flavor::Down, 1e-3, 4), 0.447382);
		flag &= double_comparison(pdf.xf_evaluate(Flavor::Down, 1e-2, 4), 0.395838);
		flag &= double_comparison(pdf.xf_evaluate(Flavor::Down, 1e-1, 4), 0.438871);
		flag &= double_comparison(pdf.xf_evaluate(Flavor::Down, 0.2, 4), 0.369744);
		flag &= double_comparison(pdf.xf_evaluate(Flavor::Down, 0.5, 4), 0.0831107);
		flag &= double_comparison(pdf.xf_evaluate(Flavor::Down, 0.9, 4), 0.000204438);

		return flag;
	}

	/// Checks LHAInterface for using the evaluate() first and then call xf() usage pattern against the simpler xf_evaluate().
	bool pdf_comparison_tests() {
		bool flag = true;
		LHAInterface pdf1("CT18NLO");
		LHAInterface pdf2("CT18NLO");

		const double x = 0.1;
		const double Q2 = 200;
		const double tol = 1e-15;
		pdf1.evaluate(x, Q2);

		flag &= double_comparison(pdf1.xf(Flavor::Up), pdf2.xf_evaluate(Flavor::Up, x, Q2), tol);
		flag &= double_comparison(pdf1.xf(Flavor::Down), pdf2.xf_evaluate(Flavor::Down, x, Q2), tol);
		flag &= double_comparison(pdf1.xf(Flavor::Charm), pdf2.xf_evaluate(Flavor::Charm, x, Q2), tol);
		flag &= double_comparison(pdf1.xf(Flavor::Strange), pdf2.xf_evaluate(Flavor::Strange, x, Q2), tol);
		flag &= double_comparison(pdf1.xf(Flavor::Top), pdf2.xf_evaluate(Flavor::Top, x, Q2), tol);
		flag &= double_comparison(pdf1.xf(Flavor::Bottom), pdf2.xf_evaluate(Flavor::Bottom, x, Q2), tol);
		flag &= double_comparison(pdf1.xg(), pdf2.xf_evaluate(Flavor::Gluon, x, Q2), tol);

		return flag;
	}

	/// Checks that the function CKM::squared works properly by comparing the results directly to the various CKM matrix elements found in Constants namespace.
	bool ckm_tests() {
		bool flag = true;
		const double tol = 1e-12;

		flag &= double_comparison(CKM::squared(Flavor::Up, Flavor::Down), Constants::V_ud, tol);
		flag &= double_comparison(CKM::squared(Flavor::Down, Flavor::Up), Constants::V_ud, tol);
		
		flag &= double_comparison(CKM::squared(Flavor::Up, Flavor::Strange), Constants::V_us, tol);
		flag &= double_comparison(CKM::squared(Flavor::Strange, Flavor::Up), Constants::V_us, tol);
		
		flag &= double_comparison(CKM::squared(Flavor::Up, Flavor::Bottom), Constants::V_ub, tol);
		flag &= double_comparison(CKM::squared(Flavor::Bottom, Flavor::Up), Constants::V_ub, tol);


		flag &= double_comparison(CKM::squared(Flavor::Charm, Flavor::Down), Constants::V_cd, tol);
		flag &= double_comparison(CKM::squared(Flavor::Down, Flavor::Charm), Constants::V_cd, tol);
		
		flag &= double_comparison(CKM::squared(Flavor::Charm, Flavor::Strange), Constants::V_cs, tol);
		flag &= double_comparison(CKM::squared(Flavor::Strange, Flavor::Charm), Constants::V_cs, tol);
		
		flag &= double_comparison(CKM::squared(Flavor::Charm, Flavor::Bottom), Constants::V_cb, tol);
		flag &= double_comparison(CKM::squared(Flavor::Bottom, Flavor::Charm), Constants::V_cb, tol);


		flag &= double_comparison(CKM::squared(Flavor::Top, Flavor::Down), Constants::V_td, tol);
		flag &= double_comparison(CKM::squared(Flavor::Down, Flavor::Top), Constants::V_td, tol);
		
		flag &= double_comparison(CKM::squared(Flavor::Top, Flavor::Strange), Constants::V_ts, tol);
		flag &= double_comparison(CKM::squared(Flavor::Strange, Flavor::Top), Constants::V_ts, tol);
		
		flag &= double_comparison(CKM::squared(Flavor::Top, Flavor::Bottom), Constants::V_tb, tol);
		flag &= double_comparison(CKM::squared(Flavor::Bottom, Flavor::Top), Constants::V_tb, tol);

		const FlavorVector upper = {Flavor::Up, Flavor::Charm, Flavor::Top};
		for (size_t i = 0; i < upper.size(); i++) {
			for (size_t j = 0; j < upper.size(); j++) {
				flag &= double_comparison(CKM::squared(upper[i], upper[j]), 0, tol);
			}
		}
		const FlavorVector lower = {Flavor::Down, Flavor::Strange, Flavor::Bottom};
		for (size_t i = 0; i < lower.size(); i++) {
			for (size_t j = 0; j < lower.size(); j++) {
				flag &= double_comparison(CKM::squared(lower[i], lower[j]), 0, tol);
			}
		}

		return flag;
	}

	bool dis_cross_section_tests() {
		bool flag = true;

		DIS dis1(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			LHAInterface("CT18ANLO"),
			20'000,
			Process {Process::Type::NeutrinoToLepton, Particle(), Particle()}
		);

		dis1.max_chi_squared_deviation = 0.2;
		dis1.max_relative_error = 1e-3;
		dis1.iter_max = 10;
		dis1.compute_differential_cross_section_directly = true;

		DIS dis2(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			LHAInterface("CT18ANLO"),
			20'000,
			Process {Process::Type::NeutrinoToLepton, Particle(), Particle()}
		);

		dis2.max_chi_squared_deviation = 0.2;
		dis2.max_relative_error = 1e-3;
		dis2.iter_max = 10;
		dis2.compute_differential_cross_section_directly = false;

		const double sqrt_s = 318.0;
		const std::vector<double> x_values = {0.002, 0.1, 0.2, 0.5};
		const std::vector<double> Q2_values = {50.0, 500.0, 1000.0, 10'000.0};

		const std::vector<std::vector<double>> comparison_values = {
			{0.0000000384639, 0, 0, 0}, 
			{0.000000000240341, 0.000000000203083, 0.000000000171673, 0.0000000000271489}, 
			{0.0000000000765554, 0.0000000000596697, 0.0000000000500484, 0.00000000000854444}, 
			{0.00000000000522944, 0.00000000000328736, 0.00000000000262646, 0.000000000000416598}
		};

		#pragma omp parallel for collapse(2)
		for (size_t i = 0; i < x_values.size(); i++) {
			for (size_t j = 0; j < Q2_values.size(); j++) {
				const double x = x_values[i];
				const double Q2 = Q2_values[j];

				const double comparison_value = comparison_values[i][j];

				const TRFKinematics kinematics = TRFKinematics::Q2_sqrt_s(x, Q2, sqrt_s, 0.0, 0.0);

				const double computed_value_1 = dis1.differential_cross_section_xQ2(kinematics).nlo;
				const double computed_value_2 = dis2.differential_cross_section_xQ2(kinematics).nlo;

				#pragma omp critical 
				{
					std::cout << "Running DIS test at x = " << x << ", Q^2 = " << Q2 << std::endl;
					flag &= double_comparison_rel(computed_value_1, comparison_value, 1e-3);
					flag &= double_comparison_rel(computed_value_2, comparison_value, 1e-3);
					std::cout << std::endl;
				}
			}
		}
		return flag;
	}

	/// Checks SIDIS cross section values against analytically computed values (with Mathematica) for several values of x and z.
	/// To obtain analytical results, simple functional forms are chosen for PDFs and FFs. These forms include flavor-dependence.
	/// Here, only quarks are active (gluon PDFs and FFs = 0).
	bool sidis_structure_function_only_quarks_tests() {
		bool flag = true;

		SIDIS sidis(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			FunctionalFormInterface([](const FlavorType flavor, const double x, [[maybe_unused]] const double Q2) {
				if (Flavor::is_gluon(flavor)) { return 0.0; }
				const double f = double(flavor);
				return (2.0 * f * f - f + 1.0) * x * (1 - x);
			}),
			FunctionalFormInterface([](const FlavorType flavor, const double x, [[maybe_unused]] const double Q2) {
				if (Flavor::is_gluon(flavor)) { return 0.0; }
				const double f = double(flavor);
				return (2.0 * f * f + f + 1.0) * x * x * (1 - x) * (1 - x);
			}),
			100'000,
			Process {Process::Type::NeutrinoToLepton, Particle(), Particle()}
		);
		sidis.max_chi_squared_deviation = 0.2;
		sidis.iter_max = 10;

		const std::vector<double> xz_values = {0.001, 0.01, 0.1, 0.2, 0.5, 0.9, 0.99};

		const std::vector<std::vector<double>> F2_comparison_values = {
			{1.1424, 1.56154, 2.07818, 2.08919, 1.2504, 0.0832725, 0.00149219}, 
			{3.68039, 7.4315, 10.0706, 9.51852, 5.6935, 0.611975, 0.0152504}, 
			{-32.1797, -0.887595, 14.693, 15.4749, 23.3157, 6.17641, 0.172884}, 
			{-90.3869, -36.6069, -9.0347, 1.07716, 41.3219, 12.7353, 0.344892}, 
			{-215.198, -135.629, -79.0762, -27.2644, 107.161, 28.2993, 0.68805}, 
			{-112.376, -83.0067, -29.6816, 26.6385, 115.398, 19.3673, 0.390172}, 
			{-16.3319, -12.3907, 1.97406, 15.5945, 30.0035, 3.85295, 0.067476}
		};

		const std::vector<std::vector<double>> FL_comparison_values = {
			{0.0108889, 0.062903, 0.200796, 0.19703, 0.0685839, 0.000653087, 0.000000669252}, 
			{0.103604, 0.598505, 1.91052, 1.87469, 0.652557, 0.00621394, 0.00000636775}, 
			{0.735085, 4.24645, 13.5553, 13.3011, 4.62996, 0.0440886, 0.0000451799}, 
			{1.04952, 6.06288, 19.3536, 18.9907, 6.61044, 0.0629476, 0.0000645056}, 
			{0.841977, 4.86395, 15.5265, 15.2353, 5.30323, 0.0504998, 0.0000517497}, 
			{0.0511244, 0.295336, 0.942758, 0.925077, 0.322009, 0.00306632, 0.00000314221}, 
			{0.000545115, 0.00314903, 0.0100522, 0.00986365, 0.00343343, 0.0000326947, 0.0000000335039}
		};

		const std::vector<std::vector<double>> xF3_comparison_values = {
			{-0.0000000000000000277556, -0.00000000000000000346945, 0.0, 0.000000000000000111022, 0.0, 0.0, 0.0}, 
			{-0.000000000000000444089, 0.0, 0.0, 0.0, 0.0, 0.0000000000000000555112, -0.000000000000000000867362}, 
			{0.00000000000000355271, 0.00000000000000177636, 0.0, 0.0, 0.0, 0.0, 0.0}, 
			{0.0, -0.00000000000000355271, 0.0, 0.00000000000000177636, 0.0, 0.0, 0.0}, 
			{-0.0000000000000142109, 0.0, 0.0, 0.00000000000000355271, 0.00000000000000710543, 0.00000000000000177636, 0.0000000000000000555112}, 
			{-0.00000000000000710543, 0.00000000000000710543, 0.00000000000000177636, 0.0, 0.0, 0.0, 0.0000000000000000277556}, 
			{0.00000000000000177636, -0.000000000000000888178, 0.0, 0.0, -0.00000000000000177636, -0.000000000000000222045, 0.0}
		};

		#pragma omp parallel for collapse(2)
		for (size_t i = 0; i < xz_values.size(); i++) {
			for (size_t j = 0; j < xz_values.size(); j++) {
				const double x = xz_values[i];
				const double z = xz_values[j];

				const double comparison_value_F2 = F2_comparison_values[i][j];
				const double comparison_value_FL = FL_comparison_values[i][j];
				const double comparison_value_xF3 = xF3_comparison_values[i][j];

				const double computed_value_F2 = sidis.F2(x, z, 1).nlo;
				const double computed_value_FL = sidis.FL(x, z, 1).nlo;
				const double computed_value_xF3 = sidis.xF3(x, z, 1).nlo;

				#pragma omp critical 
				{
					std::cout << "Running SIDIS test at x = " << x << ", z = " << z << " (only quarks):" << std::endl;
					flag &= double_comparison(computed_value_F2, comparison_value_F2, 1e-1);
					flag &= double_comparison(computed_value_FL, comparison_value_FL, 1e-1);
					flag &= double_comparison(computed_value_xF3, comparison_value_xF3, 1e-1);
					std::cout << std::endl;
				}
			}
		}
		return flag;
	}

	/// Checks SIDIS cross section values against analytically computed values (with Mathematica) for several values of x and z.
	/// To obtain analytical results, simple functional forms are chosen for PDFs and FFs. These forms include flavor-dependence.
	/// Here, PDF quarks = 0.
	bool sidis_structure_function_only_gluons_tests_1() {
		bool flag = true;

		SIDIS sidis(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			FunctionalFormInterface([](const FlavorType flavor, const double x, [[maybe_unused]] const double Q2) {
				if (Flavor::is_gluon(flavor)) { return x * x * (1 - x); }
				return 0.0;
			}),
			FunctionalFormInterface([](const FlavorType flavor, const double x, [[maybe_unused]] const double Q2) {
				if (Flavor::is_gluon(flavor)) { return x * x * (1 - x) * (1 - x); }
				const double f = double(flavor);
				return (2.0 * f * f + f + 1.0) * x * (1 - x) * (1 - x);
			}),
			100'000,
			Process {Process::Type::NeutrinoToLepton, Particle(), Particle()}
		);
		sidis.max_chi_squared_deviation = 0.2;
		sidis.iter_max = 10;

		const std::vector<double> xz_values = {0.001, 0.01, 0.1, 0.2, 0.5, 0.9, 0.99};

		const std::vector<std::vector<double>> F2_comparison_values = {
			{1.71144, 0.178639, 0.0283324, 0.0171812, 0.00485811, 0.0000894449, -0.000000323709}, 
			{15.2803, 1.55749, 0.192539, 0.100062, 0.0202343, -0.000171179, -0.0000126463}, 
			{92.5022, 10.185, 1.17213, 0.473529, 0.0139155, -0.00639044, -0.000130204}, 
			{135.067, 15.1416, 1.65274, 0.583375, -0.0474123, -0.0124968, -0.000221879}, 
			{150.496, 15.1138, 0.824673, -0.040859, -0.275829, -0.0220655, -0.000329451}, 
			{14.3418, 1.16663, -0.0746113, -0.106349, -0.0585699, -0.00326477, -0.0000431373}, 
			{0.165916, 0.0120124, -0.00202054, -0.00211138, -0.00100764, -0.0000507824, -0.000000630181}
		};

		const std::vector<std::vector<double>> FL_comparison_values = {
			{0.00223694, 0.00129224, 0.000412503, 0.000202383, 0.0000281789, 0.000000149074, 0.000000000138876}, 
			{0.121533, 0.0702074, 0.0224113, 0.0109955, 0.00153096, 0.00000809917, 0.00000000754512}, 
			{3.33423, 1.92612, 0.614848, 0.301658, 0.0420016, 0.000222199, 0.000000206999}, 
			{6.02974, 3.48328, 1.11191, 0.54553, 0.0759572, 0.000401833, 0.000000374344}, 
			{4.51795, 2.60994, 0.833132, 0.408754, 0.0569131, 0.000301084, 0.000000280488}, 
			{0.06817, 0.0393806, 0.0125709, 0.00616756, 0.000858744, 0.00000454297, 0.0000000042322}, 
			{0.0000750695, 0.0000433663, 0.0000138432, 0.00000679178, 0.000000945658, 0.00000000500277, 0.00000000000465931}
		};

		const std::vector<std::vector<double>> xF3_comparison_values = {
			{-0.749496, -0.0558365, 0.00597865, 0.00600071, 0.00211888, 0.000040437, -0.000000146823}, 
			{-6.7341, -0.541772, 0.0175123, 0.0252338, 0.00784727, -0.0000816872, -0.00000573676}, 
			{-40.1092, -3.45333, -0.100211, -0.0117224, -0.016474, -0.00300032, -0.0000591232}, 
			{-58.3411, -5.15704, -0.266749, -0.112897, -0.061358, -0.00585123, -0.00010076}, 
			{-66.6961, -6.20683, -0.586031, -0.351738, -0.157024, -0.0101441, -0.000149487}, 
			{-6.61831, -0.658788, -0.0968756, -0.0655799, -0.0275504, -0.00148256, -0.0000195586}, 
			{-0.0780417, -0.00826805, -0.00159366, -0.00113085, -0.000464368, -0.0000230295, -0.0000002857}
		};

		#pragma omp parallel for collapse(2)
		for (size_t i = 0; i < xz_values.size(); i++) {
			for (size_t j = 0; j < xz_values.size(); j++) {
				const double x = xz_values[i];
				const double z = xz_values[j];

				const double comparison_value_F2 = F2_comparison_values[i][j];
				const double comparison_value_FL = FL_comparison_values[i][j];
				const double comparison_value_xF3 = xF3_comparison_values[i][j];

				const double computed_value_F2 = sidis.F2(x, z, 1).nlo;
				const double computed_value_FL = sidis.FL(x, z, 1).nlo;
				const double computed_value_xF3 = sidis.xF3(x, z, 1).nlo;

				#pragma omp critical 
				{
					std::cout << "Running SIDIS test at x = " << x << ", z = " << z << " (only gluons on PDF side):" << std::endl;
					flag &= double_comparison(computed_value_F2, comparison_value_F2, 1e-1);
					flag &= double_comparison(computed_value_FL, comparison_value_FL, 1e-1);
					flag &= double_comparison(computed_value_xF3, comparison_value_xF3, 1e-1);
					std::cout << std::endl;
				}
			}
		}
		return flag;
	}
	/// Checks SIDIS cross section values against analytically computed values (with Mathematica) for several values of x and z.
	/// To obtain analytical results, simple functional forms are chosen for PDFs and FFs. These forms include flavor-dependence.
	/// Here, FF quarks = 0.
	bool sidis_structure_function_only_gluons_tests_2() {
		bool flag = true;

		SIDIS sidis(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			FunctionalFormInterface([](const FlavorType flavor, const double x, [[maybe_unused]] const double Q2) {
				if (Flavor::is_gluon(flavor)) { return x * x * (1 - x); }
				const double f = double(flavor);
				return (2.0 * f * f - f + 1.0) * x * (1 - x);
			}),
			FunctionalFormInterface([](const FlavorType flavor, const double x, [[maybe_unused]] const double Q2) {
				if (Flavor::is_gluon(flavor)) { return x * x * (1 - x) * (1 - x); }
				return 0.0;
			}),
			100'000,
			Process {Process::Type::NeutrinoToLepton, Particle(), Particle()}
		);
		sidis.max_chi_squared_deviation = 0.2;
		sidis.iter_max = 10;

		const std::vector<double> xz_values = {0.001, 0.01, 0.1, 0.2, 0.5, 0.9, 0.99};

		const std::vector<std::vector<double>> F2_comparison_values = {
			{-15.6057, -0.411576, 0.0928191, 0.058628, 0.0117651, 0.0000724169, 0.0000000503215}, 
			{-208.34, -9.47338, 0.344285, 0.298162, 0.06181, 0.000297243, 0.0000000790771}, 
			{-2418.76, -139.129, -2.52197, 0.00291112, 0.0681222, -0.000931589, -0.0000028851}, 
			{-4647.68, -282.533, -8.12488, -1.6551, -0.150944, -0.00350813, -0.00000694352}, 
			{-8346.65, -550.237, -23.0991, -7.04779, -0.893015, -0.00964476, -0.0000148829}, 
			{-3803.33, -276.333, -14.7207, -5.03349, -0.647906, -0.00545929, -0.00000727767}, 
			{-526.674, -40.8833, -2.4102, -0.847224, -0.108147, -0.000829393, -0.00000102315}
		};

		const std::vector<std::vector<double>> FL_comparison_values = {
			{0.0363701, 0.032505, 0.0159355, 0.00809662, 0.000844724, 0.000000986697, 0.0000000000932677}, 
			{0.346052, 0.309276, 0.151622, 0.0770371, 0.00803732, 0.00000938815, 0.000000000887416}, 
			{2.45527, 2.19435, 1.07577, 0.546587, 0.0570256, 0.0000666099, 0.00000000629632}, 
			{3.50552, 3.13299, 1.53594, 0.78039, 0.0814185, 0.0000951024, 0.00000000898958}, 
			{2.81231, 2.51344, 1.23221, 0.626069, 0.065318, 0.000076296, 0.0000000072119}, 
			{0.170762, 0.152615, 0.074819, 0.0380145, 0.00396607, 0.00000463265, 0.000000000437902}, 
			{0.00182075, 0.00162726, 0.000797759, 0.000405331, 0.0000422883, 0.0000000493957, 0.00000000000466718}
		};

		const std::vector<std::vector<double>> xF3_comparison_values = {
			{7.09212, 0.205233, -0.0223883, -0.0106752, -0.000692474, 0.00000816677, 0.0000000187826}, 
			{94.6129, 4.45626, -0.0196099, -0.0337837, -0.0012425, 0.0000898082, 0.000000190314}, 
			{1097.69, 64.1283, 1.8169, 0.428825, 0.0584417, 0.00105695, 0.00000193021}, 
			{2108.66, 129.562, 4.54784, 1.26903, 0.162736, 0.00218007, 0.00000371202}, 
			{3785.3, 250.609, 11.0811, 3.52842, 0.451669, 0.00457091, 0.00000691839}, 
			{1724.34, 125.347, 6.70807, 2.29962, 0.295675, 0.00247849, 0.00000330099}, 
			{238.772, 18.5355, 1.09305, 0.38428, 0.0490486, 0.000376036, 0.000000463858}
		};

		#pragma omp parallel for collapse(2)
		for (size_t i = 0; i < xz_values.size(); i++) {
			for (size_t j = 0; j < xz_values.size(); j++) {
				const double x = xz_values[i];
				const double z = xz_values[j];

				const double comparison_value_F2 = F2_comparison_values[i][j];
				const double comparison_value_FL = FL_comparison_values[i][j];
				const double comparison_value_xF3 = xF3_comparison_values[i][j];

				const double computed_value_F2 = sidis.F2(x, z, 1).nlo;
				const double computed_value_FL = sidis.FL(x, z, 1).nlo;
				const double computed_value_xF3 = sidis.xF3(x, z, 1).nlo;

				#pragma omp critical 
				{
					std::cout << "Running SIDIS test at x = " << x << ", z = " << z << " (only gluons on FF side):" << std::endl;
					flag &= double_comparison(computed_value_F2, comparison_value_F2, 1e-1);
					flag &= double_comparison(computed_value_FL, comparison_value_FL, 1e-1);
					flag &= double_comparison(computed_value_xF3, comparison_value_xF3, 1e-1);
					std::cout << std::endl;
				}
			}
		}
		return flag;
	}

	/// Checks SIDIS cross section values against analytically computed values (with Mathematica) for several values of x and z.
	/// To obtain analytical results, simple functional forms are chosen for PDFs and FFs. These forms include flavor-dependence.
	/// Here, all quarks and gluons are enabled.
	bool sidis_structure_function_quarks_gluons_tests() {
		bool flag = true;

		SIDIS sidis(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			FunctionalFormInterface([](const FlavorType flavor, const double x, [[maybe_unused]] const double Q2) {
				if (Flavor::is_gluon(flavor)) { return x * x * (1 - x); }
				const double f = double(flavor);
				return (2.0 * f * f - f + 1.0) * x * (1 - x);
			}),
			FunctionalFormInterface([](const FlavorType flavor, const double x, [[maybe_unused]] const double Q2) {
				if (Flavor::is_gluon(flavor)) { return x * x * (1 - x) * (1 - x); }
				const double f = double(flavor);
				return (2.0 * f * f + f + 1.0) * x * (1 - x) * (1 - x);
			}),
			100'000,
			Process {Process::Type::NeutrinoToLepton, Particle(), Particle()}
		);
		sidis.max_chi_squared_deviation = 0.2;
		sidis.iter_max = 10;

		const std::vector<double> xz_values = {0.001, 0.01, 0.1, 0.2, 0.5, 0.9, 0.99};

		const std::vector<std::vector<double>> F2_comparison_values = {
			{19.5802, 24.359, 12.6699, 8.25673, 2.42879, 0.0932363, 0.00150834}, 
			{-40.2698, 110.152, 59.0232, 37.7309, 11.613, 0.68987, 0.0154091}, 
			{-2288.32, 50.8941, 130.886, 102.382, 58.1581, 6.98666, 0.174693}, 
			{-5060.91, -331.163, 95.6154, 121.536, 109.483, 14.396, 0.34852}, 
			{-10202.3, -1265.1, 56.332, 242.241, 275.266, 31.9082, 0.695319}, 
			{-4661.11, -446.285, 309.216, 384.455, 265.526, 21.7525, 0.39436}, 
			{-572.835, 16.5798, 121.449, 119.423, 65.4672, 4.31532, 0.0681965}
		};

		const std::vector<std::vector<double>> FL_comparison_values = {
			{2.02362, 1.86103, 1.0821, 0.644802, 0.115307, 0.00070778, 0.000000674547}, 
			{19.3544, 17.7651, 10.3144, 6.14418, 1.09838, 0.00674101, 0.00000642435}, 
			{139.793, 127.473, 73.6375, 43.8173, 7.82426, 0.0479929, 0.0000457349}, 
			{200.86, 182.733, 105.37, 62.6751, 11.1871, 0.0686066, 0.0000653769}, 
			{160.82, 146.413, 84.4743, 50.2522, 8.97083, 0.0550184, 0.0000524288}, 
			{9.55876, 8.77104, 5.09122, 3.03263, 0.542107, 0.00332695, 0.00000317065}, 
			{0.101269, 0.0931449, 0.054165, 0.0322766, 0.00577202, 0.0000354302, 0.0000000337666}
		};

		const std::vector<std::vector<double>> xF3_comparison_values = {
			{6.34263, 0.149397, -0.0164097, -0.00467449, 0.0014264, 0.0000486038, -0.00000012804}, 
			{87.8788, 3.91449, -0.00209762, -0.00854993, 0.00660478, 0.000008121, -0.00000554645}, 
			{1057.58, 60.675, 1.71668, 0.417102, 0.0419678, -0.00194336, -0.000057193}, 
			{2050.32, 124.405, 4.28109, 1.15613, 0.101378, -0.00367116, -0.0000970484}, 
			{3718.6, 244.403, 10.4951, 3.17668, 0.294645, -0.00557314, -0.000142568}, 
			{1717.73, 124.688, 6.6112, 2.23404, 0.268125, 0.000995922, -0.0000162576}, 
			{238.694, 18.5272, 1.09145, 0.383149, 0.0485842, 0.000353006, 0.000000178159}
		};

		#pragma omp parallel for collapse(2)
		for (size_t i = 0; i < xz_values.size(); i++) {
			for (size_t j = 0; j < xz_values.size(); j++) {
				const double x = xz_values[i];
				const double z = xz_values[j];

				const double comparison_value_F2 = F2_comparison_values[i][j];
				const double comparison_value_FL = FL_comparison_values[i][j];
				const double comparison_value_xF3 = xF3_comparison_values[i][j];

				const double computed_value_F2 = sidis.F2(x, z, 1).nlo;
				const double computed_value_FL = sidis.FL(x, z, 1).nlo;
				const double computed_value_xF3 = sidis.xF3(x, z, 1).nlo;

				#pragma omp critical 
				{
					std::cout << "Running SIDIS test at x = " << x << ", z = " << z << " (quarks and gluons):" << std::endl;
					flag &= double_comparison(computed_value_F2, comparison_value_F2, 1e-1);
					flag &= double_comparison(computed_value_FL, comparison_value_FL, 1e-1);
					flag &= double_comparison(computed_value_xF3, comparison_value_xF3, 1e-1);
					std::cout << std::endl;
				}
			}
		}
		return flag;
	}

	bool sidis_structure_function_quarks_gluons_scale_logs_tests() {
		bool flag = true;

		SIDIS sidis(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			FunctionalFormInterface([](const FlavorType flavor, const double x, [[maybe_unused]] const double Q2) {
				if (Flavor::is_gluon(flavor)) { return x * x * (1 - x); }
				const double f = double(flavor);
				return (2.0 * f * f - f + 1.0) * x * (1 - x);
			}),
			FunctionalFormInterface([](const FlavorType flavor, const double x, [[maybe_unused]] const double Q2) {
				if (Flavor::is_gluon(flavor)) { return x * x * (1 - x) * (1 - x); }
				const double f = double(flavor);
				return (2.0 * f * f + f + 1.0) * x * (1 - x) * (1 - x);
			}),
			200'000,
			Process {Process::Type::NeutrinoToLepton, Particle(), Particle()},
			ScaleDependence::constant(5.0),
			ScaleDependence::constant(30.0)
		);
		sidis.max_chi_squared_deviation = 0.2;
		sidis.iter_max = 10;

		const std::vector<double> xz_values = {0.001, 0.01, 0.1, 0.2, 0.5, 0.9, 0.99};

		const std::vector<std::vector<double>> F2_comparison_values = {
			{20.9822, 27.874, 16.1389, 11.1343, 3.6576, 0.150557, 0.0021788}, 
			{-42.5053, 129.146, 80.3101, 55.9039, 19.75, 1.0963, 0.0204371}, 
			{-2466.3, 68.7253, 196.439, 166.484, 92.6354, 9.10165, 0.204603}, 
			{-5481.73, -401.989, 127.422, 168.546, 144.624, 17.1099, 0.391234}, 
			{-11185.8, -1695.81, -158.473, 106.705, 248.538, 32.8832, 0.729404}, 
			{-5255.72, -837.606, 36.6268, 181.382, 195.638, 19.6929, 0.382525}, 
			{-670.904, -58.5411, 64.9554, 76.1398, 49.5978, 3.7615, 0.0636219}
		};

		const std::vector<std::vector<double>> FL_comparison_values = {
			{2.02362, 1.86103, 1.0821, 0.644802, 0.115307, 0.00070778, 0.000000674547}, 
			{19.3544, 17.7651, 10.3144, 6.14418, 1.09838, 0.00674101, 0.00000642435}, 
			{139.793, 127.473, 73.6375, 43.8173, 7.82426, 0.0479929, 0.0000457349}, 
			{200.86, 182.733, 105.37, 62.6751, 11.1871, 0.0686066, 0.0000653769},
			{160.82, 146.413, 84.4743, 50.2522, 8.97083, 0.0550184, 0.0000524288}, 
			{9.55876, 8.77104, 5.09122, 3.03263, 0.542107, 0.00332695, 0.00000317065}, 
			{0.101269, 0.0931449, 0.054165, 0.0322766, 0.00577202, 0.0000354302, 0.0000000337666}
		};

		const std::vector<std::vector<double>> xF3_comparison_values = {
			{7.20047, 0.235145, -0.00764879, -0.000344872, 0.00251706, 0.0000826894, 0.000000196879}, 
			{96.3767, 4.76113, 0.082146, 0.0323208, 0.0166181, 0.000314103, -0.00000264456}, 
			{1134.74, 68.2821, 2.40841, 0.730088, 0.11012, -0.0000768031, -0.0000399631}, 
			{2187.44, 137.874, 5.46561, 1.67682, 0.208582, -0.000911088, -0.0000720002}, 
			{3932.75, 265.345, 12.2606, 3.92284, 0.435814, -0.00231401, -0.000113965}, 
			{1794.75, 132.155, 7.1871, 2.45551, 0.300528, 0.00143252, -0.0000133274}, 
			{247.164, 19.3459, 1.15259, 0.405764, 0.0514666, 0.000373753, 0.000000227695}
		};

		#pragma omp parallel for collapse(2)
		for (size_t i = 0; i < xz_values.size(); i++) {
			for (size_t j = 0; j < xz_values.size(); j++) {
				const double x = xz_values[i];
				const double z = xz_values[j];

				const double comparison_value_F2 = F2_comparison_values[i][j];
				const double comparison_value_FL = FL_comparison_values[i][j];
				const double comparison_value_xF3 = xF3_comparison_values[i][j];

				const double computed_value_F2 = sidis.F2(x, z, 20.0).nlo;
				const double computed_value_FL = sidis.FL(x, z, 20.0).nlo;
				const double computed_value_xF3 = sidis.xF3(x, z, 20.0).nlo;

				#pragma omp critical
				{
					std::cout << "Running SIDIS test at x = " << x << ", z = " << z << " (quarks and gluons + scale logs):" << std::endl;
					flag &= double_comparison_rel(computed_value_F2, comparison_value_F2, 1e-3);
					flag &= double_comparison_rel(computed_value_FL, comparison_value_FL, 1e-3);
					flag &= double_comparison_rel(computed_value_xF3, comparison_value_xF3, 1e-3);
					std::cout << std::endl;
				}
			}
		}
		return flag;
	}

	/// Tests the evaluation of DecayFunctions::decay_function against a set of values computed with Mathematica.
	bool decay_function_tests_1() {
		bool flag = true;
		
		const double z_min = 0.1;
		const DecayParametrization param(1.0, 1.4, 2.3, 2.0);

		const Particle resonance = Particle(1.8, 1.0);
		const Particle target = Particle(1.0);

		flag &= double_comparison(DecayFunctions::decay_function(0.1, 0.3, 10, z_min, param, resonance, target), 2 * M_PI * 0.00100292);
		flag &= double_comparison(DecayFunctions::decay_function(0.2, 0.3, 10, z_min, param, resonance, target), 2 * M_PI * 0.000995329);
		flag &= double_comparison(DecayFunctions::decay_function(0.3, 0.3, 10, z_min, param, resonance, target), 2 * M_PI * 0.00098196);
		flag &= double_comparison(DecayFunctions::decay_function(0.4, 0.5, 10, z_min, param, resonance, target), 2 * M_PI * 0.00180983);
		flag &= double_comparison(DecayFunctions::decay_function(0.5, 0.5, 10, z_min, param, resonance, target), 2 * M_PI * 0.00181612);
		flag &= double_comparison(DecayFunctions::decay_function(0.7, 0.9, 10, z_min, param, resonance, target), 2 * M_PI * 0.00247266);
		flag &= double_comparison(DecayFunctions::decay_function(0.99, 0.9, 10, z_min, param, resonance, target), 2 * M_PI * 0.00251066);

		return flag;
	}

	/// Tests the analytically obtained decay function DecayFunctions::decay_function against numerically integrated version 
	/// starting from the integrand DecayFunctions::decay_function_integrand.
	bool decay_function_tests_2() {
		bool flag = true;
		
		const double z_min = 0.1;
		const DecayParametrization param(1.0, 1.4, 2.3, 2.0);

		const Particle resonance = Particle(1.8, 1.0);
		const Particle target = Particle(1.0);

		const std::vector<double> xz_values = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
		const std::vector<double> Q2_values = {10.0, 20.0, 50.0, 100.0};

		#pragma omp parallel for collapse(3)
		for (const double x : xz_values) {
			for (const double z : xz_values) {
				for (const double Q2 : Q2_values) {
					std::vector<double> params = {x, z, Q2, 1.0, 1.8, 1.0, 1.4, 2.3, 2.0, 1.0};
					Integrator integrator(&DecayFunctions::decay_function_integrand, {0.1 / params[1], -1}, {1, 1}, 2'000'000, &params, 0.2, 1e-3, 10);
					const auto result = integrator.integrate();
					if (result.value == 0) {
						std::cout << "Skipping value due to non-applicable region" << std::endl;
						continue;
					}
					flag &= double_comparison_rel(DecayFunctions::decay_function(x, z, Q2, z_min, param, resonance, target), result.value, 1e-5);
				}
			}
		}

		return flag;
	}

	bool lepton_pair_lo_cross_section_integration_tests() {
		bool flag = true;

		const double sqrt_s = 120.0;
		const double x = 0.3;
		const double Q2 = 10.0;

		DecayParametrization parametrization(7.365, 1.4, 2.276, 2.0);

		// Since lambdas have a unique type, have to first convert them to e.g. std::function and then pass in a vector, otherwise the template DecayFunction cannot be uniquely defined.
		auto lambda1 = []([[maybe_unused]] const double x, [[maybe_unused]] const double z, [[maybe_unused]] const double Q2, [[maybe_unused]] const double z_min, [[maybe_unused]] const DecayParametrization &decay, [[maybe_unused]] const Particle &resonance, [[maybe_unused]] const Particle &hadron) { 
			return 3 * Constants::Particles::D0.lifetime;
		};
		Decay decay1(parametrization, Constants::Particles::D0, Constants::Particles::Proton, std::function(lambda1), 5.0);

		auto lambda2 = []([[maybe_unused]] const double x, [[maybe_unused]] const double z, [[maybe_unused]] const double Q2, [[maybe_unused]] const double z_min, [[maybe_unused]] const DecayParametrization &decay, [[maybe_unused]] const Particle &resonance, [[maybe_unused]] const Particle &hadron) { 
			return 5 * Constants::Particles::Dp.lifetime; 
		};
		Decay decay2(parametrization, Constants::Particles::D0, Constants::Particles::Proton, std::function(lambda2), 5.0);

		SIDIS sidis(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			LHAInterface("EPPS21nlo_CT18Anlo_Fe56"),
			FragmentationConfiguration(
				{
					LHAInterface("kkks08_opal_d0___mas"), 
					LHAInterface("kkks08_opal_d+___mas")
				},
				{
					decay1,
					decay2
				}
			),
			100'000,
			Process {Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino}
		);
		sidis.global_sqrt_s = sqrt_s;
		sidis.max_chi_squared_deviation = 0.2;
		sidis.max_relative_error = 1e-3;
		sidis.iter_max = 10;

		Decay decay(parametrization, Particle(), Particle(), DecayFunctions::trivial, 5.0);

		TRFKinematics kinematics = TRFKinematics::Q2_sqrt_s(x, Q2, sqrt_s, Constants::Particles::Proton.mass, 0.0);

		const double z_min = SIDISFunctions::compute_z_min(kinematics, decay);

		const PerturbativeQuantity result = sidis.lepton_pair_cross_section_xQ2(kinematics);
		std::cout << "Lepton-pair cross section value = " << result.lo << std::endl;

		SIDIS sidis1(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			LHAInterface("EPPS21nlo_CT18Anlo_Fe56"),
			LHAInterface("kkks08_opal_d0___mas"),
			10'000,
			Process {Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino}
		);
		sidis1.global_sqrt_s = sqrt_s;
		sidis1.max_chi_squared_deviation = 0.2;
		sidis1.max_relative_error = 1e-3;
		sidis1.iter_max = 10;

		SIDIS sidis2(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			LHAInterface("EPPS21nlo_CT18Anlo_Fe56"),
			LHAInterface("kkks08_opal_d+___mas"),
			10'000,
			Process {Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino}
		);
		sidis2.global_sqrt_s = sqrt_s;
		sidis2.max_chi_squared_deviation = 0.2;
		sidis2.max_relative_error = 1e-3;
		sidis2.iter_max = 10;

		Integrator integrator([&](const double input[], [[maybe_unused]] const size_t dim, [[maybe_unused]] void *params) {
			const double z = input[0];

			const double differential_cs_1 = sidis1.differential_cross_section_xQ2(z, kinematics).lo;
			const double differential_cs_2 = sidis2.differential_cross_section_xQ2(z, kinematics).lo;

			const double result = 3 * Constants::Particles::D0.lifetime * differential_cs_1 + 5 * Constants::Particles::Dp.lifetime * differential_cs_2;
			return result;
		}, {z_min}, {1}, 100, nullptr, 0.5, 1e-2, 10);
		integrator.verbose = true;
		const auto result2 = integrator.integrate();
		std::cout << "Lepton-pair cross section value = " << result2 << std::endl;
		flag &= double_comparison_rel(result.lo, result2.value, 1e-1);

		return flag;
	}

	bool lepton_pair_nlo_cross_section_integration_tests() {
		bool flag = true;

		const double sqrt_s = 120.0;
		const double x = 0.3;
		const double Q2 = 10.0;

		DecayParametrization parametrization(7.365, 1.4, 2.276, 2.0);

		// Since lambdas have a unique type, have to first convert them to e.g. std::function and then pass in a vector, otherwise the template DecayFunction cannot be uniquely defined.
		auto lambda1 = []([[maybe_unused]] const double x, [[maybe_unused]] const double z, [[maybe_unused]] const double Q2, [[maybe_unused]] const double z_min, [[maybe_unused]] const DecayParametrization &decay, [[maybe_unused]] const Particle &resonance, [[maybe_unused]] const Particle &hadron) { 
			return 3 * Constants::Particles::D0.lifetime;
		};
		Decay decay1(parametrization, Constants::Particles::D0, Constants::Particles::Proton, std::function(lambda1), 5.0);

		auto lambda2 = []([[maybe_unused]] const double x, [[maybe_unused]] const double z, [[maybe_unused]] const double Q2, [[maybe_unused]] const double z_min, [[maybe_unused]] const DecayParametrization &decay, [[maybe_unused]] const Particle &resonance, [[maybe_unused]] const Particle &hadron) { 
			return 5 * Constants::Particles::Dp.lifetime; 
		};
		Decay decay2(parametrization, Constants::Particles::D0, Constants::Particles::Proton, std::function(lambda2), 5.0);

		SIDIS sidis(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			LHAInterface("EPPS21nlo_CT18Anlo_Fe56"),
			FragmentationConfiguration(
				{
					LHAInterface("kkks08_opal_d0___mas"), 
					LHAInterface("kkks08_opal_d+___mas")
				},
				{
					decay1,
					decay2
				}
			),
			100'000,
			Process {Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino}
		);
		sidis.global_sqrt_s = sqrt_s;
		sidis.max_chi_squared_deviation = 0.2;
		sidis.max_relative_error = 1e-3;
		sidis.iter_max = 10;

		Decay decay(parametrization, Particle(), Particle(), DecayFunctions::trivial, 5.0);

		TRFKinematics kinematics = TRFKinematics::Q2_sqrt_s(x, Q2, sqrt_s, Constants::Particles::Proton.mass, 0.0);

		const double z_min = SIDISFunctions::compute_z_min(kinematics, decay);

		const PerturbativeQuantity result = sidis.lepton_pair_cross_section_xQ2(kinematics);
		std::cout << "Lepton-pair cross section value = " << result.nlo << std::endl;

		SIDIS sidis1(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			LHAInterface("EPPS21nlo_CT18Anlo_Fe56"),
			LHAInterface("kkks08_opal_d0___mas"),
			100'000,
			Process {Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino}
		);
		sidis1.global_sqrt_s = sqrt_s;
		sidis1.max_chi_squared_deviation = 0.2;
		sidis1.max_relative_error = 1e-3;
		sidis1.iter_max = 10;

		SIDIS sidis2(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			LHAInterface("EPPS21nlo_CT18Anlo_Fe56"),
			LHAInterface("kkks08_opal_d+___mas"),
			100'000,
			Process {Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino}
		);
		sidis2.global_sqrt_s = sqrt_s;
		sidis2.max_chi_squared_deviation = 0.2;
		sidis2.max_relative_error = 1e-3;
		sidis2.iter_max = 10;

		Integrator integrator([&](const double input[], [[maybe_unused]] const size_t dim, [[maybe_unused]] void *params) {
			const double z = input[0];

			const double differential_cs_1 = sidis1.differential_cross_section_xQ2(z, kinematics).nlo;
			const double differential_cs_2 = sidis2.differential_cross_section_xQ2(z, kinematics).nlo;

			const double result = 3 * Constants::Particles::D0.lifetime * differential_cs_1 + 5 * Constants::Particles::Dp.lifetime * differential_cs_2;
			return result;
		}, {z_min}, {1}, 100, nullptr, 0.5, 1e-2, 10);
		integrator.verbose = true;
		const auto result2 = integrator.integrate();
		std::cout << "Lepton-pair cross section value = " << result2 << std::endl;
		flag &= double_comparison_rel(result.nlo, result2.value, 1e-1);

		return flag;
	}

	bool sidis_lo_cross_section_integration_test() {
		bool flag = true;

		const double sqrt_s = 120.0;
		const double x = 0.3;
		const double Q2 = 10.0;
		const auto pdf = LHAInterface("EPPS21nlo_CT18Anlo_Fe56");

		DIS dis(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			pdf,
			10'000,
			Process {Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino}
		);

		dis.global_sqrt_s = sqrt_s;
		dis.max_chi_squared_deviation = 0.2;
		dis.max_relative_error = 1e-3;
		dis.iter_max = 10;

		SIDIS sidis(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			pdf,
			FunctionalFormInterface([]([[maybe_unused]] const FlavorType flavor, const double z, [[maybe_unused]] const double Q2) {
				return z;
			}),
			100,
			Process {Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino}
		);
		sidis.global_sqrt_s = sqrt_s;
		sidis.max_chi_squared_deviation = 0.2;
		sidis.max_relative_error = 1e-3;
		sidis.iter_max = 10;

		SIDIS sidis2(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			pdf,
			FunctionalFormInterface([]([[maybe_unused]] const FlavorType flavor, const double z, [[maybe_unused]] const double Q2) {
				return z;
			}),
			2,
			Process {Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino}
		);
		sidis2.global_sqrt_s = sqrt_s;
		sidis2.max_chi_squared_deviation = 0.2;
		sidis2.max_relative_error = 1e-3;
		sidis2.iter_max = 0;

		TRFKinematics kinematics = TRFKinematics::Q2_sqrt_s(x, Q2, sqrt_s, Constants::Particles::Proton.mass, 0);
		const PerturbativeQuantity result_dis = dis.differential_cross_section_xQ2(kinematics);
		const PerturbativeQuantity result_sidis = sidis.lepton_pair_cross_section_xQ2(kinematics);
		std::cout << "Lepton-pair cross section value (DIS) = " << result_dis.lo << std::endl;
		std::cout << "Lepton-pair cross section value (SIDIS integrated) = " << result_sidis.lo << std::endl;

		Integrator integrator([&](const double input[], [[maybe_unused]] const size_t dim, [[maybe_unused]] void *params) {
			const double z = input[0];
			return sidis2.differential_cross_section_xQ2(z, kinematics).lo;
		}, {0}, {1}, 100, nullptr, 0.5, 1e-2, 10);
		integrator.verbose = true;
		const auto result2 = integrator.integrate();
		std::cout << "Lepton-pair cross section value (SIDIS integrated manually) = " << result2 << std::endl;
		flag &= double_comparison_rel(result_dis.lo, result2.value, 1e-1);

		return flag;
	}

	bool run_tests() {
		bool flag = true;
		flag &= pdf_evaluation_tests();
		flag &= pdf_comparison_tests();
		flag &= ckm_tests();
		flag &= dis_cross_section_tests();
		flag &= sidis_structure_function_only_quarks_tests();
		flag &= sidis_structure_function_only_gluons_tests_1();
		flag &= sidis_structure_function_only_gluons_tests_2();
		flag &= sidis_structure_function_quarks_gluons_tests();
		flag &= sidis_structure_function_quarks_gluons_scale_logs_tests();
		flag &= decay_function_tests_1();
		flag &= decay_function_tests_2();

		return flag;
	}

	bool run_heavy_tests() {
		bool flag = true;
		flag &= lepton_pair_lo_cross_section_integration_tests();
		flag &= lepton_pair_nlo_cross_section_integration_tests();
		flag &= sidis_lo_cross_section_integration_test();

		return flag;
	}

	template <typename PDFInterface>
	void evaluate_pdf(PDFInterface pdf, const std::vector<double> Q2_list, const std::string filename, const bool divide_by_x = false) {
		std::ofstream file(filename);
		for (auto const Q2 : Q2_list) {	
			for (size_t i = 0; i < 999; i++) {
				const double x = double(i + 1) * 0.001;
				
				file << x << ", " << Q2;

				for (auto const flavor : pdf.available_flavors) {
					double pdf_value = pdf.xf_evaluate(flavor, x, Q2);
					if (divide_by_x) {
						pdf_value /= x;
					}
					file << ", " << pdf_value;
				}

				file << std::endl;
			}
		}

		file.close();
	}
}

#endif