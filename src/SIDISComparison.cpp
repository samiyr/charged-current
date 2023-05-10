#ifndef SIDIS_COMPARISON_H
#define SIDIS_COMPARISON_H

#include "SIDIS.cpp"
#include "Integrator.cpp"
#include "LHAInterface.cpp"
#include "DeltaInterface.cpp"
#include "TrivialInterface.cpp"
#include "FunctionalFormInterface.cpp"

namespace SIDISComparison {
	template <typename PDFInterface, typename FFInterface>
	struct Params {
		SIDIS<PDFInterface, FFInterface> &sidis;
		const double x;
		const double Q2;
	};
	template <typename PDFInterface, typename FFInterface>
	static double F2_integrand_gsl(double input[], size_t dim, void *params_in) {
		const double z = input[0];
		const struct Params<PDFInterface, FFInterface> *params = (struct Params<PDFInterface, FFInterface> *)params_in;
		SIDIS<PDFInterface, FFInterface> &sidis = params->sidis;
		const double x = params->x;
		const double Q2 = params->Q2;

		const double value = sidis.F2(x, z, Q2).nlo;
		return value;
	}

	// static double F2_z_integral(const double x, const double Q2, const double z_min) {
	// 	SIDIS sidis(
	// 		318,
	// 		{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
	// 		// LHAInterface("CT18ANLO"),
	// 		// LHAInterface("JAM20-SIDIS_FF_pion_nlo"),
			// FunctionalFormInterface([](const FlavorType flavor, const double x, const double Q2) {
			// 	return 2 * x * (1 - x);
			// }),
			// FunctionalFormInterface([](const FlavorType flavor, const double x, const double Q2) {
			// 	return x * (1 - x);
			// }),
	// 		// TrivialInterface(),
	// 		// DeltaInterface(1, 0.0001),
	// 		2'000'000,
	// 		Process {Process::Type::NeutrinoToLepton}
	// 	);
	// 	sidis.max_chi_squared_deviation = 0.2;
	// 	sidis.iter_max = 20;

	// 	return sidis.F2_z_integrated(x, Q2, z_min).nlo;
	// }
	// static double F2_z_integral(const double x, const double Q2, const double z_min) {
	// 	SIDIS sidis(
	// 		0,
	// 		{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
	// 		// LHAInterface("CT18ANLO"),
	// 		// TrivialInterface(),
	// 		// TrivialInterface(),
	// 		// DeltaInterface(1, 0.001),
	// 		FunctionalFormInterface([](const FlavorType flavor, const double x, const double Q2) {
	// 			return x * (1 - x);
	// 		}),
	// 		FunctionalFormInterface([](const FlavorType flavor, const double x, const double Q2) {
	// 			return x;
	// 		}),
	// 		// FunctionalFormInterface([](const FlavorType flavor, const double x, const double Q2) {
	// 		// 	return x * (1 - x);
	// 		// }),
	// 		100'000,
	// 		Process {Process::Type::NeutrinoToLepton}
	// 	);
	// 	sidis.max_chi_squared_deviation = 0.2;
	// 	sidis.iter_max = 10;

	// 	Params<FunctionalFormInterface, FunctionalFormInterface> params {sidis, x, Q2};
	// 	Integrator integrator(&F2_integrand_gsl<FunctionalFormInterface, FunctionalFormInterface>, {z_min}, {1}, 200, &params, 0.5, 1e-2, 10);
	// 	integrator.verbose = true;
	// 	Integrator::Result result = integrator.integrate();

	// 	return result.value;
	// }
}

#endif