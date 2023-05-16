#ifndef SIDIS_COMPUTATION_H
#define SIDIS_COMPUTATION_H

#include <vector>
#include "Flavor.cpp"
#include "PerturbativeResult.cpp"
#include "PDFCommon.cpp"
#include "SIDISFunctions.cpp"
#include "Integrator.cpp"

template <typename PDFInterface, typename FFInterface>
class SIDISComputation {
	public:
	double sqrt_s;
	double s;
	
	FlavorInfo flavors;

	PDFInterface pdf1;
	FFInterface ff1;

	PDFInterface pdf2;
	FFInterface ff2;

	const size_t points;
	const double max_chi_squared_deviation;
	const double max_relative_error;
	const unsigned int iter_max;

	const Process process;

	SIDISComputation (
		const double _sqrt_s, 
		const std::vector<FlavorType> _active_flavors, 
		const PDFInterface _pdf,
		const FFInterface _ff,
		const size_t _points,
		const double _max_chi_squared_deviation, 
		const double _max_relative_error,
		const unsigned int _iter_max,
		const Process _process
	) : sqrt_s(_sqrt_s), 
	s(_sqrt_s * _sqrt_s), 
	flavors(_active_flavors),
	pdf1(_pdf), 
	ff1(_ff),
	pdf2(_pdf), 
	ff2(_ff),
	points(_points), 
	max_chi_squared_deviation(_max_chi_squared_deviation), 
	max_relative_error(_max_relative_error),
	iter_max(_iter_max),
	process(_process) { }

	PerturbativeResult F2(const double x, const double z, const double Q2) {
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		pdf1.evaluate(x, Q2);
		ff1.evaluate(z, Q2);

		const double xq_zq = PDFCommon::xq_zq_sum(pdf1, ff1, flavors, false, process);

		const double lo = 2 * xq_zq / z;

		SIDISFunctions::CommonParams<PDFInterface, FFInterface> common {
			pdf1, ff1, pdf2, ff2,
			flavors,
			Q2, nlo_coefficient, s,
			process
		};
		SIDISFunctions::UnintegratedParams<PDFInterface, FFInterface> params {common, x, z, xq_zq};

		Integrator xi_integrator(&SIDISFunctions::F2_xi_integrand_gsl<PDFInterface, FFInterface>, {x}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		Integrator::Result xi_result = xi_integrator.integrate();
		const double xi_integral = xi_result.value;
		
		Integrator xip_integrator(&SIDISFunctions::F2_xip_integrand_gsl<PDFInterface, FFInterface>, {z}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		Integrator::Result xip_result = xip_integrator.integrate();
		const double xip_integral = xip_result.value;

		Integrator xi_xip_integrator(&SIDISFunctions::F2_xi_xip_integrand_gsl<PDFInterface, FFInterface>, {x, z}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		Integrator::Result xi_xip_result = xi_xip_integrator.integrate();
		const double xi_xip_integral = xi_xip_result.value;

		const double nlo1 = xq_zq * SIDISFunctions::delta_contribution(x, z);
		const double nlo2 = xi_integral + xip_integral + xi_xip_integral;

		const double nlo = 2 * nlo_coefficient * (nlo1 + nlo2) / z;

		return PerturbativeResult {lo, lo + nlo};
	}
	PerturbativeResult FL(const double x, const double z, const double Q2) {
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		SIDISFunctions::CommonParams<PDFInterface, FFInterface> common {
			pdf1, ff1, pdf2, ff2,
			flavors,
			Q2, nlo_coefficient, s,
			process
		};
		SIDISFunctions::UnintegratedParams<PDFInterface, FFInterface> params {common, x, z, 0};

		Integrator xi_xip_integrator(&SIDISFunctions::FL_xi_xip_integrand_gsl<PDFInterface, FFInterface>, {x, z}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		Integrator::Result xi_xip_result = xi_xip_integrator.integrate();
		const double xi_xip_integral = xi_xip_result.value;

		const double nlo = 2 * nlo_coefficient * xi_xip_integral / z;

		return PerturbativeResult {0, nlo};
	}
	PerturbativeResult xF3(const double x, const double z, const double Q2) {
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		pdf1.evaluate(x, Q2);
		ff1.evaluate(z, Q2);

		const double xq_zq = PDFCommon::xq_zq_sum(pdf1, ff1, flavors, true, process);

		const double lo = 2 * xq_zq / z;

		SIDISFunctions::CommonParams<PDFInterface, FFInterface> common {
			pdf1, ff1, pdf2, ff2,
			flavors,
			Q2, nlo_coefficient, s,
			process
		};
		SIDISFunctions::UnintegratedParams<PDFInterface, FFInterface> params {common, x, z, xq_zq};

		Integrator xi_integrator(&SIDISFunctions::F3_xi_integrand_gsl<PDFInterface, FFInterface>, {x}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		Integrator::Result xi_result = xi_integrator.integrate();
		const double xi_integral = xi_result.value;
		
		Integrator xip_integrator(&SIDISFunctions::F3_xip_integrand_gsl<PDFInterface, FFInterface>, {z}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		Integrator::Result xip_result = xip_integrator.integrate();
		const double xip_integral = xip_result.value;

		Integrator xi_xip_integrator(&SIDISFunctions::F3_xi_xip_integrand_gsl<PDFInterface, FFInterface>, {x, z}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		Integrator::Result xi_xip_result = xi_xip_integrator.integrate();
		const double xi_xip_integral = xi_xip_result.value;

		const double nlo1 = xq_zq * SIDISFunctions::delta_contribution(x, z);
		const double nlo2 = xi_integral + xip_integral + xi_xip_integral;

		const double nlo = 2 * nlo_coefficient * (nlo1 + nlo2) / z;

		return PerturbativeResult {lo, lo + nlo};
	}

	constexpr PerturbativeResult structure_function(const StructureFunction F, const double x, const double z, const double Q2) {
		switch (F) {
		case StructureFunction::F2: return F2(x, z, Q2);
		case StructureFunction::FL: return FL(x, z, Q2);
		case StructureFunction::xF3: return xF3(x, z, Q2);
		}
	}

	PerturbativeResult differential_cross_section(const double x, const double z, const double Q2) {
		const double prefactor = CommonFunctions::cross_section_prefactor(Q2);

		const PerturbativeResult f2 = F2(x, z, Q2);
		const PerturbativeResult fL = FL(x, z, Q2);
		const PerturbativeResult xf3 = xF3(x, z, Q2);
		// const PerturbativeResult xf3 = PerturbativeResult{0.0, 0.0};

		const std::optional<double> y_opt = CommonFunctions::compute_y(x, Q2, s);
		if (!y_opt.has_value()) {
			return {0, 0};
		}
		const double y = *y_opt;

		const double term1 = 1 - y + 0.5 * y * y;
		const double term2 = - 0.5 * y * y;
		const double term3 = y * (1 - 0.5 * y);

		const PerturbativeResult result = prefactor * (term1 * f2 + term2 * fL + double(process.W_sign()) * term3 * xf3) / x;
		return result;
	}

};

#endif