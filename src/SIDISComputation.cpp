#ifndef SIDIS_COMPUTATION_H
#define SIDIS_COMPUTATION_H

#include <vector>
#include "Flavor.cpp"
#include "PerturbativeResult.cpp"
#include "PDFCommon.cpp"
#include "SIDISFunctions.cpp"
#include "Integrator.cpp"
#include "Decay.cpp"

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

		SIDISFunctions::Parameters<PDFInterface, FFInterface> params {
			pdf1, ff1, pdf2, ff2,
			flavors,
			Q2, nlo_coefficient, s,
			process,
			x, z, xq_zq, 0.0
		};

		Integrator xi_integrator([](double input[], size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::evaluate_gsl_sidis_integrand<PDFInterface, FFInterface>(input, dim, params_in, SIDISFunctions::Integrands::F2_xi_integrand, true, false, false);
		}, {x}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		auto xi_result = xi_integrator.integrate();
		const double xi_integral = xi_result.value;
		
		Integrator xip_integrator([](double input[], size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::evaluate_gsl_sidis_integrand<PDFInterface, FFInterface>(input, dim, params_in, SIDISFunctions::Integrands::F2_xip_integrand, false, true, false);
		}, {z}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		auto xip_result = xip_integrator.integrate();
		const double xip_integral = xip_result.value;

		Integrator xi_xip_integrator([](double input[], size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::evaluate_gsl_sidis_integrand<PDFInterface, FFInterface>(input, dim, params_in, SIDISFunctions::Integrands::F2_xi_xip_integrand, true, true, false);
		}, {x, z}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		auto xi_xip_result = xi_xip_integrator.integrate();
		const double xi_xip_integral = xi_xip_result.value;

		const double nlo1 = xq_zq * SIDISFunctions::delta_contribution(x, z);
		const double nlo2 = xi_integral + xip_integral + xi_xip_integral;

		const double nlo = nlo_coefficient * (nlo1 + nlo2);

		return PerturbativeResult {lo, lo + nlo};
	}
	PerturbativeResult FL(const double x, const double z, const double Q2) {
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		SIDISFunctions::Parameters<PDFInterface, FFInterface> params {
			pdf1, ff1, pdf2, ff2,
			flavors,
			Q2, nlo_coefficient, s,
			process,
			x, z, 0.0, 0.0
		};

		Integrator xi_xip_integrator([](double input[], size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::evaluate_gsl_sidis_integrand<PDFInterface, FFInterface>(input, dim, params_in, SIDISFunctions::Integrands::FL_xi_xip_integrand, true, true, false);
		}, {x, z}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		auto xi_xip_result = xi_xip_integrator.integrate();
		const double xi_xip_integral = xi_xip_result.value;

		const double nlo = nlo_coefficient * xi_xip_integral;

		return PerturbativeResult {0, nlo};
	}
	PerturbativeResult xF3(const double x, const double z, const double Q2) {
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		pdf1.evaluate(x, Q2);
		ff1.evaluate(z, Q2);

		const double xq_zq = PDFCommon::xq_zq_sum(pdf1, ff1, flavors, true, process);

		const double lo = 2 * xq_zq / z;

		SIDISFunctions::Parameters<PDFInterface, FFInterface> params {
			pdf1, ff1, pdf2, ff2,
			flavors,
			Q2, nlo_coefficient, s,
			process,
			x, z, xq_zq, 0.0
		};

		Integrator xi_integrator([](double input[], size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::evaluate_gsl_sidis_integrand<PDFInterface, FFInterface>(input, dim, params_in, SIDISFunctions::Integrands::F3_xi_integrand, true, false, true);
		}, {x}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto xi_result = xi_integrator.integrate();
		const double xi_integral = xi_result.value;
		
		Integrator xip_integrator([](double input[], size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::evaluate_gsl_sidis_integrand<PDFInterface, FFInterface>(input, dim, params_in, SIDISFunctions::Integrands::F3_xip_integrand, false, true, true);
		}, {z}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto xip_result = xip_integrator.integrate();
		const double xip_integral = xip_result.value;

		Integrator xi_xip_integrator([](double input[], size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::evaluate_gsl_sidis_integrand<PDFInterface, FFInterface>(input, dim, params_in, SIDISFunctions::Integrands::F3_xi_xip_integrand, true, true, true);
		}, {x, z}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto xi_xip_result = xi_xip_integrator.integrate();
		const double xi_xip_integral = xi_xip_result.value;

		const double nlo1 = xq_zq * SIDISFunctions::delta_contribution(x, z);
		const double nlo2 = xi_integral + xip_integral + xi_xip_integral;

		const double nlo = nlo_coefficient * (nlo1 + nlo2);

		return PerturbativeResult {lo, lo + nlo};
	}

	constexpr PerturbativeResult structure_function(const StructureFunction F, const double x, const double z, const double Q2) {
		switch (F) {
		case StructureFunction::F2: return F2(x, z, Q2);
		case StructureFunction::FL: return FL(x, z, Q2);
		case StructureFunction::xF3: return xF3(x, z, Q2);
		}
	}

	PerturbativeResult differential_cross_section_direct(const double x, const double z, const double Q2) {
		const double prefactor = CommonFunctions::cross_section_prefactor(Q2);

		const PerturbativeResult f2 = F2(x, z, Q2);
		const PerturbativeResult fL = FL(x, z, Q2);
		const PerturbativeResult xf3 = xF3(x, z, Q2);

		const std::optional<double> y_opt = CommonFunctions::compute_y(x, Q2, s);
		if (!y_opt.has_value()) {
			return {0, 0};
		}
		const double y = *y_opt;

		const double term1 = 1 - y + 0.5 * y * y;
		const double term2 = - 0.5 * y * y;
		const double term3 = y * (1 - 0.5 * y);

		const PerturbativeResult result = (term1 * f2 + term2 * fL + double(process.W_sign()) * term3 * xf3) / x;
		return prefactor * result;
	}
	PerturbativeResult differential_cross_section_indirect(const double x, const double z, const double Q2) {
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		pdf1.evaluate(x, Q2);
		ff1.evaluate(z, Q2);

		const double xq_zq = PDFCommon::xq_zq_sum(pdf1, ff1, flavors, false, process);
		const double xq_zq_3 = PDFCommon::xq_zq_sum(pdf1, ff1, flavors, true, process);

		const double lo_1 = 2 * xq_zq / z;
		const double lo_3 = 2 * xq_zq_3 / z;
		const double lo = CommonFunctions::make_cross_section_variable(x, Q2, s, process, lo_1, 0.0, lo_3);

		SIDISFunctions::Parameters<PDFInterface, FFInterface> params {
			pdf1, ff1, pdf2, ff2,
			flavors,
			Q2, nlo_coefficient, s,
			process,
			x, z, xq_zq, xq_zq_3
		};

		Integrator xi_integrator([](double input[], size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::evaluate_gsl_sidis_cross_section_integrand<PDFInterface, FFInterface>(input, dim, params_in, 
				SIDISFunctions::Integrands::F2_xi_integrand, SIDISFunctions::Integrands::FL_xi_integrand, SIDISFunctions::Integrands::F3_xi_integrand,
				true, false, false
			);
		}, {x}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto xi_result = xi_integrator.integrate();
		const double xi_integral = xi_result.value;
		
		Integrator xip_integrator([](double input[], size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::evaluate_gsl_sidis_cross_section_integrand<PDFInterface, FFInterface>(input, dim, params_in, 
				SIDISFunctions::Integrands::F2_xip_integrand, SIDISFunctions::Integrands::FL_xip_integrand, SIDISFunctions::Integrands::F3_xip_integrand,
				false, true, false
			);
		}, {z}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto xip_result = xip_integrator.integrate();
		const double xip_integral = xip_result.value;

		Integrator xi_xip_integrator([](double input[], size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::evaluate_gsl_sidis_cross_section_integrand<PDFInterface, FFInterface>(input, dim, params_in, 
				SIDISFunctions::Integrands::F2_xi_xip_integrand, SIDISFunctions::Integrands::FL_xi_xip_integrand, SIDISFunctions::Integrands::F3_xi_xip_integrand,
				true, true, false
			);
		}, {x, z}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto xi_xip_result = xi_xip_integrator.integrate();
		const double xi_xip_integral = xi_xip_result.value;

		const double delta_contribution = SIDISFunctions::delta_contribution(x, z);
		const double delta_1 = xq_zq * delta_contribution;
		const double delta_3 = xq_zq_3 * delta_contribution;
		const double nlo1 = CommonFunctions::make_cross_section_variable(x, Q2, s, process, delta_1, 0.0, delta_3);
		const double nlo2 = xi_integral + xip_integral + xi_xip_integral;

		const double nlo = nlo_coefficient * (nlo1 + nlo2);

		const PerturbativeResult result = PerturbativeResult {lo, lo + nlo};
		const double prefactor = CommonFunctions::cross_section_prefactor(Q2);
		return prefactor * result;
	}

	template <typename DecayFunction>
	PerturbativeResult lepton_pair_cross_section(const double x, const double Q2, const Decay<DecayFunction> decay) {
		const double z_min = decay.parametrization.z_min;

		const double alpha_s = pdf1.alpha_s(Q2);
		const double nlo_coefficient = alpha_s / (2 * M_PI);

		pdf1.evaluate(x, Q2);

		SIDISFunctions::Parameters<PDFInterface, FFInterface> params {
			pdf1, ff1, pdf2, ff2,
			flavors,
			Q2, nlo_coefficient, s,
			process,
			x, 0.0, 0.0, 0.0
		};

		Integrator lo_integrator([&](double input[], size_t dim, void* params_in) {
			return SIDISFunctions::Evaluation::evaluate_gsl_sidis_decay_cross_section_integrand<PDFInterface, FFInterface>(input, dim, params_in, 
				SIDISFunctions::Integrands::F2_lo_integrand, SIDISFunctions::Integrands::FL_lo_integrand, SIDISFunctions::Integrands::F3_lo_integrand,
				decay,
				false, false, true
			);
		}, {z_min}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto lo_result = lo_integrator.integrate();
		const double lo = lo_result.value;

		Integrator xi_integrator([&](double input[], size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::evaluate_gsl_sidis_decay_cross_section_integrand<PDFInterface, FFInterface>(input, dim, params_in, 
				SIDISFunctions::Integrands::F2_xi_integrand, SIDISFunctions::Integrands::FL_xi_integrand, SIDISFunctions::Integrands::F3_xi_integrand,
				decay,
				true, false, true
			);
		}, {x, z_min}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto xi_result = xi_integrator.integrate();
		const double xi_integral = xi_result.value;
		
		Integrator xip_integrator([&](double input[], size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::evaluate_gsl_sidis_decay_cross_section_integrand<PDFInterface, FFInterface>(input, dim, params_in, 
				SIDISFunctions::Integrands::F2_xip_integrand, SIDISFunctions::Integrands::FL_xip_integrand, SIDISFunctions::Integrands::F3_xip_integrand,
				decay,
				false, true, true
			);
		}, {z_min, z_min}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto xip_result = xip_integrator.integrate();
		const double xip_integral = xip_result.value;

		Integrator xi_xip_integrator([&](double input[], size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::evaluate_gsl_sidis_decay_cross_section_integrand<PDFInterface, FFInterface>(input, dim, params_in, 
				SIDISFunctions::Integrands::F2_xi_xip_integrand, SIDISFunctions::Integrands::FL_xi_xip_integrand, SIDISFunctions::Integrands::F3_xi_xip_integrand,
				decay,
				true, true, true
			);
		}, {x, z_min, z_min}, {1, 1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto xi_xip_result = xi_xip_integrator.integrate();
		const double xi_xip_integral = xi_xip_result.value;

		// Could combine with z_integrator, but then lose separation between LO and NLO
		Integrator delta_integrator([&](double input[], size_t dim, void* params_in) {
			return SIDISFunctions::Evaluation::evaluate_gsl_sidis_decay_cross_section_integrand<PDFInterface, FFInterface>(input, dim, params_in, 
				SIDISFunctions::Integrands::F2_delta_integrand, SIDISFunctions::Integrands::FL_delta_integrand, SIDISFunctions::Integrands::F3_delta_integrand,
				decay,
				false, false, true
			);
		}, {z_min}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto delta_result = delta_integrator.integrate();
		const double delta_integral = delta_result.value;

		const double nlo = nlo_coefficient * (xi_integral + xip_integral + xi_xip_integral + delta_integral);

		const PerturbativeResult result = PerturbativeResult {lo, lo + nlo};
		const double prefactor = CommonFunctions::cross_section_prefactor(Q2);
		return prefactor * result;
	}
};

#endif