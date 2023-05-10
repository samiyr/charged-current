#ifndef DIS_COMPUTATION_H
#define DIS_COMPUTATION_H

#include "Integrator.cpp"
#include "PDFInterface.cpp"
#include "Coefficients.cpp"
#include "Constants.cpp"
#include "Utility.cpp"
#include "PerturbativeResult.cpp"
#include "Process.cpp"
#include "StructureFunction.cpp"
#include "DISFunctions.cpp"
#include "CommonFunctions.cpp"
#include "PDFCommon.cpp"

// double F2_integrand_gsl(double input[], size_t dim, void *params_in);
// double FL_integrand_gsl(double input[], size_t dim, void *params_in);
// double F3_integrand_gsl(double input[], size_t dim, void *params_in);
// double x_integrated_cross_section_gsl_single_integrand(double input[], size_t dim, void *params_in);
// double x_integrated_cross_section_gsl_double_integrand(double input[], size_t dim, void *params_in);

// constexpr double delta_contribution(const double x);

template <typename PDFInterface>
class DISComputation {
	public:
	double sqrt_s;
	double s;
	double y_max = 1.0;

	FlavorInfo flavors;

	PDFInterface pdf;

	const size_t points;
	const double max_chi_squared_deviation;
	const double max_relative_error;
	const unsigned int iter_max;

	const Process process;

	DISComputation (
		const double _sqrt_s, 
		const FlavorVector _active_flavors, 
		const PDFInterface _pdf, 
		const size_t _points,
		const double _max_chi_squared_deviation, 
		const double _max_relative_error,
		const unsigned int _iter_max,
		const Process _process
	) : sqrt_s(_sqrt_s), 
	s(_sqrt_s * _sqrt_s), 
	flavors(_active_flavors), 
	pdf(_pdf), 
	points(_points), 
	max_chi_squared_deviation(_max_chi_squared_deviation), 
	max_relative_error(_max_relative_error),
	iter_max(_iter_max),
	process(_process) { }

	PerturbativeResult F2(const double x, const double Q2) {	
		double alpha_s = pdf.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		pdf.evaluate(x, Q2);

		const double xq = PDFCommon::xq_sum(pdf, flavors, false, process);

		DISFunctions::CommonParams<PDFInterface> common {pdf, flavors, Q2, nlo_coefficient, s, y_max, process};
		DISFunctions::UnintegratedParams<PDFInterface> params {common, x, xq};

		Integrator integrator(&DISFunctions::F2_integrand_gsl<PDFInterface>, {x}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		Integrator::Result integral_result = integrator.integrate();
		const double integral = integral_result.value;

		const double lo = 2 * xq;

		const double nlo1 = xq * DISFunctions::delta_contribution(x);
		const double nlo2 = integral;
		const double nlo = 2 * nlo_coefficient * (-nlo1 + nlo2);

		return PerturbativeResult {lo, lo + nlo};
	}
	PerturbativeResult FL(const double x, const double Q2) {
		double alpha_s = pdf.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		pdf.evaluate(x, Q2);

		DISFunctions::CommonParams<PDFInterface> common {pdf, flavors, Q2, nlo_coefficient, s, y_max, process};
		DISFunctions::UnintegratedParams<PDFInterface> params {common, x, 0};

		Integrator integrator(&DISFunctions::FL_integrand_gsl<PDFInterface>, {x}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		Integrator::Result integral_result = integrator.integrate();
		const double integral = integral_result.value;

		const double nlo = 2 * nlo_coefficient * integral;

		return PerturbativeResult {0, nlo};
	}
	PerturbativeResult xF3(const double x, const double Q2) {
		double alpha_s = pdf.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		pdf.evaluate(x, Q2);

		const double xq = PDFCommon::xq_sum(pdf, flavors, true, process);

		DISFunctions::CommonParams<PDFInterface> common {pdf, flavors, Q2, nlo_coefficient, s, y_max, process};
		DISFunctions::UnintegratedParams<PDFInterface> params {common, x, xq};

		Integrator integrator(&DISFunctions::F3_integrand_gsl<PDFInterface>, {x}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		Integrator::Result integral_result = integrator.integrate();
		const double integral = integral_result.value;

		const double lo = 2 * xq;

		const double nlo1 = xq * DISFunctions::delta_contribution(x);
		const double nlo2 = integral;
		const double nlo = 2 * nlo_coefficient * (-nlo1 + nlo2);

		return PerturbativeResult {lo, lo + nlo};
	}
	constexpr PerturbativeResult structure_function(StructureFunction F, double x, double Q2) {
		switch (F) {
		case StructureFunction::F2: return F2(x, Q2);
		case StructureFunction::FL: return FL(x, Q2);
		case StructureFunction::xF3: return xF3(x, Q2);
		}
	}

	PerturbativeResult differential_cross_section(const double x, const double Q2) {
		const double prefactor = CommonFunctions::cross_section_prefactor(Q2);

		const PerturbativeResult f2 = F2(x, Q2);
		const PerturbativeResult fL = FL(x, Q2);
		const PerturbativeResult xf3 = xF3(x, Q2);

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

// template <typename F>
// double evaluate_gsl_x_single_integrand(double input[], size_t dim, void *params_in, F integrand) {
// 	const struct DISComputation::UnintegratedParams *params = (struct DISComputation::UnintegratedParams *)params_in;
// 	const double x = input[0];
// 	const DISComputation::CommonParams &common = params->common;
// 	const double Q2 = common.Q2;
// 	const double s = common.s;
// 	const double y_max = common.y_max;
// 	PDFInterface &pdf = common.pdf;

// 	pdf.evaluate(x, Q2);

// 	const FlavorVector &active_flavors = common.active_flavors;
// 	const FlavorVector &active_antiflavors = common.active_antiflavors;
// 	const FlavorVector &active_upper_flavors = common.active_upper_flavors;
// 	const FlavorVector &active_lower_flavors = common.active_lower_flavors;
// 	const FlavorVector &active_upper_antiflavors = common.active_upper_antiflavors;
// 	const FlavorVector &active_lower_antiflavors = common.active_lower_antiflavors;
// 	// const double xq_sum = pdf.xf_sum(active_lower_flavors, active_upper_flavors);
// 	// const double anti_xq_sum = pdf.xf_sum(active_upper_antiflavors, active_lower_antiflavors);
// 	// const double xq = xq_sum + anti_xq_sum;
// 	// const double xq3 = xq_sum - anti_xq_sum;
// 	const double xq = pdf.xq_sum(active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, false, common.process);
// 	const double xq3 = pdf.xq_sum(active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, true, common.process);

// 	const double xg = pdf.xg_sum(active_upper_flavors, active_lower_flavors);
	
// 	const double y = compute_y(x, Q2, s);
// 	if (y < 0 || y > y_max) { return 0; }
// 	const double alpha_s = pdf.alpha_s(Q2);
// 	const double nlo_coefficient = alpha_s / (2 * M_PI);

// 	return integrand(x, y, xq, xq3, nlo_coefficient);
// }

// template <typename F>
// double evaluate_gsl_x_double_integrand(double input[], size_t dim, void *params_in, F integrand) {
// 	const struct DISComputation::XIntegratedParams *params = (struct DISComputation::XIntegratedParams *)params_in;
// 	const double x = input[0];
// 	const double xi = input[1];
// 	if (xi <= x) { return 0; }
// 	const double x_hat = x / xi;
// 	const DISComputation::CommonParams &common = params->common;
// 	const double Q2 = common.Q2;
// 	const double s = common.s;
// 	const double y_max = common.y_max;

// 	const FlavorVector &active_flavors = common.active_flavors;
// 	const FlavorVector &active_antiflavors = common.active_antiflavors;
// 	const FlavorVector &active_upper_flavors = common.active_upper_flavors;
// 	const FlavorVector &active_lower_flavors = common.active_lower_flavors;
// 	const FlavorVector &active_upper_antiflavors = common.active_upper_antiflavors;
// 	const FlavorVector &active_lower_antiflavors = common.active_lower_antiflavors;

// 	PDFInterface &pdf = common.pdf;

// 	pdf.evaluate(x, Q2);

// 	// const double xq_sum = pdf.xf_sum(active_lower_flavors, active_upper_flavors);
// 	// const double anti_xq_sum = pdf.xf_sum(active_upper_antiflavors, active_lower_antiflavors);
// 	// const double xq = xq_sum + anti_xq_sum;
// 	// const double xq3 = xq_sum - anti_xq_sum;
// 	const double xq = pdf.xq_sum(active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, false, common.process);
// 	const double xq3 = pdf.xq_sum(active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, true, common.process);

// 	const double xg = pdf.xg_sum(active_upper_flavors, active_lower_flavors);

// 	pdf.evaluate(x_hat, Q2);

// 	// const double xq_hat_sum = pdf.xf_sum(active_lower_flavors, active_upper_flavors);
// 	// const double anti_xq_hat_sum = pdf.xf_sum(active_upper_antiflavors, active_lower_antiflavors);
// 	// const double xq_hat = xq_hat_sum + anti_xq_hat_sum;
// 	// const double xq3_hat = xq_hat_sum - anti_xq_hat_sum;
// 	const double xq_hat = pdf.xq_sum(active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, false, common.process);
// 	const double xq3_hat = pdf.xq_sum(active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, true, common.process);

// 	const double xg_hat = pdf.xg_sum(active_upper_flavors, active_lower_flavors);

// 	const double y = compute_y(x, Q2, s);
// 	if (y < 0 || y > y_max) { return 0; }
// 	const double alpha_s = pdf.alpha_s(Q2);
// 	const double nlo_coefficient = alpha_s / (2 * M_PI);

// 	return integrand(x, xi, y, xq, xq_hat, xq3, xq3_hat, xg_hat, nlo_coefficient);
// }

// constexpr double x_integrated_cross_section_single_integrand(const double x, const double y, const double xq, const double xq3, const double nlo_coefficient) {
// 	const double multiplier = 1 - nlo_coefficient * delta_contribution(x);

// 	const double term1 = (1 + std::pow(1 - y, 2)) * xq;
// 	const double term2 = y * (2 - y) * xq3;

// 	return multiplier * (term1 + term2) / x;
// }
// constexpr double x_integrated_cross_section_double_integrand(
// 	const double x, 
// 	const double xi, 
// 	const double y, 
// 	const double xq, 
// 	const double xq_hat, 
// 	const double xq3, 
// 	const double xq3_hat, 
// 	const double xg_hat, 
// 	const double nlo_coefficient) {
// 		const double f2 = F2_integrand(xi, xq_hat, xq, xg_hat);
// 		const double fL = FL_integrand(xi, xq_hat, xq, xg_hat);
// 		const double f3 = F3_integrand(xi, xq3_hat, xq3, xg_hat);

// 		const double term2 = 1 + std::pow(1 - y, 2);
// 		const double termL = - y * y;
// 		const double term3 = y * (2 - y);

// 		return (term2 * f2 + termL * fL + term3 * f3) / x;
// }

// double x_integrated_cross_section_gsl_single_integrand(double input[], size_t dim, void *params_in) {
// 	return evaluate_gsl_x_single_integrand(input, dim, params_in, x_integrated_cross_section_single_integrand);
// }
// double x_integrated_cross_section_gsl_double_integrand(double input[], size_t dim, void *params_in) {
// 	return evaluate_gsl_x_double_integrand(input, dim, params_in, x_integrated_cross_section_double_integrand);
// }

#endif