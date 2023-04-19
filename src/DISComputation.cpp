#ifndef DIS_COMPUTATION_H
#define DIS_COMPUTATION_H

#include "Integrator.cpp"
#include "PDFEvaluation.cpp"
#include "Coefficients.cpp"
#include "Constants.cpp"
#include "Utility.cpp"
#include "PerturbativeResult.cpp"
#include "Process.cpp"

double F2_integrand_gsl(double input[], size_t dim, void *params_in);
double FL_integrand_gsl(double input[], size_t dim, void *params_in);
double F3_integrand_gsl(double input[], size_t dim, void *params_in);
double x_integrated_cross_section_gsl_single_integrand(double input[], size_t dim, void *params_in);
double x_integrated_cross_section_gsl_double_integrand(double input[], size_t dim, void *params_in);

constexpr double compute_y(const double x, const double Q2, const double s);
constexpr double delta_contribution(const double x);

class DISComputation {
	public:

	PDFEvaluation pdf;

	std::vector<FlavorType> active_flavors;
	std::vector<FlavorType> active_antiflavors;

	double sqrt_s;
	double s;
	double y_max = 1.0;

	const size_t points;
	const double max_chi_squared_deviation;
	const double max_relative_error;
	const unsigned int iter_max;

	const Process process;

	DISComputation (
		const double _sqrt_s, 
		const std::vector<FlavorType> _active_flavors, 
		const std::vector<FlavorType> _active_antiflavors, 
		const std::string _pdf_set_name, 
		const size_t _points,
		const double _max_chi_squared_deviation, 
		const double _max_relative_error,
		const unsigned int _iter_max,
		const Process _process
	) : sqrt_s(_sqrt_s), 
	s(_sqrt_s * _sqrt_s), 
	active_flavors(_active_flavors), 
	active_antiflavors(_active_antiflavors), 
	pdf(_pdf_set_name), 
	points(_points), 
	max_chi_squared_deviation(_max_chi_squared_deviation), 
	max_relative_error(_max_relative_error),
	iter_max(_iter_max),
	process(_process) { }

	struct CommonParams {
		PDFEvaluation &pdf;

		const std::vector<FlavorType> &active_flavors;
		const std::vector<FlavorType> &active_antiflavors;

		const std::vector<FlavorType> &active_upper_flavors;
		const std::vector<FlavorType> &active_lower_flavors;
		const std::vector<FlavorType> &active_upper_antiflavors;
		const std::vector<FlavorType> &active_lower_antiflavors;

		const double Q2;
		const double nlo_coefficient;
		const double s;
		const double y_max;

		const Process &process;
	};

	struct UnintegratedParams {
		const CommonParams &common;

		const double unintegrated_parameter_value;
		const double xq;
	};

	struct XIntegratedParams {
		const CommonParams &common;
	};

	enum class StructureFunction {
		F2, FL, xF3
	};

	double x_integrated_cross_section(const double Q2, const double x_min) {
		const double prefactor = cross_section_prefactor(Q2);
		const std::vector<FlavorType> active_upper_flavors = upper_flavors(active_flavors);
		const std::vector<FlavorType> active_lower_flavors = lower_flavors(active_flavors);
		const std::vector<FlavorType> active_upper_antiflavors = upper_flavors(active_antiflavors);
		const std::vector<FlavorType> active_lower_antiflavors = lower_flavors(active_antiflavors);
	
		double alpha_s = pdf.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		CommonParams common {pdf, active_flavors, active_antiflavors, active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, Q2, nlo_coefficient, s, y_max, process};

		UnintegratedParams unintegrated_params {common, 0, 0};
		Integrator integrator1(&x_integrated_cross_section_gsl_single_integrand, {x_min}, {1}, points, &unintegrated_params, max_chi_squared_deviation, max_relative_error, iter_max);
		Integrator::Result integral_result1 = integrator1.integrate();
		const double integral1 = integral_result1.value;

		XIntegratedParams x_integrated_params {common};
		Integrator integrator2(&x_integrated_cross_section_gsl_double_integrand, {x_min, 0}, {1, 1}, points, &x_integrated_params, max_chi_squared_deviation, max_relative_error, iter_max);
		Integrator::Result integral_result2 = integrator2.integrate();
		const double integral2 = integral_result2.value;

		return prefactor * (integral1 + nlo_coefficient * integral2);
	}

	PerturbativeResult F2(const double x, const double Q2) {
		const std::vector<FlavorType> active_upper_flavors = upper_flavors(active_flavors);
		const std::vector<FlavorType> active_lower_flavors = lower_flavors(active_flavors);
		const std::vector<FlavorType> active_upper_antiflavors = upper_flavors(active_antiflavors);
		const std::vector<FlavorType> active_lower_antiflavors = lower_flavors(active_antiflavors);
	
		double alpha_s = pdf.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		pdf.evaluate(x, Q2);

		// const double xq_sum = pdf.xf_sum(active_lower_flavors, active_upper_flavors);
		// const double anti_xq_sum = pdf.xf_sum(active_upper_antiflavors, active_lower_antiflavors);
		// const double xq = xq_sum + anti_xq_sum;
		const double xq = pdf.xq_sum(active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, false, process);

		CommonParams common {pdf, active_flavors, active_antiflavors, active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, Q2, nlo_coefficient, s, y_max, process};
		UnintegratedParams params {common, x, xq};

		Integrator integrator(&F2_integrand_gsl, {x}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		Integrator::Result integral_result = integrator.integrate();
		const double integral = integral_result.value;

		const double lo = 2 * xq;

		const double nlo1 = xq * delta_contribution(x);
		const double nlo2 = integral;
		const double nlo = 2 * nlo_coefficient * (-nlo1 + nlo2);

		return PerturbativeResult {lo, lo + nlo};
	}
	PerturbativeResult FL(const double x, const double Q2) {
		const std::vector<FlavorType> active_upper_flavors = upper_flavors(active_flavors);
		const std::vector<FlavorType> active_lower_flavors = lower_flavors(active_flavors);
		const std::vector<FlavorType> active_upper_antiflavors = upper_flavors(active_antiflavors);
		const std::vector<FlavorType> active_lower_antiflavors = lower_flavors(active_antiflavors);
	
		double alpha_s = pdf.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		pdf.evaluate(x, Q2);

		CommonParams common {pdf, active_flavors, active_antiflavors, active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, Q2, nlo_coefficient, s, y_max, process};
		UnintegratedParams params {common, x, 0};

		Integrator integrator(&FL_integrand_gsl, {x}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		Integrator::Result integral_result = integrator.integrate();
		const double integral = integral_result.value;

		const double nlo = 2 * nlo_coefficient * integral;

		return PerturbativeResult {0, nlo};
	}
	PerturbativeResult xF3(const double x, const double Q2) {
		const std::vector<FlavorType> active_upper_flavors = upper_flavors(active_flavors);
		const std::vector<FlavorType> active_lower_flavors = lower_flavors(active_flavors);
		const std::vector<FlavorType> active_upper_antiflavors = upper_flavors(active_antiflavors);
		const std::vector<FlavorType> active_lower_antiflavors = lower_flavors(active_antiflavors);
	
		double alpha_s = pdf.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		pdf.evaluate(x, Q2);

		// const double xq_sum = pdf.xf_sum(active_lower_flavors, active_upper_flavors);
		// const double anti_xq_sum = pdf.xf_sum(active_upper_antiflavors, active_lower_antiflavors);
		// const double xq = xq_sum - anti_xq_sum;
		const double xq = pdf.xq_sum(active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, true, process);

		CommonParams common {pdf, active_flavors, active_antiflavors, active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, Q2, nlo_coefficient, s, y_max, process};
		UnintegratedParams params {common, x, xq};

		Integrator integrator(&F3_integrand_gsl, {x}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		Integrator::Result integral_result = integrator.integrate();
		const double integral = integral_result.value;

		const double lo = 2 * xq;

		const double nlo1 = xq * delta_contribution(x);
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
		const double prefactor = cross_section_prefactor(Q2);

		const PerturbativeResult f2 = F2(x, Q2);
		const PerturbativeResult fL = FL(x, Q2);
		const PerturbativeResult xf3 = xF3(x, Q2);

		const double y = compute_y(x, Q2, s);

		const double term1 = 1 - y + 0.5 * y * y;
		const double term2 = - 0.5 * y * y;
		const double term3 = y * (1 - 0.5 * y);

		const PerturbativeResult result = prefactor * (term1 * f2 + term2 * fL + double(process.W_sign()) * term3 * xf3) / x;
		return result;
	}
	private:

	constexpr double cross_section_prefactor(const double Q2) const {
		const double numerator = POW2(Constants::fermi_coupling) * POW4(Constants::mass_W);
		const double denominator = 2 * M_PI * POW2(Q2 + POW2(Constants::mass_W));

		return numerator / denominator;
	}
};

constexpr double delta_contribution(const double x) {
	const double log = std::log(1 - x);
	return Constants::C_F * (9.0 / 2.0 + (M_PI * M_PI) / 3.0 - log * log + 1.5 * log);
}

constexpr double compute_y(const double x, const double Q2, const double s) {
	return Q2 / (x * s);
}

template <typename F>
double evaluate_gsl_xi_integrand(double input[], size_t dim, void *params_in, F integrand, bool quark_minus = false) {
	const struct DISComputation::UnintegratedParams *params = (struct DISComputation::UnintegratedParams *)params_in;
	const double xi = input[0];
	if (std::abs(xi - 1) < 1e-15) { return 0; }
	const double x = params->unintegrated_parameter_value;
	const DISComputation::CommonParams &common = params->common;
	const double Q2 = common.Q2;
	const double x_hat = x / xi;
	const double xq = params->xq;
	PDFEvaluation &pdf = common.pdf;

	pdf.evaluate(x_hat, Q2);

	const std::vector<FlavorType> &active_flavors = common.active_flavors;
	const std::vector<FlavorType> &active_antiflavors = common.active_antiflavors;
	const std::vector<FlavorType> &active_upper_flavors = common.active_upper_flavors;
	const std::vector<FlavorType> &active_lower_flavors = common.active_lower_flavors;
	const std::vector<FlavorType> &active_upper_antiflavors = common.active_upper_antiflavors;
	const std::vector<FlavorType> &active_lower_antiflavors = common.active_lower_antiflavors;
	// const double xq_hat_sum = pdf.xf_sum(active_lower_flavors, active_upper_flavors);
	// const double anti_xq_hat_sum = pdf.xf_sum(active_upper_antiflavors, active_lower_antiflavors);
	// const double xq_hat = quark_minus ? xq_hat_sum - anti_xq_hat_sum : xq_hat_sum + anti_xq_hat_sum;

	const double xq_hat = pdf.xq_sum(active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, quark_minus, common.process);
	const double xg_hat = pdf.xg_sum(active_upper_flavors, active_lower_flavors);
	
	const double integrand_value = integrand(xi, xq_hat, xq, xg_hat);
	return integrand_value;
}

template <typename F>
double evaluate_gsl_x_single_integrand(double input[], size_t dim, void *params_in, F integrand) {
	const struct DISComputation::UnintegratedParams *params = (struct DISComputation::UnintegratedParams *)params_in;
	const double x = input[0];
	const DISComputation::CommonParams &common = params->common;
	const double Q2 = common.Q2;
	const double s = common.s;
	const double y_max = common.y_max;
	PDFEvaluation &pdf = common.pdf;

	pdf.evaluate(x, Q2);

	const std::vector<FlavorType> &active_flavors = common.active_flavors;
	const std::vector<FlavorType> &active_antiflavors = common.active_antiflavors;
	const std::vector<FlavorType> &active_upper_flavors = common.active_upper_flavors;
	const std::vector<FlavorType> &active_lower_flavors = common.active_lower_flavors;
	const std::vector<FlavorType> &active_upper_antiflavors = common.active_upper_antiflavors;
	const std::vector<FlavorType> &active_lower_antiflavors = common.active_lower_antiflavors;
	// const double xq_sum = pdf.xf_sum(active_lower_flavors, active_upper_flavors);
	// const double anti_xq_sum = pdf.xf_sum(active_upper_antiflavors, active_lower_antiflavors);
	// const double xq = xq_sum + anti_xq_sum;
	// const double xq3 = xq_sum - anti_xq_sum;
	const double xq = pdf.xq_sum(active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, false, common.process);
	const double xq3 = pdf.xq_sum(active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, true, common.process);

	const double xg = pdf.xg_sum(active_upper_flavors, active_lower_flavors);
	
	const double y = compute_y(x, Q2, s);
	if (y < 0 || y > y_max) { return 0; }
	const double alpha_s = pdf.alpha_s(Q2);
	const double nlo_coefficient = alpha_s / (2 * M_PI);

	return integrand(x, y, xq, xq3, nlo_coefficient);
}

template <typename F>
double evaluate_gsl_x_double_integrand(double input[], size_t dim, void *params_in, F integrand) {
	const struct DISComputation::XIntegratedParams *params = (struct DISComputation::XIntegratedParams *)params_in;
	const double x = input[0];
	const double xi = input[1];
	if (xi <= x) { return 0; }
	const double x_hat = x / xi;
	const DISComputation::CommonParams &common = params->common;
	const double Q2 = common.Q2;
	const double s = common.s;
	const double y_max = common.y_max;

	const std::vector<FlavorType> &active_flavors = common.active_flavors;
	const std::vector<FlavorType> &active_antiflavors = common.active_antiflavors;
	const std::vector<FlavorType> &active_upper_flavors = common.active_upper_flavors;
	const std::vector<FlavorType> &active_lower_flavors = common.active_lower_flavors;
	const std::vector<FlavorType> &active_upper_antiflavors = common.active_upper_antiflavors;
	const std::vector<FlavorType> &active_lower_antiflavors = common.active_lower_antiflavors;

	PDFEvaluation &pdf = common.pdf;

	pdf.evaluate(x, Q2);

	// const double xq_sum = pdf.xf_sum(active_lower_flavors, active_upper_flavors);
	// const double anti_xq_sum = pdf.xf_sum(active_upper_antiflavors, active_lower_antiflavors);
	// const double xq = xq_sum + anti_xq_sum;
	// const double xq3 = xq_sum - anti_xq_sum;
	const double xq = pdf.xq_sum(active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, false, common.process);
	const double xq3 = pdf.xq_sum(active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, true, common.process);

	const double xg = pdf.xg_sum(active_upper_flavors, active_lower_flavors);

	pdf.evaluate(x_hat, Q2);

	// const double xq_hat_sum = pdf.xf_sum(active_lower_flavors, active_upper_flavors);
	// const double anti_xq_hat_sum = pdf.xf_sum(active_upper_antiflavors, active_lower_antiflavors);
	// const double xq_hat = xq_hat_sum + anti_xq_hat_sum;
	// const double xq3_hat = xq_hat_sum - anti_xq_hat_sum;
	const double xq_hat = pdf.xq_sum(active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, false, common.process);
	const double xq3_hat = pdf.xq_sum(active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, true, common.process);

	const double xg_hat = pdf.xg_sum(active_upper_flavors, active_lower_flavors);

	const double y = compute_y(x, Q2, s);
	if (y < 0 || y > y_max) { return 0; }
	const double alpha_s = pdf.alpha_s(Q2);
	const double nlo_coefficient = alpha_s / (2 * M_PI);

	return integrand(x, xi, y, xq, xq_hat, xq3, xq3_hat, xg_hat, nlo_coefficient);
}

constexpr double F2_integrand(const double xi, const double xq_hat, const double xq, const double xg_hat) {
	const double term1 = std::log(1 - xi) * ((1 + xi * xi) * xq_hat - 2 * xq) / (1 - xi);
	const double term2 = (xq_hat - xq) / (1 - xi);
	const double term3 = xq_hat * (- (1 + xi * xi) * std::log(xi) / (1 - xi) + 3 + 2 * xi);

	const double quark_contribution = Constants::C_F * (term1 - 1.5 * term2 + term3);

	const double term4 = (xi * xi + (1 - xi) * (1 - xi)) * std::log((1 - xi) / xi);
	const double term5 = -1 + 8 * xi * (1 - xi);

	const double gluon_contribution = 0.5 * xg_hat * (term4 + term5);

	return quark_contribution + gluon_contribution;
}
double F2_integrand_gsl(double input[], size_t dim, void *params_in) {
	return evaluate_gsl_xi_integrand(input, dim, params_in, F2_integrand);
}

constexpr double FL_integrand(const double xi, const double xq_hat, const double xq, const double xg_hat) {
	const double quark_contribution = Constants::C_F * xq_hat * 2 * xi;
	const double gluon_contribution = 2 * xi * (1 - xi) * xg_hat;

	return quark_contribution + gluon_contribution;
}
double FL_integrand_gsl(double input[], size_t dim, void *params_in) {
	return evaluate_gsl_xi_integrand(input, dim, params_in, FL_integrand);
}

constexpr double F3_integrand(const double xi, const double xq_hat, const double xq, const double xg_hat) {
	const double term1 = std::log(1 - xi) * ((1 + xi * xi) * xq_hat - 2 * xq) / (1 - xi);
	const double term2 = (xq_hat - xq) / (1 - xi);
	const double term3 = xq_hat * (- (1 + xi * xi) * std::log(xi) / (1 - xi) + 2 + xi);

	const double quark_contribution = Constants::C_F * (term1 - 1.5 * term2 + term3);

	return quark_contribution;
}
double F3_integrand_gsl(double input[], size_t dim, void *params_in) {
	return evaluate_gsl_xi_integrand(input, dim, params_in, F3_integrand, true);
}

constexpr double x_integrated_cross_section_single_integrand(const double x, const double y, const double xq, const double xq3, const double nlo_coefficient) {
	const double multiplier = 1 - nlo_coefficient * delta_contribution(x);

	const double term1 = (1 + std::pow(1 - y, 2)) * xq;
	const double term2 = y * (2 - y) * xq3;

	return multiplier * (term1 + term2) / x;
}
constexpr double x_integrated_cross_section_double_integrand(
	const double x, 
	const double xi, 
	const double y, 
	const double xq, 
	const double xq_hat, 
	const double xq3, 
	const double xq3_hat, 
	const double xg_hat, 
	const double nlo_coefficient) {
		const double f2 = F2_integrand(xi, xq_hat, xq, xg_hat);
		const double fL = FL_integrand(xi, xq_hat, xq, xg_hat);
		const double f3 = F3_integrand(xi, xq3_hat, xq3, xg_hat);

		const double term2 = 1 + std::pow(1 - y, 2);
		const double termL = - y * y;
		const double term3 = y * (2 - y);

		return (term2 * f2 + termL * fL + term3 * f3) / x;
}

double x_integrated_cross_section_gsl_single_integrand(double input[], size_t dim, void *params_in) {
	return evaluate_gsl_x_single_integrand(input, dim, params_in, x_integrated_cross_section_single_integrand);
}
double x_integrated_cross_section_gsl_double_integrand(double input[], size_t dim, void *params_in) {
	return evaluate_gsl_x_double_integrand(input, dim, params_in, x_integrated_cross_section_double_integrand);
}

#endif