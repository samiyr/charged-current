#ifndef DIS_FUNCTIONS_H
#define DIS_FUNCTIONS_H

#include "Constants.cpp"
#include <cmath>
#include "Flavor.cpp"
#include "PDFInterface.cpp"
#include "FFInterface.cpp"

namespace DISFunctions {
	struct CommonParams {
		PDFInterface &pdf;

		const FlavorVector &active_flavors;
		const FlavorVector &active_antiflavors;

		const FlavorVector &active_upper_flavors;
		const FlavorVector &active_lower_flavors;
		const FlavorVector &active_upper_antiflavors;
		const FlavorVector &active_lower_antiflavors;

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

	template <typename F>
	double evaluate_gsl_xi_integrand(double input[], size_t dim, void *params_in, F integrand, bool quark_minus = false) {
		const struct UnintegratedParams *params = (struct UnintegratedParams *)params_in;
		const double xi = input[0];
		if (std::abs(xi - 1) < 1e-15) { return 0; }
		const double x = params->unintegrated_parameter_value;
		const CommonParams &common = params->common;
		const double Q2 = common.Q2;
		const double x_hat = x / xi;
		const double xq = params->xq;
		PDFInterface &pdf = common.pdf;

		pdf.evaluate(x_hat, Q2);

		const FlavorVector &active_upper_flavors = common.active_upper_flavors;
		const FlavorVector &active_lower_flavors = common.active_lower_flavors;
		const FlavorVector &active_upper_antiflavors = common.active_upper_antiflavors;
		const FlavorVector &active_lower_antiflavors = common.active_lower_antiflavors;

		const double xq_hat = pdf.xq_sum(active_upper_flavors, active_lower_flavors, active_upper_antiflavors, active_lower_antiflavors, quark_minus, common.process);
		const double xg_hat = pdf.xg_sum(active_upper_flavors, active_lower_flavors);
		
		const double integrand_value = integrand(xi, xq_hat, xq, xg_hat);
		return integrand_value;
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

	constexpr double delta_contribution(const double x) {
		const double log = std::log(1 - x);
		return Constants::C_F * (9.0 / 2.0 + (M_PI * M_PI) / 3.0 - log * log + 1.5 * log);
	}
} 


#endif