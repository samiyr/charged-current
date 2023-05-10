#ifndef DIS_FUNCTIONS_H
#define DIS_FUNCTIONS_H

#include "Constants.cpp"
#include <cmath>
#include "Flavor.cpp"
#include <functional>
#include "PDFCommon.cpp"

namespace DISFunctions {
	template <typename PDFInterface>
	struct CommonParams {
		PDFInterface &pdf;

		const FlavorInfo &flavors;

		const double Q2;
		const double nlo_coefficient;
		const double s;
		const double y_max;

		const Process &process;
	};

	template <typename PDFInterface>
	struct UnintegratedParams {
		const CommonParams<PDFInterface> &common;

		const double unintegrated_parameter_value;
		const double xq;
	};

	template <typename PDFInterface>
	static double evaluate_gsl_xi_integrand(double input[], size_t dim, void *params_in, std::function<double(double, double, double, double)> integrand, bool quark_minus = false) {
		const struct UnintegratedParams<PDFInterface> *params = (struct UnintegratedParams<PDFInterface> *)params_in;
		const double xi = input[0];
		if (std::abs(xi - 1) < 1e-15) { return 0; }
		const double x = params->unintegrated_parameter_value;
		const CommonParams<PDFInterface> &common = params->common;
		const double Q2 = common.Q2;
		const double x_hat = x / xi;
		const double xq = params->xq;
		PDFInterface &pdf = common.pdf;

		pdf.evaluate(x_hat, Q2);

		const FlavorInfo &flavors = common.flavors;

		const double xq_hat = PDFCommon::xq_sum(pdf, flavors, quark_minus, common.process);
		const double xg_hat = PDFCommon::xg_sum(pdf, flavors);
		
		const double integrand_value = integrand(xi, xq_hat, xq, xg_hat);
		return integrand_value;
	}

	static constexpr double F2_integrand(const double xi, const double xq_hat, const double xq, const double xg_hat) {
		const double term1 = std::log(1 - xi) * ((1 + xi * xi) * xq_hat - 2 * xq) / (1 - xi);
		const double term2 = (xq_hat - xq) / (1 - xi);
		const double term3 = xq_hat * (- (1 + xi * xi) * std::log(xi) / (1 - xi) + 3 + 2 * xi);

		const double quark_contribution = Constants::C_F * (term1 - 1.5 * term2 + term3);

		const double term4 = (xi * xi + (1 - xi) * (1 - xi)) * std::log((1 - xi) / xi);
		const double term5 = -1 + 8 * xi * (1 - xi);

		const double gluon_contribution = 0.5 * xg_hat * (term4 + term5);

		return quark_contribution + gluon_contribution;
	}
	template <typename PDFInterface>
	static double F2_integrand_gsl(double input[], size_t dim, void *params_in) {
		return evaluate_gsl_xi_integrand<PDFInterface>(input, dim, params_in, F2_integrand);
	}

	static constexpr double FL_integrand(const double xi, const double xq_hat, const double xq, const double xg_hat) {
		const double quark_contribution = Constants::C_F * xq_hat * 2 * xi;
		const double gluon_contribution = 2 * xi * (1 - xi) * xg_hat;

		return quark_contribution + gluon_contribution;
	}
	template <typename PDFInterface>
	static double FL_integrand_gsl(double input[], size_t dim, void *params_in) {
		return evaluate_gsl_xi_integrand<PDFInterface>(input, dim, params_in, FL_integrand);
	}

	static constexpr double F3_integrand(const double xi, const double xq_hat, const double xq, const double xg_hat) {
		const double term1 = std::log(1 - xi) * ((1 + xi * xi) * xq_hat - 2 * xq) / (1 - xi);
		const double term2 = (xq_hat - xq) / (1 - xi);
		const double term3 = xq_hat * (- (1 + xi * xi) * std::log(xi) / (1 - xi) + 2 + xi);

		const double quark_contribution = Constants::C_F * (term1 - 1.5 * term2 + term3);

		return quark_contribution;
	}
	template <typename PDFInterface>
	static double F3_integrand_gsl(double input[], size_t dim, void *params_in) {
		return evaluate_gsl_xi_integrand<PDFInterface>(input, dim, params_in, F3_integrand, true);
	}

	static constexpr double delta_contribution(const double x) {
		const double log = std::log(1 - x);
		return Constants::C_F * (9.0 / 2.0 + (M_PI * M_PI) / 3.0 - log * log + 1.5 * log);
	}
}

#endif