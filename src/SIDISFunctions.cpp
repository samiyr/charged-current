#ifndef SIDIS_FUNCTIONS_H
#define SIDIS_FUNCTIONS_H

#include "Constants.cpp"
#include <cmath>
#include "Flavor.cpp"
#include <functional>
#include "PDFCommon.cpp"

namespace SIDISFunctions {
	template <typename PDFInterface, typename FFInterface>
	struct CommonParams {
		PDFInterface &pdf1;
		FFInterface &ff1;
		PDFInterface &pdf2;
		FFInterface &ff2;

		const FlavorInfo &flavors;

		const double Q2;
		const double nlo_coefficient;
		const double s;

		const Process &process;
	};

	template <typename PDFInterface, typename FFInterface>
	struct UnintegratedParams {
		const CommonParams<PDFInterface, FFInterface> &common;

		const double x;
		const double z;

		const double xq_zq;
	};

	template <typename PDFInterface, typename FFInterface>
	static double evaluate_gsl_xi_integrand(double input[], size_t dim, void *params_in, std::function<double(double, double, double, double, double, double)> integrand, bool quark_minus = false) {
		const struct UnintegratedParams<PDFInterface, FFInterface> *params = (struct UnintegratedParams<PDFInterface, FFInterface> *)params_in;
		const double xi = input[0];
		if (std::abs(xi - 1) < 1e-15) { return 0; }
		const double x = params->x;
		const double z = params->z;
		const CommonParams<PDFInterface, FFInterface> &common = params->common;
		const double Q2 = common.Q2;
		const double x_hat = x / xi;
		const double xq_zq = params->xq_zq;

		PDFInterface &pdf2 = common.pdf2;
		FFInterface &ff1 = common.ff1;

		pdf2.evaluate(x_hat, Q2);

		const FlavorInfo &flavors = common.flavors;

		const double xq_hat_zq = PDFCommon::xq_zq_sum(pdf2, ff1, flavors, quark_minus, common.process);
		const double xg_hat_zq = PDFCommon::xg_zq_sum(pdf2, ff1, flavors, quark_minus, common.process);

		const double integrand_value = integrand(
			xi, x, z, 
			xq_zq, xq_hat_zq, 
			xg_hat_zq
		);
		return integrand_value;
	}
	template <typename PDFInterface, typename FFInterface>
	static double evaluate_gsl_xip_integrand(double input[], size_t dim, void *params_in, std::function<double(double, double, double, double, double, double)> integrand, bool quark_minus = false) {
		const struct UnintegratedParams<PDFInterface, FFInterface> *params = (struct UnintegratedParams<PDFInterface, FFInterface> *)params_in;
		const double xip = input[0];
		if (std::abs(xip - 1) < 1e-15) { return 0; }
		const double x = params->x;
		const double z = params->z;
		const CommonParams<PDFInterface, FFInterface> &common = params->common;
		const double Q2 = common.Q2;
		const double z_hat = z / xip;
		const double xq_zq = params->xq_zq;

		PDFInterface &pdf1 = common.pdf1;
		FFInterface &ff2 = common.ff2;

		ff2.evaluate(z_hat, Q2);

		const FlavorInfo &flavors = common.flavors;

		const double xq_zq_hat = PDFCommon::xq_zq_sum(pdf1, ff2, flavors, quark_minus, common.process);
		const double xq_zg_hat = PDFCommon::xq_zg_sum(pdf1, ff2, flavors, quark_minus, common.process);

		const double integrand_value = integrand(
			xip, x, z,
			xq_zq, xq_zq_hat, 
			xq_zg_hat
		);
		return integrand_value;
	}
	template <typename PDFInterface, typename FFInterface>
	static double evaluate_gsl_xi_xip_integrand(double input[], size_t dim, void *params_in, std::function<double(double, double, double, double, double, double, double, double, double, double, double, double)> integrand, bool quark_minus = false) {
		const struct UnintegratedParams<PDFInterface, FFInterface> *params = (struct UnintegratedParams<PDFInterface, FFInterface> *)params_in;
		const double xi = input[0];
		const double xip = input[1];
		if (std::abs(xi - 1) < 1e-15) { return 0; }
		if (std::abs(xip - 1) < 1e-15) { return 0; }
		const double x = params->x;
		const double z = params->z;
		const CommonParams<PDFInterface, FFInterface> &common = params->common;
		const double Q2 = common.Q2;
		const double x_hat = x / xi;
		const double z_hat = z / xip;
		const double xq_zq = params->xq_zq;

		PDFInterface &pdf1 = common.pdf1;
		FFInterface &ff1 = common.ff1;
		PDFInterface &pdf2 = common.pdf2;
		FFInterface &ff2 = common.ff2;

		pdf2.evaluate(x_hat, Q2);
		ff2.evaluate(z_hat, Q2);

		const FlavorInfo &flavors = common.flavors;

		const double xq_hat_zq = PDFCommon::xq_zq_sum(pdf2, ff1, flavors, quark_minus, common.process);
		const double xq_zq_hat = PDFCommon::xq_zq_sum(pdf1, ff2, flavors, quark_minus, common.process);
		const double xq_hat_zq_hat = PDFCommon::xq_zq_sum(pdf2, ff2, flavors, quark_minus, common.process);

		const double xq_zg_hat = PDFCommon::xq_zg_sum(pdf1, ff2, flavors, quark_minus, common.process);
		const double xg_hat_zq = PDFCommon::xg_zq_sum(pdf2, ff1, flavors, quark_minus, common.process);
		const double xq_hat_zg_hat = PDFCommon::xq_zg_sum(pdf2, ff2, flavors, quark_minus, common.process);
		const double xg_hat_zq_hat = PDFCommon::xg_zq_sum(pdf2, ff2, flavors, quark_minus, common.process);

		const double integrand_value = integrand(
			xi, xip, x, z, 
			xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, 
			xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat
		);
		return integrand_value;
	}
	template <typename PDFInterface, typename FFInterface>
	static double evaluate_gsl_xi_xip_integrand_FL(double input[], size_t dim, void *params_in, std::function<double(double, double, double, double, double, double, double)> integrand, bool quark_minus = false) {
		const struct UnintegratedParams<PDFInterface, FFInterface> *params = (struct UnintegratedParams<PDFInterface, FFInterface> *)params_in;
		const double xi = input[0];
		const double xip = input[1];
		if (std::abs(xi - 1) < 1e-15) { return 0; }
		if (std::abs(xip - 1) < 1e-15) { return 0; }
		const double x = params->x;
		const double z = params->z;
		const CommonParams<PDFInterface, FFInterface> &common = params->common;
		const double Q2 = common.Q2;
		const double x_hat = x / xi;
		const double z_hat = z / xip;

		PDFInterface &pdf2 = common.pdf2;
		FFInterface &ff2 = common.ff2;

		pdf2.evaluate(x_hat, Q2);
		ff2.evaluate(z_hat, Q2);

		const FlavorInfo &flavors = common.flavors;

		const double xq_hat_zq_hat = PDFCommon::xq_zq_sum(pdf2, ff2, flavors, quark_minus, common.process);
		const double xq_hat_zg_hat = PDFCommon::xq_zg_sum(pdf2, ff2, flavors, quark_minus, common.process);
		const double xg_hat_zq_hat = PDFCommon::xg_zq_sum(pdf2, ff2, flavors, quark_minus, common.process);

		const double integrand_value = integrand(
			xi, xip, x, z, 
			xq_hat_zq_hat, xq_hat_zg_hat, xg_hat_zq_hat
		);
		return integrand_value;
	}

	static constexpr double F2_xi_integrand(
		const double xi,
		const double x,
		const double z, 
		const double xq_zq, 
		const double xq_hat_zq, 
		const double xg_hat_zq) {
		const double term1 = (1 - xi) * (1 + std::log((1 - xi) / xi) + std::log(1 - z)) - (2 * xi * std::log(xi)) / (1 - xi);
		const double term2 = 2 * (std::log(1 - xi) + std::log(1 - z));
		const double term3 = (xi * xq_hat_zq - xq_zq)  / (1 - xi);

		const double quark_contribution = Constants::C_F * (xq_hat_zq * term1 + term2 * term3);

		const double term4 = 1 - (std::pow(xi, 2) + std::pow(1 - xi, 2)) * (1 - std::log((1 - xi) / xi));
		const double term5 = std::log(1 - z) * (1 - 2 * xi * (1 - xi));

		const double gluon_contribution = Constants::T_R * xg_hat_zq * (term4 + term5);

		return quark_contribution + gluon_contribution;
	}
	static constexpr double F3_xi_integrand(
		const double xi,
		const double x,
		const double z, 
		const double xq_zq, 
		const double xq_hat_zq, 
		const double xg_hat_zq) {
		return F2_xi_integrand(xi, x, z, xq_zq, xq_hat_zq, xg_hat_zq);
	}

	static constexpr double F2_xip_integrand(
		const double xip, 
		const double x,
		const double z,
		const double xq_zq, 
		const double xq_zq_hat, 
		const double xq_zg_hat) {
		const double term1 = (1 - xip) * (1 + std::log(xip * (1 - xip)) + std::log(1 - x)) + (2 * xip * std::log(xip)) / (1 - xip);
		const double term2 = 2 * (std::log(1 - xip) + std::log(1 - x));
		const double term3 = (xip * xq_zq_hat - xq_zq)  / (1 - xip);

		const double quark_contribution = Constants::C_F * (xq_zq_hat * term1 + term2 * term3);

		const double term4 = xip + std::log(xip * (1 - xip)) * (1 + std::pow(1 - xip, 2)) / xip;
		const double term5 = std::log(1 - x) * (xip + 2 * (1 - xip) / xip);

		const double gluon_contribution = Constants::C_F * xq_zg_hat * (term4 + term5);

		return quark_contribution + gluon_contribution;
	}
	static constexpr double F3_xip_integrand(
		const double xip, 
		const double x,
		const double z,
		const double xq_zq, 
		const double xq_zq_hat, 
		const double xq_zg_hat) {
		return F2_xip_integrand(xip, x, z, xq_zq, xq_zq_hat, xq_zg_hat);
	}

	static constexpr double F2_xi_xip_integrand(const double xi, 
	const double xip, 
	const double x, 
	const double z, 
	const double xq_zq,
	const double xq_hat_zq, 
	const double xq_zq_hat,
	const double xq_hat_zq_hat, 
	const double xq_zg_hat,
	const double xg_hat_zq,
	const double xq_hat_zg_hat,
	const double xg_hat_zq_hat) {
		const double term1 = (1 - xi) / (1 - xip) * (xq_hat_zq_hat - xq_hat_zq);
		const double term2 = (1 - xip) / (1 - xi) * (xq_hat_zq_hat - xq_zq_hat);
		const double term3 = 2 * (xi * xip * xq_hat_zq_hat - xi * xq_hat_zq - xip * xq_zq_hat + xq_zq) / ((1 - xi) * (1 - xip));
		const double term4 = 6 * xi * xip * xq_hat_zq_hat;

		const double quark_contribution = Constants::C_F * (term1 + term2 + term3 + term4);

		const double term5 = xq_hat_zg_hat * (6 * xi * (1 - xip) + (1 - xi) / xip);
		const double term6 = (xq_hat_zg_hat * (xip + 2 * xi * (1 - xip) / xip) - xq_zg_hat * (xip + 2 * (1 - xip) / xip)) / (1 - xi);
		const double gluon_contribution_1 = Constants::C_F * (term5 + term6);

		const double term7 = xg_hat_zq_hat * (12 * xi * (1 - xi) + (1 - xip) / xip);
		const double term8 = (xg_hat_zq_hat * (xip - 2 * xi * (1 - xi) / xip) - xg_hat_zq * (1 - 2 * xi * (1 - xi))) / (1 - xip);

		const double gluon_contribution_2 = Constants::T_R * (term7 + term8);

		return quark_contribution + gluon_contribution_1 + gluon_contribution_2;
	}

	static constexpr double FL_xi_xip_integrand(const double xi, 
	const double xip, 
	const double x, 
	const double z, 
	const double xq_hat_zq_hat, 
	const double xq_hat_zg_hat,
	const double xg_hat_zq_hat) {
		const double term1 = xq_hat_zq_hat * Constants::C_F * 4 * xi * xip;
		const double term2 = xq_hat_zg_hat * Constants::C_F * 4 * xi * (1 - xip);
		const double term3 = xg_hat_zq_hat * Constants::T_R * 8 * xi * (1 - xi);
		
		const double value = term1 + term2 + term3;
		return value;
	}

	static constexpr double F3_xi_xip_integrand(const double xi, 
	const double xip, 
	const double x, 
	const double z, 
	const double xq_zq,
	const double xq_hat_zq, 
	const double xq_zq_hat,
	const double xq_hat_zq_hat, 
	const double xq_zg_hat,
	const double xg_hat_zq,
	const double xq_hat_zg_hat,
	const double xg_hat_zq_hat) {
		const double F2_value = F2_xi_xip_integrand(xi, xip, x, z, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
		const double term1 = Constants::C_F * xq_hat_zq_hat * (6 * xi * xip + 2 * (1 - xi - xip));
		const double term2 = Constants::C_F * xq_hat_zg_hat * (4 * xi * (1 - xip) + 2 * (1 - xi) * xip);
		const double term3 = Constants::T_R * xg_hat_zq_hat * (12 * xi * (1 - xi) + 2 * (1 - 2 * xi * (1 - xi)) / xip - 2);

		const double value = F2_value - (term1 + term2 + term3);

		return value;
	}

	template <typename PDFInterface, typename FFInterface>
	static double F2_xi_integrand_gsl(double input[], size_t dim, void *params_in) {
		return evaluate_gsl_xi_integrand<PDFInterface, FFInterface>(input, dim, params_in, F2_xi_integrand);
	}
	template <typename PDFInterface, typename FFInterface>
	static double F2_xip_integrand_gsl(double input[], size_t dim, void *params_in) {
		return evaluate_gsl_xip_integrand<PDFInterface, FFInterface>(input, dim, params_in, F2_xip_integrand);
	}
	template <typename PDFInterface, typename FFInterface>
	static double F2_xi_xip_integrand_gsl(double input[], size_t dim, void *params_in) {
		return evaluate_gsl_xi_xip_integrand<PDFInterface, FFInterface>(input, dim, params_in, F2_xi_xip_integrand);
	}
	template <typename PDFInterface, typename FFInterface>
	static double FL_xi_xip_integrand_gsl(double input[], size_t dim, void *params_in) {
		return evaluate_gsl_xi_xip_integrand_FL<PDFInterface, FFInterface>(input, dim, params_in, FL_xi_xip_integrand);
	}
	template <typename PDFInterface, typename FFInterface>
	static double F3_xi_integrand_gsl(double input[], size_t dim, void *params_in) {
		return evaluate_gsl_xi_integrand<PDFInterface, FFInterface>(input, dim, params_in, F3_xi_integrand, true);
	}
	template <typename PDFInterface, typename FFInterface>
	static double F3_xip_integrand_gsl(double input[], size_t dim, void *params_in) {
		return evaluate_gsl_xip_integrand<PDFInterface, FFInterface>(input, dim, params_in, F3_xip_integrand, true);
	}
	template <typename PDFInterface, typename FFInterface>
	static double F3_xi_xip_integrand_gsl(double input[], size_t dim, void *params_in) {
		return evaluate_gsl_xi_xip_integrand<PDFInterface, FFInterface>(input, dim, params_in, F3_xi_xip_integrand, true);
	}

	constexpr double delta_contribution(const double x, const double z) {
		return Constants::C_F * (std::pow(std::log(1 - x) + std::log(1 - z), 2) - 8);
	}
}

#endif