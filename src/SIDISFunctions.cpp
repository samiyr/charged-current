#ifndef SIDIS_FUNCTIONS_H
#define SIDIS_FUNCTIONS_H

#include "Constants.cpp"
#include <cmath>
#include "Flavor.cpp"
#include <functional>
#include "PDFCommon.cpp"

namespace SIDISFunctions {
	template <typename PDFInterface, typename FFInterface>
	struct Parameters {
		PDFInterface &pdf1;
		FFInterface &ff1;
		PDFInterface &pdf2;
		FFInterface &ff2;

		const FlavorInfo &flavors;

		const double Q2;
		const double nlo_coefficient;
		const double s;

		const Process &process;

		const double x;
		const double z;

		const double xq_zq;
	};

	namespace Evaluation {
		template <typename PDFInterface, typename FFInterface>
		static double evaluate_gsl_sidis_integrand(double input[], size_t dim, void *params_in, std::function<double(double, double, double, double, double, double, double, double, double, double, double, double)> integrand, bool xi_int, bool xip_int, bool z_int, bool quark_minus) {
			const struct Parameters<PDFInterface, FFInterface> *params = static_cast<Parameters<PDFInterface, FFInterface> *>(params_in);

			const bool xi_xip_int = xi_int && xip_int;

			const size_t xi_index = 0;
			const size_t xip_index = xi_int ? 1 : 0;

			const double xi = xi_int ? input[xi_index] : 1.0;
			const double xip = xip_int ? input[xip_index] : 1.0;

			const double x = params->x;
			const double z = z_int ? input[2] : params->z;
			const double Q2 = params->Q2;

			const double x_hat = x / xi;
			const double z_hat = z / xip;
			const double xq_zq = params->xq_zq;

			PDFInterface &pdf1 = params->pdf1;
			FFInterface &ff1 = params->ff1;
			PDFInterface &pdf2 = params->pdf2;
			FFInterface &ff2 = params->ff2;

			if (xi_int && std::abs(xi - 1) < 1e-15) { return 0; }
			if (xip_int && std::abs(xip - 1) < 1e-15) { return 0; }

			if (xi_int) { pdf2.evaluate(x_hat, Q2); }
			if (xip_int) { ff2.evaluate(z_hat, Q2); }

			const FlavorInfo &flavors = params->flavors;

			const double xq_hat_zq = xi_int ? PDFCommon::xq_zq_sum(pdf2, ff1, flavors, quark_minus, params->process) : 0.0;
			const double xq_zq_hat = xip_int ? PDFCommon::xq_zq_sum(pdf1, ff2, flavors, quark_minus, params->process) : 0.0;
			const double xq_hat_zq_hat = xi_xip_int ? PDFCommon::xq_zq_sum(pdf2, ff2, flavors, quark_minus, params->process) : 0.0;

			const double xq_zg_hat = xip_int ? PDFCommon::xq_zg_sum(pdf1, ff2, flavors, quark_minus, params->process) : 0.0;
			const double xg_hat_zq = xi_int ? PDFCommon::xg_zq_sum(pdf2, ff1, flavors, quark_minus, params->process) : 0.0;
			const double xq_hat_zg_hat = xi_xip_int ? PDFCommon::xq_zg_sum(pdf2, ff2, flavors, quark_minus, params->process) : 0.0;
			const double xg_hat_zq_hat = xi_xip_int ? PDFCommon::xg_zq_sum(pdf2, ff2, flavors, quark_minus, params->process) : 0.0;

			const double integrand_value = integrand(
				xi, xip, x, z, 
				xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, 
				xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat
			);
			return integrand_value;
		}
	}

	namespace Integrands {
		static constexpr double F2_xi_integrand(
			const double xi, 
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

			const double log_term = std::log((1 - xi) / xi);

			const double term1 = (1 - xi) * (1 + log_term + Utility::logm1(z)) - (2 * xi * std::log(xi)) / (1 - xi);
			const double term2 = 2 * (Utility::logm1(xi) + Utility::logm1(z));
			const double term3 = (xi * xq_hat_zq - xq_zq)  / (1 - xi);

			const double quark_contribution = Constants::C_F * (xq_hat_zq * term1 + term2 * term3);

			const double term4 = 1 - (std::pow(xi, 2) + std::pow(1 - xi, 2)) * (1 - log_term);
			const double term5 = Utility::logm1(z) * (1 - 2 * xi * (1 - xi));

			const double gluon_contribution = Constants::T_R * xg_hat_zq * (term4 + term5);

			return quark_contribution + gluon_contribution;
		}
		static constexpr double F3_xi_integrand(
			const double xi, 
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
			return F2_xi_integrand(xi, xip, x, z, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
		}

		static constexpr double F2_xip_integrand(
			const double xi, 
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

			const double log_term = std::log(xip * (1 - xip));
			
			const double term1 = (1 - xip) * (1 + log_term + Utility::logm1(x)) + (2 * xip * std::log(xip)) / (1 - xip);
			const double term2 = 2 * (Utility::logm1(xip) + Utility::logm1(x));
			const double term3 = (xip * xq_zq_hat - xq_zq)  / (1 - xip);

			const double quark_contribution = Constants::C_F * (xq_zq_hat * term1 + term2 * term3);

			const double term4 = xip + log_term * (1 + std::pow(1 - xip, 2)) / xip;
			const double term5 = Utility::logm1(x) * (xip + 2 * (1 - xip) / xip);

			const double gluon_contribution = Constants::C_F * xq_zg_hat * (term4 + term5);

			return quark_contribution + gluon_contribution;
		}
		static constexpr double F3_xip_integrand(
			const double xi, 
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
			return F2_xip_integrand(xi, xip, x, z, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
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
		const double xq_zq,
		const double xq_hat_zq, 
		const double xq_zq_hat,
		const double xq_hat_zq_hat, 
		const double xq_zg_hat,
		const double xg_hat_zq,
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
	}

	template <typename PDFInterface, typename FFInterface>
	static double F2_xi_integrand_gsl(double input[], size_t dim, void *params_in) {
		return Evaluation::evaluate_gsl_sidis_integrand<PDFInterface, FFInterface>(input, dim, params_in, Integrands::F2_xi_integrand, true, false, false);
	}
	template <typename PDFInterface, typename FFInterface>
	static double F2_xip_integrand_gsl(double input[], size_t dim, void *params_in) {
		return Evaluation::evaluate_gsl_sidis_integrand<PDFInterface, FFInterface>(input, dim, params_in, Integrands::F2_xip_integrand, false, true, false);
	}
	template <typename PDFInterface, typename FFInterface>
	static double F2_xi_xip_integrand_gsl(double input[], size_t dim, void *params_in) {
		return Evaluation::evaluate_gsl_sidis_integrand<PDFInterface, FFInterface>(input, dim, params_in, Integrands::F2_xi_xip_integrand, true, true, false);
	}
	template <typename PDFInterface, typename FFInterface>
	static double FL_xi_xip_integrand_gsl(double input[], size_t dim, void *params_in) {
		return Evaluation::evaluate_gsl_sidis_integrand<PDFInterface, FFInterface>(input, dim, params_in, Integrands::FL_xi_xip_integrand, true, true, false);
	}
	template <typename PDFInterface, typename FFInterface>
	static double F3_xi_integrand_gsl(double input[], size_t dim, void *params_in) {
		return Evaluation::evaluate_gsl_sidis_integrand<PDFInterface, FFInterface>(input, dim, params_in, Integrands::F3_xi_integrand, true, false, true);
	}
	template <typename PDFInterface, typename FFInterface>
	static double F3_xip_integrand_gsl(double input[], size_t dim, void *params_in) {
		return Evaluation::evaluate_gsl_sidis_integrand<PDFInterface, FFInterface>(input, dim, params_in, Integrands::F3_xip_integrand, false, true, true);
	}
	template <typename PDFInterface, typename FFInterface>
	static double F3_xi_xip_integrand_gsl(double input[], size_t dim, void *params_in) {
		return Evaluation::evaluate_gsl_sidis_integrand<PDFInterface, FFInterface>(input, dim, params_in, Integrands::F3_xi_xip_integrand, true, true, true);
	}

	constexpr double delta_contribution(const double x, const double z) {
		return Constants::C_F * (std::pow(std::log(1 - x) + std::log(1 - z), 2) - 8);
	}
}

#endif