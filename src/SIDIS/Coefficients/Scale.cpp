#ifndef SIDIS_FUNCTIONS_SCALE_H
#define SIDIS_FUNCTIONS_SCALE_H

#include "Common/Constants.cpp"

namespace SIDISFunctions::Scale {
	static constexpr double delta_integrand(
		[[maybe_unused]] const double xi,
		[[maybe_unused]] const double xip,
		[[maybe_unused]] const double x, 
		const double z, 
		const double factorization_scale_log,
		const double fragmentation_scale_log,
		const double xq_zq,
		[[maybe_unused]] const double xq_hat_zq, 
		[[maybe_unused]] const double xq_zq_hat,
		[[maybe_unused]] const double xq_hat_zq_hat, 
		[[maybe_unused]] const double xq_zg_hat,
		[[maybe_unused]] const double xg_hat_zq,
		[[maybe_unused]] const double xq_hat_zg_hat,
		[[maybe_unused]] const double xg_hat_zq_hat,
		const double log1mx,
		const double log1mz,
		[[maybe_unused]] const double logxi,
		[[maybe_unused]] const double logxip,
		[[maybe_unused]] const double log1mxi,
		[[maybe_unused]] const double log1mxip,
		[[maybe_unused]] const double m2,
		[[maybe_unused]] const double Q2) {
		return 2 * Constants::C_F * xq_zq * ((2 * log1mz + 1.5) * fragmentation_scale_log + (2 * log1mx + 1.5) * factorization_scale_log) / (z);
	}
	static constexpr double xi_integrand(
		const double xi,
		[[maybe_unused]] const double xip,
		[[maybe_unused]] const double x, 
		const double z, 
		const double factorization_scale_log,
		[[maybe_unused]] const double fragmentation_scale_log,
		const double xq_zq,
		const double xq_hat_zq, 
		[[maybe_unused]] const double xq_zq_hat,
		[[maybe_unused]] const double xq_hat_zq_hat, 
		[[maybe_unused]] const double xq_zg_hat,
		const double xg_hat_zq,
		[[maybe_unused]] const double xq_hat_zg_hat,
		[[maybe_unused]] const double xg_hat_zq_hat,
		[[maybe_unused]] const double log1mx,
		[[maybe_unused]] const double log1mz,
		[[maybe_unused]] const double logxi,
		[[maybe_unused]] const double logxip,
		[[maybe_unused]] const double log1mxi,
		[[maybe_unused]] const double log1mxip,
		[[maybe_unused]] const double m2,
		[[maybe_unused]] const double Q2) {

		const double term1 = Constants::C_F * ((1 + xi * xi) * xq_hat_zq - 2 * xq_zq) / (1 - xi);
		const double term2 = Constants::T_R * xg_hat_zq * (std::pow(xi, 2) + std::pow(1 - xi, 2));
		const double result = 2 * factorization_scale_log * (term1 + term2) / (z);
		return result;
	}
	static constexpr double xip_integrand(
		[[maybe_unused]] const double xi,
		const double xip,
		[[maybe_unused]] const double x, 
		const double z, 
		[[maybe_unused]] const double factorization_scale_log,
		const double fragmentation_scale_log,
		const double xq_zq,
		[[maybe_unused]] const double xq_hat_zq, 
		const double xq_zq_hat,
		[[maybe_unused]] const double xq_hat_zq_hat, 
		const double xq_zg_hat,
		[[maybe_unused]] const double xg_hat_zq,
		[[maybe_unused]] const double xq_hat_zg_hat,
		[[maybe_unused]] const double xg_hat_zq_hat,
		[[maybe_unused]] const double log1mx,
		[[maybe_unused]] const double log1mz,
		[[maybe_unused]] const double logxi,
		[[maybe_unused]] const double logxip,
		[[maybe_unused]] const double log1mxi,
		[[maybe_unused]] const double log1mxip,
		[[maybe_unused]] const double m2,
		[[maybe_unused]] const double Q2) {

		const double term1 = Constants::C_F * ((1 + std::pow(xip, 2)) * xq_zq_hat - 2 * xq_zq) / (1 - xip);
		const double term2 = Constants::C_F * xq_zg_hat * (1 + std::pow(1 - xip, 2)) / xip;
		const double result = 2 * fragmentation_scale_log * (term1 + term2) / (z);
		return result;
	}
}

#endif