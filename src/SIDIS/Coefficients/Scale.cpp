#ifndef SIDIS_FUNCTIONS_SCALE_H
#define SIDIS_FUNCTIONS_SCALE_H

#include "Common/Constants.cpp"

#include "SIDIS/Coefficients/EvaluationParameters.cpp"

namespace SIDISFunctions::Scale {
	static constexpr double delta_integrand(const EvaluationParameters &p) {
		return 2.0 * Constants::C_F * p.xq_zq * ((2.0 * p.log1mz + 1.5) * p.fragmentation_scale_log + (2.0 * p.log1mx + 1.5) * p.factorization_scale_log) / p.z;
	}
	static constexpr double xi_integrand(const EvaluationParameters &p) {
		const double term1 = Constants::C_F * ((1.0 + p.xi * p.xi) * p.xq_hat_zq - 2.0 * p.xq_zq) / (1.0 - p.xi);
		const double term2 = Constants::T_R * p.xg_hat_zq * (std::pow(p.xi, 2) + std::pow(1.0 - p.xi, 2));
		const double result = 2.0 * p.factorization_scale_log * (term1 + term2) / p.z;
		return result;
	}
	static constexpr double xip_integrand(const EvaluationParameters &p) {
		const double term1 = Constants::C_F * ((1.0 + std::pow(p.xip, 2)) * p.xq_zq_hat - 2.0 * p.xq_zq) / (1.0 - p.xip);
		const double term2 = Constants::C_F * p.xq_zg_hat * (1.0 + std::pow(1.0 - p.xip, 2)) / p.xip;
		const double result = 2.0 * p.fragmentation_scale_log * (term1 + term2) / p.z;
		return result;
	}
}

#endif