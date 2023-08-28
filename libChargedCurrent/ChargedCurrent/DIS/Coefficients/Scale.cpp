#ifndef DIS_FUNCTIONS_SCALE_H
#define DIS_FUNCTIONS_SCALE_H

#include "Common/Constants.cpp"

namespace DISFunctions::Scale {
	constexpr double delta(
		[[maybe_unused]] const double xi, 
		const double x, 
		const double factorization_scale_log, 
		const double xq, 
		[[maybe_unused]] const double xq_hat, 
		[[maybe_unused]] const double xg_hat,
		[[maybe_unused]] const double m2,
		[[maybe_unused]] const double Q2) {
			
		return 2 * Constants::C_F * factorization_scale_log * xq * (2.0 * std::log(1.0 - x) + 1.5);
	}
	constexpr double integrand(
		const double xi, 
		[[maybe_unused]] const double x, 
		const double factorization_scale_log, 
		const double xq, 
		const double xq_hat, 
		[[maybe_unused]] const double xg_hat,
		[[maybe_unused]] const double m2,
		[[maybe_unused]] const double Q2) {

		return 2 * Constants::C_F * factorization_scale_log * (xq_hat * (1.0 + std::pow(xi, 2)) - 2.0 * xq) / (1.0 - xi);
	}
}

#endif