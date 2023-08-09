#ifndef DIS_FUNCTIONS_SCALE_H
#define DIS_FUNCTIONS_SCALE_H

#include "Common/Constants.cpp"

namespace DISFunctions::Scale {
	static constexpr double delta(
		const double xi, 
		const double x, 
		const double factorization_scale_log, 
		const double xq, 
		const double xq_hat, 
		const double xg_hat,
		const double m2,
		const double Q2) {
			
		return 2 * Constants::C_F * factorization_scale_log * xq * (2.0 * std::log(1.0 - x) + 1.5);
	}
	static constexpr double integrand(
		const double xi, 
		const double x, 
		const double factorization_scale_log, 
		const double xq, 
		const double xq_hat, 
		const double xg_hat,
		const double m2,
		const double Q2) {

		return 2 * Constants::C_F * factorization_scale_log * (xq_hat * (1.0 + std::pow(xi, 2)) - 2.0 * xq) / (1.0 - xi);
	}
}

#endif