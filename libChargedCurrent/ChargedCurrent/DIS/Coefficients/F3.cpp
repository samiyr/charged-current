#ifndef DIS_FUNCTIONS_F3
#define DIS_FUNCTIONS_F3

#include "Common/Constants.cpp"

#include "DIS/Coefficients/F2.cpp"

namespace DISFunctions::F3 {
	namespace LO {
		constexpr double integrand(
			[[maybe_unused]] const double xi, 
			const double x, 
			[[maybe_unused]] const double factorization_scale_log, 
			const double xq, 
			[[maybe_unused]] const double xq_hat, 
			[[maybe_unused]] const double xg_hat,
			[[maybe_unused]] const double m2,
			[[maybe_unused]] const double Q2) {

			return 2 * xq / x;
		}
	}
	
	namespace NLO {
		constexpr double delta(
			const double xi, 
			const double x, 
			const double factorization_scale_log, 
			const double xq, 
			const double xq_hat, 
			const double xg_hat,
			const double m2,
			const double Q2) {
				
			return F2::NLO::delta(xi, x, factorization_scale_log, xq, xq_hat, xg_hat, m2, Q2) / x;
		}
		constexpr double integrand(
			const double xi, 
			const double x, 
			const double factorization_scale_log, 
			const double xq, 
			const double xq_hat, 
			const double xg_hat,
			const double m2,
			const double Q2) {

			const double term1 = std::log(1 - xi) * ((1 + xi * xi) * xq_hat - 2 * xq) / (1 - xi);
			const double term2 = (xq_hat - xq) / (1 - xi);
			const double term3 = xq_hat * (- (1 + xi * xi) * std::log(xi) / (1 - xi) + 2 + xi);

			const double quark_contribution = 2 * Constants::C_F * (term1 - 1.5 * term2 + term3) / x;

			const double scale_contribution = factorization_scale_log == 0 ? 0 : Scale::integrand(
				xi, x, factorization_scale_log, xq, xq_hat, xg_hat, m2, Q2
			);

			return quark_contribution + scale_contribution;
		}
	}
}

#endif