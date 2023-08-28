#ifndef DIS_FUNCTIONS_F2
#define DIS_FUNCTIONS_F2

#include "DIS/Coefficients/Helper.cpp"
#include "DIS/Coefficients/Scale.cpp"

namespace DISFunctions::F2 {
	namespace LO {
		constexpr double integrand(
			[[maybe_unused]] const double xi, 
			[[maybe_unused]] const double x, 
			[[maybe_unused]] const double factorization_scale_log, 
			const double xq, 
			[[maybe_unused]] const double xq_hat, 
			[[maybe_unused]] const double xg_hat,
			[[maybe_unused]] const double m2,
			[[maybe_unused]] const double Q2) {

			return 2 * xq;
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

			const double term1 = Helper::delta_contribution(x) * xq;

			const double term2 = factorization_scale_log == 0 ? 0 : Scale::delta(
				xi, x, factorization_scale_log, xq, xq_hat, xg_hat, m2, Q2
			);
			return term1 + term2;
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
			const double term3 = xq_hat * (- (1 + xi * xi) * std::log(xi) / (1 - xi) + 3 + 2 * xi);

			const double quark_contribution = Constants::C_F * (term1 - 1.5 * term2 + term3);

			const double term4 = (xi * xi + (1 - xi) * (1 - xi)) * std::log((1 - xi) / xi);
			const double term5 = -1 + 8 * xi * (1 - xi);

			const double gluon_contribution = 0.5 * xg_hat * (term4 + term5);

			const double coefficient_contribution = 2 * (quark_contribution + gluon_contribution);

			const double scale_contribution = factorization_scale_log == 0 ? 0 : Scale::integrand(
				xi, x, factorization_scale_log, xq, xq_hat, xg_hat, m2, Q2
			);

			return coefficient_contribution + scale_contribution;
		}
	}
}

#endif