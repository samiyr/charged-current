#ifndef DIS_FUNCTIONS_F2
#define DIS_FUNCTIONS_F2

#include "DIS/Coefficients/Utility.cpp"

namespace DISFunctions::F2 {
	namespace LO {
		static constexpr double integrand([[maybe_unused]] const double xi, [[maybe_unused]] const double x, [[maybe_unused]] const double factorization_scale_log, const double xq, [[maybe_unused]] const double xq_hat, [[maybe_unused]] const double xg_hat) {
			return 2 * xq;
		}
	}

	namespace NLO {
		static constexpr double delta([[maybe_unused]] const double xi, const double x, const double factorization_scale_log, const double xq, [[maybe_unused]] const double xq_hat, [[maybe_unused]] const double xg_hat) {
			const double term1 = DISFunctions::delta_contribution(x) * xq;
			return term1;
		}

		static constexpr double integrand(const double xi, [[maybe_unused]] const double x, const double factorization_scale_log, const double xq, const double xq_hat, const double xg_hat) {
			const double term1 = std::log(1 - xi) * ((1 + xi * xi) * xq_hat - 2 * xq) / (1 - xi);
			const double term2 = (xq_hat - xq) / (1 - xi);
			const double term3 = xq_hat * (- (1 + xi * xi) * std::log(xi) / (1 - xi) + 3 + 2 * xi);

			const double quark_contribution = Constants::C_F * (term1 - 1.5 * term2 + term3);

			const double term4 = (xi * xi + (1 - xi) * (1 - xi)) * std::log((1 - xi) / xi);
			const double term5 = -1 + 8 * xi * (1 - xi);

			const double gluon_contribution = 0.5 * xg_hat * (term4 + term5);

			const double total_contribution = 2 * (quark_contribution + gluon_contribution);
			return total_contribution;
		}
	}
}

#endif