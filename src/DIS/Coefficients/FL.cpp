#ifndef DIS_FUNCTIONS_FL
#define DIS_FUNCTIONS_FL

#include "Common/Constants.cpp"

namespace DISFunctions::FL {
	namespace LO {
		static constexpr double integrand([[maybe_unused]] const double xi, [[maybe_unused]] const double x, [[maybe_unused]] const double factorization_scale_log, [[maybe_unused]] const double xq, [[maybe_unused]] const double xq_hat, [[maybe_unused]] const double xg_hat) {
			return 0.0;
		}
	}

	namespace NLO {
		static constexpr double delta([[maybe_unused]] const double xi, [[maybe_unused]] const double x, [[maybe_unused]] const double factorization_scale_log, [[maybe_unused]] const double xq, [[maybe_unused]] const double xq_hat, [[maybe_unused]] const double xg_hat) {
			return 0.0;
		}

		static constexpr double integrand(const double xi, [[maybe_unused]] const double x, [[maybe_unused]] const double factorization_scale_log, [[maybe_unused]] const double xq, const double xq_hat, const double xg_hat) {
			const double quark_contribution = Constants::C_F * xq_hat * 2 * xi;
			const double gluon_contribution = 2 * xi * (1 - xi) * xg_hat;

			const double total_contribution = 2 * (quark_contribution + gluon_contribution);
			return total_contribution;
		}
	}
}

#endif