#ifndef DIS_FUNCTIONS_FL
#define DIS_FUNCTIONS_FL

#include "Common/Constants.cpp"

#include "DIS/Coefficients/F2.cpp"

namespace DISFunctions::FL {
	namespace LO {
		constexpr double integrand(
			const double xi, 
			const double x, 
			const double factorization_scale_log, 
			const double xq, 
			const double xq_hat, 
			const double xg_hat,
			const double m2,
			const double Q2) {

			const double prefactor = m2 / (m2 + Q2);
			return prefactor * DISFunctions::F2::LO::integrand(
				xi, x,
				factorization_scale_log,
				xq, xq_hat, xg_hat,
				m2, Q2
			);
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

			const double prefactor = m2 / (m2 + Q2);
			return prefactor * DISFunctions::F2::NLO::delta(
				xi, x,
				factorization_scale_log,
				xq, xq_hat, xg_hat,
				m2, Q2
			);
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

			const double quark_contribution = Constants::C_F * xq_hat * 2 * xi;
			const double gluon_contribution = 2 * xi * (1 - xi) * xg_hat;

			const double longitudinal_prefactor = Q2 / (m2 + Q2);
			const double longitudinal_contribution = 2 * longitudinal_prefactor * (quark_contribution + gluon_contribution);

			const double mass_prefactor = m2 / (m2 + Q2);
			const double mass_contribution = mass_prefactor * DISFunctions::F2::NLO::integrand(
				xi, x,
				factorization_scale_log,
				xq, xq_hat, xg_hat,
				m2, Q2
			);

			return longitudinal_contribution + mass_contribution;
		}
	}
}

#endif