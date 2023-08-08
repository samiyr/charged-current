#ifndef SIDIS_FUNCTIONS_F3_H
#define SIDIS_FUNCTIONS_F3_H

#include "SIDIS/Coefficients/F2.cpp"
#include "SIDIS/Coefficients/Helper.cpp"

namespace SIDISFunctions::F3 {
	namespace LO {
		static constexpr double integrand(
			[[maybe_unused]] const double xi, 
			[[maybe_unused]] const double xip, 
			const double x, 
			const double z, 
			[[maybe_unused]] const double factorization_scale_log,
			[[maybe_unused]] const double fragmentation_scale_log,
			[[maybe_unused]] const double xq_zq,
			[[maybe_unused]] const double xq_hat_zq, 
			[[maybe_unused]] const double xq_zq_hat,
			[[maybe_unused]] const double xq_hat_zq_hat, 
			[[maybe_unused]] const double xq_zg_hat,
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
			return 2 * xq_zq / (x * z);
		}
	}

	namespace NLO {
		static constexpr double delta_integrand(
			const double xi, 
			const double xip, 
			const double x, 
			const double z, 
			const double factorization_scale_log,
			const double fragmentation_scale_log,
			const double xq_zq,
			const double xq_hat_zq, 
			const double xq_zq_hat,
			const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			const double xg_hat_zq,
			const double xq_hat_zg_hat,
			const double xg_hat_zq_hat,
			const double log1mx,
			const double log1mz,
			const double logxi,
			const double logxip,
			const double log1mxi,
			const double log1mxip,
			const double m2,
			const double Q2) {

			return F2::NLO::delta_integrand(
				xi, xip, x, z, 
				factorization_scale_log, fragmentation_scale_log, 
				xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
				xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
				log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
				m2, Q2
			) / x;
		}

		static constexpr double xi_integrand(
			const double xi, 
			const double xip, 
			const double x, 
			const double z, 
			const double factorization_scale_log,
			const double fragmentation_scale_log,
			const double xq_zq,
			const double xq_hat_zq, 
			const double xq_zq_hat,
			const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			const double xg_hat_zq,
			const double xq_hat_zg_hat,
			const double xg_hat_zq_hat,
			const double log1mx,
			const double log1mz,
			const double logxi,
			const double logxip,
			const double log1mxi,
			const double log1mxip,
			const double m2,
			const double Q2) {
			return F2::NLO::xi_integrand(
				xi, xip, x, z, 
				factorization_scale_log, fragmentation_scale_log, 
				xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
				xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
				log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
				m2, Q2) / x;
		}

		static constexpr double xip_integrand(
			const double xi, 
			const double xip, 
			const double x, 
			const double z, 
			const double factorization_scale_log,
			const double fragmentation_scale_log,
			const double xq_zq,
			const double xq_hat_zq, 
			const double xq_zq_hat,
			const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			const double xg_hat_zq,
			const double xq_hat_zg_hat,
			const double xg_hat_zq_hat,
			const double log1mx,
			const double log1mz,
			const double logxi,
			const double logxip,
			const double log1mxi,
			const double log1mxip,
			const double m2,
			const double Q2) {
			return F2::NLO::xip_integrand(
				xi, xip, x, z, 
				factorization_scale_log, fragmentation_scale_log, 
				xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
				xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
				log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
				m2, Q2) / x;
		}

		static constexpr double xi_xip_integrand(
			const double xi, 
			const double xip, 
			const double x, 
			const double z, 
			const double factorization_scale_log,
			const double fragmentation_scale_log,
			const double xq_zq,
			const double xq_hat_zq, 
			const double xq_zq_hat,
			const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			const double xg_hat_zq,
			const double xq_hat_zg_hat,
			const double xg_hat_zq_hat,
			const double log1mx,
			const double log1mz,
			const double logxi,
			const double logxip,
			const double log1mxi,
			const double log1mxip,
			const double m2,
			const double Q2) {

			const double F2_value = F2::NLO::xi_xip_integrand(
										xi, xip, x, z, 
										factorization_scale_log, fragmentation_scale_log, 
										xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
										xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
										log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
										m2, Q2);
			const double term1 = Constants::C_F * xq_hat_zq_hat * (6 * xi * xip + 2 * (1 - xi - xip));
			const double term2 = Constants::C_F * xq_hat_zg_hat * (4 * xi * (1 - xip) + 2 * (1 - xi) * xip);
			const double term3 = Constants::T_R * xg_hat_zq_hat * (12 * xi * (1 - xi) + 2 * (1 - 2 * xi * (1 - xi)) / xip - 2);

			const double value = F2_value / x - 2 * (term1 + term2 + term3) / (x * z);

			return value;
		}

		static constexpr double total_integrand(
			const double xi,
			const double xip, 
			const double x, 
			const double z, 
			const double factorization_scale_log,
			const double fragmentation_scale_log,
			const double xq_zq,
			const double xq_hat_zq, 
			const double xq_zq_hat,
			const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			const double xg_hat_zq,
			const double xq_hat_zg_hat,
			const double xg_hat_zq_hat,
			const double log1mx,
			const double log1mz,
			const double logxi,
			const double logxip,
			const double log1mxi,
			const double log1mxip,
			const double m2,
			const double Q2) {
			
			return Helper::make_nlo_integrand(
				F3::NLO::delta_integrand, F3::NLO::xi_integrand, F3::NLO::xip_integrand, F3::NLO::xi_xip_integrand, 
				xi, xip, x, z, factorization_scale_log, fragmentation_scale_log, 
				xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
				xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
				log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
				m2, Q2
			);
		}

	}
}

#endif