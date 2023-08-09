#ifndef SIDIS_FUNCTIONS_FL_H
#define SIDIS_FUNCTIONS_FL_H

#include "SIDIS/Coefficients/F2.cpp"
#include "SIDIS/Coefficients/Helper.cpp"

namespace SIDISFunctions::FL {
	namespace LO {
		static constexpr double integrand(
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

			const double prefactor = m2 / (m2 + Q2);
			return prefactor * SIDISFunctions::F2::LO::integrand(
				xi, xip, x, z, 
				factorization_scale_log, fragmentation_scale_log, 
				xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
				xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
				log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
				m2, Q2
			);
		}
	}

	namespace NLO {
		static constexpr double delta_integrand(
			const double xi, 
			const double xip, 
			const double x, 
			const double z, 
			const double fragmentation_scale_log,
			const double factorization_scale_log,
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

			const double prefactor = m2 / (m2 + Q2);
			return prefactor * SIDISFunctions::F2::NLO::delta_integrand(
				xi, xip, x, z, 
				factorization_scale_log, fragmentation_scale_log, 
				xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
				xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
				log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
				m2, Q2
			);
		}

		static constexpr double xi_integrand(
			const double xi, 
			const double xip, 
			const double x, 
			const double z, 
			const double fragmentation_scale_log,
			const double factorization_scale_log,
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

			const double prefactor = m2 / (m2 + Q2);
			return prefactor * SIDISFunctions::F2::NLO::xi_integrand(
				xi, xip, x, z, 
				factorization_scale_log, fragmentation_scale_log, 
				xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
				xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
				log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
				m2, Q2
			);
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

			const double prefactor = m2 / (m2 + Q2);
			return prefactor * SIDISFunctions::F2::NLO::xip_integrand(
				xi, xip, x, z, 
				factorization_scale_log, fragmentation_scale_log, 
				xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
				xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
				log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
				m2, Q2
			);
		}
		static constexpr double xi_xip_integrand(const double xi, 
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

			const double term1 = xq_hat_zq_hat * Constants::C_F * 4 * xi * xip;
			const double term2 = xq_hat_zg_hat * Constants::C_F * 4 * xi * (1 - xip);
			const double term3 = xg_hat_zq_hat * Constants::T_R * 8 * xi * (1 - xi);

			const double value = term1 + term2 + term3;
			const double longitudinal_prefactor = Q2 / (m2 + Q2);
			const double longitudinal_contribution = 2 * longitudinal_prefactor * value / z;

			const double mass_prefactor = m2 / (m2 + Q2);
			const double mass_contribution = mass_prefactor * SIDISFunctions::F2::NLO::xi_xip_integrand(
				xi, xip, x, z, 
				factorization_scale_log, fragmentation_scale_log, 
				xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
				xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
				log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
				m2, Q2
			);

			return longitudinal_contribution + mass_contribution;
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
				FL::NLO::delta_integrand, FL::NLO::xi_integrand, FL::NLO::xip_integrand, FL::NLO::xi_xip_integrand, 
				xi, xip, x, z, 
				factorization_scale_log, fragmentation_scale_log, 
				xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
				xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
				log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
				m2, Q2
			);
		}
	}
}

#endif