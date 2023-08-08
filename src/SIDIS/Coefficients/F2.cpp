#ifndef SIDIS_FUNCTIONS_F2_H
#define SIDIS_FUNCTIONS_F2_H

#include "SIDIS/Coefficients/Scale.cpp"
#include "SIDIS/Coefficients/Helper.cpp"

namespace SIDISFunctions::F2 {
	namespace LO {
		static constexpr double integrand(
			[[maybe_unused]] const double xi, 
			[[maybe_unused]] const double xip, 
			[[maybe_unused]] const double x, 
			const double z,
			[[maybe_unused]] const double factorization_scale_log,
			[[maybe_unused]] const double fragmentation_scale_log,
			const double xq_zq,
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

			return 2 * xq_zq / (z);
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

			const double term1 = Helper::delta_contribution(x, z, log1mx, log1mz) * xq_zq;
			const double term2 = (factorization_scale_log == 0 && fragmentation_scale_log == 0) 
									? 0
									: Scale::delta_integrand(xi, xip, x, z, 
										factorization_scale_log, fragmentation_scale_log, 
										xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
										xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
										log1mx, log1mz, logxi, logxip, log1mxi, log1mxip, 
										m2, Q2);
			return term1 + term2;
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

			const double log_term = log1mxi - logxi; // std::log((1 - xi) / xi);

			const double term1 = (1 - xi) * (1 + log_term + log1mz) - (2 * xi * logxi) / (1 - xi);
			const double term2 = 2 * (log1mxi + log1mz);
			const double term3 = (xi * xq_hat_zq - xq_zq)  / (1 - xi);

			const double quark_contribution = Constants::C_F * (xq_hat_zq * term1 + term2 * term3);

			const double term4 = 1 - (std::pow(xi, 2) + std::pow(1 - xi, 2)) * (1 - log_term);
			const double term5 = log1mz * (1 - 2 * xi * (1 - xi));

			const double gluon_contribution = Constants::T_R * xg_hat_zq * (term4 + term5);

			const double coefficient_contribution = 2 * (quark_contribution + gluon_contribution) / (z);
			const double scale_contribution = (factorization_scale_log == 0 && fragmentation_scale_log == 0) 
												? 0 
												: Scale::xi_integrand(xi, xip, x, z, 
													factorization_scale_log, fragmentation_scale_log, 
													xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
													xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
													log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
													m2, Q2);

			return coefficient_contribution + scale_contribution;
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

			const double log_term = logxip + log1mxip; // std::log(xip * (1 - xip));
			
			const double term1 = (1 - xip) * (1 + log_term + log1mx) + (2 * xip * logxip) / (1 - xip);
			const double term2 = 2 * (log1mxip + log1mx);
			const double term3 = (xip * xq_zq_hat - xq_zq)  / (1 - xip);

			const double quark_contribution = Constants::C_F * (xq_zq_hat * term1 + term2 * term3);

			const double term4 = xip + log_term * (1 + std::pow(1 - xip, 2)) / xip;
			const double term5 = log1mx * (xip + 2 * (1 - xip) / xip);

			const double gluon_contribution = Constants::C_F * xq_zg_hat * (term4 + term5);

			const double coefficient_contribution = 2 * (quark_contribution + gluon_contribution) / (z);
			const double scale_contribution = (factorization_scale_log == 0 && fragmentation_scale_log == 0) 
												? 0 
												: Scale::xip_integrand(xi, xip, x, z, 
													factorization_scale_log, fragmentation_scale_log, 
													xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
													xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
													log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
													m2, Q2);

			return coefficient_contribution + scale_contribution;
		}

		static constexpr double xi_xip_integrand(const double xi, 
			const double xip, 
			[[maybe_unused]] const double x, 
			const double z, 
			[[maybe_unused]] const double factorization_scale_log,
			[[maybe_unused]] const double fragmentation_scale_log,
			const double xq_zq,
			const double xq_hat_zq, 
			const double xq_zq_hat,
			const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			const double xg_hat_zq,
			const double xq_hat_zg_hat,
			const double xg_hat_zq_hat,
			[[maybe_unused]] const double log1mx,
			[[maybe_unused]] const double log1mz,
			[[maybe_unused]] const double logxi,
			[[maybe_unused]] const double logxip,
			[[maybe_unused]] const double log1mxi,
			[[maybe_unused]] const double log1mxip,
			[[maybe_unused]] const double m2,
			[[maybe_unused]] const double Q2) {

			const double term1 = (1 - xi) / (1 - xip) * (xq_hat_zq_hat - xq_hat_zq);
			const double term2 = (1 - xip) / (1 - xi) * (xq_hat_zq_hat - xq_zq_hat);
			const double term3 = 2 * (xi * xip * xq_hat_zq_hat - xi * xq_hat_zq - xip * xq_zq_hat + xq_zq) / ((1 - xi) * (1 - xip));
			const double term4 = 6 * xi * xip * xq_hat_zq_hat;

			const double quark_contribution = Constants::C_F * (term1 + term2 + term3 + term4);

			const double term5 = xq_hat_zg_hat * (6 * xi * (1 - xip) + (1 - xi) / xip);
			const double term6 = (xq_hat_zg_hat * (xip + 2 * xi * (1 - xip) / xip) - xq_zg_hat * (xip + 2 * (1 - xip) / xip)) / (1 - xi);
			const double gluon_contribution_1 = Constants::C_F * (term5 + term6);

			const double term7 = xg_hat_zq_hat * (12 * xi * (1 - xi) + (1 - xip) / xip);
			const double term8 = (xg_hat_zq_hat * (xip - 2 * xi * (1 - xi) / xip) - xg_hat_zq * (1 - 2 * xi * (1 - xi))) / (1 - xip);

			const double gluon_contribution_2 = Constants::T_R * (term7 + term8);

			return 2 * (quark_contribution + gluon_contribution_1 + gluon_contribution_2) / (z);
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
				F2::NLO::delta_integrand, F2::NLO::xi_integrand, F2::NLO::xip_integrand, F2::NLO::xi_xip_integrand, 
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