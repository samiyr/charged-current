#ifndef SIDIS_FUNCTIONS_HELPER_H
#define SIDIS_FUNCTIONS_HELPER_H

#include "Decay/Decay.cpp"
#include "Common/Constants.cpp"
#include <cmath>

namespace SIDISFunctions::Helper {
	template <is_decay_function DecayFunction, typename Kinematics>
	constexpr double compute_z_min(const Kinematics &kinematics, const Decay<DecayFunction> &decay) {
		return std::max({
			decay.lepton_momentum_min / (kinematics.y * kinematics.E_beam), 
			2 * kinematics.x * decay.resonance.mass * decay.hadron.mass / kinematics.Q2,
			decay.z_min_cutoff
		});
	}

	constexpr double delta_contribution([[maybe_unused]] const double x, const double z, const double log1mx, const double log1mz) {
		return 2 * Constants::C_F * (std::pow(log1mx + log1mz, 2) - 8) / (z);
	}

	template <typename Signature>
	static constexpr double make_nlo_integrand(
		const Signature delta_integrand_f,
		const Signature xi_integrand_f,
		const Signature xip_integrand_f,
		const Signature xi_xip_integrand_f,
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
		
		const double delta_integrand = delta_integrand_f(
										xi, xip, x, z, 
										factorization_scale_log, fragmentation_scale_log, 
										xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
										xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
										log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
										m2, Q2);
		const double xi_integrand = xi_integrand_f(
										xi, xip, x, z, 
										factorization_scale_log, fragmentation_scale_log, 
										xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
										xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
										log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
										m2, Q2);
		const double xip_integrand = xip_integrand_f(
										xi, xip, x, z, 
										factorization_scale_log, fragmentation_scale_log, 
										xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
										xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
										log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
										m2, Q2);
		const double xi_xip_integrand = xi_xip_integrand_f(
										xi, xip, x, z, 
										factorization_scale_log, fragmentation_scale_log, 
										xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
										xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
										log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
										m2, Q2);

		const double integrand = delta_integrand / ((1.0 - x) * (1.0 - z)) + xi_integrand / (1.0 - z) + xip_integrand / (1.0 - x) + xi_xip_integrand;
		// const double integrand = delta_integrand + (1.0 - x) * xi_integrand + (1.0 - z) * xip_integrand + (1.0 - x) * (1.0 - z) * xi_xip_integrand;
		return integrand;
	}
}

#endif