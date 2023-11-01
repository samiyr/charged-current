#ifndef SIDIS_FUNCTIONS_HELPER_H
#define SIDIS_FUNCTIONS_HELPER_H

#include <cmath>

#include "Decay/Decay.cpp"

#include "Common/Constants.cpp"

#include "SIDIS/Coefficients/EvaluationParameters.cpp"

namespace SIDISFunctions::Helper {
	template <is_decay_function DecayFunction, typename Kinematics>
	constexpr double compute_z_min(const Kinematics &kinematics, const Decay<DecayFunction> &decay) {
		return std::max({
			decay.lepton_momentum_min / (kinematics.y * kinematics.E_beam), 
			decay.resonance.mass / (kinematics.y * kinematics.E_beam),
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
		const EvaluationParameters &p) {
		
		const double delta_integrand = delta_integrand_f(p);
		const double xi_integrand = xi_integrand_f(p);
		const double xip_integrand = xip_integrand_f(p);
		const double xi_xip_integrand = xi_xip_integrand_f(p);

		const double integrand = delta_integrand / ((1.0 - p.x) * (1.0 - p.z)) + xi_integrand / (1.0 - p.z) + xip_integrand / (1.0 - p.x) + xi_xip_integrand;
		// const double integrand = delta_integrand + (1.0 - x) * xi_integrand + (1.0 - z) * xip_integrand + (1.0 - x) * (1.0 - z) * xi_xip_integrand;
		return integrand;
	}
}

#endif