#ifndef DECAY_FUNCTIONS_H
#define DECAY_FUNCTIONS_H

#include <concepts>
#include <numbers>

#include "Common/Particle.cpp"

#include "Utility/Math.cpp"

#include "Decay/DecayParametrization.cpp"

#include "Integration/Integrator.cpp"

template <typename T>
concept is_decay_function = requires(
	T decay_function, const double x, const double z, const double Q2, const double z_min, 
	const DecayParametrization &parametrization, const Particle &resonance, const Particle &hadron, const Particle &lepton) {
	{ decay_function(x, z, Q2, z_min, parametrization, resonance, hadron, lepton) } -> std::same_as<double>;
};

namespace DecayFunctions {
	inline auto trivial = [](
		[[maybe_unused]] const double x, 
		[[maybe_unused]] const double z, 
		[[maybe_unused]] const double Q2, 
		[[maybe_unused]] const double z_min, 
		[[maybe_unused]] const DecayParametrization &decay, 
		[[maybe_unused]] const Particle &resonance, 
		[[maybe_unused]] const Particle &hadron,
		[[maybe_unused]] const Particle &lepton) noexcept { return 1.0; };

	constexpr double decay_function(
		const double x, 
		const double z, 
		const double Q2, 
		const double z_min, 
		const DecayParametrization &parametrization, 
		const Particle &resonance, 
		const Particle &hadron,
		[[maybe_unused]] const Particle &lepton) {

		if (z < z_min) { return 0.0; }
		const double alpha = parametrization.alpha;
		const double beta = parametrization.beta;
		const double gamma = parametrization.gamma;
		const double N = parametrization.N;
		const double mD = resonance.mass;
		const double rho_min = z_min / z;

		const double h0 = (z * Q2) / (2 * x * hadron.mass);
		const double hv = std::sqrt(h0 * h0 - mD * mD);

		const double a = h0 * h0 / (mD * mD);
		const double b = h0 * hv / (mD * mD);

		const double beta_arg_min = (a - b) * gamma * rho_min;
		if (beta_arg_min > 1) { 
			return 0.0;
		}
		
		const double prefactor = 2 * std::numbers::pi * (N * mD * parametrization.gamma_prefactor_term * h0) / (2 * hv);
		const double beta_1 = -rho_min * parametrization.beta_term_1p_alpha_beta;
		const double beta_2 = rho_min * Math::incomplete_beta(beta_arg_min, 1 + alpha, 1 + beta);
		const double beta_3 = 1.0 / (gamma * (a - b)) * (parametrization.beta_term_2p_alpha_beta - Math::incomplete_beta(beta_arg_min, 2 + alpha, 1 + beta));

		const double result = prefactor * resonance.lifetime * (beta_1 + beta_2 + beta_3);

		return result;
	}

	// double decay_function(
	// 	const double x, 
	// 	const double z, 
	// 	const double Q2, 
	// 	const double E_min, 
	// 	const DecayParametrization &parametrization, 
	// 	const Particle &resonance,
	// 	const Particle &target,
	// 	const Particle &lepton) {

	// 	Integrator integrator([&](double input[], size_t, void *) {
	// 		const double rho = input[0];
	// 		const double c = input[1];

	// 		const double h0 = z * Q2 / (2.0 * x * target.mass);
	// 		if (h0 < resonance.mass) { return 0.0; }
	// 		const double pp0 = rho * h0;
	// 		if (pp0 < E_min) { return 0.0; }

	// 		const double mu = lepton.mass / h0;
	// 		const double reduced_rho = sqrt(std::pow(rho, 2) - std::pow(mu, 2));

	// 		const double a = std::pow(h0, 2) / std::pow(resonance.mass, 2);
	// 		const double b = h0 * sqrt(std::pow(h0, 2) - std::pow(resonance.mass, 2)) / std::pow(resonance.mass, 2);

	// 		const double cos_min = (parametrization.gamma * a * rho - 1.0) / (parametrization.gamma * b * reduced_rho);
	// 		const double cos = cos_min + (1.0 - cos_min) * c;

	// 		return (1.0 - cos_min) * DecayFunctions::differential_decay_function(cos, rho, z, x, Q2, E_min, parametrization, resonance, target, lepton);
	// 	}, {0.0, 0.0}, {1.0, 1.0}, nullptr, IntegrationMethod::CubaSuave);
	// 	const auto result = integrator.integrate();
	// 	return result.value;
	// }
}

#endif