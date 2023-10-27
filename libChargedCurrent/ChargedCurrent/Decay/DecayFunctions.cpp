#ifndef DECAY_FUNCTIONS_H
#define DECAY_FUNCTIONS_H

#include <concepts>
#include <numbers>

#include "Common/Particle.cpp"

#include "Utility/Math.cpp"

#include "Decay/DecayParametrization.cpp"

template <typename T>
concept is_decay_function = requires(
	T decay_function, const double cos, const double rho, const double z, const double x, const double Q2,
	const DecayParametrization &parametrization,
	const Particle &resonance, const Particle &target, const Particle &lepton) {
	{ decay_function(cos, rho, z, x, Q2, parametrization, resonance, target, lepton) } -> std::same_as<double>;
};

namespace DecayFunctions {
	inline auto trivial = [](
		const double,
		const double,
		const double,
		const double, 
		const double,
		const DecayParametrization&,
		const Particle&, const Particle&, const Particle&) noexcept { return 1.0; };

	constexpr double decay_function(
		const double cos, const double rho, const double z, const double x, const double Q2,
		const DecayParametrization &parametrization,
		const Particle &resonance, const Particle &target, const Particle &lepton
	) {
		const double h0 = z * Q2 / (2.0 * x * target.mass);
		if (h0 < resonance.mass) { return 0.0; }

		const double mu = lepton.mass / h0;
		const double reduced_rho = sqrt(std::pow(rho, 2) - std::pow(mu, 2));

		const double a = std::pow(h0, 2) / std::pow(resonance.mass, 2);
		const double b = h0 * sqrt(std::pow(h0, 2) - std::pow(resonance.mass, 2)) / std::pow(resonance.mass, 2);
		const double w = a * rho - b * reduced_rho * cos;

		if (w < 0.0 || w > 1.0 / parametrization.gamma) { return 0.0; }

		const double decay_value = parametrization.N * std::pow(w, parametrization.alpha) * std::pow(1.0 - parametrization.gamma * w, parametrization.beta);
		const double integrand = 2.0 * std::numbers::pi * resonance.lifetime * decay_value * reduced_rho * std::pow(h0, 2) / (2.0 * resonance.mass);

		return integrand;
	}

	// constexpr double decay_function_integrand(double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
	// 	const double rho = input[0];
	// 	const double cos = input[1];

	// 	std::vector<double> &params = *static_cast<std::vector<double> *>(params_in);
	// 	const double x = params[0];
	// 	const double z = params[1];
	// 	const double Q2 = params[2];
	// 	const double mN = params[3];
	// 	const double mD = params[4];
	// 	const double n = params[5];
	// 	const double alpha = params[6];
	// 	const double beta = params[7];
	// 	const double gamma = params[8];
	// 	const double gamma_tot = params[9];

	// 	const double h0 = z * Q2 / (2 * x * mN);
	// 	if (h0 < mD) { return 0.0; }
	// 	const double a = (h0 * h0) / (mD * mD);
	// 	const double b = h0 * std::sqrt(h0 * h0 - mD * mD) / (mD * mD);
	// 	const double w = a * rho - b * rho * cos;

	// 	if (w < 0 || w > 1.0 / gamma) { return 0.0; }

	// 	const double decay = n * std::pow(w, alpha) * std::pow(1.0 - gamma * w, beta);
	// 	const double integrand = 2 * std::numbers::pi * decay * rho * h0 * h0 / (gamma_tot * 2 * mD);

	// 	return integrand;
	// };
}

#endif