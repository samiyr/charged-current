#ifndef DECAY_FUNCTIONS_H
#define DECAY_FUNCTIONS_H

#include <concepts>
#include <numbers>
#include <filesystem>

#include "Common/Particle.cpp"

#include "Utility/Math.cpp"

#include "Decay/DecayParametrization.cpp"

#include "Interpolation/InterpolatingFunction.cpp"

template <typename T>
concept is_decay_function = requires(
	T decay_function, const double x, const double z, const double Q2, const double min, 
	const DecayParametrization &parametrization, const Particle &resonance, const Particle &target) {
	{ decay_function(x, z, Q2, min, parametrization, resonance, target) } -> std::same_as<double>;
};

namespace DecayFunctions {
	inline auto trivial = [](
		[[maybe_unused]] const double x, 
		[[maybe_unused]] const double z, 
		[[maybe_unused]] const double Q2, 
		[[maybe_unused]] const double min, 
		[[maybe_unused]] const DecayParametrization &decay, 
		[[maybe_unused]] const Particle &resonance, 
		[[maybe_unused]] const Particle &target) noexcept { return 1.0; };

	constexpr double differential_decay_function_integrand(
		const double rho, const double reduced_rho, const double cos, const double a, const double b, const double h0,
		const DecayParametrization &parametrization, const Particle &resonance
	) {
		const double w = a * rho - b * reduced_rho * cos;

		if (w < 0.0 || w > 1.0 / parametrization.gamma) { return 0.0; }

		const double decay_value = parametrization.N * std::pow(w, parametrization.alpha) * std::pow(1.0 - parametrization.gamma * w, parametrization.beta);
		const double integrand = 2.0 * std::numbers::pi * resonance.lifetime * decay_value * reduced_rho * std::pow(h0, 2) / (2.0 * resonance.mass);

		return integrand;
	}
	
	constexpr double differential_decay_function(
		const double cos, const double rho, const double zyE, const double E_min,
		const DecayParametrization &parametrization,
		const Particle &resonance, [[maybe_unused]] const Particle &target, const Particle &lepton
	) {
		const double h0 = zyE;
		if (h0 < resonance.mass) { return 0.0; }
		const double pp0 = rho * h0;
		if (pp0 < E_min) { return 0.0; }

		const double mu = lepton.mass / h0;
		const double reduced_rho = sqrt(std::pow(rho, 2) - std::pow(mu, 2));

		const double a = std::pow(h0, 2) / std::pow(resonance.mass, 2);
		const double b = h0 * sqrt(std::pow(h0, 2) - std::pow(resonance.mass, 2)) / std::pow(resonance.mass, 2);

		return differential_decay_function_integrand(rho, reduced_rho, cos, a, b, h0, parametrization, resonance);
	}

	constexpr double differential_decay_function(
		const double cos, const double rho, const double z, const double x, const double Q2, const double E_min,
		const DecayParametrization &parametrization,
		const Particle &resonance, const Particle &target, const Particle &lepton
	) {
		const double zyE = z * Q2 / (2.0 * x * target.mass);
		return differential_decay_function(cos, rho, zyE, E_min, parametrization, resonance, target, lepton);
	}

	constexpr double decay_function(
		const double x, 
		const double z, 
		const double Q2, 
		const double z_min, 
		const DecayParametrization &parametrization, 
		const Particle &resonance, 
		const Particle &target) {

		if (z < z_min) { return 0.0; }
		const double alpha = parametrization.alpha;
		const double beta = parametrization.beta;
		const double gamma = parametrization.gamma;
		const double N = parametrization.N;
		const double mD = resonance.mass;
		const double rho_min = z_min / z;

		const double h0 = (z * Q2) / (2 * x * target.mass);
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

	constexpr double decay_function_integrand(double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
		const double rho = input[0];
		const double cos = input[1];

		std::vector<double> &params = *static_cast<std::vector<double> *>(params_in);
		const double x = params[0];
		const double z = params[1];
		const double Q2 = params[2];
		const double mN = params[3];
		const double mD = params[4];
		const double n = params[5];
		const double alpha = params[6];
		const double beta = params[7];
		const double gamma = params[8];
		const double gamma_tot = params[9];

		const double h0 = z * Q2 / (2 * x * mN);
		if (h0 < mD) { return 0.0; }
		const double a = (h0 * h0) / (mD * mD);
		const double b = h0 * std::sqrt(h0 * h0 - mD * mD) / (mD * mD);
		const double w = a * rho - b * rho * cos;

		if (w < 0 || w > 1.0 / gamma) { return 0.0; }

		const double decay = n * std::pow(w, alpha) * std::pow(1.0 - gamma * w, beta);
		const double integrand = 2 * std::numbers::pi * decay * rho * h0 * h0 / (gamma_tot * 2 * mD);

		return integrand;
	};

	struct DecayGrid {
		InterpolatingFunction function;

		DecayGrid(const std::filesystem::path grid_path) : function(grid_path) {}

		double operator()(
			const double x, 
			const double z, 
			const double Q2, 
			[[maybe_unused]] const double min, 
			[[maybe_unused]] const DecayParametrization &parametrization, 
			[[maybe_unused]] const Particle &resonance, 
			const Particle &target
		) const {
			const double zyE = z * Q2 / (2.0 * x * target.mass);
			return function(zyE);
		}
		double operator()(const double zyE) const {
			return function(zyE);
		}
	};

	static inline DecayGrid decay_grid(const std::filesystem::path grid_path) {
		return DecayGrid(grid_path);
	}
}

#endif