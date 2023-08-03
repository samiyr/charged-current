#ifndef DECAY_FUNCTIONS_H
#define DECAY_FUNCTIONS_H

#include "Utility/Math.cpp"
#include "Decay/Decay.cpp"
#include <any>
#include <concepts>

namespace DecayFunctions {
	template <typename T>
	concept Concept = requires(T decay_function, const double x, const double z, const double Q2, const double z_min, const DecayParametrization &parametrization, const Particle &resonance, const Particle &hadron) {
		{ decay_function(x, z, Q2, z_min, parametrization, resonance, hadron) } -> std::same_as<double>;
	};

	static auto trivial = [](
		[[maybe_unused]] const double x, 
		[[maybe_unused]] const double z, 
		[[maybe_unused]] const double Q2, 
		[[maybe_unused]] const double z_min, 
		[[maybe_unused]] const DecayParametrization &decay, 
		[[maybe_unused]] const Particle &resonance, 
		[[maybe_unused]] const Particle &hadron) noexcept { return 1.0; };

	constexpr double decay_function(
		const double x, 
		const double z, 
		const double Q2, 
		const double z_min, 
		const DecayParametrization &parametrization, 
		const Particle &resonance, 
		const Particle &hadron) {

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
}

#endif