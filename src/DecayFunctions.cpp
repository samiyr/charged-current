#ifndef DECAY_FUNCTIONS_H
#define DECAY_FUNCTIONS_H

#include "Utility.cpp"
#include "DecayParametrization.cpp"

namespace DecayFunctions {
	constexpr double decay_function(const double x, const double z, const double Q2, const double z_min, const DecayParametrization &decay) {
		const double alpha = decay.alpha;
		const double beta = decay.beta;
		const double gamma = decay.gamma;
		const double N = decay.N;
		const double mD = decay.resonance_mass;
		const double rho_min = z_min / z;

		const double h0 = (z * Q2) / (2 * x * decay.hadron_mass);
		const double hv = std::sqrt(h0 * h0 - mD * mD);

		const double a = h0 * h0 / (mD * mD);
		const double b = h0 * hv / (mD * mD);

		const double beta_arg_min = (a - b) * gamma * rho_min;
		if (beta_arg_min > 1) { 
			return 0.0;
		}
		
		const double prefactor = (N * mD * decay.gamma_prefactor_term * h0) / (2 * decay.total_decay_width * hv * 8 * M_PI * M_PI);
		const double beta_1 = -rho_min * decay.beta_term_1p_alpha_beta;
		const double beta_2 = rho_min * Utility::incomplete_beta(beta_arg_min, 1 + alpha, 1 + beta);
		const double beta_3 = 1.0 / (gamma * (a - b)) * (decay.beta_term_2p_alpha_beta - Utility::incomplete_beta(beta_arg_min, 2 + alpha, 1 + beta));

		const double result = prefactor * (beta_1 + beta_2 + beta_3);

		return result;
	}

	constexpr double decay_function_integrand(double input[], size_t dim, void *params_in) {
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
		const double integrand = decay * rho * h0 * h0 / (gamma_tot * 2 * mD * 8 * M_PI * M_PI);

		return integrand;
	};
}

#endif