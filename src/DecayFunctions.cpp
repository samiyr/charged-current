#ifndef DECAY_FUNCTIONS_H
#define DECAY_FUNCTIONS_H

#include "Utility.cpp"
#include "DecayParametrization.cpp"

namespace DecayFunctions {
	double decay_function(const double x, const double z, const double Q2, const DecayParametrization &decay) {
		const double alpha = decay.alpha;
		const double beta = decay.beta;
		const double gamma = decay.gamma;
		const double N = decay.N;
		const double m = decay.resonance_mass;

		const double resonance_energy = (z * Q2) / (2 * x * decay.hadron_mass);

		const double pmin = decay.pmin;
		const double pmax = resonance_energy;

		const double a_tilde = resonance_energy  / (m * m);
		const double b_tilde = std::sqrt(resonance_energy * resonance_energy - m * m) / (m * m);

		const double prefactor = (N * std::pow(gamma, -1 - alpha)) / (8 * M_PI * M_PI * decay.total_decay_width * 2 * m);
		const double gamma_term = gamma * (a_tilde - b_tilde);

		const double beta_1_max = Utility::incomplete_beta(gamma_term * pmax, 2 + alpha, 1 + beta);
		const double beta_1_min = Utility::incomplete_beta(gamma_term * pmin, 2 + alpha, 1 + beta);
		const double beta_2_max = Utility::incomplete_beta(gamma_term * pmax, 1 + alpha, 1 + beta);
		const double beta_2_min = Utility::incomplete_beta(gamma_term * pmin, 1 + alpha, 1 + beta);

		const double reflection_term = decay.reflection * (pmax - pmin);

		const double result = prefactor * ((beta_1_max - beta_1_min) / gamma_term - beta_2_max + beta_2_min - reflection_term) / b_tilde;
		return result;
	}
}

#endif