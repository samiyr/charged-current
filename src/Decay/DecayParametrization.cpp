#ifndef DECAY_PARAMETRIZATION_H
#define DECAY_PARAMETRIZATION_H

#include "Utility/Math.cpp"

struct DecayParametrization {
	static DecayParametrization fit1() {
		return DecayParametrization(7.365, 1.4, 2.276, 2.04);
	}
	const double N;
	const double alpha;
	const double beta;
	const double gamma;

	const double beta_term_1p_alpha_beta;
	const double beta_term_2p_alpha_beta;

	const double gamma_prefactor_term;

	constexpr DecayParametrization() : N(0.0), alpha(0.0), beta(0.0), gamma(0.0), beta_term_1p_alpha_beta(0.0), beta_term_2p_alpha_beta(0.0), gamma_prefactor_term(0.0) {}

	constexpr DecayParametrization(const double _N,
	const double _alpha,
	const double _beta,
	const double _gamma)
	: N(_N),
	alpha(_alpha),
	beta(_beta),
	gamma(_gamma),
	beta_term_1p_alpha_beta(Math::beta(1 + alpha, 1 + beta)),
	beta_term_2p_alpha_beta(Math::beta(2 + alpha, 1 + beta)),
	gamma_prefactor_term(std::pow(gamma, -1 - alpha))
	{ }
};


#endif