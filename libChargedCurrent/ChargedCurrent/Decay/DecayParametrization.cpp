#ifndef DECAY_PARAMETRIZATION_H
#define DECAY_PARAMETRIZATION_H

#include "Utility/Math.cpp"
#include "Utility/Utility.cpp"

struct DecayParametrization {
	constexpr static DecayParametrization fit1() noexcept {
		return DecayParametrization(7.365, 1.4, 2.276, 2.04);
	}
	constexpr static DecayParametrization fit2() noexcept {
		return DecayParametrization(4.62698, 1.17383, 2.06030, 2.05650, 1.61876, 0.149058, 0.273986, 0.0580984);
	}

	double N;
	double alpha;
	double beta;
	double gamma;

	double N_sigma;
	double alpha_sigma;
	double beta_sigma;
	double gamma_sigma;

	double beta_term_1p_alpha_beta;
	double beta_term_2p_alpha_beta;

	double gamma_prefactor_term;

	constexpr DecayParametrization() noexcept
	: N(0.0), alpha(0.0), beta(0.0), gamma(0.0), 
	N_sigma(0.0), alpha_sigma(0.0), beta_sigma(0.0), gamma_sigma(0.0),
	beta_term_1p_alpha_beta(0.0), beta_term_2p_alpha_beta(0.0), gamma_prefactor_term(0.0) {}

	constexpr DecayParametrization(
		const double N,
		const double alpha,
		const double beta,
		const double gamma
	) noexcept : DecayParametrization(N, alpha, beta, gamma, 0.0, 0.0, 0.0, 0.0) { }

	constexpr DecayParametrization(
		const double _N,
		const double _alpha,
		const double _beta,
		const double _gamma,
		const double N_sigma,
		const double alpha_sigma,
		const double beta_sigma,
		const double gamma_sigma
	) noexcept
	: N(_N),
	alpha(_alpha),
	beta(_beta),
	gamma(_gamma),
	N_sigma(N_sigma),
	alpha_sigma(alpha_sigma),
	beta_sigma(beta_sigma),
	gamma_sigma(gamma_sigma),
	beta_term_1p_alpha_beta(Math::beta(1 + alpha, 1 + beta)),
	beta_term_2p_alpha_beta(Math::beta(2 + alpha, 1 + beta)),
	gamma_prefactor_term(std::pow(gamma, -1 - alpha))
	{ }

};


#endif