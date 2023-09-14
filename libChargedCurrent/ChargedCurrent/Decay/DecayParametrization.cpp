#ifndef DECAY_PARAMETRIZATION_H
#define DECAY_PARAMETRIZATION_H

#include "Utility/Math.cpp"
#include "Utility/Utility.cpp"


struct DecayParametrization {
	using Set = std::vector<DecayParametrization>;

	constexpr static DecayParametrization fit1() noexcept {
		return DecayParametrization(7.365, 1.4, 2.276, 2.04);
	}
	constexpr static DecayParametrization fit2() noexcept {
		return DecayParametrization(4.62698, 1.17383, 2.06030, 2.05650);
	}

	consteval static DecayParametrization::Set fits() noexcept {
		const double N 		= 4.62698;
		const double alpha 	= 1.17383;
		const double beta 	= 2.06030;
		const double gamma 	= 2.05650;

		const double N_sigma 		= 1.61876;
		const double alpha_sigma 	= 0.149058;
		const double beta_sigma 	= 0.273986;
		const double gamma_sigma 	= 0.0580984;

		const std::vector<std::vector<int>> multiplicative_factors{
			{-1, -1, -1, -1}, {-1, -1, -1, 0}, {-1, -1, -1, 1}, 
			{-1, -1, 0, -1}, {-1, -1, 0, 0}, {-1, -1, 0, 1}, 
			{-1, -1, 1, -1}, {-1, -1, 1, 0}, {-1, -1, 1, 1}, 
			{-1, 0, -1, -1}, {-1, 0, -1, 0}, {-1, 0, -1, 1}, 
			{-1, 0, 0, -1}, {-1, 0, 0, 0}, {-1, 0, 0, 1}, 
			{-1, 0, 1, -1}, {-1, 0, 1, 0}, {-1, 0, 1, 1}, 
			{-1, 1, -1, -1}, {-1, 1, -1, 0}, {-1, 1, -1, 1}, 
			{-1, 1, 0, -1}, {-1, 1, 0, 0}, {-1, 1, 0, 1}, 
			{-1, 1, 1, -1}, {-1, 1, 1, 0}, {-1, 1, 1, 1}, 
			{0, -1, -1, -1}, {0, -1, -1, 0}, {0, -1, -1, 1}, 
			{0, -1, 0, -1}, {0, -1, 0, 0}, {0, -1, 0, 1}, 
			{0, -1, 1, -1}, {0, -1, 1, 0}, {0, -1, 1, 1}, 
			{0, 0, -1, -1}, {0, 0, -1, 0}, {0, 0, -1, 1}, 
			{0, 0, 0, -1}, {0, 0, 0, 0}, {0, 0, 0, 1}, 
			{0, 0, 1, -1}, {0, 0, 1, 0}, {0, 0, 1, 1}, 
			{0, 1, -1, -1}, {0, 1, -1, 0}, {0, 1, -1, 1}, 
			{0, 1, 0, -1}, {0, 1, 0, 0}, {0, 1, 0, 1}, 
			{0, 1, 1, -1}, {0, 1, 1, 0}, {0, 1, 1, 1}, 
			{1, -1, -1, -1}, {1, -1, -1, 0}, {1, -1, -1, 1}, 
			{1, -1, 0, -1}, {1, -1, 0, 0}, {1, -1, 0, 1}, 
			{1, -1, 1, -1}, {1, -1, 1, 0}, {1, -1, 1, 1}, 
			{1, 0, -1, -1}, {1, 0, -1, 0}, {1, 0, -1, 1}, 
			{1, 0, 0, -1}, {1, 0, 0, 0}, {1, 0, 0, 1}, 
			{1, 0, 1, -1}, {1, 0, 1, 0}, {1, 0, 1, 1}, 
			{1, 1, -1, -1}, {1, 1, -1, 0}, {1, 1, -1, 1}, 
			{1, 1, 0, -1}, {1, 1, 0, 0}, {1, 1, 0, 1}, 
			{1, 1, 1, -1}, {1, 1, 1, 0}, {1, 1, 1, 1}
		};

		std::vector<DecayParametrization> parametrizations;

		for (const std::vector<int> &set : multiplicative_factors) {
			parametrizations.emplace_back(N + set[0] * N_sigma, alpha + set[1] * alpha_sigma, beta + set[2] * beta_sigma, gamma + set[3] * gamma_sigma);
		}

		return parametrizations;
	}

	const double N;
	const double alpha;
	const double beta;
	const double gamma;

	const double beta_term_1p_alpha_beta;
	const double beta_term_2p_alpha_beta;

	const double gamma_prefactor_term;

	constexpr DecayParametrization() noexcept
		: N(0.0), alpha(0.0), beta(0.0), gamma(0.0), beta_term_1p_alpha_beta(0.0), beta_term_2p_alpha_beta(0.0), gamma_prefactor_term(0.0) {}

	constexpr DecayParametrization(const double _N,
	const double _alpha,
	const double _beta,
	const double _gamma) noexcept
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