#ifndef DECAY_PARAMETRIZATION_H
#define DECAY_PARAMETRIZATION_H

struct DecayParametrization {
	const double N;
	const double alpha;
	const double beta;
	const double gamma;

	const double resonance_mass;
	const double hadron_mass;
	const double lifetime;
	
	const double lepton_momentum_min;
	const double z_min_cutoff;

	const double beta_term_1p_alpha_beta;
	const double beta_term_2p_alpha_beta;

	const double gamma_prefactor_term;

	DecayParametrization() : N(0.0), alpha(0.0), beta(0.0), gamma(0.0), resonance_mass(0.0), hadron_mass(0.0), lifetime(0.0), lepton_momentum_min(0.0), z_min_cutoff(0.0), beta_term_1p_alpha_beta(0.0), beta_term_2p_alpha_beta(0.0), gamma_prefactor_term(0.0) {}

	DecayParametrization(const double _N,
	const double _alpha,
	const double _beta,
	const double _gamma,
	const double _resonance_mass,
	const double _hadron_mass,
	const double _lifetime,
	const double _lepton_momentum_min,
	const double _z_min_cutoff)
	: N(_N),
	alpha(_alpha),
	beta(_beta),
	gamma(_gamma),
	resonance_mass(_resonance_mass),
	hadron_mass(_hadron_mass),
	lifetime(_lifetime),
	lepton_momentum_min(_lepton_momentum_min),
	z_min_cutoff(_z_min_cutoff),
	beta_term_1p_alpha_beta(Utility::beta(1 + alpha, 1 + beta)),
	beta_term_2p_alpha_beta(Utility::beta(2 + alpha, 1 + beta)),
	gamma_prefactor_term(std::pow(gamma, -1 - alpha))
	{ }
};


#endif