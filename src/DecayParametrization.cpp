#ifndef DECAY_PARAMETRIZATION_H
#define DECAY_PARAMETRIZATION_H

struct DecayParametrization {
	const double N;
	const double alpha;
	const double beta;
	const double gamma;

	const double resonance_mass;
	const double hadron_mass;
	const double total_decay_width;
	
	const double z_min;

	DecayParametrization(const double _N,
	const double _alpha,
	const double _beta,
	const double _gamma,
	const double _resonance_mass,
	const double _hadron_mass,
	const double _total_decay_width,
	const double _z_min)
	: N(_N),
	alpha(_alpha),
	beta(_beta),
	gamma(_gamma),
	resonance_mass(_resonance_mass),
	hadron_mass(_hadron_mass),
	total_decay_width(_total_decay_width),
	z_min(_z_min)
	{ }
};


#endif