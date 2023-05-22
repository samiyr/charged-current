#ifndef DECAY_H
#define DECAY_H

#include "DecayParametrization.cpp"

template <typename DecayFunction>
struct Decay {
	const DecayParametrization parametrization;
	const DecayFunction decay_function;

	Decay(const DecayParametrization _parametrization, const DecayFunction _decay_function) : parametrization(_parametrization), decay_function(_decay_function) {}

	double operator()(const double x, const double z, const double Q2) {
		return decay_function(x, z, Q2, parametrization);
	}
};

#endif