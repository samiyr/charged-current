#ifndef DECAY_H
#define DECAY_H

#include "DecayParametrization.cpp"
#include "DecayFunctions.cpp"
#include "Particle.cpp"

template <typename DecayFunction>
struct Decay {
	const DecayParametrization parametrization;
	const Particle resonance;
	const Particle hadron;
	const DecayFunction decay_function;

	const double lepton_momentum_min;
	const double z_min_cutoff;

	Decay(const DecayParametrization _parametrization, const Particle _resonance, const Particle _hadron, const DecayFunction _decay_function, const double _lepton_momentum_min, const double _z_min_cutoff = 0.0)
	: parametrization(_parametrization), resonance(_resonance), hadron(_hadron), decay_function(_decay_function), lepton_momentum_min(_lepton_momentum_min), z_min_cutoff(_z_min_cutoff) {}

	// Decay(const DecayParametrization _parametrization, const DecayFunction _decay_function) : parametrization(_parametrization), decay_function(_decay_function) {}
	// Decay(const DecayFunction _decay_function) : parametrization(), decay_function(_decay_function) {}

	constexpr double operator()(const double x, const double z, const double Q2, const double z_min) const {
		if (x != prev_x || z != prev_z || Q2 != prev_Q2 || z_min != prev_z_min) {
			prev_x = x;
			prev_z = z;
			prev_Q2 = Q2;
			prev_z_min = z_min;

			// std::cout << x << ", " << z << ", " << Q2 << ", " << z_min << std::endl;

			prev_value = decay_function(x, z, Q2, z_min, parametrization, resonance, hadron);
		}
		return prev_value;
	}

	private:
	mutable double prev_x;
	mutable double prev_z;
	mutable double prev_Q2;
	mutable double prev_z_min;
	mutable double prev_value;
};

const static auto TrivialDecay = Decay(DecayParametrization(), Particle(), Particle(), DecayFunctions::trivial, 0.0);

#endif