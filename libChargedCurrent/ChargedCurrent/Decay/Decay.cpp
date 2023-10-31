#ifndef DECAY_H
#define DECAY_H

#include "Decay/DecayParametrization.cpp"
#include "Decay/DecayFunctions.cpp"

#include "Common/Particle.cpp"

template <is_decay_function DecayFunction = decltype(DecayFunctions::trivial)>
struct Decay {
	DecayParametrization parametrization;
	const Particle resonance;
	const Particle hadron;
	const DecayFunction decay_function;

	const double lepton_momentum_min;
	const double z_min_cutoff;

	constexpr Decay(
		const DecayParametrization _parametrization, 
		const Particle _resonance, 
		const Particle _hadron, 
		const DecayFunction _decay_function, 
		const double _lepton_momentum_min, 
		const double _z_min_cutoff = 0.0) noexcept
	: parametrization(_parametrization), 
	resonance(_resonance), 
	hadron(_hadron), 
	decay_function(_decay_function), 
	lepton_momentum_min(_lepton_momentum_min), 
	z_min_cutoff(_z_min_cutoff) {}

	constexpr Decay() noexcept : Decay(DecayParametrization(), Particle(), Particle(), DecayFunctions::trivial, 0.0) { }

	constexpr double operator()(const double x, const double z, const double Q2, const double z_min) const noexcept(noexcept(decay_function)) {
		if (x != prev_x || z != prev_z || Q2 != prev_Q2 || z_min != prev_z_min) {
			prev_x = x;
			prev_z = z;
			prev_Q2 = Q2;
			prev_z_min = z_min;
			
			// FIXME: rename lepton_momentum_min to min energy, implement z_min in DecayFunctions::analytical_decay_function
			prev_value = decay_function(x, z, Q2, lepton_momentum_min, parametrization, resonance, hadron, Particle());
		}
		return prev_value;
	}

	private:
	mutable double prev_x = -1.0;
	mutable double prev_z = -1.0;
	mutable double prev_Q2 = -1.0;
	mutable double prev_z_min = -1.0;
	mutable double prev_value = 0.0;
};

// constexpr static auto TrivialDecay = Decay(DecayParametrization(), Particle(), Particle(), DecayFunctions::trivial, 0.0);

#endif