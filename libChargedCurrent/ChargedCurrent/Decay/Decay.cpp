#ifndef DECAY_H
#define DECAY_H

#include "Decay/DecayParametrization.cpp"
#include "Decay/DecayFunctions.cpp"

#include "Common/Particle.cpp"

template <is_decay_function DecayFunction = decltype(DecayFunctions::trivial)>
struct Decay {
	DecayParametrization parametrization;
	const Particle resonance;
	const Particle target;
	const Particle lepton;
	const DecayFunction decay_function;

	const double minimum_lepton_energy;
	const double z_min_cutoff;

	constexpr Decay(
		const DecayParametrization parametrization, 
		const Particle resonance, 
		const Particle target,
		const Particle lepton, 
		const DecayFunction decay_function, 
		const double minimum_lepton_energy, 
		const double z_min_cutoff = 0.0) noexcept
	: parametrization(parametrization), 
	resonance(resonance), 
	target(target),
	lepton(lepton), 
	decay_function(decay_function), 
	minimum_lepton_energy(minimum_lepton_energy), 
	z_min_cutoff(z_min_cutoff) {}

	constexpr Decay() noexcept : Decay(DecayParametrization(), Particle(), Particle(), Particle(), DecayFunctions::trivial, 0.0) { }

	constexpr double operator()(const double x, const double z, const double Q2) const noexcept(noexcept(decay_function)) {
		if (x != prev_x || z != prev_z || Q2 != prev_Q2) {
			prev_x = x;
			prev_z = z;
			prev_Q2 = Q2;

			prev_value = decay_function(x, z, Q2, minimum_lepton_energy, parametrization, resonance, target, lepton);
		}
		return prev_value;
	}

	private:
	mutable double prev_x = -1.0;
	mutable double prev_z = -1.0;
	mutable double prev_Q2 = -1.0;
	mutable double prev_value = 0.0;
};

// constexpr static auto TrivialDecay = Decay(DecayParametrization(), Particle(), Particle(), DecayFunctions::trivial, 0.0);

#endif