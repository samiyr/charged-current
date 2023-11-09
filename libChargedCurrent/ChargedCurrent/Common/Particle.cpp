#ifndef PARTICLE_H
#define PARTICLE_H

#include <limits>
#include <iostream>

// Container for the following particle data: mass, lifetime.
struct Particle {
	// Mass of the particle, in units of GeV.
	double mass;
	// Lifetime of the particle, in units of 1e-12 seconds.
	double lifetime;

	// Constructs an empty particle, with zero mass and zero lifetime.
	constexpr Particle() noexcept : mass(0.0), lifetime(0.0) {}
	/// @brief Constructs a particle with the given mass and lifetime.
	/// @param mass Mass of the particle in GeV.
	/// @param lifetime Lifetime of the particle, in 1e-12 seconds. Defaults to infinity, i.e. a stable particle.
	constexpr Particle(const double _mass, const double _lifetime = std::numeric_limits<double>::infinity()) noexcept : mass(_mass), lifetime(_lifetime) {}

	friend std::ostream& operator<<(std::ostream &os, const Particle &o) {
		return os << "m = " << o.mass << " | " << "tau = " << o.lifetime;
	}
};

#endif