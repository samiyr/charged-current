#ifndef PARTICLE_H
#define PARTICLE_H

#include <limits>

struct Particle {
	const double mass;
	const double lifetime;

	constexpr Particle() : mass(0.0), lifetime(0.0) {}
	constexpr Particle(const double _mass, const double _lifetime = std::numeric_limits<double>::infinity()) : mass(_mass), lifetime(_lifetime) {}
};

#endif