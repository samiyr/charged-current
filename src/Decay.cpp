#ifndef DECAY_H
#define DECAY_H

#include "DecayParametrization.cpp"
#include "DecayFunctions.cpp"

template <typename DecayFunction>
struct Decay {
	const DecayParametrization parametrization;
	const DecayFunction decay_function;

	Decay(const DecayParametrization _parametrization, const DecayFunction _decay_function) : parametrization(_parametrization), decay_function(_decay_function) {}
	Decay(const DecayFunction _decay_function) : parametrization(), decay_function(_decay_function) {}

	constexpr double operator()(const double x, const double z, const double Q2, const double z_min) {
		if (x != prev_x || z != prev_z || Q2 != prev_Q2 || z_min != prev_z_min) {
			prev_x = x;
			prev_z = z;
			prev_Q2 = Q2;
			prev_z_min = z_min;

			// std::cout << x << ", " << z << ", " << Q2 << ", " << z_min << std::endl;

			prev_value = decay_function(x, z, Q2, z_min, parametrization);
		}
		return prev_value;
	}

	private:
	double prev_x;
	double prev_z;
	double prev_Q2;
	double prev_z_min;
	double prev_value;
};

const static auto TrivialDecay = Decay(DecayParametrization(), DecayFunctions::trivial);

#endif