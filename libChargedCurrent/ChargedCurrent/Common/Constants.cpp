#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <limits>

#include "Common/Particle.cpp"

#include "Utility/Utility.cpp"
#include "Utility/Math.cpp"

// Provides access to certain physical constants.
namespace Constants {
	// QCD color factor (fundamental representation).
	constexpr inline double C_F = 4.0 / 3.0;
	// QCD color factor (adjoint representation).
	constexpr inline double C_A = 3.0;
	// Number of quark flavors.
	constexpr inline double N_f = 6.0;
	// QCD normalization of tr(t^a t^b) = T_R δ^ab in the fundamental representation.
	constexpr inline double T_R = 0.5;
	// Fermi coupling constant.
	constexpr inline double fermi_coupling = 1.1663787e-5;
	constexpr inline double lambda_QCD = 0.226;

	// The squared CKM matrix elements.

	constexpr inline double V_ud = Math::pow2(0.97435);
	constexpr inline double V_us = Math::pow2(0.22500);
	constexpr inline double V_ub = Math::pow2(0.00369);

	constexpr inline double V_cd = Math::pow2(0.22486);
	constexpr inline double V_cs = Math::pow2(0.97349);
	constexpr inline double V_cb = Math::pow2(0.04182);

	constexpr inline double V_td = Math::pow2(0.00857);
	constexpr inline double V_ts = Math::pow2(0.04110);
	constexpr inline double V_tb = Math::pow2(0.999112);

	namespace Particles {
		constexpr inline Particle Proton = Particle(0.938272);
		constexpr inline Particle D0 = Particle(1.86484, 0.4101);
		constexpr inline Particle Dp = Particle(1.86962, 1.04);
		constexpr inline Particle Ds = Particle(1.96847, 0.5);
		constexpr inline Particle LambdaC = Particle(2.28646, 0.2);
		constexpr inline Particle W = Particle(80.377, 4.7961631e-9);
		constexpr inline Particle Neutrino = Particle(0.0);
		constexpr inline Particle Muon = Particle(0.105658372);
	}

	namespace Charm {
		constexpr inline double Mass = 1.3;
	}

	// Value of the Riemann zeta function ζ(3).
	constexpr inline double zeta3 = 1.20205690315959428539973816151;
} 

#endif