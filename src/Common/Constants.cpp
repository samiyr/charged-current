#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "Utility/Utility.cpp"
#include "Common/Particle.cpp"
#include <limits>
#include "Utility/Math.cpp"

namespace Constants {
	constexpr double C_F = 4.0 / 3.0;
	constexpr double T_R = 0.5;
	constexpr double fermi_coupling = 1.1663787e-5;
	constexpr double lambda_QCD = 0.226;

	constexpr static double V_ud = Math::pow2(0.97435);
	constexpr static double V_us = Math::pow2(0.22500);
	constexpr static double V_ub = Math::pow2(0.00369);

	constexpr static double V_cd = Math::pow2(0.22486);
	constexpr static double V_cs = Math::pow2(0.97349);
	constexpr static double V_cb = Math::pow2(0.04182);

	constexpr static double V_td = Math::pow2(0.00857);
	constexpr static double V_ts = Math::pow2(0.04110);
	constexpr static double V_tb = Math::pow2(0.999112);

	namespace Particles {
		constexpr static Particle Proton = Particle(0.938272);
		constexpr static Particle D0 = Particle(1.86484, 0.4101);
		constexpr static Particle Dp = Particle(1.86962, 1.04);
		constexpr static Particle Ds = Particle(1.96847, 0.5);
		constexpr static Particle LambdaC = Particle(2.28646, 0.2);
		constexpr static Particle W = Particle(80.377, 4.7961631e-9);
		constexpr static Particle Neutrino = Particle(0.0);
	}

	namespace Charm {
		constexpr static double Mass = 1.3;
	}
} 

#endif