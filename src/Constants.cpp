#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "Utility.cpp"
#include "Particle.cpp"
#include <limits>

namespace Constants {
	constexpr double C_F = 4.0 / 3.0;
	constexpr double T_R = 0.5;
	constexpr double fermi_coupling = 1.1663787e-5;
	constexpr double lambda_QCD = 0.226;

	constexpr static double V_ud = POW2(0.97435);
	constexpr static double V_us = POW2(0.22500);
	constexpr static double V_ub = POW2(0.00369);

	constexpr static double V_cd = POW2(0.22486);
	constexpr static double V_cs = POW2(0.97349);
	constexpr static double V_cb = POW2(0.04182);

	// const static double V_td = POW2(0);
	// const static double V_ts = POW2(0);
	// const static double V_tb = POW2(0);

	constexpr static double V_td = POW2(0.00857);
	constexpr static double V_ts = POW2(0.04110);
	constexpr static double V_tb = POW2(0.999112);

	// constexpr double ckm_squared_matrix_elements[3][3] = {
	// 	/*
	// 	  	V_ud				V_us				V_ub
	// 	  	V_cd				V_cs				V_cb
	// 	  	V_td				V_ts				V_tb
	// 	*/
	// 	{POW2(0.97435), 	POW2(0.22500), 		POW2(0.00369)},
	// 	{POW2(0.22486),		POW2(0.97349),		POW2(0.04182)},
	// 	{POW2(0.00857),		POW2(0.04110),		POW2(0.999112)}
	// };
	// constexpr double ckm_squared_matrix_elements[3][3] = {
	// 	/*
	// 	  	V_ud				V_us				V_ub
	// 	  	V_cd				V_cs				V_cb
	// 	  	V_td				V_ts				V_tb
	// 	*/
	// 	{POW2(0.97435), 	0, 		0},
	// 	{0,		POW2(0.97349),		0},
	// 	{0,		0,		0}
	// };

	namespace Particles {
		constexpr static Particle Proton = Particle { .mass = 0.938272, .lifetime = std::numeric_limits<double>::infinity() };
		constexpr static Particle D0 = Particle { .mass = 1.86484, .lifetime = 0.4101 };
		constexpr static Particle Dp = Particle { .mass = 1.86962, .lifetime = 1.04 };
		constexpr static Particle Ds = Particle { .mass = 1.96847, .lifetime = 0.5 };
		constexpr static Particle LambdaC = Particle { .mass = 2.28646, .lifetime = 0.2 };
		constexpr static Particle W = Particle { .mass = 80.377, .lifetime = 4.7961631e-9 };
	}

	namespace Charm {
		constexpr static double Mass = 1.5;
	}
	namespace Bottom {
		constexpr static double Mass = 5.0;
	}
} 

#endif