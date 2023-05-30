#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "Utility.cpp"

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

	namespace D0 {
		constexpr static double Mass = 1.86484;
		// constexpr static double DecayWidth = 1.605e-12;
		constexpr static double Lifetime = 0.4101;
	}

	namespace Dp {
		constexpr static double Mass = 1.86962;
		// constexpr static double DecayWidth = 6.32896e-13;
		constexpr static double Lifetime = 1.04;
	}
	namespace Ds {
		constexpr static double Mass = 1.96847;
		// constexpr static double DecayWidth = 1.3164e-12;
		constexpr static double Lifetime = 0.5;
	}

	namespace Proton {
		constexpr static double Mass = 0.938272;
	}

	namespace WBoson {
		constexpr static double Mass = 80.377;
	}

	namespace Fe56 {
		constexpr static double Mass = 52.1031;
	}

	namespace Charm {
		constexpr static double Mass = 1.5;
	}
} 

#endif