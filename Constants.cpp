#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "Utility.cpp"

namespace Constants {
	constexpr double C_F = 4.0 / 3.0;
	constexpr double fermi_coupling = 1.1663787e-5;
	constexpr double mass_W = 80.377;
	constexpr double lambda_QCD = 0.226;

	constexpr double ckm_squared_matrix_elements[3][3] = {
		/*
		  	V_ud				V_us				V_ub
		  	V_cd				V_cs				V_cb
		  	V_td				V_ts				V_tb
		*/
		{POW2(0.97435), 	POW2(0.22500), 		POW2(0.00369)},
		{POW2(0.22486),		POW2(0.97349),		POW2(0.04182)},
		{POW2(0.00857),		POW2(0.04110),		POW2(0.999112)}
	};
} 

#endif