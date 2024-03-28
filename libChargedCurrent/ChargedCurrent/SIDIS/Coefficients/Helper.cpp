#ifndef SIDIS_FUNCTIONS_HELPER_H
#define SIDIS_FUNCTIONS_HELPER_H

#include <cmath>

#include "Decay/Decay.cpp"

#include "Common/Constants.cpp"

#include "SIDIS/Coefficients/EvaluationParameters.cpp"

namespace SIDISFunctions::Helper {
	template <is_decay_function DecayFunction, typename Kinematics>
	constexpr double compute_z_min(const Kinematics &kinematics, const Decay<DecayFunction> &decay) {
		return std::max({
			decay.lepton_momentum_min / (kinematics.y * kinematics.E_beam), 
			decay.resonance.mass / (kinematics.y * kinematics.E_beam),
			decay.z_min_cutoff
		});
	}
}

#endif