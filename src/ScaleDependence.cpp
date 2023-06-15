#ifndef SCALE_DEPENDENCE_H
#define SCALE_DEPENDENCE_H

#include "TRFKinematics.cpp"
#include <optional>

namespace ScaleDependence {
	static auto Q2 = std::make_optional([](const TRFKinematics &kinematics) { return kinematics.Q2; });
	static auto trivial = [](const TRFKinematics &kinematics) { return kinematics.Q2; };
	static auto constant(const double value) {
		return std::make_optional([value]([[maybe_unused]] const TRFKinematics &kinematics) { return value; });
	}
	static auto multiplicative(const double factor) {
		return std::make_optional([factor](const TRFKinematics &kinematics) { return kinematics.Q2 * factor; });
	}
}

#endif