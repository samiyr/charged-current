#ifndef SCALE_DEPENDENCE_H
#define SCALE_DEPENDENCE_H

#include "TRFKinematics.cpp"
#include <optional>
#include <concepts>

namespace ScaleDependence {
	template <typename T>
	concept Concept = requires(T scale_function, const TRFKinematics &kinematics) {
		{ scale_function(kinematics) } -> std::same_as<double>;
	};

	[[maybe_unused]] static auto Q2 = std::make_optional([](const TRFKinematics &kinematics) { return kinematics.Q2; });
	static auto trivial = [](const TRFKinematics &kinematics) { return kinematics.Q2; };
	[[maybe_unused]] static auto constant(const double value) {
		return std::make_optional([value]([[maybe_unused]] const TRFKinematics &kinematics) { return value; });
	}
	[[maybe_unused]] static auto multiplicative(const double factor) {
		return std::make_optional([factor](const TRFKinematics &kinematics) { return kinematics.Q2 * factor; });
	}
}

#endif