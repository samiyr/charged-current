#ifndef DIS_FUNCTIONS_HELPER_H
#define DIS_FUNCTIONS_HELPER_H

#include <numbers>

#include "Common/Constants.cpp"

namespace DISFunctions::Helper {
	constexpr double delta_contribution(const double x) {
		const double log = std::log(1 - x);
		return -2.0 * Constants::C_F * (9.0 / 2.0 + (std::numbers::pi * std::numbers::pi) / 3.0 - log * log + 1.5 * log);
	}
}

#endif