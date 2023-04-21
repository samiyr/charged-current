#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H

#include "Constants.cpp"
#include <optional>

namespace CommonFunctions {
	constexpr std::optional<double> compute_y(const double x, const double Q2, const double s) {
		const double y = Q2 / (x * s);
		if (y < 0 || y > 1) { return std::nullopt; }
		return y;
	}
	constexpr double cross_section_prefactor(const double Q2) {
		const double numerator = POW2(Constants::fermi_coupling) * POW4(Constants::mass_W);
		const double denominator = 2 * M_PI * POW2(Q2 + POW2(Constants::mass_W));

		return numerator / denominator;
	}
}


#endif