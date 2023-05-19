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

	double make_cross_section_variable(const double x, const double Q2, const double s, const Process process, const double f2, const double fL, const double xf3) {
		const std::optional<double> y = CommonFunctions::compute_y(x, Q2, s);

		const double term1 = y.has_value() ? 1 - *y + 0.5 * *y * *y : 0.0;
		const double term2 = y.has_value() ? - 0.5 * *y * *y : 0.0;
		const double term3 = y.has_value() ? *y * (1 - 0.5 * *y) : 0.0;

		const double result = (term1 * f2 + term2 * fL + double(process.W_sign()) * term3 * xf3) / x;

		return result;
	}
}


#endif