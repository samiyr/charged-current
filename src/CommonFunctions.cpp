#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H

#include "Constants.cpp"
#include <optional>

namespace CommonFunctions {
	constexpr std::optional<double> compute_y(const double x, const double Q2, const double s, const double target_mass, const double projectile_mass = 0.0) {
		const double y = Q2 / (x * (s - std::pow(target_mass, 2) - std::pow(projectile_mass, 2)));
		if (y < 0 || y > 1) { return std::nullopt; }
		return y;
	}
	constexpr double cross_section_prefactor(const double Q2) {
		constexpr double numerator = POW2(Constants::fermi_coupling) * POW4(Constants::Particles::W.mass);
		const double denominator = 2 * M_PI * POW2(Q2 + POW2(Constants::Particles::W.mass));

		return numerator / denominator;
	}

	constexpr double make_cross_section_variable(const double x, const double Q2, const double s, const Process process, const double f2, const double fL, const double xf3) {
		const std::optional<double> y_opt = CommonFunctions::compute_y(x, Q2, s, process.target.mass, process.projectile.mass);

		if (!y_opt.has_value()) { return 0.0; }

		const double y = *y_opt;
		const double M2 = std::pow(process.target.mass, 2);

		const double term1 = 1 - y + 0.5 * y * y - (x * y * M2) / (s - M2);
		const double term2 = - 0.5 * y * y;
		const double term3 = y * (1 - 0.5 * y);

		const double result = (term1 * f2 + term2 * fL + double(process.W_sign()) * term3 * xf3) / x;

		return result;
	}

	constexpr double compute_momentum_fraction_mass_correction(const double x_0, const double Q2, const double mass_scale, const double target_mass) {
		return x_0 * (1 + std::pow(mass_scale, 2) / Q2) * (1 - std::pow(x_0 * target_mass, 2) / Q2);
	}
}


#endif