#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H

#include <optional>
#include <numbers>

#include "Common/Constants.cpp"
#include "Common/Process.cpp"
#include "Common/Flavor.cpp"
#include "Common/TRFKinematics.cpp"

namespace CommonFunctions {
	// constexpr std::optional<double> compute_y(const double x, const double Q2, const double s, const double target_mass, const double projectile_mass = 0.0) {
	// 	const double y = Q2 / (x * (s - std::pow(target_mass, 2) - std::pow(projectile_mass, 2)));
	// 	if (y < 0 || y > 1) { return std::nullopt; }
	// 	return y;
	// }

	// Computes the prefactor G_W^2 M_W^4 / 2π (Q^2 + M_W^2)^2.
	template <typename Kinematics>
	constexpr double cross_section_prefactor(const Kinematics &kinematics) noexcept {
		constexpr double numerator = Math::pow2(Constants::fermi_coupling) * Math::pow4(Constants::Particles::W.mass);
		const double denominator = 2.0 * std::numbers::pi * Math::pow2(kinematics.Q2 + Math::pow2(Constants::Particles::W.mass));

		return numerator / denominator;
	}
	// Computes the prefactor 100π M_W^4 / 2π (Q^2 + M_W^2)^2 M E.
	template <typename Kinematics>
	constexpr double cross_section_modified_prefactor(const Kinematics &kinematics) {
		constexpr double numerator = 50.0 * Math::pow4(Constants::Particles::W.mass);
		const double denominator = Math::pow2(kinematics.Q2 + Math::pow2(Constants::Particles::W.mass)) * kinematics.target_mass * kinematics.E_beam;

		return numerator / denominator;
	}
	// Computes the mass correction to the momentum fraction x, using x_mass = x (1 + m^2 / Q^2) (1 - M^2 / Q2), where M is the target mass.
	constexpr double compute_momentum_fraction_mass_correction(const double x_0, const double Q2, const double mass_scale, const double target_mass) {
		return x_0 * (1 + std::pow(mass_scale, 2) / Q2) * (1 - std::pow(x_0 * target_mass, 2) / Q2);
	}
	// Computes the Jacobian of the transformation (x, Q^2) -> (x, y).
	template <typename Kinematics>
	constexpr static double xy_jacobian(const Kinematics &kinematics, const Process &process) {
		return (kinematics.s - std::pow(process.target.mass, 2) - std::pow(process.projectile.mass, 2)) * kinematics.x;
	}
	// Combines the F2, FL and F3 structure function values into the appropriate cross section value. Takes the sign of neutrino vs antineutrino into account.
	template <typename Kinematics, typename T>
	constexpr T make_cross_section_variable(
		const Kinematics &kinematics, 
		const Process &process, 
		const T f2, 
		const T fL, 
		const T f3) {

		const double x = kinematics.x;
		const double y = kinematics.y;
		const double s = kinematics.s;
		const double M2 = std::pow(process.target.mass, 2);

		const double term1 = 1 - y + 0.5 * y * y - (x * y * M2) / (s - M2);
		const double term2 = - 0.5 * y * y;
		const double term3 = y * (1 - 0.5 * y);

		const T result = (term1 * f2 + term2 * fL) / x + double(process.W_sign()) * term3 * f3;

		return result;
	}
}

#endif