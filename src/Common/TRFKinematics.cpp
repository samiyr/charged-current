#ifndef TRF_KINEMATICS_H
#define TRF_KINEMATICS_H

#include <cmath>

#include "Common/Constants.cpp"

// A structure encapsulating target rest-frame kinematics.
struct TRFKinematics {
	// Björken-x (momentum fraction), in the range [0, 1].
	double x;
	// Inelasticity, in the range [0, 1].
	double y;
	// Beam energy in units of GeV.
	double E_beam;

	// Target (i.e. proton or nucleus) mass in units of GeV.
	double target_mass;
	// Projectile (i.e. neutrino) mass in units of GeV.
	double projectile_mass;

	// Mandelstam s, in units of GeV^2.
	double s;
	// Virtuality, in units of GeV^2.
	double Q2;

	/// @brief Constructs the kinematics from x, y and the beam energy.
	/// @param x Björken-x, must be in the range [0, 1].
	/// @param y Inelasticity, must be in the range [0, 1].
	/// @param E_beam Beam energy in units of GeV.
	/// @param target_mass Target (i.e. proton or nucleus) mass in units of GeV.
	/// @param projectile_mass Projectile (i.e. neutrino) mass in units of GeV.
	constexpr static TRFKinematics y_E_beam(
		const double x, const double y, const double E_beam, 
		const double target_mass, const double projectile_mass) {

		const double s = TRFKinematics::s_from_beam_energy(E_beam, target_mass, projectile_mass);
		const double Q2 = TRFKinematics::Q2_from_y(x, y, E_beam, target_mass);

		return TRFKinematics {x, y, E_beam, target_mass, projectile_mass, s, Q2};
	}
	/// @brief Constructs the kinematics from x, Q^2 and s.
	/// @param x Björken-x, must be in the range [0, 1].
	/// @param Q2 Virtuality in units of GeV^2.
	/// @param s Mandelstam s in units of GeV^2.
	/// @param target_mass Target (i.e. proton or nucleus) mass in units of GeV.
	/// @param projectile_mass Projectile (i.e. neutrino) mass in units of GeV.
	constexpr static TRFKinematics Q2_s(
		const double x, const double Q2, const double s, 
		const double target_mass, const double projectile_mass) {

		const double y = y_from_Q2(x, Q2, s, target_mass, projectile_mass);
		const double E_beam = beam_energy_from_s(s, target_mass, projectile_mass);

		return TRFKinematics {x, y, E_beam, target_mass, projectile_mass, s, Q2};
	}
	/// @brief Constructs the kinematics from x, Q^2 and sqrt(s).
	/// @param x Björken-x, must be in the range [0, 1].
	/// @param Q2 Virtuality in units of GeV^2.
	/// @param sqrt_s Square root of Mandelstam s in units of GeV.
	/// @param target_mass Target (i.e. proton or nucleus) mass in units of GeV.
	/// @param projectile_mass Projectile (i.e. neutrino) mass in units of GeV.
	constexpr static TRFKinematics Q2_sqrt_s(
		const double x, const double Q2, const double sqrt_s, 
		const double target_mass, const double projectile_mass) {
		return TRFKinematics::Q2_s(x, Q2, sqrt_s * sqrt_s, target_mass, projectile_mass);
	}

	// Checks whether the kinematical variables are in their physically valid ranges.
	constexpr bool is_valid() const noexcept {
		return (y >= 0.0 && y <= 1.0 && x >= 0.0 && x <= 1.0 && E_beam >= 0.0 && target_mass >= 0.0 && projectile_mass >= 0.0 && s >= 0.0 && Q2 >= 0.0);
	}

	private:
	// Conversion from beam energy to Mandelstam s.
	constexpr static double s_from_beam_energy(const double E_beam, const double target_mass, const double projectile_mass) {
		return std::pow(target_mass, 2) + std::pow(projectile_mass, 2) + 2 * target_mass * E_beam;
	}
	// Conversion from Mandelstam s to beam energy.
	constexpr static double beam_energy_from_s(const double s, const double target_mass, const double projectile_mass) {
		return (s - std::pow(target_mass, 2) - std::pow(projectile_mass, 2)) / (2 * target_mass);
	}
	// Conversion from y to Q^2.
	constexpr static double Q2_from_y(const double x, const double y, const double E_beam, const double target_mass) {
		return 2 * target_mass * E_beam * x * y;
	}
	// Conversion from Q^2 to y.
	constexpr static double y_from_Q2(const double x, const double Q2, const double s, const double target_mass, const double projectile_mass) {
		return Q2 / ((s - std::pow(target_mass, 2) - std::pow(projectile_mass, 2)) * x);
	}
};


#endif