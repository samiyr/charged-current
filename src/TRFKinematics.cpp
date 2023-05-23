#ifndef TRF_KINEMATICS_H
#define TRF_KINEMATICS_H

#include <cmath>

struct TRFKinematics {
	double x;
	double y;
	double E_beam;

	double target_mass;
	double projectile_mass;

	double s;
	double Q2;

	constexpr static TRFKinematics y_E_beam(const double x, const double y, const double E_beam, const double target_mass, const double projectile_mass) {
		const double s = TRFKinematics::s_from_beam_energy(E_beam, target_mass, projectile_mass);
		const double Q2 = TRFKinematics::Q2_from_y(x, y, E_beam, target_mass);

		return TRFKinematics {x, y, E_beam, target_mass, projectile_mass, s, Q2};
	}
	constexpr static TRFKinematics Q2_s(const double x, const double Q2, const double s, const double target_mass, const double projectile_mass) {
		const double y = y_from_Q2(x, Q2, s, target_mass, projectile_mass);
		const double E_beam = beam_energy_from_s(s, target_mass, projectile_mass);

		return TRFKinematics {x, y, E_beam, target_mass, projectile_mass, s, Q2};
	}
	constexpr static TRFKinematics Q2_sqrt_s(const double x, const double Q2, const double sqrt_s, const double target_mass, const double projectile_mass) {
		return TRFKinematics::Q2_s(x, Q2, sqrt_s * sqrt_s, target_mass, projectile_mass);
	}

	private:
	constexpr static double s_from_beam_energy(const double E_beam, const double target_mass, const double projectile_mass) {
		return std::pow(target_mass, 2) + std::pow(projectile_mass, 2) + 2 * target_mass * E_beam;
	}
	constexpr static double beam_energy_from_s(const double s, const double target_mass, const double projectile_mass) {
		return (s - std::pow(target_mass, 2) - std::pow(projectile_mass, 2)) / (2 * target_mass);
	}
	constexpr static double Q2_from_y(const double x, const double y, const double E_beam, const double target_mass) {
		return 2 * target_mass * E_beam * x * y;
	}
	constexpr static double y_from_Q2(const double x, const double Q2, const double s, const double target_mass, const double projectile_mass) {
		return Q2 / ((s - std::pow(target_mass, 2) - std::pow(projectile_mass, 2)) * x);
	}
};


#endif