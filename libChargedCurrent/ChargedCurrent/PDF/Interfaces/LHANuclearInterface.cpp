#ifndef LHA_NUCLEAR_INTERFACE_H
#define LHA_NUCLEAR_INTERFACE_H

#include <string>

#include "LHAPDF/LHAPDF.h"

#include "PDF/Interfaces/LHAInterface.cpp"

template <std::derived_from<LHAPDF::Extrapolator> Extrapolator = ZeroExtrapolator>
class LHANuclearInterface {
	public:
	const std::string proton_set_name;
	const double Z;
	const double A;

	LHANuclearInterface(
		const std::string proton_set_name,
		const double Z,
		const double A
	) noexcept : proton_set_name(proton_set_name), Z(Z), A(A), proton(proton_set_name), flavor_values(TOTAL_FLAVORS, 0.0) { }

	void evaluate(const double x, const double Q2) const {
		if (x == prev_x && Q2 == prev_Q2) {	return;	}

		proton.evaluate(x, Q2);

		const double Z_factor = Z / A;
		const double N_factor = (A - Z) / A;

		const double xtbar = proton.xf(Flavor::AntiTop);
		const double xbbar = proton.xf(Flavor::AntiBottom);
		const double xcbar = proton.xf(Flavor::AntiCharm);
		const double xsbar = proton.xf(Flavor::AntiStrange);
		const double xubar = proton.xf(Flavor::AntiUp);
		const double xdbar = proton.xf(Flavor::AntiDown);

		const double xg = proton.xg();
		
		const double xd = proton.xf(Flavor::Down);
		const double xu = proton.xf(Flavor::Up);
		const double xs = proton.xf(Flavor::Strange);
		const double xc = proton.xf(Flavor::Charm);
		const double xb = proton.xf(Flavor::Bottom);
		const double xt = proton.xf(Flavor::Top);

		set_flavor(Flavor::AntiTop, xtbar);
		set_flavor(Flavor::AntiBottom, xbbar);
		set_flavor(Flavor::AntiCharm, xcbar);
		set_flavor(Flavor::AntiStrange, xsbar);

		set_flavor(Flavor::AntiUp, Z_factor * xubar + N_factor * xdbar);
		set_flavor(Flavor::AntiDown, Z_factor * xdbar + N_factor * xubar);

		set_flavor(Flavor::Gluon, xg);

		set_flavor(Flavor::Down, Z_factor * xd + N_factor * xu);
		set_flavor(Flavor::Up, Z_factor * xu + N_factor * xd);

		set_flavor(Flavor::Strange, xs);
		set_flavor(Flavor::Charm, xc);
		set_flavor(Flavor::Bottom, xb);
		set_flavor(Flavor::Top, xt);

		prev_x = x;
		prev_Q2 = Q2;
	}

	double xf_evaluate(const FlavorType flavor, const double x, const double Q2) const {
		const std::vector<double> current_flavor_values = flavor_values;

		evaluate(x, Q2);
		const double output_value = flavor_values[static_cast<std::size_t>(flavor + 6)];

		flavor_values = current_flavor_values;
		return output_value;
	}

	constexpr double xf(const FlavorType flavor) const {
		return flavor_values[static_cast<std::size_t>(flavor + 6)];
	}
	constexpr double xg() const {
		return xf(Flavor::Gluon);
	}
	double alpha_s(const double Q2) const {
		return proton.alpha_s(Q2);
	}

	private:
	const LHAInterface<Extrapolator> proton;
	mutable std::vector<double> flavor_values;

	mutable double prev_x = -1.0;
	mutable double prev_Q2 = -1.0;

	void set_flavor(const FlavorType flavor, const double value) const {
		flavor_values[static_cast<std::size_t>(flavor + 6)] = value;
	}
}

#endif