#ifndef LHA_NUCLEAR_INTERFACE_H
#define LHA_NUCLEAR_INTERFACE_H

#include <string>

#include "LHAPDF/LHAPDF.h"

#include "PDF/Interfaces/LHAInterface.cpp"

template <std::derived_from<LHAPDF::Extrapolator> Extrapolator = ZeroExtrapolator>
class LHANuclearInterface {
	public:
	const std::string set_name;
	const double Z;
	const double A;

	LHANuclearInterface(
		const std::string proton_set_name,
		const double Z,
		const double A
	) noexcept : set_name(proton_set_name), Z(Z), A(A), base(proton_set_name), flavor_values(TOTAL_FLAVORS, 0.0) { }

	void evaluate(const double x, const double Q2) const {
		if (x == prev_x && Q2 == prev_Q2) {	return;	}

		base.evaluate(x, Q2);

		const double Z_factor = Z / A;
		const double N_factor = (A - Z) / A;

		const double xtbar = base.xf(Flavor::AntiTop);
		const double xbbar = base.xf(Flavor::AntiBottom);
		const double xcbar = base.xf(Flavor::AntiCharm);
		const double xsbar = base.xf(Flavor::AntiStrange);
		const double xubar = base.xf(Flavor::AntiUp);
		const double xdbar = base.xf(Flavor::AntiDown);

		const double xg = base.xg();
		
		const double xd = base.xf(Flavor::Down);
		const double xu = base.xf(Flavor::Up);
		const double xs = base.xf(Flavor::Strange);
		const double xc = base.xf(Flavor::Charm);
		const double xb = base.xf(Flavor::Bottom);
		const double xt = base.xf(Flavor::Top);

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
		return base.alpha_s(Q2);
	}

	private:
	const LHAInterface<Extrapolator> base;
	mutable std::vector<double> flavor_values;

	mutable double prev_x = -1.0;
	mutable double prev_Q2 = -1.0;

	void set_flavor(const FlavorType flavor, const double value) const {
		flavor_values[static_cast<std::size_t>(flavor + 6)] = value;
	}
};

#endif