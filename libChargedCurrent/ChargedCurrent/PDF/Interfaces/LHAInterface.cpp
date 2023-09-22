#ifndef LHA_INTERFACE_H
#define LHA_INTERFACE_H

#include <string>
#include <stdexcept>

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"

#include "Common/Flavor.cpp"
#include "Common/CKM.cpp"
#include "Common/Process.cpp"

#include "PDF/ZeroExtrapolator.cpp"

#include "Utility/Globals.cpp"

template <typename explicit_isospin = std::false_type, std::derived_from<LHAPDF::Extrapolator> Extrapolator = ZeroExtrapolator>
class LHAInterface {
	public:
	std::string set_name;
	int set_member_number;

	double Z = 1.0;
	double A = 1.0;

	LHAInterface(
		const std::string _set_name, 
		const int _set_member_number, 
		const bool _use_multipliers, 
		const std::vector<double> _multipliers,
		const bool _use_global_multiplier = false,
		const double _global_multiplier = 1.0,
		const double Z = 1.0,
		const double A = 1.0) noexcept
	: set_name(_set_name), 
	set_member_number(_set_member_number), 
	Z(Z),
	A(A),
	use_multipliers(_use_multipliers), 
	multipliers(_multipliers),
	use_global_multiplier(_use_global_multiplier),
	global_multiplier(_global_multiplier),
	flavor_values(TOTAL_FLAVORS, 0.0),
	prev_x(-1.0),
	prev_Q2(-1.0),
	total_hits(0),
	cache_hits(0) { }

	LHAInterface(std::string _set_name, const std::vector<double> _multipliers, int _set_member_number = 0) noexcept
	: LHAInterface(_set_name, _set_member_number, _multipliers.size() == TOTAL_FLAVORS, _multipliers) { }

	LHAInterface(std::string _set_name, int _set_member_number = 0) noexcept
	: LHAInterface(_set_name, _set_member_number, false, {}) { }

	LHAInterface(const LHAInterface &o) noexcept : LHAInterface(
		o.set_name, o.set_member_number, o.use_multipliers, o.multipliers, o.use_global_multiplier, o.global_multiplier, o.Z, o.A
	) { }
	LHAInterface(LHAInterface &&o) = default;

	void enable_global_multiplier(const double multiplier) {
		use_global_multiplier = true;
		global_multiplier = multiplier;
	}
	void disable_global_multiplier() {
		use_global_multiplier = false;
		global_multiplier = 1.0;
	}

	void activate() const {
		if (activated) [[likely]] { return; }

		pdf = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(set_name, set_member_number));

		LHAPDF::GridPDF *grid_pdf = static_cast<LHAPDF::GridPDF *>(pdf.get());
		LHAPDF::Extrapolator *extrapolator = new Extrapolator();
		// Ownership of extrapolator is passed to LHAPDF::GridPDF, it handles deletion.
		grid_pdf->setExtrapolator(extrapolator);

		activated = true;
	}

	static void disable_verbosity() {
		LHAPDF::Info &config = LHAPDF::getConfig();
		config.set_entry("Verbosity", 0);
	}

	void evaluate(const double x, const double Q2) const {
		activate();

		if constexpr (Globals::LHAInterfaceCacheStats) {
			std::cout << "Cache hit ratio: " << 100 * double(cache_hits) / double(total_hits);
			std::cout << " (cache hits: " << cache_hits << ", total hits: " << total_hits << ")" << IO::endl;
			total_hits++;
		}
		if (x == prev_x && Q2 == prev_Q2) { 
			if constexpr (Globals::LHAInterfaceCacheStats) {
				cache_hits++; 
			}
			return;
		}

		pdf->xfxQ2(x, Q2, flavor_values);

		prev_x = x;
		prev_Q2 = Q2;

		if constexpr (explicit_isospin::value) {
			const double Z_factor = Z / A;
			const double N_factor = (A - Z) / A;

			const double xtbar = xf(Flavor::AntiTop);
			const double xbbar = xf(Flavor::AntiBottom);
			const double xcbar = xf(Flavor::AntiCharm);
			const double xsbar = xf(Flavor::AntiStrange);
			const double xubar = xf(Flavor::AntiUp);
			const double xdbar = xf(Flavor::AntiDown);

			const double xg = xf(Flavor::Gluon);
			
			const double xd = xf(Flavor::Down);
			const double xu = xf(Flavor::Up);
			const double xs = xf(Flavor::Strange);
			const double xc = xf(Flavor::Charm);
			const double xb = xf(Flavor::Bottom);
			const double xt = xf(Flavor::Top);

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
		}

		if (use_multipliers) {
			for (const auto &[value, multiplier] : std::views::zip(flavor_values, multipliers)) {
				value *= multiplier * global_multiplier;
			}
		} else if (use_global_multiplier) {
			for (auto &value : flavor_values) {
				value *= global_multiplier;
			}
		}
	}
	
	double xf_evaluate(const FlavorType flavor, const double x, const double Q2) const {
		activate();
		return pdf->xfxQ2(flavor, x, Q2);
	}
	constexpr double xf(const FlavorType flavor) const {
		return flavor_values[std::size_t(flavor + 6)];
	}
	constexpr double xg() const {
		return xf(Flavor::Gluon);
	}
	double alpha_s(const double Q2) const {
		activate();
		return pdf->alphasQ2(Q2);
	}

	private:
	const bool use_multipliers;
	const std::vector<double> multipliers;

	bool use_global_multiplier = false;
	double global_multiplier = 1.0;

	mutable std::unique_ptr<LHAPDF::PDF> pdf;
	mutable std::vector<double> flavor_values;

	mutable bool activated = false;

	mutable double prev_x = -1.0;
	mutable double prev_Q2 = -1.0;

	mutable std::size_t total_hits = 0;
	mutable std::size_t cache_hits = 0;

	void set_flavor(const FlavorType flavor, const double value) const {
		flavor_values[static_cast<std::size_t>(flavor + 6)] = value;
	}
};

template <typename explicit_isospin, std::derived_from<LHAPDF::Extrapolator> Extrapolator> 
const LHAInterface<explicit_isospin, Extrapolator> operator*(const double lhs, const LHAInterface<explicit_isospin, Extrapolator> &rhs) {
	LHAInterface<explicit_isospin, Extrapolator> copy(rhs);
	copy.enable_global_multiplier(lhs);
	return copy;
}

template <typename explicit_isospin, std::derived_from<LHAPDF::Extrapolator> Extrapolator> 
const LHAInterface<explicit_isospin, Extrapolator> operator*(const LHAInterface<explicit_isospin, Extrapolator> &lhs, const double rhs) {
	return rhs * lhs;
}

#endif