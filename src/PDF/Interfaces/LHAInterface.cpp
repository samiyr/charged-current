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

class LHAInterface {
	public:
	std::string set_name;
	int set_member_number;

	private:
	LHAInterface(const std::string _set_name, const int _set_member_number, const bool _use_multipliers, const std::vector<double> _multipliers)
	: set_name(_set_name), 
	set_member_number(_set_member_number), 
	use_multipliers(_use_multipliers), 
	multipliers(_multipliers),
	flavor_values(TOTAL_FLAVORS, 0.0),
	prev_x(-1.0),
	prev_Q2(-1.0) {
		if constexpr (Globals::LHAInterfaceCacheStats) {
			total_hits = 0;
			cache_hits = 0;
		}
	}

	public:
	LHAInterface(std::string _set_name, const std::vector<double> _multipliers, int _set_member_number = 0)
	: LHAInterface(_set_name, _set_member_number, _multipliers.size() == TOTAL_FLAVORS, _multipliers) { }

	LHAInterface(std::string _set_name, int _set_member_number = 0)
	: LHAInterface(_set_name, _set_member_number, false, {}) { }

	void activate() const {
		if (activated) { return; }

		pdf = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(set_name, set_member_number));

		LHAPDF::GridPDF *grid_pdf = static_cast<LHAPDF::GridPDF *>(pdf.get());
		LHAPDF::Extrapolator *zero_extrapolator = new ZeroExtrapolator();
		grid_pdf->setExtrapolator(zero_extrapolator);

		activated = true;
	}

	LHAInterface(const LHAInterface &o) : LHAInterface(o.set_name, o.set_member_number, o.use_multipliers, o.multipliers) { }

	static void disable_verbosity() {
		LHAPDF::Info &cfg = LHAPDF::getConfig();
		cfg.set_entry("Verbosity", 0);
	}

	void evaluate(const double x, const double Q2) const {
		activate();

		if constexpr (Globals::LHAInterfaceCacheStats) {
			std::cout << "Cache hit ratio: " << 100 * double(cache_hits) / double(total_hits) << " (cache hits: " << cache_hits << ", total hits: " << total_hits << ")" << IO::endl;
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

		if (use_multipliers) {
			Collections::multiply(flavor_values, multipliers);
		}
	}
	
	double xf_evaluate(const FlavorType flavor, const double x, const double Q2) const {
		activate();
		return pdf->xfxQ2(flavor, x, Q2);
	}
	constexpr double xf(const FlavorType flavor) const {
		return flavor_values[size_t(flavor + 6)];
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

	mutable std::unique_ptr<LHAPDF::PDF> pdf;
	mutable std::vector<double> flavor_values;

	mutable bool activated = false;

	mutable double prev_x = -1.0;
	mutable double prev_Q2 = -1.0;

	mutable size_t total_hits = 0;
	mutable size_t cache_hits = 0;
};

#endif