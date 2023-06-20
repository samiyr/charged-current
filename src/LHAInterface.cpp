#ifndef LHA_INTERFACE_H
#define LHA_INTERFACE_H

#include <string>
#include <stdexcept>
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"
#include "Flavor.cpp"
#include "CKM.cpp"
#include "Process.cpp"
#include "ZeroExtrapolator.cpp"

#define CACHE_STATS false

class LHAInterface {
	public:
	std::string set_name;
	int set_member_number;

	FlavorVector available_flavors;

	LHAInterface(std::string _set_name, int _set_member_number = 0)
	: set_name(_set_name), 
	set_member_number(_set_member_number), 
	flavor_values(TOTAL_FLAVORS, 0.0),
	use_multipliers(false),
	multipliers({}) {
		initialize();
		available_flavors = _pdf->flavors();
	}

	LHAInterface(std::string _set_name, const std::vector<double> _multipliers, int _set_member_number = 0)
	: set_name(_set_name), 
	set_member_number(_set_member_number), 
	flavor_values(TOTAL_FLAVORS, 0.0),
	use_multipliers(true),
	multipliers(_multipliers) {
		initialize();
		available_flavors = _pdf->flavors();
	}

	static void disable_verbosity() {
		LHAPDF::Info &cfg = LHAPDF::getConfig();
		cfg.set_entry("Verbosity", 0);
	}

	void evaluate(const double x, const double Q2) const {
		#if CACHE_STATS
		std::cout << "Cache hit ratio: " << 100 * double(cache_hits) / double(total_hits) << " (cache hits: " << cache_hits << ", total hits: " << total_hits << ")" << std::endl;
		total_hits++;
		#endif
		if (x == prev_x && Q2 == prev_Q2) { 
			#if CACHE_STATS
			cache_hits++; 
			#endif
			return; 
		}

		_pdf->xfxQ2(x, Q2, flavor_values);

		prev_x = x;
		prev_Q2 = Q2;

		if (use_multipliers) {
			Utility::multiply(flavor_values, multipliers);
		}
	}
	
	double xf_evaluate(const FlavorType flavor, const double x, const double Q2) const {
		return _pdf->xfxQ2(flavor, x, Q2);
	}
	constexpr double xf(const FlavorType flavor) const {
		return flavor_values[size_t(flavor + 6)];
	}
	constexpr double xg() const {
		return xf(Flavor::Gluon);
	}
	double alpha_s(const double Q2) const {
		return _pdf->alphasQ2(Q2);
	}

	LHAInterface(const LHAInterface &o) : set_name(o.set_name), set_member_number(o.set_member_number), use_multipliers(o.use_multipliers), multipliers(o.multipliers) {
		initialize();
	}

	private:
	mutable std::vector<double> flavor_values;
	mutable std::unique_ptr<LHAPDF::PDF> _pdf;

	const bool use_multipliers;
	const std::vector<double> multipliers;
	
	void initialize() {
		_pdf = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(set_name, set_member_number));
		available_flavors = _pdf->flavors();
		flavor_values = std::vector<double>(TOTAL_FLAVORS, 0.0);
		prev_x = -1.0;
		prev_Q2 = -1.0;

		LHAPDF::GridPDF *grid_pdf = static_cast<LHAPDF::GridPDF *>(_pdf.get());
		LHAPDF::Extrapolator *zero_extrapolator = new ZeroExtrapolator();
		grid_pdf->setExtrapolator(zero_extrapolator);
		
		#if CACHE_STATS
		total_hits = 0;
		cache_hits = 0;
		#endif
	}

	mutable double prev_x = -1.0;
	mutable double prev_Q2 = -1.0;

	#if CACHE_STATS
	mutable size_t total_hits = 0;
	mutable size_t cache_hits = 0;
	#endif
};

#endif