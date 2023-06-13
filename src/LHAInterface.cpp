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

#define TOTAL_FLAVORS 13

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
	_pdf(LHAPDF::mkPDF(set_name, set_member_number)) {
		std::shared_ptr<LHAPDF::GridPDF> grid_pdf = std::static_pointer_cast<LHAPDF::GridPDF>(_pdf);
		LHAPDF::Extrapolator *zero_extrapolator = new ZeroExtrapolator();
		grid_pdf->setExtrapolator(zero_extrapolator);
		available_flavors = _pdf->flavors();
	}

	static void disable_verbosity() {
		LHAPDF::Info &cfg = LHAPDF::getConfig();
		cfg.set_entry("Verbosity", 0);
	}

	void evaluate(const double x, const double Q2) {
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
	}
	
	double xf_evaluate(const FlavorType flavor, const double x, const double Q2) {
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

	private:
	std::vector<double> flavor_values;
	std::shared_ptr<LHAPDF::PDF> _pdf;

	double prev_x = -1.0;
	double prev_Q2 = -1.0;

	#if CACHE_STATS
	size_t total_hits = 0;
	size_t cache_hits = 0;
	#endif
};

#endif