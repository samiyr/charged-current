#ifndef LHA_INTERFACE_H
#define LHA_INTERFACE_H

#include <string>
#include <stdexcept>
#include "LHAPDF/LHAPDF.h"
#include "Flavor.cpp"
#include "CKM.cpp"
#include "Process.cpp"
// #include <boost/compute/detail/lru_cache.hpp>

#define TOTAL_FLAVORS 13

#define CACHE_STATS false

class LHAInterface {
	// using Cache = boost::compute::detail::lru_cache<std::pair<double, double>, std::vector<double>>;
	public:
	std::string set_name;
	int set_member_number;

	FlavorVector available_flavors;

	LHAInterface(std::string _set_name, int _set_member_number = 0)
	: set_name(_set_name), set_member_number(_set_member_number), flavor_values(TOTAL_FLAVORS, 0.0) {
		LHAPDF::Info &cfg = LHAPDF::getConfig();
		cfg.set_entry("Verbosity", 0);
		initialize();
	}

	void evaluate(const double x, const double Q2) {
		#if CACHE_STATS
		std::cout << "Cache hit ratio: " << 100 * double(cache_hits) / double(total_hits) << " (cache hits: " << cache_hits << ", total hits: " << total_hits << ")" << std::endl;
		total_hits++;
		#endif
		// if (const boost::optional<std::vector<double>> values = cache.get({x, Q2})) {
		// 	flavor_values = *values;
		// 	return;
		// }
		if (x == prev_x && Q2 == prev_Q2) { 
			#if CACHE_STATS
			cache_hits++; 
			#endif
			return; 
		}

		if (_pdf->inRangeXQ2(x, Q2)) {
			_pdf->xfxQ2(x, Q2, flavor_values);

			// cache.insert({x, Q2}, flavor_values);
			prev_x = x;
			prev_Q2 = Q2;
		} else {
			std::fill(flavor_values.begin(), flavor_values.end(), 0.0);
		}
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

	LHAInterface(const LHAInterface &o) {
		set_name = o.set_name;
		set_member_number = o.set_member_number;
		initialize();
	}

	private:
	std::unique_ptr<LHAPDF::PDF> _pdf;
	std::vector<double> flavor_values;
	
	void initialize() {
		_pdf = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(set_name, set_member_number));
		available_flavors = _pdf->flavors();
		flavor_values = std::vector<double>(TOTAL_FLAVORS, 0.0);
		// cache = Cache(TOTAL_FLAVORS);
		prev_x = -1.0;
		prev_Q2 = -1.0;
		#if CACHE_STATS
		total_hits = 0;
		cache_hits = 0;
		#endif
	}

	// Cache cache{TOTAL_FLAVORS};

	double prev_x;
	double prev_Q2;

	#if CACHE_STATS
	size_t total_hits = 0;
	size_t cache_hits = 0;
	#endif
};

#endif