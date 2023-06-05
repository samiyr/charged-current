#ifndef LHA_INTERFACE_H
#define LHA_INTERFACE_H

#include <string>
#include <stdexcept>
#include "LHAPDF/LHAPDF.h"
#include "Flavor.cpp"
#include "CKM.cpp"
#include "Process.cpp"

#define TOTAL_FLAVORS 13

class LHAInterface {
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
		if (x == prev_x && Q2 == prev_Q2) { return; }

		if (_pdf->inRangeXQ2(x, Q2)) {
			_pdf->xfxQ2(x, Q2, flavor_values);

			prev_x = x;
			prev_Q2 = Q2;

			// if (set_name.find("kkks08") != std::string::npos) {
			// 	const FlavorType isolated = Flavor::Gluon;
			// 	const double isolated_value = flavor_values[size_t(isolated + 6)];
			// 	flavor_values.assign(available_flavors.size(), 0.0);
			// 	flavor_values[size_t(isolated + 6)] = isolated_value;
			// }

			// for (size_t i = 0; i < flavor_values.size(); i++) {
			// 	flavor_values[i] = std::max(flavor_values[i], 0.0);
			// }
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
		prev_x = -1.0;
		prev_Q2 = -1.0;
	}

	double prev_x;
	double prev_Q2;

};

#endif