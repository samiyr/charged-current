#ifndef LHA_INTERFACE_H
#define LHA_INTERFACE_H

#include <string>
#include <stdexcept>
#include "LHAPDF/LHAPDF.h"
#include "Flavor.cpp"
#include "CKM.cpp"
#include "Process.cpp"

class LHAInterface {
	public:
	std::string set_name;
	int set_member_number;

	FlavorVector available_flavors;

	LHAInterface(std::string _set_name, int _set_member_number = 0)
	: set_name(_set_name), set_member_number(_set_member_number) {
		LHAPDF::Info &cfg = LHAPDF::getConfig();
		cfg.set_entry("Verbosity", 0);
		_pdf = LHAPDF::mkPDF(set_name, set_member_number);
		available_flavors = _pdf->flavors();
	}

	void evaluate(const double x, const double Q2) {
		if (_pdf->inRangeXQ2(x, Q2)) {
			_pdf->xfxQ2(x, Q2, flavor_values);

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
			flavor_values.assign(available_flavors.size(), 0.0);
		}
	}
	double xf_evaluate(const FlavorType flavor, const double x, const double Q2) {
		return _pdf->xfxQ2(flavor, x, Q2);
	}
	double xf(const FlavorType flavor) const {
		return flavor_values[size_t(flavor + 6)];
	}
	double xg() const {
		return xf(Flavor::Gluon);
	}
	double alpha_s(const double Q2) const {
		return _pdf->alphasQ2(Q2);
	}
	private:
	LHAPDF::PDF *_pdf;
	std::vector<double> flavor_values;
};

#endif