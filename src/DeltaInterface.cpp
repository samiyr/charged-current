#ifndef DELTA_INTERFACE_H
#define DELTA_INTERFACE_H

#include "Flavor.cpp"

class DeltaInterface {
	public:
	DeltaInterface(const double location, const double width) : _peak_location(location), _peak_width(width) {}

	void evaluate(const double x, [[maybe_unused]] const double Q2) {
		_current_value = x * _get_value(x);
	}
	double xf_evaluate(const FlavorType flavor, const double x, const double Q2) {
		evaluate(x, Q2);
		return xf(flavor);
	}
	double xf([[maybe_unused]] const FlavorType flavor) const {
		return _current_value;
	}
	double xg() const {
		return xf(Flavor::Gluon);
	}
	double alpha_s([[maybe_unused]] const double Q2) const {
		return 1;
	}

	private:
	const double _peak_location;
	const double _peak_width;

	double _current_value;

	double _get_value(const double x) {
		return std::abs(x - _peak_location) < _peak_width ? 1.0 / _peak_width : 0.0;
	}
};

#endif