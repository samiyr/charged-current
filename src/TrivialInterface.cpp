#ifndef TRIVIAL_INTERFACE_H
#define TRIVIAL_INTERFACE_H

#include "Flavor.cpp"

class TrivialInterface {
	public:
	TrivialInterface(const double constant = 1) : _constant(constant) {}

	void evaluate(const double x, [[maybe_unused]] const double Q2) {
		_x = x;
	}
	double xf_evaluate([[maybe_unused]] const FlavorType flavor, const double x, [[maybe_unused]] const double Q2) {
		return _constant * x;
	}
	double xf([[maybe_unused]] const FlavorType flavor) const {
		return _constant * _x;
	}
	double xg() const {
		return xf(Flavor::Gluon);
	}
	double alpha_s([[maybe_unused]] const double Q2) const {
		return 1;
	}

	private:
	const double _constant;
	double _x;
};

#endif