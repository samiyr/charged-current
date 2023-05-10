#ifndef FUNCTIONAL_FORM_H
#define FUNCTIONAL_FORM_H

#include "Flavor.cpp"
#include <functional>

class FunctionalFormInterface {
	using Signature = std::function<double(const FlavorType, const double, const double)>;
	public:
	const Signature function;
	const FlavorVector available_flavors;

	FunctionalFormInterface(const Signature _function, const FlavorVector _available_flavors = Flavor::all_flavors) : function(_function), available_flavors(_available_flavors), flavor_values(available_flavors.size(), 0) {

	}

	void evaluate(const double x, const double Q2) {
		for (auto const flavor : available_flavors) {
			const double function_value = function(flavor, x, Q2);
			flavor_values[size_t(flavor + 6)] = function_value;
		}
	}
	double xf_evaluate(const FlavorType flavor, const double x, const double Q2) {
		return function(flavor, x, Q2);
	}
	double xf(const FlavorType flavor) const {
		return flavor_values[size_t(flavor + 6)];
	}
	double xg() const {
		return xf(Flavor::Gluon);
	}
	double alpha_s(const double Q2) const {
		// fix this!
		return 1;
	}
	private:
	std::vector<double> flavor_values;
};

#endif