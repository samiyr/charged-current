#ifndef FUNCTIONAL_FORM_H
#define FUNCTIONAL_FORM_H

#include <functional>

#include "Common/Flavor.cpp"

template <typename Signature>
class FunctionalFormInterface {
	public:
	const Signature function;
	const FlavorVector available_flavors;

	const std::string set_name = "functional_form";
	const int set_member_number = 0;

	constexpr FunctionalFormInterface(const Signature _function, const FlavorVector _available_flavors = Flavor::all_flavors) noexcept
	: function(_function), available_flavors(_available_flavors), flavor_values(available_flavors.size(), 0) { }

	void evaluate(const double x, const double Q2) const {
		for (auto const flavor : available_flavors) {
			const double function_value = function(flavor, x, Q2);
			flavor_values[std::size_t(flavor + 6)] = function_value;
		}
	}
	constexpr double xf_evaluate(const FlavorType flavor, const double x, const double Q2) const {
		return function(flavor, x, Q2);
	}
	constexpr double xf(const FlavorType flavor) const {
		return flavor_values[std::size_t(flavor + 6)];
	}
	constexpr double xg() const {
		return xf(Flavor::Gluon);
	}
	constexpr double alpha_s([[maybe_unused]] const double Q2) const {
		return 1.0;
	}
	constexpr double quark_mass(const FlavorType) const {
		return 0.0;
	}
	auto size() const { return 1; }
	using size_type = int;

	auto operator[](int) const { return *this; }

	private:
	mutable std::vector<double> flavor_values;
};

#endif