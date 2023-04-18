#ifndef CKM_H
#define CKM_H

#include "Flavor.cpp"
#include "Constants.cpp"

struct CKM {
	constexpr static double squared(const FlavorType flavor1, FlavorType flavor2) {
		const FlavorType f1 = std::abs(flavor1);
		const FlavorType f2 = std::abs(flavor2);
		if (is_upper_flavor(f1) && is_lower_flavor(f2)) {
			return _matrix_element(f1, f2);
		} else if (is_lower_flavor(f1) && is_upper_flavor(f2)) {
			return _matrix_element(f2, f1);
		}
		throw std::runtime_error("CKM matrix elements are defined only for a pair of upper and lower flavors");
	}
	constexpr static double sum_of_squares(std::vector<FlavorType> flavors, FlavorType fixed) {
		double sum = 0.0;
		for (auto const &flavor : flavors) {
			const double square = CKM::squared(flavor, fixed);
			sum += square;
		}
		return sum;
	}
	private:
	constexpr static double _matrix_element(FlavorType upper, FlavorType lower) {
		const int upper_index = upper / 2 - 1;
		const int lower_index = lower / 2;

		return Constants::ckm_squared_matrix_elements[upper_index][lower_index];
	}
};


#endif