#ifndef FRAGMENTATION_CONFIGURATION_H
#define FRAGMENTATION_CONFIGURATION_H

#include <vector>

template <typename Interface>
struct FragmentationConfiguration {
	std::vector<Interface> interfaces;
	const std::vector<double> coefficients;

	FragmentationConfiguration(const std::vector<Interface> _interfaces) : interfaces(_interfaces), coefficients(interfaces.size(), 1.0) {}
	FragmentationConfiguration(const std::initializer_list<Interface> _interfaces) : interfaces(_interfaces), coefficients(interfaces.size(), 1.0) {}
	FragmentationConfiguration(const std::vector<Interface> _interfaces, const std::vector<double> _coefficients) : interfaces(_interfaces), coefficients(_coefficients) {}
	FragmentationConfiguration(const std::initializer_list<Interface> _interfaces, const std::vector<double> _coefficients) : interfaces(_interfaces), coefficients(_coefficients) {}

	void evaluate(const double x, const double Q2) {
		for (auto &interface : interfaces) {
			interface.evaluate(x, Q2);
		}
	}

	std::vector<Interface>::iterator begin() {
		return interfaces.begin();
    }
    std::vector<Interface>::iterator end() {
    	return interfaces.end();
    }
    std::vector<Interface>::const_iterator begin() const {
    	return interfaces.begin();
    }
    std::vector<Interface>::const_iterator end() const {
    	return interfaces.end();
    }

	std::vector<Interface>::size_type size() const {
		return interfaces.size();
	}

	Interface &operator[](const std::vector<Interface>::size_type index) {
		return interfaces[index];
	}
};


#endif