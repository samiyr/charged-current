#ifndef FRAGMENTATION_CONFIGURATION_H
#define FRAGMENTATION_CONFIGURATION_H

#include <vector>
#include "Decay.cpp"
#include "DecayFunctions.cpp"
#include "PDFConcept.cpp"

template <PDFConcept Interface, DecayFunctions::Concept DecayFunction>
struct FragmentationConfiguration {
	const std::vector<Interface> interfaces;
	const std::vector<Decay<DecayFunction>> decays;

	FragmentationConfiguration(const std::vector<Interface> _interfaces) : interfaces(_interfaces), decays(interfaces.size(), TrivialDecay) {}
	FragmentationConfiguration(const std::initializer_list<Interface> _interfaces) : interfaces(_interfaces), decays(interfaces.size(), TrivialDecay) {}
	FragmentationConfiguration(const std::vector<Interface> _interfaces, const std::vector<Decay<DecayFunction>> _decays) : interfaces(_interfaces), decays(_decays) {}
	FragmentationConfiguration(const std::initializer_list<Interface> _interfaces, const std::initializer_list<Decay<DecayFunction>> _decays) : interfaces(_interfaces), decays(_decays) {}

	void evaluate(const double x, const double Q2) const {
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

	const Interface &operator[](const std::vector<Interface>::size_type index) const {
		return interfaces[index];
	}
};


#endif