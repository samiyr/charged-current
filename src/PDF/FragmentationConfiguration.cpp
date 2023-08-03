#ifndef FRAGMENTATION_CONFIGURATION_H
#define FRAGMENTATION_CONFIGURATION_H

#include <vector>
#include "Decay/Decay.cpp"
#include "Decay/DecayFunctions.cpp"
#include "PDFConcept.cpp"

template <PDFConcept Interface, DecayFunctions::Concept DecayFunction>
struct FragmentationConfiguration {
	const std::vector<Interface> interfaces;
	const std::vector<Decay<DecayFunction>> decays;

	FragmentationConfiguration(const std::vector<Interface> _interfaces) noexcept 
	: interfaces(_interfaces), decays(interfaces.size(), TrivialDecay) {}

	FragmentationConfiguration(const std::initializer_list<Interface> _interfaces) noexcept 
	: interfaces(_interfaces), decays(interfaces.size(), TrivialDecay) {}

	FragmentationConfiguration(const std::vector<Interface> _interfaces, const std::vector<Decay<DecayFunction>> _decays) noexcept
	: interfaces(_interfaces), decays(_decays) {}

	FragmentationConfiguration(const std::initializer_list<Interface> _interfaces, const std::initializer_list<Decay<DecayFunction>> _decays) noexcept
	: interfaces(_interfaces), decays(_decays) {}

	void evaluate(const double x, const double Q2) const {
		for (auto &interface : interfaces) {
			interface.evaluate(x, Q2);
		}
	}

	void activate() const {
		if constexpr (PDFConceptActivation<Interface>) {
			for (auto &interface : interfaces) {
				interface.activate();
			}
		}
	}

	typename std::vector<Interface>::iterator begin() {
		return interfaces.begin();
    }
    typename std::vector<Interface>::iterator end() {
    	return interfaces.end();
    }
    typename std::vector<Interface>::const_iterator cbegin() const {
    	return interfaces.cbegin();
    }
    typename std::vector<Interface>::const_iterator cend() const {
    	return interfaces.cend();
    }

	typename std::vector<Interface>::size_type size() const {
		return interfaces.size();
	}

	const Interface &operator[](const typename std::vector<Interface>::size_type index) const {
		return interfaces[index];
	}
};


#endif