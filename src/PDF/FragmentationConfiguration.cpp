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
	: interfaces(_interfaces), decays(interfaces.size(), Decay()) {}

	FragmentationConfiguration(const std::initializer_list<Interface> _interfaces) noexcept 
	: interfaces(_interfaces), decays(interfaces.size(), Decay()) {}

	FragmentationConfiguration(const std::vector<Interface> _interfaces, const std::vector<Decay<DecayFunction>> _decays) noexcept
	: interfaces(_interfaces), decays(_decays) {}

	FragmentationConfiguration(const std::initializer_list<Interface> _interfaces, const std::initializer_list<Decay<DecayFunction>> _decays) noexcept
	: interfaces(_interfaces), decays(_decays) {}

	void activate() const {
		if constexpr (PDFConceptActivation<Interface>) {
			for (const auto &interface : interfaces) {
				interface.activate();
			}
		}
	}
};


#endif