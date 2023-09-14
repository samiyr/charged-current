#ifndef FRAGMENTATION_CONFIGURATION_H
#define FRAGMENTATION_CONFIGURATION_H

#include <vector>

#include "Decay/Decay.cpp"
#include "Decay/DecayFunctions.cpp"

#include "PDFConcept.cpp"

template <is_pdf_interface Interface, is_decay_function DecayFunction>
struct FragmentationConfiguration {
	std::vector<Interface> interfaces;
	std::vector<Decay<DecayFunction>> decays;

	FragmentationConfiguration(const std::vector<Interface> _interfaces) noexcept 
	: interfaces(_interfaces), decays(interfaces.size(), Decay()) {}

	FragmentationConfiguration(const std::initializer_list<Interface> _interfaces) noexcept 
	: interfaces(_interfaces), decays(interfaces.size(), Decay()) {}

	FragmentationConfiguration(const std::vector<Interface> _interfaces, const std::vector<Decay<DecayFunction>> _decays) noexcept
	: interfaces(_interfaces), decays(_decays) {}

	FragmentationConfiguration(const std::initializer_list<Interface> _interfaces, const std::initializer_list<Decay<DecayFunction>> _decays) noexcept
	: interfaces(_interfaces), decays(_decays) {}

	void activate() const {
		if constexpr (has_pdf_activation<Interface>) {
			for (const auto &interface : interfaces) {
				interface.activate();
			}
		}
	}

	void evaluate(const double x, const double Q2) const {
		for (const auto &interface : interfaces) {
			interface.evaluate(x, Q2);
		}
	}
};

#endif