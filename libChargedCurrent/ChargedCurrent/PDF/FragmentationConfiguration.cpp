#ifndef FRAGMENTATION_CONFIGURATION_H
#define FRAGMENTATION_CONFIGURATION_H

#include <vector>
#include <ranges>

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

	double Q2_min() const {
		const auto Q2_values = std::views::transform(interfaces, [](const Interface &interface) { return interface.Q2_min(); });
		return *std::ranges::min_element(Q2_values.begin(), Q2_values.end());
	}
	double Q2_max() const {
		const auto Q2_values = std::views::transform(interfaces, [](const Interface &interface) { return interface.Q2_max(); });
		return *std::ranges::max_element(Q2_values.begin(), Q2_values.end());
	}
};

#endif