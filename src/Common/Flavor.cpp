#ifndef FLAVOR_H
#define FLAVOR_H

#include <vector>
#include <array>
#include <string>
#include "Utility/Utility.cpp"
#include "Common/Constants.cpp"

#define TOTAL_FLAVORS 13

using FlavorType = int;
using FlavorVector = std::vector<FlavorType>;

namespace Flavor {
	constexpr static const int Up = 2;
	constexpr static const int Down = 1;
	constexpr static const int Charm = 4;
	constexpr static const int Strange = 3;
	constexpr static const int Top = 6;
	constexpr static const int Bottom = 5;
	constexpr static const int Gluon = 0;
	constexpr static const int AntiUp = -2;
	constexpr static const int AntiDown = -1;
	constexpr static const int AntiCharm = -4;
	constexpr static const int AntiStrange = -3;
	constexpr static const int AntiTop = -6;
	constexpr static const int AntiBottom = -5;

	static FlavorVector all_upper_flavors = {Flavor::Up, Flavor::Charm, Flavor::Top, Flavor::AntiUp, Flavor::AntiCharm, Flavor::AntiTop};
	static FlavorVector all_lower_flavors = {Flavor::Down, Flavor::Strange, Flavor::Bottom, Flavor::AntiDown, Flavor::AntiStrange, Flavor::AntiBottom};
	static FlavorVector all_flavors = {-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6};

	constexpr FlavorType conjugate_flavor(const FlavorType flavor) {
		return -flavor;
	}
	constexpr FlavorVector conjugate_flavors(const FlavorVector &flavors) {
		FlavorVector conjugated(flavors.size());
		std::transform(flavors.begin(), flavors.end(), conjugated.begin(), conjugate_flavor);
		// for (auto const &flavor : flavors) {
		// 	conjugated.push_back(conjugate_flavor(flavor));
		// }
		return conjugated;
	}

	constexpr FlavorVector upper_flavors(FlavorVector &flavors) {
		return Collections::vector_intersection<FlavorType>(Flavor::all_upper_flavors, flavors);
	}
	constexpr FlavorVector lower_flavors(FlavorVector &flavors) {
		return Collections::vector_intersection<FlavorType>(Flavor::all_lower_flavors, flavors);
	}
	constexpr bool is_upper_flavor(const FlavorType flavor) {
		return flavor != 0 && flavor % 2 == 0;
	}
	constexpr bool is_lower_flavor(const FlavorType flavor) {
		return flavor % 2 == 1;
	}
	constexpr bool is_antiflavor(const FlavorType flavor) {
		return flavor < 0;
	}
	constexpr bool is_gluon(const FlavorType flavor) {
		return flavor == Flavor::Gluon;
	}

	constexpr FlavorType reflect_flavor(const FlavorType flavor) {
		if (is_antiflavor(flavor)) {
			return conjugate_flavor(reflect_flavor(conjugate_flavor(flavor)));
		}
		if (is_upper_flavor(flavor)) {
			return flavor - 1;
		}
		if (is_lower_flavor(flavor)) {
			return flavor + 1;
		}
		return flavor;
	}
	void reflect_flavors(FlavorVector &flavors) {
		for (FlavorType &flavor : flavors) {
			flavor = reflect_flavor(flavor);
		}
	}
};

struct FlavorInfo {
	FlavorVector active_flavors;
	FlavorVector active_antiflavors;

	FlavorVector upper_flavors;
	FlavorVector lower_flavors;
	FlavorVector upper_antiflavors;
	FlavorVector lower_antiflavors;

	const std::array<double, TOTAL_FLAVORS> flavor_masses;

	FlavorInfo(const FlavorVector _active_flavors, const std::array<double, TOTAL_FLAVORS> _masses = {}) 
	: active_flavors(_active_flavors), active_antiflavors(Flavor::conjugate_flavors(active_flavors)), flavor_masses(_masses) {
		upper_flavors = Flavor::upper_flavors(active_flavors);
		lower_flavors = Flavor::lower_flavors(active_flavors);
		upper_antiflavors = Flavor::upper_flavors(active_antiflavors);
		lower_antiflavors = Flavor::lower_flavors(active_antiflavors);
	}

	constexpr double mass(const FlavorType flavor) const {
		return flavor_masses[size_t(flavor + 6)];
	}
};

#endif