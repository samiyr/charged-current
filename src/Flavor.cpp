#ifndef FLAVOR_H
#define FLAVOR_H

#include <vector>
#include <string>
#include "Utility.cpp"

using FlavorType = int;
using FlavorVector = std::vector<FlavorType>;

namespace Flavor {
	static const int Up = 2;
	static const int Down = 1;
	static const int Charm = 4;
	static const int Strange = 3;
	static const int Top = 6;
	static const int Bottom = 5;
	static const int Gluon = 0;
	static const int AntiUp = -2;
	static const int AntiDown = -1;
	static const int AntiCharm = -4;
	static const int AntiStrange = -3;
	static const int AntiTop = -6;
	static const int AntiBottom = -5;
};

constexpr FlavorType conjugate_flavor(const FlavorType flavor) {
	return -flavor;
}
constexpr std::vector<FlavorType> conjugate_flavors(const std::vector<FlavorType> &flavors) {
	std::vector<FlavorType> conjugated;
	for (auto const &flavor : flavors) {
		conjugated.push_back(conjugate_flavor(flavor));
	}
	return conjugated;
}

constexpr std::vector<FlavorType> upper_flavors(const std::vector<FlavorType> &flavors) {
	return vector_intersection<FlavorType>({Flavor::Up, Flavor::Charm, Flavor::Top, Flavor::AntiUp, Flavor::AntiCharm, Flavor::AntiTop}, flavors);
}
constexpr std::vector<FlavorType> lower_flavors(const std::vector<FlavorType> &flavors) {
	return vector_intersection<FlavorType>({Flavor::Down, Flavor::Strange, Flavor::Bottom, Flavor::AntiDown, Flavor::AntiStrange, Flavor::AntiBottom}, flavors);
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
	throw std::runtime_error("Cannot reflect flavor");
}
void reflect_flavors(std::vector<FlavorType> &flavors) {
	for (int i = 0; i < flavors.size(); i++) {
		flavors[i] = reflect_flavor(flavors[i]);
	}
}

#endif