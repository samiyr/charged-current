#ifndef FLAVOR_H
#define FLAVOR_H

#include <vector>
#include <array>
#include <string>

#include "Common/Constants.cpp"
#include "Utility/Utility.cpp"

// Total number of flavors, including the gluon and antiquark flavors.
static constexpr std::size_t TOTAL_FLAVORS = 13;

using FlavorType = int;
using FlavorVector = std::vector<FlavorType>;

// Provides access and convenience methods for handling quark flavors (+ gluon) using the MC particle identification scheme.
// An exception is made for the gluon, which has the ID 0 instead of 21, since this avoids exceptions for the gluon and is supported
// by LHAPDF.
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

	// Upper flavors: u, c, t, ubar, cbar, tbar.
	static FlavorVector all_upper_flavors = {Flavor::Up, Flavor::Charm, Flavor::Top, Flavor::AntiUp, Flavor::AntiCharm, Flavor::AntiTop};
	// Lower flavors: d, s, b, dbar, sbar, bbar.
	static FlavorVector all_lower_flavors = {Flavor::Down, Flavor::Strange, Flavor::Bottom, Flavor::AntiDown, Flavor::AntiStrange, Flavor::AntiBottom};
	// All quark and antiquark flavors, as well as the gluon.
	static FlavorVector all_flavors = {-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6};

	// Converts the given flavor to its conjugate flavor, i.e. turn a quark flavor to an antiquark flavor and vice versa.
	// Does nothing to a gluon.
	static constexpr FlavorType conjugate_flavor(const FlavorType flavor) noexcept {
		return -flavor;
	}
	// Conjugates a vector of flavor values; see Flavor::conjugate_flavor.
	static constexpr FlavorVector conjugate_flavors(const FlavorVector &flavors) noexcept {
		FlavorVector conjugated(flavors.size());
		std::transform(flavors.begin(), flavors.end(), conjugated.begin(), conjugate_flavor);
		return conjugated;
	}

	// Returns a flavor vector containing only upper flavors from the given flavor vector; see Flavor::all_upper_flavors. Order is not preserved.
	constexpr FlavorVector upper_flavors(FlavorVector &flavors) {
		return Collections::intersection(Flavor::all_upper_flavors, flavors);
	}
	// Returns a flavor vector containing only lower flavors from the given flavor vector; see Flavor::all_lower_flavors. Order is not preserved.
	constexpr FlavorVector lower_flavors(FlavorVector &flavors) {
		return Collections::intersection(Flavor::all_lower_flavors, flavors);
	}

	// Checks whether a given flavor is an upper flavor; see Flavor::all_upper_flavors. A gluon is not an upper flavor.
	constexpr bool is_upper_flavor(const FlavorType flavor) noexcept {
		return flavor != 0 && flavor % 2 == 0;
	}
	// Checks whether a given flavor is a lower flavor; see Flavor::all_lower_flavors. A gluon is not a lower flavor.
	constexpr bool is_lower_flavor(const FlavorType flavor) noexcept {
		return flavor % 2 == 1;
	}
	// Checks whether a given flavor is an antiflavor. A gluon is not an antiflavor.
	constexpr bool is_antiflavor(const FlavorType flavor) noexcept {
		return flavor < 0;
	}
	// Checks whether the given flavor is the gluon.
	constexpr bool is_gluon(const FlavorType flavor) noexcept {
		return flavor == Flavor::Gluon;
	}

	// Reflects a given flavor. Reflecting an upper flavor turns it into a lower flavor and vice versa, while maintaining the flavor/antiflavor signedness.
	// Does nothing to the gluon.
	static inline constexpr FlavorType reflect_flavor(const FlavorType flavor) noexcept {
		// If the flavor is an antiflavor, reflect the conjugation (i.e. non-antiflavor) and then conjugate back to an antiflavor.
		if (is_antiflavor(flavor)) {
			return conjugate_flavor(reflect_flavor(conjugate_flavor(flavor)));
		}
		// Turn an upper flavor into a lower flavor by subtracting 1 from the PID, since u = 2 --> d = 1, etc.
		if (is_upper_flavor(flavor)) {
			return flavor - 1;
		}
		// Turn a lower flavor into an upper flavor by adding 1 to the PID, since d = 1 --> u = 2, etc.
		if (is_lower_flavor(flavor)) {
			return flavor + 1;
		}
		// In case the flavor is not an antiflavor or an upper or lower flavor, return the flavor unchanged.
		// This means that the flavor is the gluon.
		return flavor;
	}
	// Reflects a given flavor vector in place, modifying the argument. See Flavor::reflect_flavor.
	[[maybe_unused]] static inline void reflect_flavors(FlavorVector &flavors) noexcept {
		for (FlavorType &flavor : flavors) {
			flavor = reflect_flavor(flavor);
		}
	}
};

// A structure encapsulating the necessary information needed for flavor summation, as well as precomputed convenience variables.
struct FlavorInfo {
	// A vector of active quark flavors. Does not contain antiflavors or the gluon.
	FlavorVector active_flavors;
	// A vector of active antiquark flavors, computed from active_flavors.
	FlavorVector active_antiflavors;

	// A vector of active upper quark flavors, computed from active_flavors.
	FlavorVector upper_flavors;
	// A vector of active lower quark flavors, computed from active_flavors.
	FlavorVector lower_flavors;
	// A vector of active upper antiquark flavors, computed from active_flavors.
	FlavorVector upper_antiflavors;
	// A vector of active lower antiquark flavors, computed from active_flavors.
	FlavorVector lower_antiflavors;

	// An array of quark masses (+ gluon), given in the order tbar, bbar, ..., gluon, d, u, ..., t.
	const std::array<double, TOTAL_FLAVORS> flavor_masses;

	/// @brief Initializes with the given active flavors and quark masses.
	/// @param active_flavors A vector of active quark flavors. Must not contain antiflavors or the gluon.
	/// @param masses An array of quark masses in the order tbar, bbar, ..., gluon, d, u, ..., t. Must contain all 13 entries, even if not all are active flavors.
	FlavorInfo(const FlavorVector _active_flavors, const std::array<double, TOTAL_FLAVORS> _masses) 
	: active_flavors(_active_flavors), active_antiflavors(Flavor::conjugate_flavors(active_flavors)), flavor_masses(_masses) {
		upper_flavors = Flavor::upper_flavors(active_flavors);
		lower_flavors = Flavor::lower_flavors(active_flavors);
		upper_antiflavors = Flavor::upper_flavors(active_antiflavors);
		lower_antiflavors = Flavor::lower_flavors(active_antiflavors);
	}

	// Returns the assigned quark mass for the given flavor. The flavor must be in the range [-6, 6], otherwise results in undefined behaviour.
	constexpr double mass(const FlavorType flavor) const {
		return flavor_masses[std::size_t(flavor + 6)];
	}
};

#endif