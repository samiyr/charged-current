#ifndef PROCESS_H
#define PROCESS_H

#include <iostream>

#include "Common/Particle.cpp"

// A structure encapsulating a reaction process.
struct Process {
	// Type of the process, either neutrino or antineutrino scattering.
	enum class Type {
		// neutrino + target --> lepton + X
		NeutrinoToLepton,
		// antineutrino + target --> antilepton + X
		AntiNeutrinoToAntiLepton, 
		// LeptonToNeutrino,
		// AntiLeptonToAntiNeutrino
	};

	// Process type, either neutrino of antineutrino scattering.
	Type type;
	// Target particle (i.e. the proton or a nucleus).
	Particle target;
	// Projectile particle (i.e. neutrino or antineutrino).
	Particle projectile;

	/// @brief Constructs a process.
	/// @param type Process type, either neutrino or antineutrino scattering.
	/// @param target Target particle (i.e. the proton or a nucleus).
	/// @param projectile Projectile particle (i.e. neutrino or antineutrino).
	constexpr Process(Type type, Particle target, Particle projectile) noexcept : type(type), target(target), projectile(projectile) { }

	// Gives the sign of the process (+1 for neutrino scattering, -1 for antineutrino scattering) used in the cross section formula;
	// see CommonFunctions::make_cross_section_variable.
	constexpr int W_sign() const noexcept {
		return type == Type::NeutrinoToLepton ? +1 : -1;
		// return (type == Type::NeutrinoToLepton || type == Type::AntiLeptonToAntiNeutrino) ? +1 : -1;
	}

	// Checks whether the sign of the process is positive or not; see Process::W_sign.
	constexpr bool positive_W() const noexcept {
		return W_sign() == +1;
	}

	friend std::ostream& operator<<(std::ostream &os, const Process &o) {
		std::string type_string;
		switch (o.type) {
		case Type::NeutrinoToLepton:
			type_string = "neutrino to lepton";
			break;
		case Type::AntiNeutrinoToAntiLepton:
			type_string = "antineutrino to antilepton";
			break;
		}

		return os << type_string << " [target = (" << o.target << ") | projectile = (" << o.projectile << ")]";
	}
};


#endif