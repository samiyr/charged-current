#ifndef PROCESS_H
#define PROCESS_H

#include "Particle.cpp"

struct Process {
	enum class Type {
		NeutrinoToLepton, 
		AntiNeutrinoToAntiLepton, 
		LeptonToNeutrino,
		AntiLeptonToAntiNeutrino
	};

	const Type type;
	const Particle target;
	const Particle projectile;

	constexpr Process(Type type, Particle target, Particle projectile) noexcept : type(type), target(target), projectile(projectile) { }

	constexpr int W_sign() const noexcept {
		return (type == Type::NeutrinoToLepton || type == Type::AntiLeptonToAntiNeutrino) ? +1 : -1;
	}

	constexpr bool positive_W() const noexcept {
		return W_sign() == +1;
	}
};


#endif