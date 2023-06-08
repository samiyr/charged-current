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

	constexpr int W_sign() const {
		return (type == Type::NeutrinoToLepton || type == Type::AntiLeptonToAntiNeutrino) ? +1 : -1;
	}

	constexpr bool positive_W() const {
		return W_sign() == +1;
	}
};


#endif