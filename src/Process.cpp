#ifndef PROCESS_H
#define PROCESS_H

struct Process {
	enum class Type {
		NeutrinoToLepton, 
		AntiNeutrinoToAntiLepton, 
		LeptonToNeutrino,
		AntiLeptonToAntiNeutrino
	};

	const Type type;

	int W_sign() const {
		return (type == Type::NeutrinoToLepton || type == Type::AntiLeptonToAntiNeutrino) ? +1 : -1;
	}

	bool positive_W() const {
		return W_sign() == +1;
	}
};


#endif