#ifndef ANALYSIS_PARAMETERS_H
#define ANALYSIS_PARAMETERS_H

#include <string>
#include "Common/Process.cpp"
#include "Common/Constants.cpp"

struct AnalysisParameters {
	std::string pdf_set = "EPPS21nlo_CT18Anlo_Fe56";

	Process process = Process(Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino);

	double charm_mass = Constants::Charm::Mass;
	double minimum_lepton_momentum = 5.0;
};

#endif