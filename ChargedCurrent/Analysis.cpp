#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <filesystem>

#include <ChargedCurrent/DIS/DIS.cpp>

#include <ChargedCurrent/SIDIS/SIDIS.cpp>

#include <ChargedCurrent/PDF/Interfaces/LHAInterface.cpp>

#include <ChargedCurrent/Common/Particle.cpp>
#include <ChargedCurrent/Common/Constants.cpp>

#include "DIS.cpp"
#include "SIDIS.cpp"
#include "Utility.cpp"
#include "Parameters.cpp"

template <
	is_scale_dependence RenormalizationScale = decltype(ScaleDependence::trivial),
	is_scale_dependence FactorizationScale = decltype(ScaleDependence::trivial),
	is_scale_dependence FragmentationScale = decltype(ScaleDependence::trivial)
>
struct Analysis {
	AnalysisParameters params;
	const RenormalizationScale renormalization;
	const FactorizationScale factorization;
	const FragmentationScale fragmentation;

	Analysis(
		const RenormalizationScale renormalization = ScaleDependence::trivial, 
		const FactorizationScale factorization = ScaleDependence::trivial, 
		const FragmentationScale fragmentation = ScaleDependence::trivial)
	: renormalization(renormalization), factorization(factorization), fragmentation(fragmentation) { }

	auto dis() {
		return DISAnalysis(params, renormalization, factorization);
	}
	auto sidis() {
		return SIDISAnalysis(params, renormalization, factorization, fragmentation);
	}
	auto utility() {
		return UtilityAnalysis(params);
	}
};

#endif