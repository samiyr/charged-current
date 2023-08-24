#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <filesystem>

#include "DIS/DIS.cpp"

#include "SIDIS/SIDIS.cpp"

#include "PDF/Interfaces/LHAInterface.cpp"

#include "Common/Particle.cpp"
#include "Common/Constants.cpp"

#include "Analysis/DIS.cpp"
#include "Analysis/SIDIS.cpp"
#include "Analysis/Parameters.cpp"

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
};

#endif