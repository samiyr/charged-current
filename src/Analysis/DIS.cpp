#ifndef ANALYSIS_DIS_H
#define ANALYSIS_DIS_H

#include <vector>
#include "Analysis/Parameters.cpp"
#include "Analysis/Constants.cpp"
#include "Common/Flavor.cpp"
#include "Common/Constants.cpp"
#include "Common/Process.cpp"
#include "PDF/Interfaces/LHAInterface.cpp"
#include "DIS/DIS.cpp"

struct DISAnalysis {
	const AnalysisParameters params;

	DISAnalysis(const AnalysisParameters params) : params(params) { }

	void charm_production(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path filename, 
		const std::string comment = "") {
		DIS dis(
			{Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			LHAInterface(params.pdf_set),
			params.process
		);

		dis.compute_differential_cross_section_directly = false;
		dis.use_modified_cross_section_prefactor = true;

		dis.charm_mass = params.charm_mass;

		dis.differential_cross_section_xy(x_bins, y_bins, E_beam_bins, filename, comment);
	}

	void charm_production(const AnalysisSet set, const std::vector<double> x_bins, const std::filesystem::path filename, const std::string comment = "") {
		switch (set) {
		case AnalysisSet::NuTeV:
			return charm_production(x_bins, AnalysisConstants::NuTeV::New::y_bins, AnalysisConstants::NuTeV::New::E_bins, filename, comment);
		case AnalysisSet::NuTeV_old:
			return charm_production(x_bins, AnalysisConstants::NuTeV::Old::y_bins, AnalysisConstants::NuTeV::Old::E_bins, filename, comment);
		case AnalysisSet::CCFR:
			return charm_production(x_bins, AnalysisConstants::CCFR::y_bins, AnalysisConstants::CCFR::E_bins, filename, comment);
		}
	}
};

#endif