#ifndef ANALYSIS_DIS_H
#define ANALYSIS_DIS_H

#include <vector>

#include "Analysis/Parameters.cpp"
#include "Analysis/Constants.cpp"

#include "Common/Flavor.cpp"
#include "Common/Constants.cpp"
#include "Common/Process.cpp"
#include "Common/ScaleDependence.cpp"

#include "PDF/Interfaces/LHAInterface.cpp"

#include "DIS/DIS.cpp"

template <
	is_scale_dependence RenormalizationScale,
	is_scale_dependence FactorizationScale
>
struct DISAnalysis {
	AnalysisParameters params;
	const RenormalizationScale renormalization;
	const FactorizationScale factorization;

	DISAnalysis(const AnalysisParameters params, const RenormalizationScale renormalization, const FactorizationScale factorization) 
	: params(params), renormalization(renormalization), factorization(factorization) { }

	template <is_scale_dependence Renormalization, is_scale_dependence Factorization>
	void charm_production(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path filename, 
		const std::string comment,
		const Renormalization renormalization_scale,
		const Factorization factorization_scale
		) {
		DIS dis(
			{Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			LHAInterface(params.pdf_set),
			params.process,
			renormalization_scale, factorization_scale
		);

		dis.compute_differential_cross_section_directly = false;
		dis.use_modified_cross_section_prefactor = true;

		dis.charm_mass = params.charm_mass;

		dis.parallelize = params.parallelize;
		dis.number_of_threads = params.number_of_threads;

		dis.differential_cross_section_xy(x_bins, y_bins, E_beam_bins, filename, comment);
	}

	void charm_production(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path filename, 
		const std::string comment) {

		charm_production(x_bins, y_bins, E_beam_bins, filename, comment, renormalization, factorization);
	}

	template <is_scale_dependence Renormalization = RenormalizationScale, is_scale_dependence Factorization = FactorizationScale>
	void charm_production(
		const AnalysisSet set, const std::vector<double> x_bins, const std::filesystem::path filename, const std::string comment,
		const Renormalization renormalization_scale,
		const Factorization factorization_scale) {

		return charm_production(
			x_bins, AnalysisConstants::get_y_bins(set, params.process), AnalysisConstants::get_E_bins(set, params.process),
			filename, comment, renormalization_scale, factorization_scale
		);
	}

	void charm_production(const AnalysisSet set, const std::vector<double> x_bins, const std::filesystem::path filename, const std::string comment = "") {
		charm_production(set, x_bins, filename, comment, renormalization, factorization);
	}

	void charm_production_mass_scaling_comparison(
		const double mass1, const double mass2,
		const AnalysisSet set, const std::vector<double> x_bins,
		const std::filesystem::path folder,
		const std::filesystem::path filename1, const std::filesystem::path filename2,
		const std::string comment1 = "", const std::string comment2 = "") {

		const double current_mass = params.charm_mass;

		auto scale1 = ScaleDependence::Function([mass1](const TRFKinematics &kinematics) { return kinematics.Q2 + std::pow(mass1, 2); });
		params.charm_mass = mass1;
		charm_production(set, x_bins, folder / filename1, comment1, scale1, scale1);

		auto scale2 = ScaleDependence::Function([mass2](const TRFKinematics &kinematics) { return kinematics.Q2 + std::pow(mass2, 2); });
		params.charm_mass = mass2;
		charm_production(set, x_bins, folder / filename2, comment2, scale2, scale2);

		params.charm_mass = current_mass;
	}
};

#endif