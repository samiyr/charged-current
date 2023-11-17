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

	template <typename PDFInterface, is_scale_dependence Renormalization, is_scale_dependence Factorization>
	void charm_production(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path filename, 
		const std::string comment,
		const PDFInterface &pdf,
		const Renormalization renormalization_scale,
		const Factorization factorization_scale
		) {
		DIS dis(
			{Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			pdf,
			params.process,
			renormalization_scale, factorization_scale
		);

		dis.use_modified_cross_section_prefactor = true;

		dis.charm_mass = params.charm_mass;

		dis.parallelize = params.parallelize;
		dis.number_of_threads = params.number_of_threads;

		dis.primary_muon_min_energy = params.primary_muon_min_energy;
		dis.hadronic_min_energy = params.hadronic_min_energy;

		dis.freeze_factorization_scale = params.freeze_factorization;

		if constexpr (is_pdf_interface<PDFInterface>) {
			dis.differential_cross_section_xy(x_bins, y_bins, E_beam_bins, filename, comment);
		} else if constexpr (is_instance<PDFInterface, LHASetInterface>) {
			dis.differential_cross_section_xy_error_sets(x_bins, y_bins, E_beam_bins, filename, comment);
		}
	}

	void charm_production(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path filename, 
		const std::string comment) {

		const double Z = params.Z;
		const double A = params.A;

		if (params.pdf_error_sets) {
			if (params.explicit_isospin) {
				LHASetInterface<std::true_type> pdf(params.pdf_set);
				pdf.Z = Z;
				pdf.A = A;
				charm_production(x_bins, y_bins, E_beam_bins, filename, comment, pdf, renormalization, factorization);
			} else {
				charm_production(x_bins, y_bins, E_beam_bins, filename, comment, LHASetInterface<std::false_type>(params.pdf_set), renormalization, factorization);
			}
		} else {
			if (params.explicit_isospin) {
				LHAInterface<std::true_type> pdf(params.pdf_set);
				pdf.Z = Z;
				pdf.A = A;
				charm_production(x_bins, y_bins, E_beam_bins, filename, comment, pdf, renormalization, factorization);
			} else {
				charm_production(x_bins, y_bins, E_beam_bins, filename, comment, LHAInterface<std::false_type>(params.pdf_set), renormalization, factorization);
			}
		}
	}

	template <is_scale_dependence Renormalization = RenormalizationScale, is_scale_dependence Factorization = FactorizationScale>
	void charm_production(
		const AnalysisSet set, const std::vector<double> x_bins, const std::filesystem::path filename, const std::string comment,
		const Renormalization renormalization_scale,
		const Factorization factorization_scale) {

		const double Z = params.Z;
		const double A = params.A;

		if (params.pdf_error_sets) {
			if (params.explicit_isospin) {
				LHASetInterface<std::true_type> pdf(params.pdf_set);
				pdf.Z = Z;
				pdf.A = A;
				charm_production(
					x_bins, AnalysisConstants::get_y_bins(set, params.process), AnalysisConstants::get_E_bins(set, params.process),
					filename, comment, pdf, renormalization_scale, factorization_scale
				);
			} else {
				charm_production(
					x_bins, AnalysisConstants::get_y_bins(set, params.process), AnalysisConstants::get_E_bins(set, params.process),
					filename, comment, LHASetInterface<std::false_type>(params.pdf_set), renormalization_scale, factorization_scale
				);
			}
		} else {
			if (params.explicit_isospin) {
				LHAInterface<std::true_type> pdf(params.pdf_set);
				pdf.Z = Z;
				pdf.A = A;
				charm_production(
					x_bins, AnalysisConstants::get_y_bins(set, params.process), AnalysisConstants::get_E_bins(set, params.process),
					filename, comment, pdf, renormalization_scale, factorization_scale
				);
			} else {
				charm_production(
					x_bins, AnalysisConstants::get_y_bins(set, params.process), AnalysisConstants::get_E_bins(set, params.process),
					filename, comment, LHAInterface<std::false_type>(params.pdf_set), renormalization_scale, factorization_scale
				);
			}
		}
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

	template <typename PDFInterface, is_scale_dependence Renormalization, is_scale_dependence Factorization>
	void integrated(
		const std::vector<double> E_beam_bins,
		const double Q2_min,
		const std::filesystem::path filename, 
		const std::string comment,
		const PDFInterface &pdf,
		const Renormalization renormalization_scale,
		const Factorization factorization_scale
		) {
		DIS dis(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			pdf,
			params.process,
			renormalization_scale, factorization_scale
		);

		dis.use_modified_cross_section_prefactor = true;

		dis.charm_mass = params.charm_mass;

		dis.parallelize = params.parallelize;
		dis.number_of_threads = params.number_of_threads;

		dis.primary_muon_min_energy = params.primary_muon_min_energy;
		dis.hadronic_min_energy = params.hadronic_min_energy;

		dis.freeze_factorization_scale = params.freeze_factorization;

		if constexpr (is_pdf_interface<PDFInterface>) {
			dis.integrated_cross_section(E_beam_bins, Q2_min, filename, comment);
		} else if constexpr (is_instance<PDFInterface, LHASetInterface>) {
			dis.integrated_cross_section_error_sets(E_beam_bins, Q2_min, filename, comment);
		}
	}

	void integrated(
		const std::vector<double> E_beam_bins,
		const double Q2_min,
		const std::filesystem::path filename, 
		const std::string comment) {

		const double Z = params.Z;
		const double A = params.A;

		if (params.pdf_error_sets) {
			if (params.explicit_isospin) {
				LHASetInterface<std::true_type> pdf(params.pdf_set);
				pdf.Z = Z;
				pdf.A = A;
				integrated(E_beam_bins, Q2_min, filename, comment, pdf, renormalization, factorization);
			} else {
				integrated(E_beam_bins, Q2_min, filename, comment, LHASetInterface<std::false_type>(params.pdf_set), renormalization, factorization);
			}
		} else {
			if (params.explicit_isospin) {
				LHAInterface<std::true_type> pdf(params.pdf_set);
				pdf.Z = Z;
				pdf.A = A;
				integrated(E_beam_bins, Q2_min, filename, comment, pdf, renormalization, factorization);
			} else {
				integrated(E_beam_bins, Q2_min, filename, comment, LHAInterface<std::false_type>(params.pdf_set), renormalization, factorization);
			}
		}
	}

	void integrated(const AnalysisSet set, const std::filesystem::path filename, const std::string comment = "") {
		integrated(AnalysisConstants::get_E_bins(set, params.process), params.Q2_min, filename, comment);
	}

	void integrated(const std::vector<double> &E_beam_bins, const std::filesystem::path filename, const std::string comment = "") {
		integrated(E_beam_bins, params.Q2_min, filename, comment);
	}

	template <typename PDFInterface, is_scale_dependence Renormalization, is_scale_dependence Factorization>
	void integrated_charm_production(
		const std::vector<double> E_beam_bins,
		const double Q2_min,
		const std::filesystem::path filename, 
		const std::string comment,
		const PDFInterface &pdf,
		const Renormalization renormalization_scale,
		const Factorization factorization_scale
		) {
		DIS dis(
			{Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			pdf,
			params.process,
			renormalization_scale, factorization_scale
		);

		dis.use_modified_cross_section_prefactor = true;

		dis.charm_mass = params.charm_mass;

		dis.parallelize = params.parallelize;
		dis.number_of_threads = params.number_of_threads;

		dis.primary_muon_min_energy = params.primary_muon_min_energy;
		dis.hadronic_min_energy = params.hadronic_min_energy;

		dis.freeze_factorization_scale = params.freeze_factorization;

		if constexpr (is_pdf_interface<PDFInterface>) {
			dis.integrated_cross_section(E_beam_bins, Q2_min, filename, comment);
		} else if constexpr (is_instance<PDFInterface, LHASetInterface>) {
			dis.integrated_cross_section_error_sets(E_beam_bins, Q2_min, filename, comment);
		}
	}

	void integrated_charm_production(
		const std::vector<double> E_beam_bins,
		const double Q2_min,
		const std::filesystem::path filename, 
		const std::string comment) {

		const double Z = params.Z;
		const double A = params.A;

		if (params.pdf_error_sets) {
			if (params.explicit_isospin) {
				LHASetInterface<std::true_type> pdf(params.pdf_set);
				pdf.Z = Z;
				pdf.A = A;
				integrated_charm_production(E_beam_bins, Q2_min, filename, comment, pdf, renormalization, factorization);
			} else {
				integrated_charm_production(E_beam_bins, Q2_min, filename, comment, LHASetInterface<std::false_type>(params.pdf_set), renormalization, factorization);
			}
		} else {
			if (params.explicit_isospin) {
				LHAInterface<std::true_type> pdf(params.pdf_set);
				pdf.Z = Z;
				pdf.A = A;
				integrated_charm_production(E_beam_bins, Q2_min, filename, comment, pdf, renormalization, factorization);
			} else {
				integrated_charm_production(E_beam_bins, Q2_min, filename, comment, LHAInterface<std::false_type>(params.pdf_set), renormalization, factorization);
			}
		}
	}

	void integrated_charm_production(const AnalysisSet set, const std::filesystem::path filename, const std::string comment = "") {
		integrated_charm_production(AnalysisConstants::get_E_bins(set, params.process), params.Q2_min, filename, comment);
	}
	
	void integrated_charm_production(const std::vector<double> &E_beam_bins, const std::filesystem::path filename, const std::string comment = "") {
		integrated_charm_production(E_beam_bins, params.Q2_min, filename, comment);
	}
};

#endif