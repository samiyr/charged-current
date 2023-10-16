#ifndef ANALYSIS_SIDIS_H
#define ANALYSIS_SIDIS_H

#include <vector>

#include "Analysis/Parameters.cpp"
#include "Analysis/Constants.cpp"

#include "Common/Flavor.cpp"
#include "Common/Constants.cpp"
#include "Common/Process.cpp"

#include "PDF/Interfaces/LHAInterface.cpp"
#include "PDF/Interfaces/LHASetInterface.cpp"

#include "SIDIS/SIDIS.cpp"

template <
	is_scale_dependence RenormalizationScale,
	is_scale_dependence FactorizationScale,
	is_scale_dependence FragmentationScale
>
struct SIDISAnalysis {
	const AnalysisParameters params;
	const RenormalizationScale renormalization;
	const FactorizationScale factorization;
	const FragmentationScale fragmentation;

	SIDISAnalysis(
		const AnalysisParameters params, 
		const RenormalizationScale renormalization, const FactorizationScale factorization, const FragmentationScale fragmentation
	) : params(params), renormalization(renormalization), factorization(factorization), fragmentation(fragmentation) { }

	template <typename PDFInterface, is_pdf_interface FFInterface, is_decay_function DecayFunction>
	void muon_pair_production(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path filename, 
		const PDFInterface &pdf,
		const FragmentationConfiguration<FFInterface, DecayFunction> &ff,
		const std::string comment = "") requires is_pdf_interface<PDFInterface> || is_instance<PDFInterface, LHASetInterface> {
		SIDIS sidis(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			pdf, ff,
			params.process,
			renormalization, factorization, fragmentation
		);

		sidis.charm_mass = params.charm_mass;

		sidis.combine_integrals = true;
		sidis.use_modified_cross_section_prefactor = true;
		sidis.scale_variation = params.scale_variation;
		sidis.order = params.order;
		sidis.use_nlp_nlo = params.use_nlp_nlo;

		sidis.integration_parameters = params.integration;

		sidis.decay_variations = params.decay_parametrization_set;
		sidis.decay_variation = params.decay_variation;

		sidis.parallelize = params.parallelize;
		sidis.number_of_threads = params.number_of_threads;

		sidis.freeze_factorization_scale = params.freeze_factorization;
		sidis.freeze_fragmentation_scale = params.freeze_fragmentation;

		if constexpr (is_pdf_interface<PDFInterface>) {
			sidis.lepton_pair_cross_section_xy(x_bins, y_bins, E_beam_bins, filename, comment);
		} else if constexpr (is_instance<PDFInterface, LHASetInterface>) {
			sidis.lepton_pair_cross_section_xy_error_sets(x_bins, y_bins, E_beam_bins, filename, comment);
		}
	}

	void muon_pair_production(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path filename, 
		const std::string comment = "") {
		
		const double minimum_lepton_momentum = params.minimum_lepton_momentum;
		const Particle target = Constants::Particles::Proton;
		const auto decay_function = DecayFunctions::decay_function;

		const DecayParametrization parametrization = params.decay_parametrization;

		const auto muon_pair_production_lambda = [&]<typename PDF>(PDF &&pdf) {
			muon_pair_production(
				x_bins, y_bins, E_beam_bins, filename,
				pdf,
				FragmentationConfiguration(
					{
						LHAInterface("kkks08_opal_d0___mas", {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
						LHAInterface("kkks08_opal_d+___mas", {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
						LHAInterface("bkk05_D3_d_s_nlo", {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
						1.14 * LHAInterface("bkk05_D3_lambda_c_nlo", {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0})
					},
					{
						Decay(parametrization, Constants::Particles::D0, target, decay_function, minimum_lepton_momentum),
						Decay(parametrization, Constants::Particles::Dp, target, decay_function, minimum_lepton_momentum),
						Decay(parametrization, Constants::Particles::Ds, target, decay_function, minimum_lepton_momentum),
						Decay(parametrization, Constants::Particles::LambdaC, target, decay_function, minimum_lepton_momentum)
					}
				),
				comment
			);
		};

		const double Z = params.Z;
		const double A = params.A;

		if (params.pdf_error_sets) {
			if (params.explicit_isospin) {
				LHASetInterface<std::true_type> pdf(params.pdf_set);
				pdf.Z = Z;
				pdf.A = A;
				muon_pair_production_lambda(pdf);
			} else {
				muon_pair_production_lambda(LHASetInterface<std::false_type>(params.pdf_set));
			}
		} else {
			if (params.explicit_isospin) {
				LHAInterface<std::true_type> pdf(params.pdf_set);
				pdf.Z = Z;
				pdf.A = A;
				muon_pair_production_lambda(pdf);
			} else {
				muon_pair_production_lambda(LHAInterface<std::false_type>(params.pdf_set));
			}
		}
	}

	void muon_pair_production(const AnalysisSet set, const std::vector<double> x_bins, const std::filesystem::path filename, const std::string comment = "") {
		return muon_pair_production(x_bins, AnalysisConstants::get_y_bins(set, params.process), AnalysisConstants::get_E_bins(set, params.process), filename, comment);
	}

	void muon_pair_production_only_D0(const AnalysisSet set, const std::vector<double> x_bins, const std::filesystem::path filename, const std::string comment = "") {
		return muon_pair_production_only_D0(x_bins, AnalysisConstants::get_y_bins(set, params.process), AnalysisConstants::get_E_bins(set, params.process), filename, comment);
	}

	/*
		INTEGRATED MUON PAIR PRODUCTION
	*/

	template <typename PDFInterface, is_pdf_interface FFInterface, is_decay_function DecayFunction>
	void integrated_muon_pair_production(
		const std::vector<double> E_beam_bins,
		const double Q2_min,
		const std::filesystem::path filename, 
		const PDFInterface &pdf,
		const FragmentationConfiguration<FFInterface, DecayFunction> &ff,
		const std::string comment = "") requires is_pdf_interface<PDFInterface> || is_instance<PDFInterface, LHASetInterface> {
		SIDIS sidis(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			pdf, ff,
			params.process,
			renormalization, factorization, fragmentation
		);

		sidis.charm_mass = params.charm_mass;

		sidis.combine_integrals = true;
		sidis.use_modified_cross_section_prefactor = true;
		sidis.scale_variation = params.scale_variation;
		sidis.order = params.order;
		sidis.use_nlp_nlo = params.use_nlp_nlo;

		sidis.integration_parameters = params.integration;

		sidis.decay_variations = params.decay_parametrization_set;
		sidis.decay_variation = params.decay_variation;

		sidis.parallelize = params.parallelize;
		sidis.number_of_threads = params.number_of_threads;

		sidis.freeze_factorization_scale = params.freeze_factorization;
		sidis.freeze_fragmentation_scale = params.freeze_fragmentation;

		if constexpr (is_pdf_interface<PDFInterface>) {
			sidis.integrated_lepton_pair_cross_section(E_beam_bins, Q2_min, filename, comment);
		} else if constexpr (is_instance<PDFInterface, LHASetInterface>) {
			sidis.integrated_lepton_pair_cross_section_error_sets(E_beam_bins, Q2_min, filename, comment);
		}
	}

	void integrated_muon_pair_production(
		const std::vector<double> E_beam_bins,
		const double Q2_min,
		const std::filesystem::path filename, 
		const std::string comment = "") {
		
		const double minimum_lepton_momentum = params.minimum_lepton_momentum;
		const Particle target = Constants::Particles::Proton;
		const auto decay_function = DecayFunctions::decay_function;

		const DecayParametrization parametrization = params.decay_parametrization;

		const auto muon_pair_production_lambda = [&]<typename PDF>(PDF &&pdf) {
			integrated_muon_pair_production(
				E_beam_bins, Q2_min, filename,
				pdf,
				FragmentationConfiguration(
					{
						LHAInterface("kkks08_opal_d0___mas", {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
						LHAInterface("kkks08_opal_d+___mas", {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
						LHAInterface("bkk05_D3_d_s_nlo", {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
						1.14 * LHAInterface("bkk05_D3_lambda_c_nlo", {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0})
					},
					{
						Decay(parametrization, Constants::Particles::D0, target, decay_function, minimum_lepton_momentum),
						Decay(parametrization, Constants::Particles::Dp, target, decay_function, minimum_lepton_momentum),
						Decay(parametrization, Constants::Particles::Ds, target, decay_function, minimum_lepton_momentum),
						Decay(parametrization, Constants::Particles::LambdaC, target, decay_function, minimum_lepton_momentum)
					}
				),
				comment
			);
		};

		const double Z = params.Z;
		const double A = params.A;

		if (params.pdf_error_sets) {
			if (params.explicit_isospin) {
				LHASetInterface<std::true_type> pdf(params.pdf_set);
				pdf.Z = Z;
				pdf.A = A;
				muon_pair_production_lambda(pdf);
			} else {
				muon_pair_production_lambda(LHASetInterface<std::false_type>(params.pdf_set));
			}
		} else {
			if (params.explicit_isospin) {
				LHAInterface<std::true_type> pdf(params.pdf_set);
				pdf.Z = Z;
				pdf.A = A;
				muon_pair_production_lambda(pdf);
			} else {
				muon_pair_production_lambda(LHAInterface<std::false_type>(params.pdf_set));
			}
		}
	}

	void integrated_muon_pair_production(const AnalysisSet set, const std::filesystem::path filename, const std::string comment = "") {
		return integrated_muon_pair_production(AnalysisConstants::get_E_bins(set, params.process), params.Q2_min, filename, comment);
	}

	/*
		FLAVOR CHANNEL DECOMPOSITION
	*/

	void muon_pair_production_separated_channels(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path filename,
		const std::vector<bool> incoming,
		const std::vector<bool> outgoing,
		const std::string comment = "") {
		
		const double minimum_lepton_momentum = params.minimum_lepton_momentum;
		const Particle target = params.process.target;
		const auto decay_function = DecayFunctions::decay_function;

		const DecayParametrization parametrization = DecayParametrization::fit1();

		const std::vector<double> incoming_multipliers(incoming.begin(), incoming.end());
		const std::vector<double> outgoing_multipliers(outgoing.begin(), outgoing.end());

		std::vector<double> multipliers(outgoing.begin(), outgoing.end());
		multipliers *= std::vector{1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0};

		muon_pair_production(
			x_bins, y_bins, E_beam_bins, filename,
			LHAInterface(params.pdf_set, incoming_multipliers),
			FragmentationConfiguration(
				{
					LHAInterface("kkks08_opal_d0___mas", multipliers), 
					LHAInterface("kkks08_opal_d+___mas", multipliers), 
					LHAInterface("bkk05_D3_d_s_nlo", multipliers), 
					1.14 * LHAInterface("bkk05_D3_lambda_c_nlo", multipliers)
				},
				{
					Decay(parametrization, Constants::Particles::D0, target, decay_function, minimum_lepton_momentum),
					Decay(parametrization, Constants::Particles::Dp, target, decay_function, minimum_lepton_momentum),
					Decay(parametrization, Constants::Particles::Ds, target, decay_function, minimum_lepton_momentum),
					Decay(parametrization, Constants::Particles::LambdaC, target, decay_function, minimum_lepton_momentum)
				}
			),
			comment
		);
	}
	void muon_pair_production_quark_to_quark(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path filename,
		const std::string comment = "") {

		muon_pair_production_separated_channels(x_bins, y_bins, E_beam_bins, filename,
		{
			true, true, true, true, true, true, false, true, true, true, true, true, true
		},
		{
			true, true, true, true, true, true, false, true, true, true, true, true, true
		}, comment);
	}
	void muon_pair_production_quark_to_gluon(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path filename,
		const std::string comment = "") {

		muon_pair_production_separated_channels(x_bins, y_bins, E_beam_bins, filename,
		{
			true, true, true, true, true, true, false, true, true, true, true, true, true
		},
		{
			false, false, false, false, false, false, true, false, false, false, false, false, false
		}, comment);
	}
	void muon_pair_production_gluon_to_quark(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path filename,
		const std::string comment = "") {

		muon_pair_production_separated_channels(x_bins, y_bins, E_beam_bins, filename,
		{
			false, false, false, false, false, false, true, false, false, false, false, false, false
		},
		{
			true, true, true, true, true, true, false, true, true, true, true, true, true
		}, comment);
	}
	void muon_pair_production_gluon_to_gluon(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path filename,
		const std::string comment = "") {

		muon_pair_production_separated_channels(x_bins, y_bins, E_beam_bins, filename,
		{
			false, false, false, false, false, false, true, false, false, false, false, false, false
		},
		{
			false, false, false, false, false, false, true, false, false, false, false, false, false
		}, comment);
	}
	void muon_pair_production_quark_gluon_channels(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::vector<std::filesystem::path> filenames) {
		
		muon_pair_production_quark_to_quark(x_bins, y_bins, E_beam_bins, filenames[0], "quark to quark");
		muon_pair_production_quark_to_gluon(x_bins, y_bins, E_beam_bins, filenames[1], "quark to gluon");
		muon_pair_production_gluon_to_quark(x_bins, y_bins, E_beam_bins, filenames[2], "gluon to quark");
		muon_pair_production_gluon_to_gluon(x_bins, y_bins, E_beam_bins, filenames[3], "gluon to gluon");
	}
	
	void muon_pair_production_flavor_to_flavor(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins,
		const FlavorType incoming,
		const FlavorType outgoing,
		const std::filesystem::path filename,
		const std::string comment = "") {

		std::vector<bool> incoming_vector(TOTAL_FLAVORS, false);
		std::vector<bool> outgoing_vector(TOTAL_FLAVORS, false);

		incoming_vector[std::size_t(incoming + 6)] = true;
		outgoing_vector[std::size_t(outgoing + 6)] = true;

		muon_pair_production_separated_channels(x_bins, y_bins, E_beam_bins, filename, incoming_vector, outgoing_vector, comment);
	}

	void muon_pair_production_flavor_decomposition_quark_to_quark(
		const std::vector<double> x_bins,
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins,
		const std::string base_filename) {

		muon_pair_production_flavor_to_flavor(x_bins, y_bins, E_beam_bins, Flavor::Down, Flavor::Up, base_filename + "_du.csv", "d -> u");
		muon_pair_production_flavor_to_flavor(x_bins, y_bins, E_beam_bins, Flavor::Down, Flavor::Charm, base_filename + "_dc.csv", "d -> c");
		muon_pair_production_flavor_to_flavor(x_bins, y_bins, E_beam_bins, Flavor::Strange, Flavor::Up, base_filename + "_su.csv", "s -> u");
		muon_pair_production_flavor_to_flavor(x_bins, y_bins, E_beam_bins, Flavor::Strange, Flavor::Charm, base_filename + "_sc.csv", "s -> c");
		muon_pair_production_flavor_to_flavor(x_bins, y_bins, E_beam_bins, Flavor::Bottom, Flavor::Up, base_filename + "_bu.csv", "b -> u");
		muon_pair_production_flavor_to_flavor(x_bins, y_bins, E_beam_bins, Flavor::Bottom, Flavor::Charm, base_filename + "_bc.csv", "b -> c");

		muon_pair_production_flavor_to_flavor(x_bins, y_bins, E_beam_bins, Flavor::AntiUp, Flavor::AntiDown, base_filename + "_ubar_dbar.csv", "ubar -> dbar");
		muon_pair_production_flavor_to_flavor(x_bins, y_bins, E_beam_bins, Flavor::AntiCharm, Flavor::AntiDown, base_filename + "_cbar_dbar.csv", "cbar -> dbar");
		muon_pair_production_flavor_to_flavor(x_bins, y_bins, E_beam_bins, Flavor::AntiUp, Flavor::AntiStrange, base_filename + "_ubar_sbar.csv", "ubar -> sbar");
		muon_pair_production_flavor_to_flavor(x_bins, y_bins, E_beam_bins, Flavor::AntiCharm, Flavor::AntiStrange, base_filename + "_cbar_sbar.csv", "cbar -> sbar");
		muon_pair_production_flavor_to_flavor(x_bins, y_bins, E_beam_bins, Flavor::AntiUp, Flavor::AntiBottom, base_filename + "_ubar_bbar.csv", "ubar -> bbar");
		muon_pair_production_flavor_to_flavor(x_bins, y_bins, E_beam_bins, Flavor::AntiCharm, Flavor::AntiBottom, base_filename + "_cbar_bbar.csv", "cbar -> bbar");
	}
	void muon_pair_production_flavor_decomposition_gluon_to_quark(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins,
		const std::string base_filename) {

		muon_pair_production_flavor_to_flavor(x_bins, y_bins, E_beam_bins, Flavor::Gluon, Flavor::Up, base_filename + "_gu.csv", "g -> u");
		muon_pair_production_flavor_to_flavor(x_bins, y_bins, E_beam_bins, Flavor::Gluon, Flavor::Charm, base_filename + "_gc.csv", "g -> c");

		muon_pair_production_flavor_to_flavor(x_bins, y_bins, E_beam_bins, Flavor::Gluon, Flavor::AntiDown, base_filename + "_g_dbar.csv", "g -> dbar");
		muon_pair_production_flavor_to_flavor(x_bins, y_bins, E_beam_bins, Flavor::Gluon, Flavor::AntiStrange, base_filename + "_g_sbar.csv", "g -> sbar");
		muon_pair_production_flavor_to_flavor(x_bins, y_bins, E_beam_bins, Flavor::Gluon, Flavor::AntiBottom, base_filename + "_g_bbar.csv", "g -> bbar");
	}

	/*
		FRAGMENTATION DECOMPOSITION
	*/

	template <is_pdf_interface FFInterface, is_decay_function DecayFunction>
	void muon_pair_production_fragmentation_channel(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins,
		const FFInterface &ff,
		const Decay<DecayFunction> decay,
		const std::filesystem::path filename, 
		const std::string comment = "") {
		
		muon_pair_production(
			x_bins, y_bins, E_beam_bins, filename,
			LHAInterface(params.pdf_set),
			FragmentationConfiguration({ ff }, { decay }),
			comment
		);
	}

	void muon_pair_production_fragmentation_decomposition(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins,
		const std::string base_filename) {

		const double minimum_lepton_momentum = params.minimum_lepton_momentum;
		const Particle target = params.process.target;
		const auto decay_function = DecayFunctions::decay_function;

		const DecayParametrization parametrization = DecayParametrization::fit1();

		muon_pair_production_fragmentation_channel(x_bins, y_bins, E_beam_bins,
			LHAInterface("kkks08_opal_d0___mas", {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}),
			Decay(parametrization, Constants::Particles::D0, target, decay_function, minimum_lepton_momentum),
			base_filename + "_d0.csv", "d0 only"
		);
		muon_pair_production_fragmentation_channel(x_bins, y_bins, E_beam_bins,
			LHAInterface("kkks08_opal_d+___mas"),
			Decay(parametrization, Constants::Particles::Dp, target, decay_function, minimum_lepton_momentum),
			base_filename + "_d+.csv", "d+ only"
		);
		muon_pair_production_fragmentation_channel(x_bins, y_bins, E_beam_bins,
			LHAInterface("bkk05_D3_d_s_nlo"),
			Decay(parametrization, Constants::Particles::Ds, target, decay_function, minimum_lepton_momentum),
			base_filename + "_d_s.csv", "d_s only"
		);
		muon_pair_production_fragmentation_channel(x_bins, y_bins, E_beam_bins,
			LHAInterface("bkk05_D3_lambda_c_nlo"),
			Decay(parametrization, Constants::Particles::LambdaC, target, decay_function, minimum_lepton_momentum),
			base_filename + "_lambda_c.csv", "lambda_c only"
		);
	}

	void muon_pair_production_only_D0(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path filename, 
		const std::string comment = "") {

		const double minimum_lepton_momentum = params.minimum_lepton_momentum;
		const Particle target = params.process.target;
		const auto decay_function = DecayFunctions::decay_function;

		const DecayParametrization parametrization = DecayParametrization::fit1();

		SIDIS sidis(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			LHAInterface(params.pdf_set),
			FragmentationConfiguration(
				{
					2.3065 * LHAInterface("kkks08_opal_d0___mas", {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
				},
				{
					Decay(parametrization, Constants::Particles::D0, target, decay_function, minimum_lepton_momentum),
				}
			),
			params.process,
			renormalization, factorization, fragmentation
		);

		sidis.charm_mass = params.charm_mass;

		sidis.combine_integrals = true;
		sidis.use_modified_cross_section_prefactor = true;

		sidis.parallelize = params.parallelize;
		sidis.number_of_threads = params.number_of_threads;

		sidis.freeze_factorization_scale = params.freeze_factorization;
		sidis.freeze_fragmentation_scale = params.freeze_fragmentation;

		sidis.lepton_pair_cross_section_xy(x_bins, y_bins, E_beam_bins, filename, comment);
	}
};

#endif