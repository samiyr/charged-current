#ifndef ANALYSIS_SIDIS_H
#define ANALYSIS_SIDIS_H

#include <vector>
#include "Analysis/Parameters.cpp"
#include "Analysis/Constants.cpp"
#include "Common/Flavor.cpp"
#include "Common/Constants.cpp"
#include "Common/Process.cpp"
#include "PDF/Interfaces/LHAInterface.cpp"
#include "SIDIS/SIDIS.cpp"

struct SIDISAnalysis {
	const AnalysisParameters params;

	SIDISAnalysis(const AnalysisParameters params) : params(params) { }

	template <PDFConcept PDFInterface, PDFConcept FFInterface, DecayFunctions::Concept DecayFunction>
	void muon_pair_production(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path filename, 
		const PDFInterface &pdf,
		const FragmentationConfiguration<FFInterface, DecayFunction> &ff,
		const std::string comment = "") {
		SIDIS sidis(
			{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
			pdf, ff,
			params.process				
		);

		sidis.charm_mass = params.charm_mass;

		sidis.combine_integrals = true;
		sidis.use_modified_cross_section_prefactor = true;

		sidis.lepton_pair_cross_section_xy(x_bins, y_bins, E_beam_bins, filename, comment);
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

		const DecayParametrization parametrization = DecayParametrization::fit1();

		muon_pair_production(
			x_bins, y_bins, E_beam_bins, filename,
			LHAInterface(params.pdf_set),
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
	}

	void muon_pair_production(const AnalysisSet set, const std::vector<double> x_bins, const std::filesystem::path filename, const std::string comment = "") {
		switch (set) {
		case AnalysisSet::NuTeV:
			return muon_pair_production(x_bins, AnalysisConstants::NuTeV::New::y_bins, AnalysisConstants::NuTeV::New::E_bins, filename, comment);
		case AnalysisSet::NuTeV_old:
			return muon_pair_production(x_bins, AnalysisConstants::NuTeV::Old::y_bins, AnalysisConstants::NuTeV::Old::E_bins, filename, comment);
		case AnalysisSet::CCFR:
			return muon_pair_production(x_bins, AnalysisConstants::CCFR::y_bins, AnalysisConstants::CCFR::E_bins, filename, comment);
		}
	}

	void muon_pair_production_only_D0(const AnalysisSet set, const std::vector<double> x_bins, const std::filesystem::path filename, const std::string comment = "") {
		switch (set) {
		case AnalysisSet::NuTeV:
			return muon_pair_production_only_D0(x_bins, AnalysisConstants::NuTeV::New::y_bins, AnalysisConstants::NuTeV::New::E_bins, filename, comment);
		case AnalysisSet::NuTeV_old:
			return muon_pair_production_only_D0(x_bins, AnalysisConstants::NuTeV::Old::y_bins, AnalysisConstants::NuTeV::Old::E_bins, filename, comment);
		case AnalysisSet::CCFR:
			return muon_pair_production_only_D0(x_bins, AnalysisConstants::CCFR::y_bins, AnalysisConstants::CCFR::E_bins, filename, comment);
		}
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

	template <PDFConcept FFInterface, DecayFunctions::Concept DecayFunction>
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
					LHAInterface("kkks08_opal_d0___mas", {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
				},
				{
					Decay(parametrization, Constants::Particles::D0, target, decay_function, minimum_lepton_momentum),
				}
			),
			params.process
		);

		sidis.charm_mass = params.charm_mass;

		sidis.combine_integrals = true;
		sidis.use_modified_cross_section_prefactor = true;

		sidis.lepton_pair_cross_section_xy(x_bins, y_bins, E_beam_bins, filename, comment);
	}
};

#endif