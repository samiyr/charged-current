#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "DIS.cpp"
#include "SIDIS.cpp"
#include "LHAInterface.cpp"
#include "Particle.cpp"

namespace Analysis {
	namespace NuTeV {
		namespace Old {
			const static std::vector<double> y_bins = {0.334, 0.573, 0.790};
			const static std::vector<double> E_bins = {90.18, 174.37, 244.72};
		}
		namespace New {
			const static std::vector<double> y_bins = {0.324, 0.558, 0.771};
			const static std::vector<double> E_bins = {88.29, 174.29, 247.0};
		}
	}
	namespace CCFR {
		const static std::vector<double> y_bins = {0.32, 0.57, 0.795};
		const static std::vector<double> E_bins = {109.46, 209.89, 332.7};
	}
	namespace Inclusive {
		static void charm_production(
			const std::vector<double> x_bins, 
			const std::vector<double> y_bins, 
			const std::vector<double> E_beam_bins, 
			const std::string filename, 
			const std::string comment = "") {
			DIS dis(
				{Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
				LHAInterface("EPPS21nlo_CT18Anlo_Fe56"),
				20'000,
				Process(Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino)
			);

			dis.max_chi_squared_deviation = 0.2;
			dis.max_relative_error = 1e-3;
			dis.iter_max = 10;
			dis.compute_differential_cross_section_directly = false;
			dis.use_modified_cross_section_prefactor = true;

			dis.charm_mass = 1.3;

			dis.differential_cross_section_xy(x_bins, y_bins, E_beam_bins, filename, comment);
		}

		static void charm_production_nutev_old(const std::vector<double> x_bins, const std::string filename, const std::string comment = "") {
			Analysis::Inclusive::charm_production(x_bins, NuTeV::Old::y_bins, NuTeV::Old::E_bins, filename, comment);
		}
		static void charm_production_nutev_new(const std::vector<double> x_bins, const std::string filename, const std::string comment = "") {
			Analysis::Inclusive::charm_production(x_bins, NuTeV::New::y_bins, NuTeV::New::E_bins, filename, comment);
		}
		static void charm_production_ccfr(const std::vector<double> x_bins, const std::string filename, const std::string comment = "") {
			Analysis::Inclusive::charm_production(x_bins, CCFR::y_bins, CCFR::E_bins, filename, comment);
		}
	}

	namespace SemiInclusive {
		template <PDFConcept PDFInterface, PDFConcept FFInterface, DecayFunctions::Concept DecayFunction>
		static void muon_pair_production(
			const std::vector<double> x_bins, 
			const std::vector<double> y_bins, 
			const std::vector<double> E_beam_bins, 
			const std::string filename, 
			const PDFInterface &pdf,
			const FragmentationConfiguration<FFInterface, DecayFunction> &ff,
			const std::string comment = "") {
			SIDIS sidis(
				{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
				pdf, ff,
				200'000,
				Process {Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino }				
			);

			sidis.charm_mass = 1.3;

			sidis.max_chi_squared_deviation = 0.5;
			sidis.max_relative_error = 1e-2;
			sidis.iter_max = 5;
			sidis.combine_integrals = true;
			sidis.use_modified_cross_section_prefactor = true;

			sidis.lepton_pair_cross_section_xy(x_bins, y_bins, E_beam_bins, filename, comment);
		}
		static void muon_pair_production(
			const std::vector<double> x_bins, 
			const std::vector<double> y_bins, 
			const std::vector<double> E_beam_bins, 
			const std::string filename, 
			const std::string comment = "") {
			
			const double minimum_lepton_momentum = 5.0;
			const Particle target = Constants::Particles::Proton;
			const auto decay_function = DecayFunctions::decay_function;

			const DecayParametrization parametrization = DecayParametrization::fit1();

			muon_pair_production(
				x_bins, y_bins, E_beam_bins, filename,
				LHAInterface("EPPS21nlo_CT18Anlo_Fe56"),
				FragmentationConfiguration(
					{
						LHAInterface("kkks08_opal_d0___mas", {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
						LHAInterface("kkks08_opal_d+___mas"), 
						LHAInterface("bkk05_D3_d_s_nlo"), 
						LHAInterface("bkk05_D3_lambda_c_nlo")
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

		namespace ChannelDecomposition {
			static void muon_pair_production_separated_channels(
				const std::vector<double> x_bins, 
				const std::vector<double> y_bins, 
				const std::vector<double> E_beam_bins, 
				const std::string filename,
				const std::vector<bool> incoming,
				const std::vector<bool> outgoing,
				const std::string comment = "") {
				
				const double minimum_lepton_momentum = 5.0;
				const Particle target = Constants::Particles::Proton;
				const auto decay_function = DecayFunctions::decay_function;

				const DecayParametrization parametrization = DecayParametrization::fit1();

				const std::vector<double> incoming_multipliers(incoming.begin(), incoming.end());
				const std::vector<double> outgoing_multipliers(outgoing.begin(), outgoing.end());

				std::vector<double> d0_multipliers(outgoing.begin(), outgoing.end());
				Utility::multiply(d0_multipliers, {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0});

				muon_pair_production(
					x_bins, y_bins, E_beam_bins, filename,
					LHAInterface("EPPS21nlo_CT18Anlo_Fe56", incoming_multipliers),
					FragmentationConfiguration(
						{
							LHAInterface("kkks08_opal_d0___mas", d0_multipliers), 
							LHAInterface("kkks08_opal_d+___mas", outgoing_multipliers), 
							LHAInterface("bkk05_D3_d_s_nlo", outgoing_multipliers), 
							LHAInterface("bkk05_D3_lambda_c_nlo", outgoing_multipliers)
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
			static void muon_pair_production_quark_to_quark(
				const std::vector<double> x_bins, 
				const std::vector<double> y_bins, 
				const std::vector<double> E_beam_bins, 
				const std::string filename,
				const std::string comment = "") {

				muon_pair_production_separated_channels(x_bins, y_bins, E_beam_bins, filename,
				{
					true, true, true, true, true, true, false, true, true, true, true, true, true
				},
				{
					true, true, true, true, true, true, false, true, true, true, true, true, true
				}, comment);
			}
			static void muon_pair_production_quark_to_gluon(
				const std::vector<double> x_bins, 
				const std::vector<double> y_bins, 
				const std::vector<double> E_beam_bins, 
				const std::string filename,
				const std::string comment = "") {

				muon_pair_production_separated_channels(x_bins, y_bins, E_beam_bins, filename,
				{
					true, true, true, true, true, true, false, true, true, true, true, true, true
				},
				{
					false, false, false, false, false, false, true, false, false, false, false, false, false
				}, comment);
			}
			static void muon_pair_production_gluon_to_quark(
				const std::vector<double> x_bins, 
				const std::vector<double> y_bins, 
				const std::vector<double> E_beam_bins, 
				const std::string filename,
				const std::string comment = "") {

				muon_pair_production_separated_channels(x_bins, y_bins, E_beam_bins, filename,
				{
					false, false, false, false, false, false, true, false, false, false, false, false, false
				},
				{
					true, true, true, true, true, true, false, true, true, true, true, true, true
				}, comment);
			}
			static void muon_pair_production_gluon_to_gluon(
				const std::vector<double> x_bins, 
				const std::vector<double> y_bins, 
				const std::vector<double> E_beam_bins, 
				const std::string filename,
				const std::string comment = "") {

				muon_pair_production_separated_channels(x_bins, y_bins, E_beam_bins, filename,
				{
					false, false, false, false, false, false, true, false, false, false, false, false, false
				},
				{
					false, false, false, false, false, false, true, false, false, false, false, false, false
				}, comment);
			}
			static void muon_pair_production_quark_gluon_channels(
				const std::vector<double> x_bins, 
				const std::vector<double> y_bins, 
				const std::vector<double> E_beam_bins, 
				const std::vector<std::string> filenames) {
				
				muon_pair_production_quark_to_quark(x_bins, y_bins, E_beam_bins, filenames[0], "quark to quark");
				muon_pair_production_quark_to_gluon(x_bins, y_bins, E_beam_bins, filenames[1], "quark to gluon");
				muon_pair_production_gluon_to_quark(x_bins, y_bins, E_beam_bins, filenames[2], "gluon to quark");
				muon_pair_production_gluon_to_gluon(x_bins, y_bins, E_beam_bins, filenames[3], "gluon to gluon");
			}
		}
		
		namespace FlavorDecomposition {
			static void muon_pair_production_flavor_to_flavor(
				const std::vector<double> x_bins, 
				const std::vector<double> y_bins, 
				const std::vector<double> E_beam_bins,
				const FlavorType incoming,
				const FlavorType outgoing,
				const std::string filename,
				const std::string comment = "") {

				std::vector<bool> incoming_vector(TOTAL_FLAVORS, false);
				std::vector<bool> outgoing_vector(TOTAL_FLAVORS, false);

				incoming_vector[size_t(incoming + 6)] = true;
				outgoing_vector[size_t(outgoing + 6)] = true;

				Analysis::SemiInclusive::ChannelDecomposition::muon_pair_production_separated_channels(x_bins, y_bins, E_beam_bins, filename, incoming_vector, outgoing_vector, comment);
			}

			static void muon_pair_production_quark_to_quark(
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
			static void muon_pair_production_gluon_to_quark(
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
		}

		namespace FragmentationDecomposition {
			template <PDFConcept FFInterface, DecayFunctions::Concept DecayFunction>
			static void muon_pair_production_fragmentation_channel(
				const std::vector<double> x_bins, 
				const std::vector<double> y_bins, 
				const std::vector<double> E_beam_bins,
				const FFInterface &ff,
				const Decay<DecayFunction> decay,
				const std::string filename, 
				const std::string comment = "") {
				
				Analysis::SemiInclusive::muon_pair_production(
					x_bins, y_bins, E_beam_bins, filename,
					LHAInterface("EPPS21nlo_CT18Anlo_Fe56"),
					FragmentationConfiguration({ ff }, { decay }),
					comment
				);
			}

			static void muon_pair_production(
				const std::vector<double> x_bins, 
				const std::vector<double> y_bins, 
				const std::vector<double> E_beam_bins,
				const std::string base_filename) {

				const double minimum_lepton_momentum = 5.0;
				const Particle target = Constants::Particles::Proton;
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
		}

		static void muon_pair_production_only_D0(
			const std::vector<double> x_bins, 
			const std::vector<double> y_bins, 
			const std::vector<double> E_beam_bins, 
			const std::string filename, 
			const std::string comment = "") {
			const double minimum_lepton_momentum = 5.0;
			const Particle target = Constants::Particles::Proton;
			const auto decay_function = DecayFunctions::decay_function;

			const DecayParametrization parametrization = DecayParametrization::fit1();

			SIDIS sidis(
				{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
				LHAInterface("EPPS21nlo_CT18Anlo_Fe56"),
				FragmentationConfiguration(
					{
						LHAInterface("kkks08_opal_d0___mas", {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
					},
					{
						Decay(parametrization, Constants::Particles::D0, target, decay_function, minimum_lepton_momentum),
					}
				),
				300'000,
				Process(Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino)
			);

			sidis.charm_mass = 1.3;

			sidis.max_chi_squared_deviation = 0.5;
			sidis.max_relative_error = 1e-2;
			sidis.iter_max = 5;
			sidis.combine_integrals = true;
			sidis.use_modified_cross_section_prefactor = true;

			sidis.lepton_pair_cross_section_xy(x_bins, y_bins, E_beam_bins, filename, comment);
		}
		static void muon_pair_production_nutev_old(const std::vector<double> x_bins, const std::string filename, const std::string comment = "") {
			Analysis::SemiInclusive::muon_pair_production(x_bins, {0.334, 0.573, 0.790}, {90.18, 174.37, 244.72}, filename, comment);
		}
		static void muon_pair_production_nutev_new(const std::vector<double> x_bins, const std::string filename, const std::string comment = "") {
			Analysis::SemiInclusive::muon_pair_production(x_bins, {0.324, 0.558, 0.771}, {88.29, 174.29, 247.0}, filename, comment);
		}
		static void muon_pair_production_ccfr(const std::vector<double> x_bins, const std::string filename, const std::string comment = "") {
			Analysis::SemiInclusive::muon_pair_production(x_bins, {0.32, 0.57, 0.795}, {109.46, 209.89, 332.7}, filename, comment);
		}

		static void muon_pair_production_nutev_old_only_D0(const std::vector<double> x_bins, const std::string filename, const std::string comment = "") {
			Analysis::SemiInclusive::muon_pair_production_only_D0(x_bins, {0.334, 0.573, 0.790}, {90.18, 174.37, 244.72}, filename, comment);
		}
		static void muon_pair_production_nutev_new_only_D0(const std::vector<double> x_bins, const std::string filename, const std::string comment = "") {
			Analysis::SemiInclusive::muon_pair_production_only_D0(x_bins, {0.324, 0.558, 0.771}, {88.29, 174.29, 247.0}, filename, comment);
		}
		static void muon_pair_production_ccfr_only_D0(const std::vector<double> x_bins, const std::string filename, const std::string comment = "") {
			Analysis::SemiInclusive::muon_pair_production_only_D0(x_bins, {0.32, 0.57, 0.795}, {109.46, 209.89, 332.7}, filename, comment);
		}
	}
}


#endif