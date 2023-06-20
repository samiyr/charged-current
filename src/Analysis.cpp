#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "DIS.cpp"
#include "SIDIS.cpp"
#include "LHAInterface.cpp"
#include "Particle.cpp"

namespace Analysis {
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
				Process {Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino}
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
			Analysis::Inclusive::charm_production(x_bins, {0.334, 0.573, 0.790}, {90.18, 174.37, 244.72}, filename, comment);
		}
		static void charm_production_nutev_new(const std::vector<double> x_bins, const std::string filename, const std::string comment = "") {
			Analysis::Inclusive::charm_production(x_bins, {0.324, 0.558, 0.771}, {88.29, 174.29, 247.0}, filename, comment);
		}
		static void charm_production_ccfr(const std::vector<double> x_bins, const std::string filename, const std::string comment = "") {
			Analysis::Inclusive::charm_production(x_bins, {0.32, 0.57, 0.795}, {109.46, 209.89, 332.7}, filename, comment);
		}
	}

	namespace SemiInclusive {
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

			SIDIS sidis(
				{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom},
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
				300'000,
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