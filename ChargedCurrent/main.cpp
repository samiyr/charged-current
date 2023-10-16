#ifndef CHARGED_CURRENT_DIS_H
#define CHARGED_CURRENT_DIS_H

#include <iostream>
#include <array>
#include <format>
#include <chrono>

#include <ChargedCurrent/DIS/DIS.cpp>
#include <ChargedCurrent/SIDIS/SIDIS.cpp>
#include <ChargedCurrent/PDF/Interfaces/LHAInterface.cpp>
#include <ChargedCurrent/Common/TRFKinematics.cpp>
#include <ChargedCurrent/Analysis/Analysis.cpp>
#include <ChargedCurrent/Decay/DecayParametrization.cpp>

int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv) {
	LHAInterface<>::disable_verbosity();

	const std::string separator = "========================================";

	std::vector<std::string> arguments(argv + 1, argv + argc);

	auto nthreads_iterator = std::find(arguments.begin(), arguments.end(), "--nthreads");
	unsigned int number_of_threads = 1;

	if (nthreads_iterator == arguments.end()) {
		std::cout << "No --nthreads argument given, defaulting to 1" << IO::endl;
	} else {
		nthreads_iterator++;
		number_of_threads = static_cast<unsigned int>(std::stoul(*nthreads_iterator));
	}

	const bool all = arguments.size() == 0 || arguments.size() == 2 || Collections::contains(arguments, "all");
	const auto run = [&](const std::string arg) { return Collections::contains(arguments, arg) || all; };

	const std::vector<double> x_bins = {
		0.02, 0.02125, 0.0225, 0.02375, 0.025, 0.02625, 0.0275, 0.02875, 
		0.03, 0.03125, 0.0325, 0.03375, 0.035, 0.03625, 0.0375, 0.03875, 
		0.04, 0.04125, 0.0425, 0.04375, 0.045, 0.04625, 0.0475, 0.04875,
		0.05, 0.0625, 0.075, 0.0875, 
		0.1, 0.125, 0.15, 0.175, 
		0.2, 0.225, 0.25, 0.275, 
		0.3, 0.325, 0.35, 0.375, 
		0.4
	};
	const std::vector<double> x_bins_all = {
		1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6, 8e-6, 9e-6,
		1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5,
		1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4,
		1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3,
		0.01, 0.01125, 0.0125, 0.01375, 0.015, 0.01625, 0.0175, 0.01875, 
		0.02, 0.02125, 0.0225, 0.02375, 0.025, 0.02625, 0.0275, 0.02875, 
		0.03, 0.03125, 0.0325, 0.03375, 0.035, 0.03625, 0.0375, 0.03875, 
		0.04, 0.04125, 0.0425, 0.04375, 0.045, 0.04625, 0.0475, 0.04875,
		0.05, 0.0625, 0.075, 0.0875, 
		0.1, 0.125, 0.15, 0.175, 
		0.2, 0.225, 0.25, 0.275, 
		0.3, 0.325, 0.35, 0.375, 
		0.4, 0.425, 0.45, 0.475, 
		0.5, 0.525, 0.55, 0.575, 
		0.6, 0.625, 0.65, 0.675, 
		0.7, 0.725, 0.75, 0.775, 
		0.8, 0.825, 0.85, 0.875, 
		0.9, 0.925, 0.95, 0.975, 
	};

	const auto bar = []<
		is_scale_dependence RenormalizationScale, is_scale_dependence FactorizationScale, is_scale_dependence FragmentationScale
	>(Analysis<RenormalizationScale, FactorizationScale, FragmentationScale> analysis) {
		analysis.params.process.type = Process::Type::AntiNeutrinoToAntiLepton;
		return analysis;
	};

	const auto measure = [](auto &&func, auto &&...params) {
		const auto &start = std::chrono::high_resolution_clock::now();

		std::forward<decltype(func)>(func)(std::forward<decltype(params)>(params)...);
		
		const auto &stop = std::chrono::high_resolution_clock::now();
		
		const std::chrono::duration<double> elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
		
		std::cout << "Took " << elapsed.count() << " seconds" << IO::endl;
    };

	const std::vector<std::string> pdfs = {
		"EPPS21nlo_CT18Anlo_Fe56",
		"nCTEQ15HQ_FullNuc_56_26",
		"nNNPDF30_nlo_as_0118_A56_Z26"
	};
	const std::vector<std::string> free_pdfs = {
		"CT18ANLO",
		"nCTEQ15HQ_FullNuc_1_1",
		"nNNPDF30_nlo_as_0118_p"
	};

	Analysis base;
	base.params.charm_mass = 1.3;
	base.params.scale_variation = ScaleVariation::None;
	base.params.integration.cuba.maximum_relative_error = 1e-2;
	base.params.order = PerturbativeOrder::NLO;
	base.params.process.type = Process::Type::NeutrinoToLepton;
	base.params.number_of_threads = number_of_threads;

	Analysis nlo = base;

	Analysis nlo_nlp = base;
	nlo_nlp.params.use_nlp_nlo = true;

	Analysis nnlo = base;
	nnlo.params.order = PerturbativeOrder::NNLO;

	Analysis scale = base;
	scale.params.scale_variation = ScaleVariation::All;
	scale.params.freeze_factorization = true;

	Analysis epps_errors = base;
	epps_errors.params.pdf_error_sets = true;

	Analysis nomad_errors = epps_errors;
	nomad_errors.params.Q2_min = 1.69;
	nomad_errors.params.minimum_lepton_momentum = 3.0;

	Analysis unity_nuclear_epps = base;
	unity_nuclear_epps.params.pdf_error_sets = true;
	unity_nuclear_epps.params.pdf_set = "CT18ANLO";
	unity_nuclear_epps.params.explicit_isospin = true;
	unity_nuclear_epps.params.Z = 26.0;
	unity_nuclear_epps.params.A = 56.0;

	Analysis unity_nuclear_ncteq = unity_nuclear_epps;
	unity_nuclear_ncteq.params.pdf_error_sets = false;
	unity_nuclear_ncteq.params.pdf_set = "nCTEQ15HQ_FullNuc_1_1";

	Analysis unity_nuclear_nnnpdf = unity_nuclear_epps;
	unity_nuclear_nnnpdf.params.pdf_set = "nNNPDF30_nlo_as_0118_p";

	if (run("charmproduction")) {
		std::cout << "========= DIS charm production =========" << IO::endl;
		
		measure([&] {
			epps_errors.dis().charm_production(AnalysisSet::NuTeV_old, x_bins, "Data/DIS/CharmProduction/nutev_old_neutrino.csv");
			epps_errors.dis().charm_production(AnalysisSet::NuTeV, x_bins, "Data/DIS/CharmProduction/nutev_new_neutrino.csv");
			epps_errors.dis().charm_production(AnalysisSet::CCFR, x_bins, "Data/DIS/CharmProduction/ccfr_neutrino.csv");

			bar(epps_errors).dis().charm_production(AnalysisSet::NuTeV_old, x_bins, "Data/DIS/CharmProduction/nutev_old_antineutrino.csv");
			bar(epps_errors).dis().charm_production(AnalysisSet::NuTeV, x_bins, "Data/DIS/CharmProduction/nutev_new_antineutrino.csv");
			bar(epps_errors).dis().charm_production(AnalysisSet::CCFR, x_bins, "Data/DIS/CharmProduction/ccfr_antineutrino.csv");
		});

		std::cout << separator << IO::endl;
	}

	if (run("disintegrated")) {
		std::cout << "======= DIS integrated inclusive =======" << IO::endl;
		
		measure([&] {
			nomad_errors.dis().integrated(AnalysisSet::NOMAD, "Data/DIS/Integrated/nomad_neutrino.csv");
		});

		std::cout << separator << IO::endl;
	}
	
	if (run("sidisintegrated")) {
		std::cout << "====== SIDIS integrated inclusive ======" << IO::endl;
		
		measure([&] {
			nomad_errors.sidis().integrated_muon_pair_production(AnalysisSet::NOMAD, "Data/SIDIS/MuonPairProduction/CharmedHadrons/Integrated/nomad_neutrino.csv");
		});

		std::cout << separator << IO::endl;
	}

	if (run("nnlo")) {
		std::cout << "============== SIDIS NNLO ==============" << IO::endl;

		measure([&] {
			nnlo.sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/nutev_old_neutrino.csv");
			nnlo.sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/nutev_new_neutrino.csv");
			nnlo.sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ccfr_neutrino.csv");

			bar(nnlo).sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/nutev_old_antineutrino.csv");
			bar(nnlo).sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/nutev_new_antineutrino.csv");
			bar(nnlo).sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ccfr_antineutrino.csv");
		});

		std::cout << separator << IO::endl;
	}

	if (run("nlonlp")) {
		std::cout << "============ SIDIS NLO NLP =============" << IO::endl;

		measure([&] {
			nlo_nlp.sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/NLO_NLP/nutev_old_neutrino.csv");
			nlo_nlp.sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/NLO_NLP/nutev_new_neutrino.csv");
			nlo_nlp.sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/NLO_NLP/ccfr_neutrino.csv");

			bar(nlo_nlp).sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/NLO_NLP/nutev_old_antineutrino.csv");
			bar(nlo_nlp).sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/NLO_NLP/nutev_new_antineutrino.csv");
			bar(nlo_nlp).sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/NLO_NLP/ccfr_antineutrino.csv");
		});

		std::cout << separator << IO::endl;
	}

	if (run("unitynuclear")) {
		std::cout << "========= SIDIS unity nuclear ==========" << IO::endl;

		measure([&] {
			unity_nuclear_epps.sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/CT18ANLO/nutev_old_neutrino.csv");
			unity_nuclear_epps.sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/CT18ANLO/nutev_new_neutrino.csv");
			unity_nuclear_epps.sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/CT18ANLO/ccfr_neutrino.csv");

			bar(unity_nuclear_epps).sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/CT18ANLO/nutev_old_antineutrino.csv");
			bar(unity_nuclear_epps).sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/CT18ANLO/nutev_new_antineutrino.csv");
			bar(unity_nuclear_epps).sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/CT18ANLO/ccfr_antineutrino.csv");
		});

		measure([&] {
			unity_nuclear_ncteq.sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/nCTEQ15HQ_FullNuc_1_1/nutev_old_neutrino.csv");
			unity_nuclear_ncteq.sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/nCTEQ15HQ_FullNuc_1_1/nutev_new_neutrino.csv");
			unity_nuclear_ncteq.sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/nCTEQ15HQ_FullNuc_1_1/ccfr_neutrino.csv");

			bar(unity_nuclear_ncteq).sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/nCTEQ15HQ_FullNuc_1_1/nutev_old_antineutrino.csv");
			bar(unity_nuclear_ncteq).sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/nCTEQ15HQ_FullNuc_1_1/nutev_new_antineutrino.csv");
			bar(unity_nuclear_ncteq).sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/nCTEQ15HQ_FullNuc_1_1/ccfr_antineutrino.csv");
		});

		measure([&] {
			unity_nuclear_nnnpdf.sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/nNNPDF30_nlo_as_0118_p/nutev_old_neutrino.csv");
			unity_nuclear_nnnpdf.sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/nNNPDF30_nlo_as_0118_p/nutev_new_neutrino.csv");
			unity_nuclear_nnnpdf.sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/nNNPDF30_nlo_as_0118_p/ccfr_neutrino.csv");

			bar(unity_nuclear_nnnpdf).sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/nNNPDF30_nlo_as_0118_p/nutev_old_antineutrino.csv");
			bar(unity_nuclear_nnnpdf).sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/nNNPDF30_nlo_as_0118_p/nutev_new_antineutrino.csv");
			bar(unity_nuclear_nnnpdf).sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/UnityNuclear/nNNPDF30_nlo_as_0118_p/ccfr_antineutrino.csv");
		});

		std::cout << separator << IO::endl;
	}

	if (run("nnlod0")) {
		std::cout << "========== SIDIS NNLO D0 only ==========" << IO::endl;

		measure([&] {
			nnlo.sidis().muon_pair_production_only_D0(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/D0/nutev_old_neutrino.csv");
			nnlo.sidis().muon_pair_production_only_D0(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/D0/nutev_new_neutrino.csv");
			nnlo.sidis().muon_pair_production_only_D0(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/D0/ccfr_neutrino.csv");

			bar(nnlo).sidis().muon_pair_production_only_D0(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/D0/nutev_old_antineutrino.csv");
			bar(nnlo).sidis().muon_pair_production_only_D0(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/D0/nutev_new_antineutrino.csv");
			bar(nnlo).sidis().muon_pair_production_only_D0(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/D0/ccfr_antineutrino.csv");
		});

		std::cout << separator << IO::endl;
	}

	if (run("scale")) {
		std::cout << "======== SIDIS scale variations ========" << IO::endl;

		measure([&] {
			scale.sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ScaleDependence/nutev_old_neutrino.csv");
			scale.sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ScaleDependence/nutev_new_neutrino.csv");
			scale.sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ScaleDependence/ccfr_neutrino.csv");

			bar(scale).sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ScaleDependence/nutev_old_antineutrino.csv");
			bar(scale).sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ScaleDependence/nutev_new_antineutrino.csv");
			bar(scale).sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/ScaleDependence/ccfr_antineutrino.csv");
		});

		std::cout << separator << IO::endl;
	}

	if (run("errors")) {
		std::cout << "=========== SIDIS error sets ===========" << IO::endl;
		for (const auto &pdf_set_name : pdfs) {
			std::cout << "PDF set: " << pdf_set_name << IO::endl;

			Analysis errors = base;
			errors.params.pdf_error_sets = true;
			errors.params.pdf_set = pdf_set_name;

			const std::string output_folder = "Data/SIDIS/MuonPairProduction/CharmedHadrons/ErrorSets/" + pdf_set_name + "/";

			measure([&] {
				errors.sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, output_folder + "nutev_old_neutrino.csv");
				errors.sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, output_folder + "nutev_new_neutrino.csv");
				errors.sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, output_folder + "ccfr_neutrino.csv");

				bar(errors).sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, output_folder + "nutev_old_antineutrino.csv");
				bar(errors).sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, output_folder + "nutev_new_antineutrino.csv");
				bar(errors).sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, output_folder + "ccfr_antineutrino.csv");
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("channels")) {
		std::cout << "======== SIDIS parton channels =========" << IO::endl;

		measure([&] {
			nlo.sidis().muon_pair_production_quark_gluon_channels(x_bins, AnalysisConstants::NuTeV::New::Neutrino::y_bins, AnalysisConstants::NuTeV::New::Neutrino::E_bins, 
				{
					"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_quark_to_quark.csv",
					"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_quark_to_gluon.csv",
					"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_gluon_to_quark.csv",
					"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_gluon_to_gluon.csv"
				}
			);
		});

		std::cout << separator << IO::endl;
	}

	if (run("flavors")) {
		std::cout << "======== SIDIS flavor channels =========" << IO::endl;

		measure([&] {
			nlo.sidis().muon_pair_production_flavor_decomposition_quark_to_quark(
				x_bins, AnalysisConstants::NuTeV::New::Neutrino::y_bins, AnalysisConstants::NuTeV::New::Neutrino::E_bins,
				"Data/SIDIS/MuonPairProduction/CharmedHadrons/FlavorDecomposition/nutev_new"
			);
			nlo.sidis().muon_pair_production_flavor_decomposition_gluon_to_quark(
				x_bins, AnalysisConstants::NuTeV::New::Neutrino::y_bins, AnalysisConstants::NuTeV::New::Neutrino::E_bins,
				"Data/SIDIS/MuonPairProduction/CharmedHadrons/FlavorDecomposition/nutev_new"
			);
		});

		std::cout << separator << IO::endl;
	}

	if (run("fragmentation")) {
		std::cout << "===== SIDIS fragmentation channels =====" << IO::endl;

		measure([&] {
			nlo.sidis().muon_pair_production_fragmentation_decomposition(
				x_bins, AnalysisConstants::NuTeV::New::Neutrino::y_bins, AnalysisConstants::NuTeV::New::Neutrino::E_bins,
				"Data/SIDIS/MuonPairProduction/CharmedHadrons/FragmentationDecomposition/nutev_new"
			);
		});

		std::cout << separator << IO::endl;
	}

	Analysis proton_nlo;
	proton_nlo.params.pdf_set = "CT14nnlo_NF3";
	nlo.params.integration.cuba.maximum_relative_error = 1e-2;
	nlo.params.order = PerturbativeOrder::NLO;
	nlo.params.number_of_threads = number_of_threads;

	if (run("massscaling")) {
		std::cout << "========== SIDIS mass scaling ==========" << IO::endl;

		measure([&] {
			proton_nlo.dis().charm_production_mass_scaling_comparison(
				1.3, 1.5, AnalysisSet::NuTeV_old, x_bins, "Data/DIS/CharmProduction/MassScaling/", "nutev_old_13.csv", "nutev_old_15.csv"
			);
			proton_nlo.dis().charm_production_mass_scaling_comparison(
				1.3, 1.5, AnalysisSet::NuTeV, x_bins, "Data/DIS/CharmProduction/MassScaling/", "nutev_new_13.csv", "nutev_new_15.csv"
			);
			proton_nlo.dis().charm_production_mass_scaling_comparison(
				1.3, 1.5, AnalysisSet::CCFR, x_bins, "Data/DIS/CharmProduction/MassScaling/", "ccfr_13.csv", "ccfr_15.csv"
			);
		});

		std::cout << separator << IO::endl;
	}

	Analysis decay1 = base;
	decay1.params.decay_parametrization_set = DecayParametrization::fit_set_1();
	decay1.params.decay_variation = true;

	Analysis decay2 = base;
	decay2.params.decay_parametrization_set = DecayParametrization::fit_set_2();
	decay2.params.decay_variation = true;

	if (run("decay")) {
		std::cout << "======== SIDIS decay variations ========" << IO::endl;

		measure([&] {
			decay1.sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet1/nutev_old_neutrino.csv");
			decay1.sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet1/nutev_new_neutrino.csv");
			decay1.sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet1/ccfr_neutrino.csv");

			bar(decay1).sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet1/nutev_old_antineutrino.csv");
			bar(decay1).sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet1/nutev_new_antineutrino.csv");
			bar(decay1).sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet1/ccfr_antineutrino.csv");


			decay2.sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet2/nutev_old_neutrino.csv");
			decay2.sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet2/nutev_new_neutrino.csv");
			decay2.sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet2/ccfr_neutrino.csv");

			bar(decay2).sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet2/nutev_old_antineutrino.csv");
			bar(decay2).sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet2/nutev_new_antineutrino.csv");
			bar(decay2).sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet2/ccfr_antineutrino.csv");
		});

		std::cout << separator << IO::endl;
	}

	if (run("pdfvalues")) {
		std::cout << "============== PDF values ==============" << IO::endl;

		std::vector<std::string> all_pdfs;
		all_pdfs.insert(all_pdfs.end(), pdfs.begin(), pdfs.end());
		all_pdfs.insert(all_pdfs.end(), free_pdfs.begin(), free_pdfs.end());
		all_pdfs.push_back("nNNPDF30_nlo_as_0118_A208_Z82");
		all_pdfs.push_back("EPPS21nlo_CT18Anlo_Pb208");

		measure([&] {
			for (const auto &pdf_set_name : all_pdfs) {
				Analysis analysis;
				analysis.params.pdf_set = pdf_set_name;

				const std::string output_folder = "Data/PDF/" + pdf_set_name + "/xf.csv";

				analysis.utility().pdf_values(x_bins_all, {10.0, 100.0}, output_folder);
			}
		});

		std::cout << separator << IO::endl;
	}

	return 0;
}

#endif