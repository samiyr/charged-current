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
#include <ChargedCurrent/Decay/DecayParametrization.cpp>
#include <ChargedCurrent/Interpolation/GridGenerator.cpp>

#include "Constants.cpp"

int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv) {
	using namespace Constants;

	LHAInterface<>::disable_verbosity();

	const std::string separator = "================================================================================";

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

	const std::vector<double> E_beam_bins = Math::linear_space(15.0, 210.0, 5.0);

	const auto measure = [](auto &&func, auto &&...params) {
		const auto &start = std::chrono::high_resolution_clock::now();

		std::forward<decltype(func)>(func)(std::forward<decltype(params)>(params)...);
		
		const auto &stop = std::chrono::high_resolution_clock::now();
		
		const std::chrono::duration<double> elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
		
		std::cout << "Took " << elapsed.count() << " seconds" << IO::endl;
    };

	const std::string epps_set = "EPPS21nlo_CT18Anlo_Fe56";
	const std::string ncteq_set = "nCTEQ15HQ_FullNuc_56_26";
	const std::string nnnpdf_set = "nNNPDF30_nlo_as_0118_A56_Z26";

	const LHAInterface epps(epps_set);
	const LHAInterface ncteq(ncteq_set);
	const LHAInterface nnnpdf(nnnpdf_set);

	const LHASetInterface epps_errors(epps_set);
	const LHASetInterface ncteq_errors(ncteq_set);
	const LHASetInterface nnnpdf_errors(nnnpdf_set);

	const std::vector pdfs = { epps, ncteq, nnnpdf };

	const std::vector pdf_errors = { epps_errors, ncteq_errors, nnnpdf_errors };
	
	const std::vector free_pdfs = {
		LHAInterface<std::true_type>("CT18ANLO"),
		LHAInterface<std::true_type>("nCTEQ15HQ_FullNuc_1_1"),
		LHAInterface<std::true_type>("nNNPDF30_nlo_as_0118_p")
	};

	const std::vector free_pdf_errors = {
		LHASetInterface<std::true_type>("CT18ANLO"),
		LHASetInterface<std::true_type>("nCTEQ15HQ_FullNuc_1_1"),
		LHASetInterface<std::true_type>("nNNPDF30_nlo_as_0118_p")
	};

	const FlavorVector flavors = {Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom};
	const FlavorVector charm_flavors = {Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom};

	const Process process(Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino);
	const Process anti_process(Process::Type::AntiNeutrinoToAntiLepton, Constants::Particles::Proton, Constants::Particles::Neutrino);

	const double charm_mass = 1.3;
	const std::vector<double> charm_masses = {0.0, 1.2, 1.3, 1.4, 1.5};

	const auto renormalization = ScaleDependence::trivial;

	const auto pdf_scale = ScaleDependence::Function([](const TRFKinematics &kinematics) {
		return std::max(1.69 + 1e-6, kinematics.Q2);
	});

	const auto ff_scale = ScaleDependence::Function([](const TRFKinematics &kinematics) {
		return std::max(2.25 + 1e-6, kinematics.Q2);
	});

	const std::filesystem::path decay_grid_folder = std::filesystem::current_path() / "DecayGrids";

	const auto construct_grid_fragmentation_configuration = [&](
		const double decay_muon_min_energy, const DecayParametrization &parametrization, 
		const Particle &target, const Particle &decay_lepton, 
		const std::array<std::string, 4> fragmentation_sets
	) {
		const std::string D0_decay_grid = GridGenerator::grid_filename(
			decay_muon_min_energy, parametrization, Constants::Particles::D0, target, decay_lepton
		);
		const std::string Dp_decay_grid = GridGenerator::grid_filename(
			decay_muon_min_energy, parametrization, Constants::Particles::Dp, target, decay_lepton
		);
		const std::string Ds_decay_grid = GridGenerator::grid_filename(
			decay_muon_min_energy, parametrization, Constants::Particles::Ds, target, decay_lepton
		);
		const std::string LambdaC_decay_grid = GridGenerator::grid_filename(
			decay_muon_min_energy, parametrization, Constants::Particles::LambdaC, target, decay_lepton
		);

		const auto D0_decay_function = DecayFunctions::decay_grid(decay_grid_folder / D0_decay_grid);
		const auto Dp_decay_function = DecayFunctions::decay_grid(decay_grid_folder / Dp_decay_grid);
		const auto Ds_decay_function = DecayFunctions::decay_grid(decay_grid_folder / Ds_decay_grid);
		const auto LambdaC_decay_function = DecayFunctions::decay_grid(decay_grid_folder / LambdaC_decay_grid);

		return FragmentationConfiguration(
			{
				LHAInterface(fragmentation_sets[0], {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
				LHAInterface(fragmentation_sets[1], {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
				LHAInterface(fragmentation_sets[2], {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
				1.14 * LHAInterface(fragmentation_sets[3], {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0})
			},
			{
				Decay(parametrization, Constants::Particles::D0, target, D0_decay_function, decay_muon_min_energy),
				Decay(parametrization, Constants::Particles::Dp, target, Dp_decay_function, decay_muon_min_energy),
				Decay(parametrization, Constants::Particles::Ds, target, Ds_decay_function, decay_muon_min_energy),
				Decay(parametrization, Constants::Particles::LambdaC, target, LambdaC_decay_function, decay_muon_min_energy)
			}
		);
	};
	const auto construct_D0_grid_fragmentation_configuration = [&](
		const double decay_muon_min_energy, const DecayParametrization &parametrization, 
		const Particle &target, const Particle &decay_lepton, 
		const std::string fragmentation_set
	) {
		const std::string D0_decay_grid = GridGenerator::grid_filename(
			decay_muon_min_energy, parametrization, Constants::Particles::D0, target, decay_lepton
		);

		const auto D0_decay_function = DecayFunctions::decay_grid(decay_grid_folder / D0_decay_grid);

		return FragmentationConfiguration(
			{
				LHAInterface(fragmentation_set, {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
			},
			{
				Decay(parametrization, Constants::Particles::D0, target, D0_decay_function, decay_muon_min_energy),
			}
		);
	};

	// const auto construct_analytic_fragmentation_configuration = [](
	// 	const double primary_muon_min_energy, const DecayParametrization &parametrization,
	// 	const Particle &target,
	// 	const std::array<std::string, 4> fragmentation_sets
	// ) {
	// 	const auto decay_function = DecayFunctions::decay_function;

	// 	return FragmentationConfiguration(
	// 		{
	// 			LHAInterface(fragmentation_sets[0], {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
	// 			LHAInterface(fragmentation_sets[1], {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
	// 			LHAInterface(fragmentation_sets[2], {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
	// 			1.14 * LHAInterface(fragmentation_sets[3], {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0})
	// 		},
	// 		{
	// 			Decay(parametrization, Constants::Particles::D0, target, decay_function, primary_muon_min_energy),
	// 			Decay(parametrization, Constants::Particles::Dp, target, decay_function, primary_muon_min_energy),
	// 			Decay(parametrization, Constants::Particles::Ds, target, decay_function, primary_muon_min_energy),
	// 			Decay(parametrization, Constants::Particles::LambdaC, target, decay_function, primary_muon_min_energy)
	// 		}
	// 	);
	// };

	const std::array<std::string, 4> opal_fragmentation = {
		"kkks08_opal_d0___mas",
		"kkks08_opal_d+___mas",
		"bkk05_D3_d_s_nlo",
		"bkk05_D3_lambda_c_nlo"
	};

	const auto grid_fragmentation_set = [&](
		const double decay_muon_min_energy, const Particle &decay_lepton, 
		const std::array<std::string, 4> fragmentation_sets, const DecayParametrization &parametrization = DecayParametrization::fit1()
	) {
		return construct_grid_fragmentation_configuration(decay_muon_min_energy, parametrization, Constants::Particles::Proton, decay_lepton, fragmentation_sets);
	};
	const auto D0_grid_fragmentation_set = [&](
		const double decay_muon_min_energy, const Particle &decay_lepton, 
		const std::string fragmentation_set, const DecayParametrization &parametrization = DecayParametrization::fit1()
	) {
		return construct_D0_grid_fragmentation_configuration(decay_muon_min_energy, parametrization, Constants::Particles::Proton, decay_lepton, fragmentation_set);
	};

	const auto grid_fragmentation = [&](
		const double decay_muon_min_energy, const Particle &decay_lepton, 
		const DecayParametrization &parametrization = DecayParametrization::fit1()
	) {
		return grid_fragmentation_set(decay_muon_min_energy, decay_lepton, opal_fragmentation, parametrization);
	};
	const auto D0_grid_fragmentation = [&](
		const double decay_muon_min_energy, const Particle &decay_lepton, 
		const DecayParametrization &parametrization = DecayParametrization::fit1()
	) {
		return D0_grid_fragmentation_set(decay_muon_min_energy, decay_lepton, opal_fragmentation[0], parametrization);
	};

	const DIS dis(flavors, process, number_of_threads);
	const DIS charm_dis(charm_flavors, process, number_of_threads);
	
	const DIS anti_dis(flavors, anti_process, number_of_threads);
	const DIS anti_charm_dis(charm_flavors, anti_process, number_of_threads);

	const SIDIS sidis(flavors, process, number_of_threads);
	const SIDIS anti_sidis(flavors, anti_process, number_of_threads);

	if (run("dis.differential.charm.errors")) {
		std::cout << "======================== dis.differential.charm.errors =========================" << IO::endl;

		for (const auto &pdf : pdf_errors) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/DIS/CharmProduction/Differential/ErrorSets/" + pdf.set_name + "/";

			measure([&] {
				charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nutev_neutrino.csv"
				);

				charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nutev_antineutrino.csv"
				);

				anti_charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "ccfr_antineutrino.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("dis.differential.total.errors")) {
		std::cout << "========================= dis.differential.total.errors ========================" << IO::endl;

		for (const auto &pdf : pdf_errors) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/DIS/TotalProduction/Differential/ErrorSets/" + pdf.set_name + "/";

			measure([&] {
				dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nutev_neutrino.csv"
				);

				dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nutev_antineutrino.csv"
				);

				anti_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "ccfr_antineutrino.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("dis.integrated.charm.errors")) {
		std::cout << "========================== dis.integrated.charm.errors =========================" << IO::endl;

		for (const auto &pdf : pdf_errors) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/DIS/CharmProduction/Integrated/ErrorSets/" + pdf.set_name + "/";

			measure([&] {
				charm_dis.integrated_errors(
					E_beam_bins, 1.00,
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nomad_neutrino_100.csv"
				);
				charm_dis.integrated_errors(
					E_beam_bins, 1.69,
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nomad_neutrino_169.csv"
				);
				charm_dis.integrated_errors(
					E_beam_bins, 2.25,
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nomad_neutrino_225.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("dis.integrated.total.x")) {
		std::cout << "============================ dis.integrated.total.x ============================" << IO::endl;

		const std::vector<double> E_beams = {20.0, 50.0, 100.0, 150.0, 200.0};
		const std::vector<double> Q2s = Math::linear_space(1.0, 50.0, 1.0);

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/DIS/TotalProduction/xIntegrated/" + pdf.set_name + "/";

			measure([&] {
				dis.x_integrated(
					E_beams, Q2s, charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "min0.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("dis.integrated.total.errors")) {
		std::cout << "========================== dis.integrated.total.errors =========================" << IO::endl;

		for (const auto &pdf : pdf_errors) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/DIS/TotalProduction/Integrated/ErrorSets/" + pdf.set_name + "/";

			measure([&] {
				dis.integrated_errors(
					E_beam_bins, 1.00,
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nomad_neutrino_100.csv"
				);
				dis.integrated_errors(
					E_beam_bins, 1.69,
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nomad_neutrino_169.csv"
				);
				dis.integrated_errors(
					E_beam_bins, 2.25,
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nomad_neutrino_225.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}


	if (run("dis.differential.charm.scales")) {
		std::cout << "======================== dis.differential.charm.scales =========================" << IO::endl;

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/DIS/CharmProduction/Differential/Scales/" + pdf.set_name + "/";

			measure([&] {
				charm_dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					charm_mass, 0.0, 0.0,
					pdf,
					out + "nutev_neutrino.csv"
				);

				charm_dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					charm_mass, 0.0, 0.0,
					pdf,
					out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_charm_dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					charm_mass, 0.0, 0.0,
					pdf,
					out + "nutev_antineutrino.csv"
				);

				anti_charm_dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					charm_mass, 0.0, 0.0,
					pdf,
					out + "ccfr_antineutrino.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("dis.differential.total.scales")) {
		std::cout << "========================= dis.differential.total.scales ========================" << IO::endl;

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/DIS/TotalProduction/Differential/Scales/" + pdf.set_name + "/";

			measure([&] {
				dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					charm_mass, 0.0, 0.0,
					pdf,
					out + "nutev_neutrino.csv"
				);

				dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					charm_mass, 0.0, 0.0,
					pdf,
					out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					charm_mass, 0.0, 0.0,
					pdf,
					out + "nutev_antineutrino.csv"
				);

				anti_dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					charm_mass, 0.0, 0.0,
					pdf,
					out + "ccfr_antineutrino.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("dis.integrated.charm.scales")) {
		std::cout << "========================== dis.integrated.charm.scales =========================" << IO::endl;

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/DIS/CharmProduction/Integrated/Scales/" + pdf.set_name + "/";

			measure([&] {
				charm_dis.integrated_scales(
					E_beam_bins, 1.00,
					charm_mass, 0.0, 0.0,
					pdf,
					out + "nomad_neutrino_100.csv"
				);
				charm_dis.integrated_scales(
					E_beam_bins, 1.69,
					charm_mass, 0.0, 0.0,
					pdf,
					out + "nomad_neutrino_169.csv"
				);
				charm_dis.integrated_scales(
					E_beam_bins, 2.25,
					charm_mass, 0.0, 0.0,
					pdf,
					out + "nomad_neutrino_225.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("dis.integrated.total.scales")) {
		std::cout << "========================== dis.integrated.total.scales =========================" << IO::endl;

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/DIS/TotalProduction/Integrated/Scales/" + pdf.set_name + "/";

			measure([&] {
				dis.integrated_scales(
					E_beam_bins, 1.00,
					charm_mass, 0.0, 0.0,
					pdf,
					out + "nomad_neutrino_100.csv"
				);
				dis.integrated_scales(
					E_beam_bins, 1.69,
					charm_mass, 0.0, 0.0,
					pdf,
					out + "nomad_neutrino_169.csv"
				);
				dis.integrated_scales(
					E_beam_bins, 2.25,
					charm_mass, 0.0, 0.0,
					pdf,
					out + "nomad_neutrino_225.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}


	if (run("sidis.differential.errors")) {
		std::cout <<"=========================== sidis.differential.errors ===========================" << IO::endl;

		const double min_E = 5.0;

		for (const auto &pdf : pdf_errors) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Differential/ErrorSets/" + pdf.set_name + "/";

			measure([&] {
				sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_neutrino.csv"
				);
				sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_antineutrino.csv"
				);
				anti_sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_antineutrino.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("sidis.differential.scales")) {
		std::cout <<"=========================== sidis.differential.scales ===========================" << IO::endl;

		const double min_E = 5.0;

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Differential/Scales/" + pdf.set_name + "/";

			measure([&] {
				sidis.lepton_pair_xy_scales(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					ScaleVariation::All,
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					out + "nutev_neutrino.csv"
				);
				sidis.lepton_pair_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					ScaleVariation::All,
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy_scales(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					ScaleVariation::All,
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					out + "nutev_antineutrino.csv"
				);
				anti_sidis.lepton_pair_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					ScaleVariation::All,
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					out + "ccfr_antineutrino.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("sidis.differential.ffs")) {
		std::cout <<"============================= sidis.differential.ffs ============================" << IO::endl;

		const double min_E = 5.0;

		const std::vector<std::array<std::string, 4>> ff_variations = {
			{"kkks08_belle_d0_mas", "kkks08_belle_d+_mas", "bkk05_D3_d_s_nlo", "bkk05_D3_lambda_c_nlo"},
			{"kkks08_global_d0_mas", "kkks08_global_d+_mas", "bkk05_D3_d_s_nlo", "bkk05_D3_lambda_c_nlo"},
		};
		const std::vector<std::string> folders = {"Belle", "Global"};

		for (const auto &pdf : pdfs) {
			for (const auto &[fragmentation_set, fragmentation_name] : std::views::zip(ff_variations, folders)) {
				std::cout << "PDF set: " << pdf.set_name << IO::endl;
				std::cout << "FF set: " << fragmentation_name << IO::endl;

				const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Differential/FragmentationVariations/" + fragmentation_name + "/" + pdf.set_name + "/";

				const auto grid = grid_fragmentation_set(min_E, Constants::Particles::Muon, fragmentation_set);

				measure([&] {
					sidis.lepton_pair_xy_scales(
						x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
						ScaleVariation::All,
						PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, grid,
						out + "nutev_neutrino.csv"
					);
					sidis.lepton_pair_xy_scales(
						x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
						ScaleVariation::All,
						PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, grid,
						out + "ccfr_neutrino.csv"
					);
				});

				measure([&] {
					anti_sidis.lepton_pair_xy_scales(
						x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
						ScaleVariation::All,
						PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, grid,
						out + "nutev_antineutrino.csv"
					);
					anti_sidis.lepton_pair_xy_scales(
						x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
						ScaleVariation::All,
						PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, grid,
						out + "ccfr_antineutrino.csv"
					);
				});
			}
		}

		std::cout << separator << IO::endl;
	}

	if (run("sidis.differential.nlp")) {
		std::cout <<"============================= sidis.differential.nlp ============================" << IO::endl;

		const double min_E = 5.0;

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Differential/NLP/" + pdf.set_name + "/";

			measure([&] {
				sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					PerturbativeOrder::NLO, true, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_neutrino.csv"
				);
				sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NLO, true, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NLO, true, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_antineutrino.csv"
				);
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NLO, true, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_antineutrino.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("sidis.differential.nnlo")) {
		std::cout <<"============================ sidis.differential.nnlo ============================" << IO::endl;

		const double min_E = 5.0;

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Differential/NNLO/" + pdf.set_name + "/";

			measure([&] {
				sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					PerturbativeOrder::NNLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_neutrino.csv"
				);
				sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NNLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NNLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_antineutrino.csv"
				);
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NNLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_antineutrino.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("sidis.integrated.x")) {
		std::cout << "============================== sidis.integrated.x ==============================" << IO::endl;

		const double min_E_massive = Constants::Particles::Muon.mass;
		const double min_E_massless = Constants::Particles::MasslessMuon.mass;

		const std::vector<double> E_beams = {20.0, 50.0, 100.0, 150.0, 200.0};
		const std::vector<double> Q2s = Math::linear_space(1.0, 50.0, 1.0);

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/xIntegrated/" + pdf.set_name + "/";

			measure([&] {
				sidis.x_integrated_lepton_pair(
					E_beams, Q2s, PerturbativeOrder::NLO, false, charm_mass, min_E_massive,
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "massive_min0.csv"
				);
			});

			const std::string out_massless = "Data/SIDIS/MuonPairProduction/CharmedHadrons/xIntegrated/Massless/Min0/" + pdf.set_name + "/";

			measure([&] {
				sidis.x_integrated_lepton_pair(
					E_beams, Q2s, PerturbativeOrder::NLO, false, charm_mass, min_E_massless,
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out_massless + "massless_min0.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("sidis.integrated.errors.min3")) {
		std::cout << "========================= sidis.integrated.errors.min3 =========================" << IO::endl;
		
		const double min_E = 3.0;

		for (const auto &pdf : pdf_errors) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Integrated/ErrorSets/Massive/Min3/" + pdf.set_name + "/";

			measure([&] {
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.00, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon), 
					renormalization, pdf_scale, ff_scale, 
					out + "nomad_neutrino_100.csv"
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.69, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon), 
					renormalization, pdf_scale, ff_scale, 
					out + "nomad_neutrino_169.csv"
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 2.25, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon), 
					renormalization, pdf_scale, ff_scale, 
					out + "nomad_neutrino_225.csv"
				);
			});

			const std::string out_massless = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Integrated/ErrorSets/Massless/Min3/" + pdf.set_name + "/";

			measure([&] {
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.00, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					renormalization, pdf_scale, ff_scale, 
					out_massless + "nomad_neutrino_100.csv"
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.69, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					renormalization, pdf_scale, ff_scale, 
					out_massless + "nomad_neutrino_169.csv"
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 2.25, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					renormalization, pdf_scale, ff_scale, 
					out_massless + "nomad_neutrino_225.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("sidis.integrated.errors.min0")) {
		std::cout << "========================= sidis.integrated.errors.min0 =========================" << IO::endl;
		
		const double min_E_massive = Constants::Particles::Muon.mass;
		const double min_E_massless = Constants::Particles::MasslessMuon.mass;

		for (const auto &pdf : pdf_errors) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Integrated/ErrorSets/Massive/Min0/" + pdf.set_name + "/";

			measure([&] {
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.00, PerturbativeOrder::NLO, false, charm_mass, min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					renormalization, pdf_scale, ff_scale, 
					out + "nomad_neutrino_100.csv"
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.69, PerturbativeOrder::NLO, false, charm_mass, min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					renormalization, pdf_scale, ff_scale, 
					out + "nomad_neutrino_169.csv"
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 2.25, PerturbativeOrder::NLO, false, charm_mass, min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					renormalization, pdf_scale, ff_scale, 
					out + "nomad_neutrino_225.csv"
				);
			});

			const std::string out_massless = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Integrated/ErrorSets/Massless/Min0/" + pdf.set_name + "/";

			measure([&] {
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.00, PerturbativeOrder::NLO, false, charm_mass, min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					renormalization, pdf_scale, ff_scale, 
					out_massless + "nomad_neutrino_100.csv"
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.69, PerturbativeOrder::NLO, false, charm_mass, min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					renormalization, pdf_scale, ff_scale, 
					out_massless + "nomad_neutrino_169.csv"
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 2.25, PerturbativeOrder::NLO, false, charm_mass, min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					renormalization, pdf_scale, ff_scale, 
					out_massless + "nomad_neutrino_225.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("sidis.integrated.scales.min3")) {
		std::cout << "========================= sidis.integrated.scales.min3 =========================" << IO::endl;
		
		const double min_E = 3.0;

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Integrated/Scales/Massive/Min3/" + pdf.set_name + "/";

			measure([&] {
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.00, ScaleVariation::All, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon), 
					out + "nomad_neutrino_100.csv"
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.69, ScaleVariation::All, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon), 
					out + "nomad_neutrino_169.csv"
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 2.25, ScaleVariation::All, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon), 
					out + "nomad_neutrino_225.csv"
				);
			});

			const std::string out_massless = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Integrated/Scales/Massless/Min3/" + pdf.set_name + "/";

			measure([&] {
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.00, ScaleVariation::All, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					out_massless + "nomad_neutrino_100.csv"
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.69, ScaleVariation::All, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					out_massless + "nomad_neutrino_169.csv"
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 2.25, ScaleVariation::All, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					out_massless + "nomad_neutrino_225.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("sidis.integrated.scales.min0")) {
		std::cout << "========================= sidis.integrated.scales.min0 =========================" << IO::endl;
		
		const double min_E_massive = Constants::Particles::Muon.mass;
		const double min_E_massless = Constants::Particles::MasslessMuon.mass;

		for (const auto &pdf : free_pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Integrated/Scales/Massive/Min0/" + pdf.set_name + "/";

			measure([&] {
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.00, ScaleVariation::All, PerturbativeOrder::NLO, false, charm_mass, min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					out + "nomad_neutrino_100.csv"
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.69, ScaleVariation::All, PerturbativeOrder::NLO, false, charm_mass, min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					out + "nomad_neutrino_169.csv"
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 2.25, ScaleVariation::All, PerturbativeOrder::NLO, false, charm_mass, min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					out + "nomad_neutrino_225.csv"
				);
			});

			const std::string out_massless = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Integrated/Scales/Massless/Min0/" + pdf.set_name + "/";

			measure([&] {
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.00, ScaleVariation::All, PerturbativeOrder::NLO, false, charm_mass, min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					out_massless + "nomad_neutrino_100.csv"
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.69, ScaleVariation::All, PerturbativeOrder::NLO, false, charm_mass, min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					out_massless + "nomad_neutrino_169.csv"
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 2.25, ScaleVariation::All, PerturbativeOrder::NLO, false, charm_mass, min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					out_massless + "nomad_neutrino_225.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("sidis.differential.isospin.errors")) {
		std::cout << "======================= sidis.differential.isospin.errors ======================" << IO::endl;
		
		const double min_E = 5.0;

		for (auto pdf : free_pdf_errors) {
			pdf.Z = 26;
			pdf.A = 56;

			const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Differential/Isospin/ErrorSets/" + pdf.set_name + "/";

			measure([&] {
				sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_neutrino.csv"
				);
				sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_antineutrino.csv"
				);
				anti_sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_antineutrino.csv"
				);
			});
		}
	}

	if (run("sidis.differential.d0")) {
		std::cout <<"============================= sidis.differential.d0 =============================" << IO::endl;

		const double min_E = 5.0;

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/SIDIS/MuonPairProduction/D0/Differential/" + pdf.set_name + "/";

			measure([&] {
				sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, D0_grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_neutrino.csv"
				);
				sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, D0_grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, D0_grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_antineutrino.csv"
				);
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, D0_grid_fragmentation(min_E, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_antineutrino.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}
	
	// if (run("channels")) {
	// 	std::cout << "======== SIDIS parton channels =========" << IO::endl;

	// 	measure([&] {
	// 		nlo.sidis().muon_pair_production_quark_gluon_channels(x_bins, NuTeV::New::Neutrino::y_bins, NuTeV::New::Neutrino::E_bins, 
	// 			{
	// 				"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_quark_to_quark.csv",
	// 				"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_quark_to_gluon.csv",
	// 				"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_gluon_to_quark.csv",
	// 				"Data/SIDIS/MuonPairProduction/CharmedHadrons/ChannelDecomposition/nutev_new_gluon_to_gluon.csv"
	// 			}
	// 		);
	// 	});

	// 	std::cout << separator << IO::endl;
	// }

	// if (run("flavors")) {
	// 	std::cout << "======== SIDIS flavor channels =========" << IO::endl;

	// 	measure([&] {
	// 		nlo.sidis().muon_pair_production_flavor_decomposition_quark_to_quark(
	// 			x_bins, NuTeV::New::Neutrino::y_bins, NuTeV::New::Neutrino::E_bins,
	// 			"Data/SIDIS/MuonPairProduction/CharmedHadrons/FlavorDecomposition/nutev_new"
	// 		);
	// 		nlo.sidis().muon_pair_production_flavor_decomposition_gluon_to_quark(
	// 			x_bins, NuTeV::New::Neutrino::y_bins, NuTeV::New::Neutrino::E_bins,
	// 			"Data/SIDIS/MuonPairProduction/CharmedHadrons/FlavorDecomposition/nutev_new"
	// 		);
	// 	});

	// 	std::cout << separator << IO::endl;
	// }

	// if (run("fragmentation")) {
	// 	std::cout << "===== SIDIS fragmentation channels =====" << IO::endl;

	// 	measure([&] {
	// 		nlo.sidis().muon_pair_production_fragmentation_decomposition(
	// 			x_bins, NuTeV::New::Neutrino::y_bins, NuTeV::New::Neutrino::E_bins,
	// 			"Data/SIDIS/MuonPairProduction/CharmedHadrons/FragmentationDecomposition/nutev_new"
	// 		);
	// 	});

	// 	std::cout << separator << IO::endl;
	// }

	// if (run("masses")) {
	// 	std::cout << "=========== SIDIS charm mass ===========" << IO::endl;

	// 	measure([&] {
	// 		for (const double charm_mass : charm_masses) {
	// 			Analysis mass = nlo;
	// 			mass.params.charm_mass = charm_mass;

	// 			const std::string folder = std::to_string(static_cast<int>(charm_mass * 10));

	// 			mass.sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/Masses/" + folder + "/nutev_old_neutrino.csv");
	// 			mass.sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/Masses/" + folder + "/nutev_new_neutrino.csv");
	// 			mass.sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/Masses/" + folder + "/ccfr_neutrino.csv");

	// 			bar(mass).sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/Masses/" + folder + "/nutev_old_antineutrino.csv");
	// 			bar(mass).sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/Masses/" + folder + "/nutev_new_antineutrino.csv");
	// 			bar(mass).sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/Masses/" + folder + "/ccfr_antineutrino.csv");
	// 		}
	// 	});

	// 	std::cout << separator << IO::endl;
	// }

	// Analysis proton_nlo;
	// proton_nlo.params.pdf_set = "CT14nnlo_NF3";
	// nlo.params.integration.cuba.maximum_relative_error = 1e-2;
	// nlo.params.order = PerturbativeOrder::NLO;
	// nlo.params.number_of_threads = number_of_threads;

	// if (run("massscaling")) {
	// 	std::cout << "========== SIDIS mass scaling ==========" << IO::endl;

	// 	measure([&] {
	// 		proton_nlo.dis().charm_production_mass_scaling_comparison(
	// 			1.3, 1.5, AnalysisSet::NuTeV_old, x_bins, "Data/DIS/CharmProduction/MassScaling/", "nutev_old_13.csv", "nutev_old_15.csv"
	// 		);
	// 		proton_nlo.dis().charm_production_mass_scaling_comparison(
	// 			1.3, 1.5, AnalysisSet::NuTeV, x_bins, "Data/DIS/CharmProduction/MassScaling/", "nutev_new_13.csv", "nutev_new_15.csv"
	// 		);
	// 		proton_nlo.dis().charm_production_mass_scaling_comparison(
	// 			1.3, 1.5, AnalysisSet::CCFR, x_bins, "Data/DIS/CharmProduction/MassScaling/", "ccfr_13.csv", "ccfr_15.csv"
	// 		);
	// 	});

	// 	std::cout << separator << IO::endl;
	// }

	// Analysis decay1 = base;
	// decay1.params.decay_parametrization_set = DecayParametrization::fit_set_1();
	// decay1.params.decay_variation = true;

	// Analysis decay2 = base;
	// decay2.params.decay_parametrization_set = DecayParametrization::fit_set_2();
	// decay2.params.decay_variation = true;

	// if (run("decay")) {
	// 	std::cout << "======== SIDIS decay variations ========" << IO::endl;

	// 	measure([&] {
	// 		decay1.sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet1/nutev_old_neutrino.csv");
	// 		decay1.sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet1/nutev_new_neutrino.csv");
	// 		decay1.sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet1/ccfr_neutrino.csv");

	// 		bar(decay1).sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet1/nutev_old_antineutrino.csv");
	// 		bar(decay1).sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet1/nutev_new_antineutrino.csv");
	// 		bar(decay1).sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet1/ccfr_antineutrino.csv");


	// 		decay2.sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet2/nutev_old_neutrino.csv");
	// 		decay2.sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet2/nutev_new_neutrino.csv");
	// 		decay2.sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet2/ccfr_neutrino.csv");

	// 		bar(decay2).sidis().muon_pair_production(AnalysisSet::NuTeV_old, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet2/nutev_old_antineutrino.csv");
	// 		bar(decay2).sidis().muon_pair_production(AnalysisSet::NuTeV, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet2/nutev_new_antineutrino.csv");
	// 		bar(decay2).sidis().muon_pair_production(AnalysisSet::CCFR, x_bins, "Data/SIDIS/MuonPairProduction/CharmedHadrons/DecayVariation/FitSet2/ccfr_antineutrino.csv");
	// 	});

	// 	std::cout << separator << IO::endl;
	// }

	// if (run("pdfvalues")) {
	// 	std::cout << "============== PDF values ==============" << IO::endl;

	// 	std::vector<std::string> all_pdfs;
	// 	all_pdfs.insert(all_pdfs.end(), pdfs.begin(), pdfs.end());
	// 	all_pdfs.insert(all_pdfs.end(), free_pdfs.begin(), free_pdfs.end());
	// 	all_pdfs.push_back("nNNPDF30_nlo_as_0118_A208_Z82");
	// 	all_pdfs.push_back("EPPS21nlo_CT18Anlo_Pb208");

	// 	measure([&] {
	// 		for (const auto &pdf_set_name : all_pdfs) {
	// 			Analysis analysis;
	// 			analysis.params.pdf_set = pdf_set_name;

	// 			const std::string output_folder = "Data/PDF/" + pdf_set_name + "/xf.csv";

	// 			analysis.utility().pdf_values(x_bins_all, {10.0, 100.0}, output_folder);
	// 		}
	// 	});

	// 	std::cout << separator << IO::endl;
	// }

	if (run("utility.decay.grid")) {
		std::cout << "============================== utility.decay.grid ==============================" << IO::endl;

		const std::vector<double> zyE_bins = Math::linear_space(1.0, 300.0, 1e-2);

		GridGenerator generator(number_of_threads);

		std::vector<DecayParametrization> parametrizations;
		parametrizations.push_back(DecayParametrization::fit1());
		parametrizations.push_back(DecayParametrization::fit2());

		generator.generate_decay_grids(
			"DecayGrids", {1.0, 5.0, 10.0, 100.0, 300.0}, {Constants::Particles::Muon.mass}, parametrizations,
			{
				Constants::Particles::D0, Constants::Particles::Dp, Constants::Particles::Ds, Constants::Particles::LambdaC
			}, Constants::Particles::Proton, Constants::Particles::Muon
		);
		generator.generate_decay_grids(
			"DecayGrids", zyE_bins, {3.0, 5.0}, parametrizations,
			{
				Constants::Particles::D0, Constants::Particles::Dp, Constants::Particles::Ds, Constants::Particles::LambdaC
			}, Constants::Particles::Proton, Constants::Particles::Muon
		);

		generator.generate_decay_grids(
			"DecayGrids", {1.0, 5.0, 10.0, 100.0, 300.0}, {Constants::Particles::MasslessMuon.mass}, parametrizations,
			{
				Constants::Particles::D0, Constants::Particles::Dp, Constants::Particles::Ds, Constants::Particles::LambdaC
			}, Constants::Particles::Proton, Constants::Particles::MasslessMuon
		);
		generator.generate_decay_grids(
			"DecayGrids", zyE_bins, {3.0, 5.0}, parametrizations,
			{
				Constants::Particles::D0, Constants::Particles::Dp, Constants::Particles::Ds, Constants::Particles::LambdaC
			}, Constants::Particles::Proton, Constants::Particles::MasslessMuon
		);

		std::cout << separator << IO::endl;
	}

	// if (run("sidisdoubledifferential")) {
	// 	std::cout << "====== SIDIS double differential =======" << IO::endl;
	// 	for (const auto &pdf_set_name : pdfs) {
	// 		std::cout << "PDF set: " << pdf_set_name << IO::endl;

	// 		Analysis analysis = base;
	// 		analysis.params.pdf_error_sets = false;
	// 		analysis.params.pdf_set = pdf_set_name;

	// 		const std::string output_folder = "Data/SIDIS/MuonPairProduction/CharmedHadrons/DoubleDifferential/" + pdf_set_name + "/";

	// 		const std::vector<double> y_bins = {1.0};
	// 		const std::vector<double> E_bins = Math::linear_space(4.0, 100.0, 2.0);

	// 		measure([&] {
	// 			analysis.sidis().muon_pair_production(x_bins, y_bins, E_bins, output_folder + "neutrino.csv");
	// 			bar(analysis).sidis().muon_pair_production(x_bins, y_bins, E_bins, output_folder + "antineutrino.csv");
	// 		});
	// 	}

	// 	std::cout << separator << IO::endl;
	// }
	// if (run("disdoubledifferential")) {
	// 	std::cout << "======= DIS double differential ========" << IO::endl;
	// 	for (const auto &pdf_set_name : pdfs) {
	// 		std::cout << "PDF set: " << pdf_set_name << IO::endl;

	// 		Analysis analysis = base;
	// 		analysis.params.pdf_error_sets = false;
	// 		analysis.params.pdf_set = pdf_set_name;

	// 		const std::string output_folder = "Data/DIS/TotalProduction/DoubleDifferential/" + pdf_set_name + "/";

	// 		const std::vector<double> y_bins = {1.0};
	// 		const std::vector<double> E_bins = Math::linear_space(4.0, 100.0, 2.0);

	// 		measure([&] {
	// 			analysis.dis().total_production(x_bins, y_bins, E_bins, output_folder + "neutrino.csv");
	// 			bar(analysis).dis().total_production(x_bins, y_bins, E_bins, output_folder + "antineutrino.csv");
	// 		});
	// 	}

	// 	std::cout << separator << IO::endl;
	// }

	return 0;
}

#endif