#ifndef CHARGED_CURRENT_DIS_H
#define CHARGED_CURRENT_DIS_H

#define BS_THREAD_POOL_ENABLE_PAUSE

#include <iostream>
#include <array>
#include <format>
#include <chrono>
#include <ranges>

#include <ChargedCurrent/DIS/DIS.cpp>
#include <ChargedCurrent/SIDIS/SIDIS.cpp>
#include <ChargedCurrent/PDF/Interfaces/LHAInterface.cpp>
#include <ChargedCurrent/Common/TRFKinematics.cpp>
#include <ChargedCurrent/Decay/DecayParametrization.cpp>
#include <ChargedCurrent/Interpolation/GridGenerator.cpp>
#include <ChargedCurrent/Threading/BS_thread_pool.hpp>
#include <ChargedCurrent/Threading/ThreadResult.cpp>

#include "Constants.cpp"

int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv) {
	using namespace Constants;

	LHAInterface<>::disable_verbosity();

	const std::string separator = "=================================================================================";

	std::vector<std::string> arguments(argv + 1, argv + argc);

	auto nthreads_iterator = std::find(arguments.begin(), arguments.end(), "--nthreads");
	unsigned int number_of_threads = 1;

	if (nthreads_iterator == arguments.end()) {
		std::cout << "No --nthreads argument given, defaulting to 1" << IO::endl;
	} else {
		nthreads_iterator++;
		number_of_threads = static_cast<unsigned int>(std::stoul(*nthreads_iterator));
	}

	BS::thread_pool thread_pool(number_of_threads);
	thread_pool.pause();

	ThreadResultCollection results;

	auto output_directory_iterator = std::find(arguments.begin(), arguments.end(), "--output");
	std::filesystem::path output_path = std::filesystem::current_path();

	if (output_directory_iterator != arguments.end()) {
		output_directory_iterator++;
		output_path = std::filesystem::path(*output_directory_iterator);
	}
	const std::string output_dir = output_path.string() + "/";
	std::cout << "Current output directory: " << output_dir << IO::endl;

	auto pdf_iterator = std::find(arguments.begin(), arguments.end(), "--pdf");
	std::optional<std::string> pdf_argument = std::nullopt;

	if (pdf_iterator != arguments.end()) {
		pdf_iterator++;
		pdf_argument = *pdf_iterator;
		std::cout << "Selected PDF: " << *pdf_argument << IO::endl;
	}

	unsigned int variation_range_count = 0;

	auto variation_start_iterator = std::find(arguments.begin(), arguments.end(), "--variation-start");
	std::size_t variation_start_index = 0;

	if (variation_start_iterator != arguments.end()) {
		variation_start_iterator++;
		variation_start_index = static_cast<std::size_t>(std::stoull(*variation_start_iterator));
		variation_range_count++;
	}
	auto variation_end_iterator = std::find(arguments.begin(), arguments.end(), "--variation-end");
	std::size_t variation_end_index = 0;

	if (variation_end_iterator != arguments.end()) {
		variation_end_iterator++;
		variation_end_index = static_cast<std::size_t>(std::stoull(*variation_end_iterator));
		variation_range_count++;
	}

	const bool custom_variation_range = variation_range_count == 2;

	if (custom_variation_range) {
		std::cout << "Selected variation range: " << variation_start_index << "..." << variation_end_index;
		std::cout << " (length: " << variation_end_index - variation_start_index + 1 << ")" << IO::endl;
	}

	const std::optional<std::pair<std::size_t, std::size_t>> variation_range = custom_variation_range 
																				? std::make_optional(std::make_pair(variation_start_index, variation_end_index))
																				: std::nullopt;

	const bool all = Collections::contains(arguments, "all");
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

	const std::vector pdfs = pdf_argument.has_value() ? std::vector{ LHAInterface(*pdf_argument) } : std::vector{ epps, ncteq, nnnpdf };

	const std::vector pdf_errors = pdf_argument.has_value() ? std::vector{ LHASetInterface(*pdf_argument) } : std::vector{ epps_errors, ncteq_errors, nnnpdf_errors };
	
	const std::vector free_pdfs = pdf_argument.has_value() ? std::vector{ LHAInterface<std::true_type>(*pdf_argument) } : std::vector{
		LHAInterface<std::true_type>("CT18ANLO"),
		LHAInterface<std::true_type>("nCTEQ15HQ_FullNuc_1_1"),
		LHAInterface<std::true_type>("nNNPDF30_nlo_as_0118_p")
	};

	const std::vector free_pdf_errors = pdf_argument.has_value() ? std::vector{ LHASetInterface<std::true_type>(*pdf_argument) } : std::vector{
		LHASetInterface<std::true_type>("CT18ANLO"),
		LHASetInterface<std::true_type>("nCTEQ15HQ_FullNuc_1_1"),
		LHASetInterface<std::true_type>("nNNPDF30_nlo_as_0118_p")
	};

	const FlavorVector flavors = {Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom};
	const FlavorVector charm_flavors = {Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom};

	const Process process(Process::Type::NeutrinoToLepton, Constants::Particles::Proton, Constants::Particles::Neutrino);
	const Process anti_process(Process::Type::AntiNeutrinoToAntiLepton, Constants::Particles::Proton, Constants::Particles::Neutrino);

	const std::vector<double> charm_masses = {0.0, 1.2, 1.3, 1.4, 1.5};


	const auto scale = [](const auto &pdf) {
		return ScaleDependence::Function([&](const TRFKinematics &kinematics) {
			return kinematics.Q2;
		}, std::pow(pdf.quark_mass(Flavor::Charm), 2));
	};

	const auto ff_scale = ScaleDependence::trivial;

	const std::filesystem::path decay_grid_folder = std::filesystem::current_path() / "DecayGrids";

	const auto construct_grid_fragmentation_configuration = [&](
		const double decay_muon_min_energy, const DecayParametrization &parametrization, 
		const Particle &target, const Particle &decay_lepton, 
		const std::array<std::string, 4> fragmentation_sets,
		const std::optional<std::vector<bool>> &outgoing = std::nullopt
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

		std::vector<double> outgoing_multipliers{1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0};
		if (outgoing) {
			const std::vector<double> converted((*outgoing).begin(), (*outgoing).end());
			outgoing_multipliers *= converted;
		}

		return FragmentationConfiguration(
			{
				LHAInterface<std::false_type, FreezeExtrapolator>(fragmentation_sets[0], outgoing_multipliers), 
				LHAInterface<std::false_type, FreezeExtrapolator>(fragmentation_sets[1], outgoing_multipliers), 
				LHAInterface<std::false_type, FreezeExtrapolator>(fragmentation_sets[2], outgoing_multipliers), 
				1.14 * LHAInterface<std::false_type, FreezeExtrapolator>(fragmentation_sets[3], outgoing_multipliers)
			},
			{
				Decay(parametrization, Constants::Particles::D0, target, D0_decay_function, decay_muon_min_energy),
				Decay(parametrization, Constants::Particles::Dp, target, Dp_decay_function, decay_muon_min_energy),
				Decay(parametrization, Constants::Particles::Ds, target, Ds_decay_function, decay_muon_min_energy),
				Decay(parametrization, Constants::Particles::LambdaC, target, LambdaC_decay_function, decay_muon_min_energy)
			}
		);
	};
	const auto construct_single_grid_fragmentation_configuration = [&](
		const double decay_muon_min_energy, const DecayParametrization &parametrization, 
		const Particle &target, const Particle &decay_lepton, const Particle &resonance,
		const std::string fragmentation_set,
		const std::optional<std::vector<bool>> &outgoing = std::nullopt
	) {
		const std::string decay_grid = GridGenerator::grid_filename(
			decay_muon_min_energy, parametrization, resonance, target, decay_lepton
		);

		const auto decay_function = DecayFunctions::decay_grid(decay_grid_folder / decay_grid);

		std::vector<double> outgoing_multipliers{1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0};
		if (outgoing) {
			const std::vector<double> converted((*outgoing).begin(), (*outgoing).end());
			outgoing_multipliers *= converted;
		}

		return FragmentationConfiguration(
			{
				LHAInterface<std::false_type, FreezeExtrapolator>(fragmentation_set, outgoing_multipliers)
			},
			{
				Decay(parametrization, resonance, target, decay_function, decay_muon_min_energy)
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
				LHAInterface<std::false_type, FreezeExtrapolator>(fragmentation_set, {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
			},
			{
				Decay(parametrization, Constants::Particles::D0, target, D0_decay_function, decay_muon_min_energy),
			}
		);
	};
	
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
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nutev_neutrino.csv",
					variation_range
				);

				charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nutev_antineutrino.csv",
					variation_range
				);

				anti_charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "ccfr_antineutrino.csv",
					variation_range
				);
			});
		}

		std::cout << separator << IO::endl;
	}
	if (run("dis.differential.charm.errors.fixed")) {
		std::cout << "===================== dis.differential.charm.errors.fixed ======================" << IO::endl;

		const auto constant_scale = ScaleDependence::constant(20.0);

		for (const auto &pdf : pdf_errors) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/DIS/CharmProduction/Differential/ErrorSets/FixedScale/" + pdf.set_name + "/";

			measure([&] {
				charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					constant_scale, constant_scale,
					output_dir + out + "nutev_neutrino.csv",
					variation_range
				);

				charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					constant_scale, constant_scale,
					output_dir + out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					constant_scale, constant_scale,
					output_dir + out + "nutev_antineutrino.csv",
					variation_range
				);

				anti_charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					constant_scale, constant_scale,
					output_dir + out + "ccfr_antineutrino.csv",
					variation_range
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
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nutev_neutrino.csv",
					variation_range
				);

				dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nutev_antineutrino.csv",
					variation_range
				);

				anti_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "ccfr_antineutrino.csv",
					variation_range
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
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nomad_neutrino_100.csv",
					variation_range
				);
				charm_dis.integrated_errors(
					E_beam_bins, 1.69,
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nomad_neutrino_169.csv",
					variation_range
				);
				charm_dis.integrated_errors(
					E_beam_bins, 2.25,
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nomad_neutrino_225.csv",
					variation_range
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("dis.integrated.total.x")) {
		std::cout << "============================ dis.integrated.total.x ============================" << IO::endl;

		const std::vector<double> E_beams = {20.0, 50.0, 100.0, 150.0, 200.0};
		const std::vector<double> Q2s_base = Math::linear_space(1.0, 50.0, 1.0);
		const std::vector<double> small_Q2s = Math::linear_space(0.1, 4.0, 0.01);
		const std::vector<double> Q2s = Collections::join(Q2s_base, small_Q2s);

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/DIS/TotalProduction/xIntegrated/" + pdf.set_name + "/";

			measure([&] {
				dis.x_integrated(
					E_beams, Q2s, pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "min0.csv"
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
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nomad_neutrino_100.csv",
					variation_range
				);
				dis.integrated_errors(
					E_beam_bins, 1.69,
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nomad_neutrino_169.csv",
					variation_range
				);
				dis.integrated_errors(
					E_beam_bins, 2.25,
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nomad_neutrino_225.csv",
					variation_range
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
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nutev_neutrino.csv",
					variation_range
				);

				charm_dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_charm_dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nutev_antineutrino.csv",
					variation_range
				);

				anti_charm_dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "ccfr_antineutrino.csv",
					variation_range
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
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nutev_neutrino.csv",
					variation_range
				);

				dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nutev_antineutrino.csv",
					variation_range
				);

				anti_dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "ccfr_antineutrino.csv",
					variation_range
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
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nomad_neutrino_100.csv",
					variation_range
				);
				charm_dis.integrated_scales(
					E_beam_bins, 1.69,
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nomad_neutrino_169.csv",
					variation_range
				);
				charm_dis.integrated_scales(
					E_beam_bins, 2.25,
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nomad_neutrino_225.csv",
					variation_range
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
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nomad_neutrino_100.csv",
					variation_range
				);
				dis.integrated_scales(
					E_beam_bins, 1.69,
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nomad_neutrino_169.csv",
					variation_range
				);
				dis.integrated_scales(
					E_beam_bins, 2.25,
					pdf.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf,
					scale(pdf), scale(pdf),
					output_dir + out + "nomad_neutrino_225.csv",
					variation_range
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("dis.differential.charm.channels")) {
		std::cout << "======================== dis.differential.charm.channels =======================" << IO::endl;

		const auto run_analysis = [&](const auto &pdf_in, const std::string out) {
			measure([&] {
				charm_dis.differential_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					pdf_in.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf_in,
					scale(pdf_in), scale(pdf_in),
					out + "nutev_neutrino.csv"
				);

				charm_dis.differential_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					pdf_in.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf_in,
					scale(pdf_in), scale(pdf_in),
					out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_charm_dis.differential_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					pdf_in.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf_in,
					scale(pdf_in), scale(pdf_in),
					out + "nutev_antineutrino.csv"
				);

				anti_charm_dis.differential_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					pdf_in.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf_in,
					scale(pdf_in), scale(pdf_in),
					out + "ccfr_antineutrino.csv"
				);
			});
		};

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/DIS/CharmProduction/Differential/Channels/" + pdf.set_name + "/";

			const auto up_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0}
			);
			const auto anti_up_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);
			
			const auto down_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);
			const auto anti_down_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);

			const auto charm_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0}
			);
			const auto anti_charm_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);
			
			const auto strange_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}
			);
			const auto anti_strange_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);
			
			const auto bottom_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0}
			);
			const auto anti_bottom_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);

			const auto gluon_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);

			run_analysis(up_pdf, output_dir + out + "u/");
			run_analysis(down_pdf, output_dir + out + "d/");
			run_analysis(charm_pdf, output_dir + out + "c/");
			run_analysis(strange_pdf, output_dir + out + "s/");
			run_analysis(bottom_pdf, output_dir + out + "b/");

			run_analysis(gluon_pdf, output_dir + out + "g/");

			run_analysis(anti_up_pdf, output_dir + out + "ubar/");
			run_analysis(anti_down_pdf, output_dir + out + "dbar/");
			run_analysis(anti_charm_pdf, output_dir + out + "cbar/");
			run_analysis(anti_strange_pdf, output_dir + out + "sbar/");
			run_analysis(anti_bottom_pdf, output_dir + out + "bbar/");
		}

		std::cout << separator << IO::endl;
	}
	if (run("dis.differential.total.channels")) {
		std::cout << "======================== dis.differential.total.channels =======================" << IO::endl;

		const auto run_analysis = [&](const auto &pdf_in, const std::string out) {
			measure([&] {
				dis.differential_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					pdf_in.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf_in,
					scale(pdf_in), scale(pdf_in),
					out + "nutev_neutrino.csv"
				);

				dis.differential_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					pdf_in.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf_in,
					scale(pdf_in), scale(pdf_in),
					out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_dis.differential_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					pdf_in.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf_in,
					scale(pdf_in), scale(pdf_in),
					out + "nutev_antineutrino.csv"
				);

				anti_dis.differential_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					pdf_in.quark_mass(Flavor::Charm), 0.0, 0.0,
					pdf_in,
					scale(pdf_in), scale(pdf_in),
					out + "ccfr_antineutrino.csv"
				);
			});
		};

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/DIS/TotalProduction/Differential/Channels/" + pdf.set_name + "/";

			const auto up_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0}
			);
			const auto anti_up_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);
			
			const auto down_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);
			const auto anti_down_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);

			const auto charm_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0}
			);
			const auto anti_charm_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);
			
			const auto strange_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}
			);
			const auto anti_strange_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);
			
			const auto bottom_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0}
			);
			const auto anti_bottom_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);

			const auto gluon_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);

			run_analysis(up_pdf, output_dir + out + "u/");
			run_analysis(down_pdf, output_dir + out + "d/");
			run_analysis(charm_pdf, output_dir + out + "c/");
			run_analysis(strange_pdf, output_dir + out + "s/");
			run_analysis(bottom_pdf, output_dir + out + "b/");

			run_analysis(gluon_pdf, output_dir + out + "g/");

			run_analysis(anti_up_pdf, output_dir + out + "ubar/");
			run_analysis(anti_down_pdf, output_dir + out + "dbar/");
			run_analysis(anti_charm_pdf, output_dir + out + "cbar/");
			run_analysis(anti_strange_pdf, output_dir + out + "sbar/");
			run_analysis(anti_bottom_pdf, output_dir + out + "bbar/");
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
				results.add(sidis.lepton_pair_xy_errors(thread_pool,
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nutev_neutrino.csv",
					variation_range
				));
				results.add(sidis.lepton_pair_xy_errors(thread_pool,
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "ccfr_neutrino.csv",
					variation_range
				));
			});

			measure([&] {
				results.add(anti_sidis.lepton_pair_xy_errors(thread_pool,
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nutev_antineutrino.csv",
					variation_range
				));
				results.add(anti_sidis.lepton_pair_xy_errors(thread_pool,
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "ccfr_antineutrino.csv",
					variation_range
				));
			});
		}

		std::cout << separator << IO::endl;
	}
	/*

	if (run("sidis.differential.errors.fixed")) {
		std::cout <<"======================== sidis.differential.errors.fixed ========================" << IO::endl;

		const double min_E = 5.0;

		const auto constant_scale = ScaleDependence::constant(20.0);

		for (const auto &pdf : pdf_errors) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Differential/ErrorSets/FixedScale/" + pdf.set_name + "/";

			measure([&] {
				sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					constant_scale, constant_scale, constant_scale,
					output_dir + out + "nutev_neutrino.csv",
					variation_range
				);
				sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					constant_scale, constant_scale, constant_scale,
					output_dir + out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					constant_scale, constant_scale, constant_scale,
					output_dir + out + "nutev_antineutrino.csv",
					variation_range
				);
				anti_sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					constant_scale, constant_scale, constant_scale,
					output_dir + out + "ccfr_antineutrino.csv",
					variation_range
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
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nutev_neutrino.csv",
					variation_range
				);
				sidis.lepton_pair_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					ScaleVariation::All,
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy_scales(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					ScaleVariation::All,
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nutev_antineutrino.csv",
					variation_range
				);
				anti_sidis.lepton_pair_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					ScaleVariation::All,
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "ccfr_antineutrino.csv",
					variation_range
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("sidis.differential.ffs")) {
		std::cout <<"============================= sidis.differential.ffs ============================" << IO::endl;

		const double min_E = 5.0;

		const std::vector<std::array<std::string, 4>> ff_variations = {
			{"kkks08_belle_d0__mas", "kkks08_belle_d+__mas", "bkk05_D3_d_s_nlo", "bkk05_D3_lambda_c_nlo"},
			{"kkks08_global_d0_mas", "kkks08_global_d+_mas", "bkk05_D3_d_s_nlo", "bkk05_D3_lambda_c_nlo"},
		};
		const std::vector<std::string> folders = {"Belle", "Global"};

		for (const auto &pdf : pdf_errors) {
			for (const auto &[fragmentation_set, fragmentation_name] : std::views::zip(ff_variations, folders)) {
				std::cout << "PDF set: " << pdf.set_name << IO::endl;
				std::cout << "FF set: " << fragmentation_name << IO::endl;

				const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Differential/FragmentationVariations/" + fragmentation_name + "/" + pdf.set_name + "/";

				const auto grid = grid_fragmentation_set(min_E, Constants::Particles::MasslessMuon, fragmentation_set);

				measure([&] {
					sidis.lepton_pair_xy_errors(
						x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
						PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
						pdf, grid,
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "nutev_neutrino.csv",
						variation_range
					);
					sidis.lepton_pair_xy_errors(
						x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
						PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
						pdf, grid,
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "ccfr_neutrino.csv",
						variation_range
					);
				});

				measure([&] {
					anti_sidis.lepton_pair_xy_errors(
						x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
						PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
						pdf, grid,
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "nutev_antineutrino.csv",
						variation_range
					);
					anti_sidis.lepton_pair_xy_errors(
						x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
						PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
						pdf, grid,
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "ccfr_antineutrino.csv",
						variation_range
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
					PerturbativeOrder::NLO, true, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nutev_neutrino.csv"
				);
				sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NLO, true, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NLO, true, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nutev_antineutrino.csv"
				);
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NLO, true, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "ccfr_antineutrino.csv"
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
					PerturbativeOrder::NNLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nutev_neutrino.csv"
				);
				sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NNLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NNLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nutev_antineutrino.csv"
				);
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NNLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "ccfr_antineutrino.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("sidis.differential.zprofile")) {
		std::cout <<"========================== sidis.differential.zprofile ==========================" << IO::endl;

		const std::vector<double> min_Es = {0.0, 3.0, 5.0};
		const std::vector<std::string> folders = {"ZeroLimit", "3_0", "5_0"};

		const std::vector<double> xs{0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4};
		const std::vector<double> z_bins = Math::linear_space(0.0, 1.0, 1e-2);

		for (const auto &[current_min_E, folder] : std::views::zip(min_Es, folders)) {
			for (const auto &pdf : pdfs) {
				std::cout << "PDF set: " << pdf.set_name << IO::endl;

				const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Differential/zprofile/" + pdf.set_name + "/" + folder + "/";
				
				// To fix the compiler error 'capturing a structured binding is not yet supported in OpenMP'
				const double min_E = current_min_E;

				measure([&] {
					sidis.lepton_pair_zxy(
						z_bins, xs, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
						PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "nutev_neutrino.csv"
					);
					sidis.lepton_pair_zxy(
						z_bins, xs, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
						PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "ccfr_neutrino.csv"
					);
				});

				measure([&] {
					anti_sidis.lepton_pair_zxy(
						z_bins, xs, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
						PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "nutev_antineutrino.csv"
					);
					anti_sidis.lepton_pair_zxy(
						z_bins, xs, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
						PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "ccfr_antineutrino.csv"
					);
				});
			}
		}

		std::cout << separator << IO::endl;
	}

	if (run("sidis.integrated.x")) {
		std::cout << "============================== sidis.integrated.x ==============================" << IO::endl;

		const double min_E_massive = Constants::Particles::Muon.mass;
		const double min_E_massless = Constants::Particles::MasslessMuon.mass;

		const std::vector<double> E_beams = {20.0, 50.0, 100.0, 150.0, 200.0};
		const std::vector<double> Q2s_base = Math::linear_space(1.0, 50.0, 1.0);
		const std::vector<double> small_Q2s = Math::linear_space(0.1, 4.0, 0.01);
		const std::vector<double> Q2s = Collections::join(Q2s_base, small_Q2s);

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/xIntegrated/" + pdf.set_name + "/";

			measure([&] {
				sidis.x_integrated_lepton_pair(
					E_beams, Q2s, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E_massive,
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "massive_min0.csv"
				);
				sidis.x_integrated_lepton_pair(
					E_beams, Q2s, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E_massless,
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "massless_min0.csv"
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
					E_beam_bins, 1.00, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon), 
					scale(pdf), scale(pdf), ff_scale, 
					output_dir + out + "nomad_neutrino_100.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.69, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon), 
					scale(pdf), scale(pdf), ff_scale, 
					output_dir + out + "nomad_neutrino_169.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 2.25, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon), 
					scale(pdf), scale(pdf), ff_scale, 
					output_dir + out + "nomad_neutrino_225.csv",
					variation_range
				);
			});

			const std::string out_massless = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Integrated/ErrorSets/Massless/Min3/" + pdf.set_name + "/";

			measure([&] {
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.00, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					scale(pdf), scale(pdf), ff_scale, 
					output_dir + out_massless + "nomad_neutrino_100.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.69, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					scale(pdf), scale(pdf), ff_scale, 
					output_dir + out_massless + "nomad_neutrino_169.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 2.25, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					scale(pdf), scale(pdf), ff_scale, 
					output_dir + out_massless + "nomad_neutrino_225.csv",
					variation_range
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
					E_beam_bins, 1.00, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					scale(pdf), scale(pdf), ff_scale, 
					output_dir + out + "nomad_neutrino_100.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.69, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					scale(pdf), scale(pdf), ff_scale, 
					output_dir + out + "nomad_neutrino_169.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 2.25, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					scale(pdf), scale(pdf), ff_scale, 
					output_dir + out + "nomad_neutrino_225.csv",
					variation_range
				);
			});

			const std::string out_massless = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Integrated/ErrorSets/Massless/Min0/" + pdf.set_name + "/";

			measure([&] {
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.00, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					scale(pdf), scale(pdf), ff_scale, 
					output_dir + out_massless + "nomad_neutrino_100.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.69, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					scale(pdf), scale(pdf), ff_scale, 
					output_dir + out_massless + "nomad_neutrino_169.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 2.25, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					scale(pdf), scale(pdf), ff_scale, 
					output_dir + out_massless + "nomad_neutrino_225.csv",
					variation_range
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
					E_beam_bins, 1.00, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon), 
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nomad_neutrino_100.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.69, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon), 
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nomad_neutrino_169.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 2.25, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon), 
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nomad_neutrino_225.csv",
					variation_range
				);
			});

			const std::string out_massless = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Integrated/Scales/Massless/Min3/" + pdf.set_name + "/";

			measure([&] {
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.00, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					scale(pdf), scale(pdf), ff_scale,
					out_massless + "nomad_neutrino_100.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.69, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					scale(pdf), scale(pdf), ff_scale,
					out_massless + "nomad_neutrino_169.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 2.25, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					scale(pdf), scale(pdf), ff_scale,
					out_massless + "nomad_neutrino_225.csv",
					variation_range
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
					E_beam_bins, 1.00, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nomad_neutrino_100.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.69, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nomad_neutrino_169.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 2.25, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nomad_neutrino_225.csv",
					variation_range
				);
			});

			const std::string out_massless = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Integrated/Scales/Massless/Min0/" + pdf.set_name + "/";

			measure([&] {
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.00, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					scale(pdf), scale(pdf), ff_scale,
					out_massless + "nomad_neutrino_100.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.69, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					scale(pdf), scale(pdf), ff_scale,
					out_massless + "nomad_neutrino_169.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 2.25, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					scale(pdf), scale(pdf), ff_scale,
					out_massless + "nomad_neutrino_225.csv",
					variation_range
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
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nutev_neutrino.csv",
					variation_range
				);
				sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nutev_antineutrino.csv",
					variation_range
				);
				anti_sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "ccfr_antineutrino.csv",
					variation_range
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
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, D0_grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nutev_neutrino.csv"
				);
				sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, D0_grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, D0_grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nutev_antineutrino.csv"
				);
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf, D0_grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "ccfr_antineutrino.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}
	
	if (run("sidis.differential.channels")) {
		std::cout << "========================== sidis.differential.channels =========================" << IO::endl;

		const double min_E = 5.0;

		const auto charm_fragmentation = construct_grid_fragmentation_configuration(
			min_E, DecayParametrization::fit1(), Constants::Particles::Proton, Constants::Particles::MasslessMuon,
			opal_fragmentation, std::make_optional(std::vector{false, false, false, false, false, false, false, false, false, false, true, false, false})
		);
		const auto anti_charm_fragmentation = construct_grid_fragmentation_configuration(
			min_E, DecayParametrization::fit1(), Constants::Particles::Proton, Constants::Particles::MasslessMuon,
			opal_fragmentation, std::make_optional(std::vector{false, false, true, false, false, false, false, false, false, false, false, false, false})
		);

		const auto run_analysis = [&](const auto &pdf_in, const auto &ff_in, const std::string out) {
			measure([&] {
				sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					PerturbativeOrder::NLO, false, pdf_in.quark_mass(Flavor::Charm), 0.0,
					pdf_in, ff_in,
					scale(pdf_in), scale(pdf_in), ff_scale,
					out + "nutev_neutrino.csv"
				);
				sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NLO, false, pdf_in.quark_mass(Flavor::Charm), 0.0,
					pdf_in, ff_in,
					scale(pdf_in), scale(pdf_in), ff_scale,
					out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NLO, false, pdf_in.quark_mass(Flavor::Charm), 0.0,
					pdf_in, ff_in,
					scale(pdf_in), scale(pdf_in), ff_scale,
					out + "nutev_antineutrino.csv"
				);
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NLO, false, pdf_in.quark_mass(Flavor::Charm), 0.0,
					pdf_in, ff_in,
					scale(pdf_in), scale(pdf_in), ff_scale,
					out + "ccfr_antineutrino.csv"
				);
			});
		};

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Differential/Channels/" + pdf.set_name + "/";

			const auto down_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);
			const auto anti_down_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);

			const auto strange_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}
			);
			const auto anti_strange_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);

			const auto gluon_pdf = LHAInterface(
				pdf.set_name,
				{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			);

			run_analysis(down_pdf, charm_fragmentation, output_dir + out + "d2c/");
			run_analysis(strange_pdf, charm_fragmentation, output_dir + out + "s2c/");
			run_analysis(gluon_pdf, charm_fragmentation, output_dir + out + "g2c/");

			run_analysis(anti_down_pdf, anti_charm_fragmentation, output_dir + out + "dbar2cbar/");
			run_analysis(anti_strange_pdf, anti_charm_fragmentation, output_dir + out + "sbar2cbar/");
			run_analysis(gluon_pdf, anti_charm_fragmentation, output_dir + out + "g2cbar/");
		}

		std::cout << separator << IO::endl;
	}

	if (run("sidis.differential.fragmentation")) {
		std::cout << "======================= sidis.differential.fragmentation =======================" << IO::endl;

		const double min_E = 5.0;
		const std::array<Particle, 4> resonances{Constants::Particles::D0, Constants::Particles::Dp, Constants::Particles::Ds, Constants::Particles::LambdaC};

		for (std::size_t i = 0; i < opal_fragmentation.size(); i++) {
			const auto fragmentation = construct_single_grid_fragmentation_configuration(
				min_E, DecayParametrization::fit1(), Constants::Particles::Proton, Constants::Particles::MasslessMuon, 
				resonances[i], opal_fragmentation[i]
			);

			for (const auto &pdf : pdfs) {
				std::cout << "PDF set: " << pdf.set_name << IO::endl;

				const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Differential/Fragmentation/" + pdf.set_name + "/" + opal_fragmentation[i] + "/";

				measure([&] {
					sidis.lepton_pair_xy(
						x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
						PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
						pdf, fragmentation,
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "nutev_neutrino.csv"
					);
					sidis.lepton_pair_xy(
						x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
						PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
						pdf, fragmentation,
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "ccfr_neutrino.csv"
					);
				});

				measure([&] {
					anti_sidis.lepton_pair_xy(
						x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
						PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
						pdf, fragmentation,
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "nutev_antineutrino.csv"
					);
					anti_sidis.lepton_pair_xy(
						x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
						PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
						pdf, fragmentation,
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "ccfr_antineutrino.csv"
					);
				});
			}
		}

		std::cout << separator << IO::endl;
	}
	
	if (run("sidis.differential.masses")) {
		std::cout <<"=========================== sidis.differential.masses ===========================" << IO::endl;

		const double min_E = 5.0;
		const std::vector<std::string> folders = {"Massless", "1_2", "1_3", "1_4", "1_5"};

		for (const auto &[current_mass, folder] : std::views::zip(charm_masses, folders)) {
			for (const auto &pdf : pdfs) {
				std::cout << "PDF set: " << pdf.set_name << IO::endl;

				const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Differential/CharmMasses/" + pdf.set_name + "/" + folder + "/";

				// To fix the compiler error 'capturing a structured binding is not yet supported in OpenMP'
				const double mass = current_mass;

				measure([&] {
					sidis.lepton_pair_xy(
						x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
						PerturbativeOrder::NLO, false, mass, 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "nutev_neutrino.csv"
					);
					sidis.lepton_pair_xy(
						x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
						PerturbativeOrder::NLO, false, mass, 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "ccfr_neutrino.csv"
					);
				});

				measure([&] {
					anti_sidis.lepton_pair_xy(
						x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
						PerturbativeOrder::NLO, false, mass, 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "nutev_antineutrino.csv"
					);
					anti_sidis.lepton_pair_xy(
						x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
						PerturbativeOrder::NLO, false, mass, 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "ccfr_antineutrino.csv"
					);
				});
			}
		}

		std::cout << separator << IO::endl;
	}

	if (run("dis.differential.charm.masses")) {
		const LHASetInterface pdf("CT14nnlo_NF3");

		const std::vector<std::string> folders = {"1_3", "1_5"};

		const auto alt_scale = [](const double charm_mass) {
			return ScaleDependence::Function([&](const TRFKinematics &kinematics) {
				return kinematics.Q2 + std::pow(charm_mass, 2);
			});
		};

		for (const auto &[current_mass, folder] : std::views::zip(std::vector{1.3, 1.5}, folders)) {
			const std::string out = "Data/DIS/CharmProduction/Differential/Masses/" + folder + "/";

			// To fix the compiler error 'capturing a structured binding is not yet supported in OpenMP'
			const double mass = current_mass;

			measure([&] {
				charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					mass, 0.0, 0.0,
					pdf,
					alt_scale(mass), alt_scale(mass),
					output_dir + out + "nutev_neutrino.csv",
					variation_range
				);

				charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					mass, 0.0, 0.0,
					pdf,
					alt_scale(mass), alt_scale(mass),
					output_dir + out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					mass, 0.0, 0.0,
					pdf,
					alt_scale(mass), alt_scale(mass),
					output_dir + out + "nutev_antineutrino.csv",
					variation_range
				);

				anti_charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					mass, 0.0, 0.0,
					pdf,
					alt_scale(mass), alt_scale(mass),
					output_dir + out + "ccfr_antineutrino.csv",
					variation_range
				);
			});
		}
	}

	if (run("sidis.differential.decays")) {
		std::cout <<"=========================== sidis.differential.decays ===========================" << IO::endl;

		const double min_E = 5.0;

		const std::vector<DecayParametrization> &fit_set_1_parametrizations = DecayParametrization::fit_set_1();
		const std::vector<DecayParametrization> &fit_set_2_parametrizations = DecayParametrization::fit_set_2();

		std::vector<FragmentationConfiguration<LHAInterface<std::false_type, FreezeExtrapolator>, DecayFunctions::DecayGrid>> fit_set_1;
		for (const DecayParametrization &parametrization : fit_set_1_parametrizations) {
			fit_set_1.push_back(
				grid_fragmentation(min_E, Constants::Particles::MasslessMuon, parametrization)
			);
		}
		
		std::vector<FragmentationConfiguration<LHAInterface<std::false_type, FreezeExtrapolator>, DecayFunctions::DecayGrid>> fit_set_2;
		for (const DecayParametrization &parametrization : fit_set_2_parametrizations) {
			fit_set_2.push_back(
				grid_fragmentation(min_E, Constants::Particles::MasslessMuon, parametrization)
			);
		}
		
		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			for (const auto &[fit_set, current_folder] : std::views::zip(std::vector{fit_set_1, fit_set_2}, std::vector{"FitSet1", "FitSet2"})) {
				// To fix the compiler error 'capturing a structured binding is not yet supported in OpenMP'
				const auto &parametrizations = fit_set;
				const std::string folder = current_folder;

				const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Differential/Decays/" + pdf.set_name + "/" + folder + "/";

				measure([&] {
					sidis.lepton_pair_xy_decays(
						x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
						parametrizations, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
						pdf,
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "nutev_neutrino.csv",
						variation_range
					);
					sidis.lepton_pair_xy_decays(
						x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
						parametrizations, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
						pdf,
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "ccfr_neutrino.csv",
						variation_range
					);
				});

				measure([&] {
					anti_sidis.lepton_pair_xy_decays(
						x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
						parametrizations, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
						pdf,
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "nutev_antineutrino.csv",
						variation_range
					);
					anti_sidis.lepton_pair_xy_decays(
						x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
						parametrizations, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
						pdf,
						scale(pdf), scale(pdf), ff_scale,
						output_dir + out + "ccfr_antineutrino.csv",
						variation_range
					);
				});
			}
		}

		std::cout << separator << IO::endl;
	}
	if (run("sidis.differential.decays.fitset3")) {
		std::cout <<"======================= sidis.differential.decays.fitset3 =======================" << IO::endl;

		const double min_E = 5.0;

		const std::vector<DecayParametrization> &fit_set_3_parametrizations = DecayParametrization::fit_set_3();
		
		std::vector<FragmentationConfiguration<LHAInterface<std::false_type, FreezeExtrapolator>, DecayFunctions::DecayGrid>> parametrizations;
		for (const DecayParametrization &parametrization : fit_set_3_parametrizations) {
			parametrizations.push_back(
				grid_fragmentation(min_E, Constants::Particles::MasslessMuon, parametrization)
			);
		}

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Differential/Decays/" + pdf.set_name + "/" + "FitSet3" + "/";

			measure([&] {
				sidis.lepton_pair_xy_decays(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					parametrizations, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf,
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nutev_neutrino.csv",
					variation_range
				);
				sidis.lepton_pair_xy_decays(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					parametrizations, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf,
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy_decays(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					parametrizations, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf,
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "nutev_antineutrino.csv",
					variation_range
				);
				anti_sidis.lepton_pair_xy_decays(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					parametrizations, PerturbativeOrder::NLO, false, pdf.quark_mass(Flavor::Charm), 0.0,
					pdf,
					scale(pdf), scale(pdf), ff_scale,
					output_dir + out + "ccfr_antineutrino.csv",
					variation_range
				);
			});
			
		}

		std::cout << separator << IO::endl;
	}

	if (run("utility.decay.grid")) {
		std::cout << "============================== utility.decay.grid ==============================" << IO::endl;

		const std::vector<double> zyE_bins = Math::linear_space(1.0, 300.0, 1.0);

		GridGenerator generator(number_of_threads);

		std::vector<DecayParametrization> parametrizations;
		parametrizations.push_back(DecayParametrization::fit1());
		parametrizations.push_back(DecayParametrization::fit2());

		const std::vector<DecayParametrization> &fit_set_1 = DecayParametrization::fit_set_1();
		const std::vector<DecayParametrization> &fit_set_2 = DecayParametrization::fit_set_2();

		parametrizations.insert(parametrizations.end(), fit_set_1.begin(), fit_set_1.end());
		parametrizations.insert(parametrizations.end(), fit_set_2.begin(), fit_set_2.end());

		const std::size_t start_index = custom_variation_range ? (*variation_range).first : 0;
		const std::size_t end_index = custom_variation_range ? (*variation_range).second : std::numeric_limits<std::size_t>::max();

		if (!custom_variation_range || Math::in_interval(start_index, end_index, static_cast<std::size_t>(0))) {
			generator.generate_decay_grids(
				output_dir + "DecayGrids", {1.0, 5.0, 10.0, 100.0, 300.0}, {Constants::Particles::MasslessMuon.mass}, {DecayParametrization::fit1()},
				{
					Constants::Particles::D0, Constants::Particles::Dp, Constants::Particles::Ds, Constants::Particles::LambdaC
				}, Constants::Particles::Proton, Constants::Particles::MasslessMuon
			);
		}

		if (!custom_variation_range || Math::in_interval(start_index, end_index, static_cast<std::size_t>(1))) {
			generator.generate_decay_grids(
				output_dir + "DecayGrids", zyE_bins, {3.0, 5.0}, parametrizations,
				{
					Constants::Particles::D0, Constants::Particles::Dp, Constants::Particles::Ds, Constants::Particles::LambdaC
				}, Constants::Particles::Proton, Constants::Particles::Muon
			);
		}

		if (!custom_variation_range || Math::in_interval(start_index, end_index, static_cast<std::size_t>(2))) {
			generator.generate_decay_grids(
				output_dir + "DecayGrids", {1.0, 5.0, 10.0, 100.0, 300.0}, {Constants::Particles::Muon.mass}, {DecayParametrization::fit1()},
				{
					Constants::Particles::D0, Constants::Particles::Dp, Constants::Particles::Ds, Constants::Particles::LambdaC
				}, Constants::Particles::Proton, Constants::Particles::Muon
			);
		}
		
		if (!custom_variation_range || Math::in_interval(start_index, end_index, static_cast<std::size_t>(3))) {
			generator.generate_decay_grids(
				output_dir + "DecayGrids", zyE_bins, {3.0, 5.0}, parametrizations,
				{
					Constants::Particles::D0, Constants::Particles::Dp, Constants::Particles::Ds, Constants::Particles::LambdaC
				}, Constants::Particles::Proton, Constants::Particles::MasslessMuon
			);
		}
		
		if (!custom_variation_range || Math::in_interval(start_index, end_index, static_cast<std::size_t>(4))) {
			generator.generate_decay_grids(
				output_dir + "DecayGrids", zyE_bins, {5.0}, DecayParametrization::fit_set_3(),
				{
					Constants::Particles::D0, Constants::Particles::Dp, Constants::Particles::Ds, Constants::Particles::LambdaC
				}, Constants::Particles::Proton, Constants::Particles::MasslessMuon
			);
		}

		std::cout << separator << IO::endl;
	}


	if (run("utility.pdf.values")) {
		const auto pdf_reader = [](const auto &pdf, const std::vector<double> &xs, const double Q2, const FlavorType flavor, const std::filesystem::path output) {
			IO::create_directory_tree(output);
			std::ofstream file(output);

			file << "x,y,E,Q2,NLO" << IO::endl;

			for (const double x : xs) {
				file << x << ", " << 1 << ", " << 1 << ", " << Q2 << ", " << pdf.xf_evaluate(flavor, x, Q2) << IO::endl;
			}

			file.close();
		};

		const std::vector<double> xs = Math::linear_space(0.01, 0.4, 1e-2);
		const std::vector<double> Q2s{5.0, 10.0, 20.0, 100.0};
		const std::vector<std::string> flavors{"xtbar", "xbbar", "xcbar", "xsbar", "xubar", "xdbar", "xg", "xd", "xu", "xs", "xc", "xb", "xt"};

		for (const auto &pdf_set : pdf_errors) {
			for (const double Q2 : Q2s) {
				const std::string out = "Data/PDF/" + pdf_set.set_name + "/" + std::to_string(Q2) + "/";

				for (unsigned int variation_index = 0; variation_index < pdf_set.size(); variation_index++) {
					for (FlavorType flavor = -6; flavor <= 6; flavor++) {
						const std::string filename = IO::leading_zeroes(variation_index, 4);
						pdf_reader(
							pdf_set[variation_index], xs, Q2, flavor,
							output_dir + out + flavors[static_cast<std::size_t>(flavor + 6)] + "/" + filename + ".csv"
						);
					}
				}
			}
		}
	}*/

	if (run("utility.pdf.testvalues")) {
		const auto pdf_reader = [](const auto &pdf, const std::vector<double> &xs, const double Q2, const FlavorType flavor, const std::filesystem::path output) {
			IO::create_directory_tree(output);
			std::ofstream file(output);

			file << "x,y,E,Q2,NLO" << IO::endl;

			for (const double x : xs) {
				pdf.evaluate(x, Q2);
				file << x << ", " << 1 << ", " << 1 << ", " << Q2 << ", " << pdf.xf(flavor) << IO::endl;
			}

			file.close();
		};

		const std::vector<double> xs = Math::log_space(1e-6, 0.999, 1e-2);
		const std::vector<double> Q2s{10.0, 100.0};
		const std::vector<std::string> flavors{"xtbar", "xbbar", "xcbar", "xsbar", "xubar", "xdbar", "xg", "xd", "xu", "xs", "xc", "xb", "xt"};

		auto epps_isospin_pdf = LHASetInterface<std::true_type>("CT18ANLO");
		epps_isospin_pdf.Z = 82.0;
		epps_isospin_pdf.A = 208.0;

		auto nnnpdf_isospin_pdf = LHASetInterface<std::true_type>("nNNPDF30_nlo_as_0118_p");
		nnnpdf_isospin_pdf.Z = 54.0;
		nnnpdf_isospin_pdf.A = 108.0;

		const auto test_sets = {
			// LHASetInterface("EPPS21nlo_CT18Anlo_Pb208"),
			// LHASetInterface("CT18ANLO"),
			// LHASetInterface("nNNPDF30_nlo_as_0118_A108_Z54"),
			// isospin_pdf
			epps_isospin_pdf
		};

		for (const auto &pdf_set : test_sets) {
			for (const double Q2 : Q2s) {
				const std::string out = "Data/PDF/" + pdf_set.set_name + "/" + std::to_string(Q2) + "/";

				for (unsigned int variation_index = 0; variation_index < pdf_set.size(); variation_index++) {
					for (FlavorType flavor = -6; flavor <= 6; flavor++) {
						const std::string filename = IO::leading_zeroes(variation_index, 4);
						pdf_reader(
							pdf_set[variation_index], xs, Q2, flavor,
							output_dir + out + flavors[static_cast<std::size_t>(flavor + 6)] + "/" + filename + ".csv"
						);
					}
				}
			}
		}
	}

	const auto repeat = [](std::function<void(void)> func, unsigned int interval) {
		std::thread([func, interval]() { 
    		while (true) { 
      			auto x = std::chrono::steady_clock::now() + std::chrono::milliseconds(interval);
				func();
				std::this_thread::sleep_until(x);
			}
		}).detach();
	};

	std::cout << thread_pool.get_tasks_queued() << " tasks queued" << IO::endl;
	repeat([&]() {
		std::cout << "Running: " << thread_pool.get_tasks_running() << "\tQueued: " << thread_pool.get_tasks_queued() << IO::endl;
	}, 1000);
	thread_pool.unpause();
	thread_pool.wait();
	std::cout << "All tasks finished, writing to file" << IO::endl;
	results.write();
	std::cout << "All results written to file" << IO::endl;

	return 0;
}

#endif