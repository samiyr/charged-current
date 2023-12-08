#ifndef CHARGED_CURRENT_DIS_H
#define CHARGED_CURRENT_DIS_H

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
				LHAInterface(fragmentation_sets[0], outgoing_multipliers), 
				LHAInterface(fragmentation_sets[1], outgoing_multipliers), 
				LHAInterface(fragmentation_sets[2], outgoing_multipliers), 
				1.14 * LHAInterface(fragmentation_sets[3], outgoing_multipliers)
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
				LHAInterface(fragmentation_set, outgoing_multipliers)
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
				LHAInterface(fragmentation_set, {1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0}), 
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
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nutev_neutrino.csv",
					variation_range
				);

				charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nutev_antineutrino.csv",
					variation_range
				);

				anti_charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "ccfr_antineutrino.csv",
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
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nutev_neutrino.csv",
					variation_range
				);

				dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nutev_antineutrino.csv",
					variation_range
				);

				anti_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "ccfr_antineutrino.csv",
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
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nomad_neutrino_100.csv",
					variation_range
				);
				charm_dis.integrated_errors(
					E_beam_bins, 1.69,
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nomad_neutrino_169.csv",
					variation_range
				);
				charm_dis.integrated_errors(
					E_beam_bins, 2.25,
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nomad_neutrino_225.csv",
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
					out + "nomad_neutrino_100.csv",
					variation_range
				);
				dis.integrated_errors(
					E_beam_bins, 1.69,
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nomad_neutrino_169.csv",
					variation_range
				);
				dis.integrated_errors(
					E_beam_bins, 2.25,
					charm_mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nomad_neutrino_225.csv",
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
					charm_mass, 0.0, 0.0,
					pdf,
					0.0, 1.69 + 1e-6,
					out + "nutev_neutrino.csv",
					variation_range
				);

				charm_dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					charm_mass, 0.0, 0.0,
					pdf,
					0.0, 1.69 + 1e-6,
					out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_charm_dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					charm_mass, 0.0, 0.0,
					pdf,
					0.0, 1.69 + 1e-6,
					out + "nutev_antineutrino.csv",
					variation_range
				);

				anti_charm_dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					charm_mass, 0.0, 0.0,
					pdf,
					0.0, 1.69 + 1e-6,
					out + "ccfr_antineutrino.csv",
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
					charm_mass, 0.0, 0.0,
					pdf,
					0.0, 1.69 + 1e-6,
					out + "nutev_neutrino.csv",
					variation_range
				);

				dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					charm_mass, 0.0, 0.0,
					pdf,
					0.0, 1.69 + 1e-6,
					out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					charm_mass, 0.0, 0.0,
					pdf,
					0.0, 1.69 + 1e-6,
					out + "nutev_antineutrino.csv",
					variation_range
				);

				anti_dis.differential_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					charm_mass, 0.0, 0.0,
					pdf,
					0.0, 1.69 + 1e-6,
					out + "ccfr_antineutrino.csv",
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
					charm_mass, 0.0, 0.0,
					pdf,
					0.0, 1.69 + 1e-6,
					out + "nomad_neutrino_100.csv",
					variation_range
				);
				charm_dis.integrated_scales(
					E_beam_bins, 1.69,
					charm_mass, 0.0, 0.0,
					pdf,
					0.0, 1.69 + 1e-6,
					out + "nomad_neutrino_169.csv",
					variation_range
				);
				charm_dis.integrated_scales(
					E_beam_bins, 2.25,
					charm_mass, 0.0, 0.0,
					pdf,
					0.0, 1.69 + 1e-6,
					out + "nomad_neutrino_225.csv",
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
					charm_mass, 0.0, 0.0,
					pdf,
					0.0, 1.69 + 1e-6,
					out + "nomad_neutrino_100.csv",
					variation_range
				);
				dis.integrated_scales(
					E_beam_bins, 1.69,
					charm_mass, 0.0, 0.0,
					pdf,
					0.0, 1.69 + 1e-6,
					out + "nomad_neutrino_169.csv",
					variation_range
				);
				dis.integrated_scales(
					E_beam_bins, 2.25,
					charm_mass, 0.0, 0.0,
					pdf,
					0.0, 1.69 + 1e-6,
					out + "nomad_neutrino_225.csv",
					variation_range
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
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_neutrino.csv",
					variation_range
				);
				sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_antineutrino.csv",
					variation_range
				);
				anti_sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_antineutrino.csv",
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
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					0.0, 1.69 + 1e-6, 2.25 + 1e-6,
					out + "nutev_neutrino.csv",
					variation_range
				);
				sidis.lepton_pair_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					ScaleVariation::All,
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					0.0, 1.69 + 1e-6, 2.25 + 1e-6,
					out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy_scales(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					ScaleVariation::All,
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					0.0, 1.69 + 1e-6, 2.25 + 1e-6,
					out + "nutev_antineutrino.csv",
					variation_range
				);
				anti_sidis.lepton_pair_xy_scales(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					ScaleVariation::All,
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					0.0, 1.69 + 1e-6, 2.25 + 1e-6,
					out + "ccfr_antineutrino.csv",
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
						PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, grid,
						renormalization, pdf_scale, ff_scale,
						out + "nutev_neutrino.csv",
						variation_range
					);
					sidis.lepton_pair_xy_errors(
						x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
						PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, grid,
						renormalization, pdf_scale, ff_scale,
						out + "ccfr_neutrino.csv",
						variation_range
					);
				});

				measure([&] {
					anti_sidis.lepton_pair_xy_errors(
						x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
						PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, grid,
						renormalization, pdf_scale, ff_scale,
						out + "nutev_antineutrino.csv",
						variation_range
					);
					anti_sidis.lepton_pair_xy_errors(
						x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
						PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, grid,
						renormalization, pdf_scale, ff_scale,
						out + "ccfr_antineutrino.csv",
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
					PerturbativeOrder::NLO, true, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_neutrino.csv"
				);
				sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NLO, true, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NLO, true, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_antineutrino.csv"
				);
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NLO, true, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
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
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_neutrino.csv"
				);
				sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NNLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NNLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_antineutrino.csv"
				);
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NNLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_antineutrino.csv"
				);
			});
		}

		std::cout << separator << IO::endl;
	}

	if (run("sidis.differential.zprofile")) {
		std::cout <<"========================== sidis.differential.zprofile ==========================" << IO::endl;

		const std::vector<double> min_Es = {0.0, Constants::Particles::Muon.mass, 3.0, 5.0};
		const std::vector<std::string> folders = {"ZeroLimit", "0_106", "3_0", "5_0"};

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
						PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						renormalization, pdf_scale, ff_scale,
						out + "nutev_neutrino.csv"
					);
					sidis.lepton_pair_zxy(
						z_bins, xs, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
						PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						renormalization, pdf_scale, ff_scale,
						out + "ccfr_neutrino.csv"
					);
				});

				measure([&] {
					anti_sidis.lepton_pair_zxy(
						z_bins, xs, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
						PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						renormalization, pdf_scale, ff_scale,
						out + "nutev_antineutrino.csv"
					);
					anti_sidis.lepton_pair_zxy(
						z_bins, xs, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
						PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						renormalization, pdf_scale, ff_scale,
						out + "ccfr_antineutrino.csv"
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
					E_beams, Q2s, PerturbativeOrder::NLO, false, charm_mass, min_E_massive,
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon),
					renormalization, pdf_scale, ff_scale,
					out + "massive_min0.csv"
				);
				sidis.x_integrated_lepton_pair(
					E_beams, Q2s, PerturbativeOrder::NLO, false, charm_mass, min_E_massless,
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "massless_min0.csv"
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
					out + "nomad_neutrino_100.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.69, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon), 
					renormalization, pdf_scale, ff_scale, 
					out + "nomad_neutrino_169.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 2.25, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon), 
					renormalization, pdf_scale, ff_scale, 
					out + "nomad_neutrino_225.csv",
					variation_range
				);
			});

			const std::string out_massless = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Integrated/ErrorSets/Massless/Min3/" + pdf.set_name + "/";

			measure([&] {
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.00, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					renormalization, pdf_scale, ff_scale, 
					out_massless + "nomad_neutrino_100.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.69, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					renormalization, pdf_scale, ff_scale, 
					out_massless + "nomad_neutrino_169.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 2.25, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					renormalization, pdf_scale, ff_scale, 
					out_massless + "nomad_neutrino_225.csv",
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
					E_beam_bins, 1.00, PerturbativeOrder::NLO, false, charm_mass, min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					renormalization, pdf_scale, ff_scale, 
					out + "nomad_neutrino_100.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.69, PerturbativeOrder::NLO, false, charm_mass, min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					renormalization, pdf_scale, ff_scale, 
					out + "nomad_neutrino_169.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 2.25, PerturbativeOrder::NLO, false, charm_mass, min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					renormalization, pdf_scale, ff_scale, 
					out + "nomad_neutrino_225.csv",
					variation_range
				);
			});

			const std::string out_massless = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Integrated/ErrorSets/Massless/Min0/" + pdf.set_name + "/";

			measure([&] {
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.00, PerturbativeOrder::NLO, false, charm_mass, min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					renormalization, pdf_scale, ff_scale, 
					out_massless + "nomad_neutrino_100.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 1.69, PerturbativeOrder::NLO, false, charm_mass, min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					renormalization, pdf_scale, ff_scale, 
					out_massless + "nomad_neutrino_169.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_errors(
					E_beam_bins, 2.25, PerturbativeOrder::NLO, false, charm_mass, min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					renormalization, pdf_scale, ff_scale, 
					out_massless + "nomad_neutrino_225.csv",
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
					E_beam_bins, 1.00, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon), 
					0.0, 1.69 + 1e-6, 2.25 + 1e-6,
					out + "nomad_neutrino_100.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.69, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon), 
					0.0, 1.69 + 1e-6, 2.25 + 1e-6,
					out + "nomad_neutrino_169.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 2.25, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::Muon), 
					0.0, 1.69 + 1e-6, 2.25 + 1e-6,
					out + "nomad_neutrino_225.csv",
					variation_range
				);
			});

			const std::string out_massless = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Integrated/Scales/Massless/Min3/" + pdf.set_name + "/";

			measure([&] {
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.00, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					0.0, 1.69 + 1e-6, 2.25 + 1e-6,
					out_massless + "nomad_neutrino_100.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.69, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					0.0, 1.69 + 1e-6, 2.25 + 1e-6,
					out_massless + "nomad_neutrino_169.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 2.25, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, charm_mass, min_E, 
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon), 
					0.0, 1.69 + 1e-6, 2.25 + 1e-6,
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
					E_beam_bins, 1.00, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, charm_mass, min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					0.0, 1.69 + 1e-6, 2.25 + 1e-6,
					out + "nomad_neutrino_100.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.69, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, charm_mass, min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					0.0, 1.69 + 1e-6, 2.25 + 1e-6,
					out + "nomad_neutrino_169.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 2.25, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, charm_mass, min_E_massive, 
					pdf, grid_fragmentation(min_E_massive, Constants::Particles::Muon), 
					0.0, 1.69 + 1e-6, 2.25 + 1e-6,
					out + "nomad_neutrino_225.csv",
					variation_range
				);
			});

			const std::string out_massless = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Integrated/Scales/Massless/Min0/" + pdf.set_name + "/";

			measure([&] {
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.00, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, charm_mass, min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					0.0, 1.69 + 1e-6, 2.25 + 1e-6,
					out_massless + "nomad_neutrino_100.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 1.69, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, charm_mass, min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					0.0, 1.69 + 1e-6, 2.25 + 1e-6,
					out_massless + "nomad_neutrino_169.csv",
					variation_range
				);
				sidis.integrated_lepton_pair_scales(
					E_beam_bins, 2.25, ScaleVariation::RenormalizationFactorization, PerturbativeOrder::NLO, false, charm_mass, min_E_massless, 
					pdf, grid_fragmentation(min_E_massless, Constants::Particles::MasslessMuon), 
					0.0, 1.69 + 1e-6, 2.25 + 1e-6,
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
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_neutrino.csv",
					variation_range
				);
				sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_antineutrino.csv",
					variation_range
				);
				anti_sidis.lepton_pair_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_antineutrino.csv",
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
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, D0_grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_neutrino.csv"
				);
				sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, D0_grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, D0_grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "nutev_antineutrino.csv"
				);
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf, D0_grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_antineutrino.csv"
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
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf_in, ff_in,
					renormalization, pdf_scale, ff_scale,
					out + "nutev_neutrino.csv"
				);
				sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf_in, ff_in,
					renormalization, pdf_scale, ff_scale,
					out + "ccfr_neutrino.csv"
				);
			});

			measure([&] {
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf_in, ff_in,
					renormalization, pdf_scale, ff_scale,
					out + "nutev_antineutrino.csv"
				);
				anti_sidis.lepton_pair_xy(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					PerturbativeOrder::NLO, false, charm_mass, 0.0,
					pdf_in, ff_in,
					renormalization, pdf_scale, ff_scale,
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

			run_analysis(down_pdf, charm_fragmentation, out + "d2c/");
			run_analysis(strange_pdf, charm_fragmentation, out + "s2c/");
			run_analysis(gluon_pdf, charm_fragmentation, out + "g2c/");

			run_analysis(anti_down_pdf, anti_charm_fragmentation, out + "dbar2cbar/");
			run_analysis(anti_strange_pdf, anti_charm_fragmentation, out + "sbar2cbar/");
			run_analysis(gluon_pdf, anti_charm_fragmentation, out + "g2cbar/");
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
						PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, fragmentation,
						renormalization, pdf_scale, ff_scale,
						out + "nutev_neutrino.csv"
					);
					sidis.lepton_pair_xy(
						x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
						PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, fragmentation,
						renormalization, pdf_scale, ff_scale,
						out + "ccfr_neutrino.csv"
					);
				});

				measure([&] {
					anti_sidis.lepton_pair_xy(
						x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
						PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, fragmentation,
						renormalization, pdf_scale, ff_scale,
						out + "nutev_antineutrino.csv"
					);
					anti_sidis.lepton_pair_xy(
						x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
						PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, fragmentation,
						renormalization, pdf_scale, ff_scale,
						out + "ccfr_antineutrino.csv"
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
						renormalization, pdf_scale, ff_scale,
						out + "nutev_neutrino.csv"
					);
					sidis.lepton_pair_xy(
						x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
						PerturbativeOrder::NLO, false, mass, 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						renormalization, pdf_scale, ff_scale,
						out + "ccfr_neutrino.csv"
					);
				});

				measure([&] {
					anti_sidis.lepton_pair_xy(
						x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
						PerturbativeOrder::NLO, false, mass, 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						renormalization, pdf_scale, ff_scale,
						out + "nutev_antineutrino.csv"
					);
					anti_sidis.lepton_pair_xy(
						x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
						PerturbativeOrder::NLO, false, mass, 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						renormalization, pdf_scale, ff_scale,
						out + "ccfr_antineutrino.csv"
					);
				});
			}
		}

		std::cout << separator << IO::endl;
	}

	if (run("dis.differential.charm.masses")) {
		const LHASetInterface pdf("CT14nnlo_NF3");

		const std::vector<std::string> folders = {"1_3", "1_5"};

		for (const auto &[current_mass, folder] : std::views::zip(std::vector{1.3, 1.5}, folders)) {
			const std::string out = "Data/DIS/CharmProduction/Differential/Masses/" + folder + "/";

			// To fix the compiler error 'capturing a structured binding is not yet supported in OpenMP'
			const double mass = current_mass;

			measure([&] {
				charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
					mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nutev_neutrino.csv",
					variation_range
				);

				charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
					mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "ccfr_neutrino.csv",
					variation_range
				);
			});

			measure([&] {
				anti_charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
					mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "nutev_antineutrino.csv",
					variation_range
				);

				anti_charm_dis.differential_xy_errors(
					x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
					mass, 0.0, 0.0,
					pdf,
					renormalization, pdf_scale,
					out + "ccfr_antineutrino.csv",
					variation_range
				);
			});
		}
	}

	if (run("sidis.differential.decays")) {
		std::cout <<"=========================== sidis.differential.decays ===========================" << IO::endl;

		const double min_E = 5.0;

		const std::vector<DecayParametrization> &fit_set_1 = DecayParametrization::fit_set_1();
		const std::vector<DecayParametrization> &fit_set_2 = DecayParametrization::fit_set_2();

		for (const auto &pdf : pdfs) {
			std::cout << "PDF set: " << pdf.set_name << IO::endl;

			for (const auto &[fit_set, current_folder] : std::views::zip(std::vector{fit_set_1, fit_set_2}, std::vector{"FitSet1", "FitSet2"})) {
				// To fix the compiler error 'capturing a structured binding is not yet supported in OpenMP'
				const std::vector<DecayParametrization> &parametrizations = fit_set;
				const std::string folder = current_folder;

				const std::string out = "Data/SIDIS/MuonPairProduction/CharmedHadrons/Differential/Decays/" + pdf.set_name + "/" + folder + "/";

				measure([&] {
					sidis.lepton_pair_xy_decays(
						x_bins, get_y_bins(AnalysisSet::NuTeV, process), get_E_bins(AnalysisSet::NuTeV, process),
						parametrizations, PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						renormalization, pdf_scale, ff_scale,
						out + "nutev_neutrino.csv",
						variation_range
					);
					sidis.lepton_pair_xy_decays(
						x_bins, get_y_bins(AnalysisSet::CCFR, process), get_E_bins(AnalysisSet::CCFR, process),
						parametrizations, PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						renormalization, pdf_scale, ff_scale,
						out + "ccfr_neutrino.csv",
						variation_range
					);
				});

				measure([&] {
					anti_sidis.lepton_pair_xy_decays(
						x_bins, get_y_bins(AnalysisSet::NuTeV, anti_process), get_E_bins(AnalysisSet::NuTeV, anti_process),
						parametrizations, PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						renormalization, pdf_scale, ff_scale,
						out + "nutev_antineutrino.csv",
						variation_range
					);
					anti_sidis.lepton_pair_xy_decays(
						x_bins, get_y_bins(AnalysisSet::CCFR, anti_process), get_E_bins(AnalysisSet::CCFR, anti_process),
						parametrizations, PerturbativeOrder::NLO, false, charm_mass, 0.0,
						pdf, grid_fragmentation(min_E, Constants::Particles::MasslessMuon),
						renormalization, pdf_scale, ff_scale,
						out + "ccfr_antineutrino.csv",
						variation_range
					);
				});
			}
		}

		std::cout << separator << IO::endl;
	}

	if (run("utility.decay.grid")) {
		std::cout << "============================== utility.decay.grid ==============================" << IO::endl;

		const std::vector<double> zyE_bins = Math::linear_space(1.0, 300.0, 1e-1);

		GridGenerator generator(number_of_threads);

		std::vector<DecayParametrization> parametrizations;
		parametrizations.push_back(DecayParametrization::fit1());
		parametrizations.push_back(DecayParametrization::fit2());

		const std::vector<DecayParametrization> &fit_set_1 = DecayParametrization::fit_set_1();
		const std::vector<DecayParametrization> &fit_set_2 = DecayParametrization::fit_set_2();

		parametrizations.insert(parametrizations.end(), fit_set_1.begin(), fit_set_1.end());
		parametrizations.insert(parametrizations.end(), fit_set_2.begin(), fit_set_2.end());

		const std::size_t start_index = custom_variation_range ? (*variation_range).first : 0;

		if (!custom_variation_range || start_index <= 0) {
			generator.generate_decay_grids(
				"DecayGrids", {1.0, 5.0, 10.0, 100.0, 300.0}, {Constants::Particles::Muon.mass}, parametrizations,
				{
					Constants::Particles::D0, Constants::Particles::Dp, Constants::Particles::Ds, Constants::Particles::LambdaC
				}, Constants::Particles::Proton, Constants::Particles::Muon
			);
		}

		if (!custom_variation_range || start_index <= 1) {
			generator.generate_decay_grids(
				"DecayGrids", zyE_bins, {3.0, 5.0}, parametrizations,
				{
					Constants::Particles::D0, Constants::Particles::Dp, Constants::Particles::Ds, Constants::Particles::LambdaC
				}, Constants::Particles::Proton, Constants::Particles::Muon
			);
		}

		if (!custom_variation_range || start_index <= 2) {
			generator.generate_decay_grids(
				"DecayGrids", {1.0, 5.0, 10.0, 100.0, 300.0}, {Constants::Particles::MasslessMuon.mass}, parametrizations,
				{
					Constants::Particles::D0, Constants::Particles::Dp, Constants::Particles::Ds, Constants::Particles::LambdaC
				}, Constants::Particles::Proton, Constants::Particles::MasslessMuon
			);
		}
		
		if (!custom_variation_range || start_index <= 3) {
			generator.generate_decay_grids(
				"DecayGrids", zyE_bins, {3.0, 5.0}, parametrizations,
				{
					Constants::Particles::D0, Constants::Particles::Dp, Constants::Particles::Ds, Constants::Particles::LambdaC
				}, Constants::Particles::Proton, Constants::Particles::MasslessMuon
			);
		}

		std::cout << separator << IO::endl;
	}

	return 0;
}

#endif