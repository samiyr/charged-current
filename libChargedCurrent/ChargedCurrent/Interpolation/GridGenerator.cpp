#ifndef GRID_GENERATOR_H
#define GRID_GENERATOR_H

#include <filesystem>
#include <vector>
#include <string>
#include <format>
#include <fstream>

#include "Decay/Decay.cpp"

#include "Common/Particle.cpp"

#include "Utility/Utility.cpp"

#include "Integration/Integrator2D.cpp"

struct GridGenerator {
	const bool parallelize;
	const unsigned int number_of_threads;
	double maximum_relative_error = 1e-4;
	int maximum_evaluations = 100'000'000;

	GridGenerator(const unsigned int number_of_threads = Utility::get_default_thread_count() / 2) : parallelize(number_of_threads > 1), number_of_threads(number_of_threads) {}

	static std::string grid_filename(const double E_min, const DecayParametrization &parametrization, const Particle &resonance, const Particle &target, const Particle &lepton) {
		return std::format(
			"dg_E_{:.2f}_r_{:.2f}_{:.2f}_t_{:.2f}_l_{:.2f}_p_{:.4f}_{:.4f}_{:.4f}_{:.4f}.dat",
			E_min, resonance.mass, resonance.lifetime, target.mass, lepton.mass, parametrization.N, parametrization.alpha, parametrization.beta, parametrization.gamma
		);
	}

	// void generate_decay_grid(
	// 	const std::filesystem::path grid_folder, const std::vector<double> &zyE_bins,
	// 	const double E_min, const DecayParametrization &parametrization, const Particle &resonance, const Particle &target, const Particle &lepton
	// ) {
	// 	const std::string filename = GridGenerator::grid_filename(E_min, parametrization, resonance, target, lepton);

	// 	const std::filesystem::path grid_path = grid_folder / filename;

	// 	IO::create_directory_tree(grid_path);
	// 	std::ofstream file(grid_path);

	// 	file << "E_min = " << E_min << ", resonance = " << resonance << ", target = " << target << ", lepton = " << lepton;
	// 	file << ", N = " << parametrization.N << ", alpha = " << parametrization.alpha << ", beta = " << parametrization.beta << ", gamma = " << parametrization.gamma;
	// 	file << IO::endl;

	// 	for (const double value : zyE_bins) {
	// 		file << std::format("{:.8f}", value) << " ";
	// 	}
	// 	file << IO::endl;
	// 	file << "-----" << IO::endl;

	// 	const std::size_t zyE_count = zyE_bins.size();
	// 	const std::size_t total_count = zyE_count;

	// 	unsigned int calculated_values = 0;

	// 	std::streamsize original_precision = std::cout.precision();

	// 	std::vector<double> grid(total_count);

	// 	#pragma omp parallel if(parallelize) num_threads(number_of_threads)
	// 	{
	// 		#pragma omp for schedule(static)
	// 		for (std::size_t i = 0; i < zyE_count; i++) {
	// 			const double zyE = zyE_bins[i];

	// 			Integrator integrator([&](double input[], size_t, void *) {
	// 				const double rho = input[0];
	// 				const double c = input[1];

	// 				const double h0 = zyE;
	// 				if (h0 < resonance.mass) { return 0.0; }
	// 				const double pp0 = rho * h0;
	// 				if (pp0 < E_min) { return 0.0; }

	// 				const double mu = lepton.mass / h0;
	// 				const double reduced_rho = sqrt(std::pow(rho, 2) - std::pow(mu, 2));

	// 				const double a = std::pow(h0, 2) / std::pow(resonance.mass, 2);
	// 				const double b = h0 * sqrt(std::pow(h0, 2) - std::pow(resonance.mass, 2)) / std::pow(resonance.mass, 2);

	// 				const double cos_min = (parametrization.gamma * a * rho - 1.0) / (parametrization.gamma * b * reduced_rho);
	// 				const double cos = cos_min + (1.0 - cos_min) * c;

	// 				return (1.0 - cos_min) * DecayFunctions::differential_decay_function_integrand(rho, reduced_rho, cos, a, b, h0, parametrization, resonance);
	// 			}, {0.0, 0.0}, {1.0, 1.0});

	// 			const auto result = integrator.integrate();
	// 			const double value = result.value;

	// 			#pragma omp critical
	// 			{
	// 				grid[i] = value;

	// 				calculated_values++;

	// 				const int base_precision = 5;
	// 				const int zyE_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(zyE)));

	// 				std::cout << std::fixed << std::setprecision(base_precision);
	// 				std::cout << "[grid] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
	// 				std::cout << " / " << total_count;
	// 				std::cout << ": " << value;
	// 				std::cout << std::setprecision(zyE_precision) << " (zyE = " << zyE << ")";
	// 				std::cout << "\r" << std::flush;
	// 			}
	// 		}
	// 	}

	// 	for (std::size_t i = 0; i < zyE_count; i++) {
	// 		file << std::format("{:.8f}", grid[i]) << IO::endl;
	// 	}

	// 	file.close();

	// 	std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
	// }

	void generate_decay_grids(
		const std::filesystem::path grid_folder, const std::vector<double> &zyE_bins,
		const std::vector<double> &E_min_bins, const std::vector<DecayParametrization> &parametrization_bins, 
		const std::vector<Particle> &resonance_bins, const Particle &target, const Particle &lepton
	) const {
		const std::size_t member_count = E_min_bins.size() * parametrization_bins.size() * resonance_bins.size();
		std::size_t member_index = 0;

		for (std::size_t parametrization_index = 0; parametrization_index < parametrization_bins.size(); parametrization_index++) {
			const DecayParametrization &parametrization = parametrization_bins[parametrization_index];

			for (std::size_t resonance_index = 0; resonance_index < resonance_bins.size(); resonance_index++) {
				const Particle &resonance = resonance_bins[resonance_index];

				for (std::size_t Emin_index = 0; Emin_index < E_min_bins.size(); Emin_index++) {
					const double E_min = E_min_bins[Emin_index];
					member_index++;

					const std::string filename = GridGenerator::grid_filename(E_min, parametrization, resonance, target, lepton);

					const std::filesystem::path grid_path = grid_folder / filename;

					IO::create_directory_tree(grid_path);
					std::ofstream file(grid_path);

					file << "E_min = " << E_min << ", resonance = " << resonance << ", target = " << target << ", lepton = " << lepton;
					file << ", N = " << parametrization.N << ", alpha = " << parametrization.alpha << ", beta = " << parametrization.beta << ", gamma = " << parametrization.gamma;
					file << IO::endl;

					for (const double value : zyE_bins) {
						file << std::format("{:.8f}", value) << " ";
					}
					file << IO::endl;
					file << "-----" << IO::endl;

					const std::size_t zyE_count = zyE_bins.size();
					const std::size_t total_count = zyE_count;

					unsigned int calculated_values = 0;

					std::streamsize original_precision = std::cout.precision();

					std::vector<double> grid(total_count);

					#pragma omp parallel if(parallelize) num_threads(number_of_threads)
					{
						#pragma omp for schedule(static)
						for (std::size_t i = 0; i < zyE_count; i++) {
							const double zyE = zyE_bins[i];

							const Integrator2D integrator;

							const auto result = integrator.integrate(
								[&](const double rho, const double cos) {
									const double h0 = zyE;
									if (h0 < resonance.mass) { return 0.0; }
									const double pp0 = rho * h0;
									if (pp0 < E_min) { return 0.0; }

									const double mu = lepton.mass / h0;
									const double reduced_rho = sqrt(std::pow(rho, 2) - std::pow(mu, 2));

									const double a = std::pow(h0, 2) / std::pow(resonance.mass, 2);
									const double b = h0 * sqrt(std::pow(h0, 2) - std::pow(resonance.mass, 2)) / std::pow(resonance.mass, 2);

									return DecayFunctions::differential_decay_function_integrand(rho, reduced_rho, cos, a, b, h0, parametrization, resonance);
								},
								0.0, 1.0,
								-1.0, 1.0,
								1e-12, 0
							);
							// Integrator integrator([&](double input[], size_t, void *) {
							// 	const double rho = input[0];
							// 	const double c = input[1];

							// 	const double h0 = zyE;
							// 	if (h0 < resonance.mass) { return 0.0; }
							// 	const double pp0 = rho * h0;
							// 	if (pp0 < E_min) { return 0.0; }

							// 	const double mu = lepton.mass / h0;
							// 	const double reduced_rho = sqrt(std::pow(rho, 2) - std::pow(mu, 2));

							// 	const double a = std::pow(h0, 2) / std::pow(resonance.mass, 2);
							// 	const double b = h0 * sqrt(std::pow(h0, 2) - std::pow(resonance.mass, 2)) / std::pow(resonance.mass, 2);

							// 	const double cos_min = (parametrization.gamma * a * rho - 1.0) / (parametrization.gamma * b * reduced_rho);
							// 	const double cos = cos_min + (1.0 - cos_min) * c;

							// 	return (1.0 - cos_min) * DecayFunctions::differential_decay_function_integrand(rho, reduced_rho, cos, a, b, h0, parametrization, resonance);
							// }, {0.0, 0.0}, {1.0, 1.0}, nullptr, IntegrationMethod::CubaSuave);
							// integrator.cuba.maximum_evaluations = maximum_evaluations;
							// integrator.cuba.maximum_relative_error = maximum_relative_error;

							// const auto result = integrator.integrate();
							const double value = result.value;

							#pragma omp critical
							{
								grid[i] = value;

								calculated_values++;

								const int base_precision = 5;
								const int zyE_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(zyE)));

								std::cout << std::fixed << std::setprecision(base_precision);
								std::cout << "[grid] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
								std::cout << " / " << total_count;
								std::cout << " [" << "grid " << IO::leading_zeroes(member_index, Math::number_of_digits(member_count));
								std::cout << " / " << member_count << "]";
								std::cout << ": " << std::format("{:.8f}", value);
								std::cout << std::setprecision(zyE_precision) << " (zyE = " << zyE << ")";
								std::cout << "\r" << std::flush;
							}
						}
					}

					for (std::size_t i = 0; i < zyE_count; i++) {
						file << std::format("{:.8f}", grid[i]) << IO::endl;
					}

					file.close();

					std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
				}
			}
		}
	}
};


#endif