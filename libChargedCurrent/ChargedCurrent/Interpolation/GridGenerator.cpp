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

#include "Integration/Integrator.cpp"

struct GridGenerator {
	const bool parallelize;
	const unsigned int number_of_threads;

	GridGenerator(const unsigned int number_of_threads = Utility::get_default_thread_count() / 2) : parallelize(number_of_threads > 1), number_of_threads(number_of_threads) {}

	void generate_decay_grid(
		const std::filesystem::path grid_folder, const std::vector<double> &z_bins, const std::vector<double> &yE_bins,
		const double E_min, const DecayParametrization &parametrization, const Particle &resonance, const Particle &target, const Particle &lepton
	) {
		const std::string filename = std::format(
			"dg_E_{:.2f}_r_{:.2f}_{:.2f}_t_{:.2f}_l_{:.2f}_p_{:.2f}_{:.2f}_{:.2f}_{:.2f}.dat",
			E_min, resonance.mass, resonance.lifetime, target.mass, lepton.mass, parametrization.N, parametrization.alpha, parametrization.beta, parametrization.gamma
		);

		const std::filesystem::path grid_path = grid_folder / filename;

		IO::create_directory_tree(grid_path);
		std::ofstream file(grid_path);

		file << "E_min = " << E_min << ", resonance = " << resonance << ", target = " << target << ", lepton = " << lepton;
		file << ", N = " << parametrization.N << ", alpha = " << parametrization.alpha << ", beta = " << parametrization.beta << ", gamma = " << parametrization.gamma;
		file << IO::endl;

		for (const auto &bin : {z_bins, yE_bins}) {
			for (const double value : bin) {
				file << value << " ";
			}
			file << IO::endl;
		}
		file << "-----" << IO::endl;

		const std::size_t z_count = z_bins.size();
		const std::size_t yE_count = yE_bins.size();
		const std::size_t total_count = z_count * yE_count;

		unsigned int calculated_values = 0;

		std::streamsize original_precision = std::cout.precision();

		std::vector<double> grid(total_count);

		#pragma omp parallel if(parallelize) num_threads(number_of_threads)
		{
			#pragma omp for collapse(2) schedule(static)
			for (std::size_t i = 0; i < z_count; i++) {
				for (std::size_t j = 0; j < yE_count; j++) {
					const double z = z_bins[i];
					const double yE = yE_bins[j];

					Integrator integrator([&](double input[], size_t, void *) {
						const double rho = input[0];
						const double c = input[1];

						const double h0 = z * yE;
						if (h0 < resonance.mass) { return 0.0; }
						const double pp0 = rho * h0;
						if (pp0 < E_min) { return 0.0; }

						const double mu = lepton.mass / h0;
						const double reduced_rho = sqrt(std::pow(rho, 2) - std::pow(mu, 2));

						const double a = std::pow(h0, 2) / std::pow(resonance.mass, 2);
						const double b = h0 * sqrt(std::pow(h0, 2) - std::pow(resonance.mass, 2)) / std::pow(resonance.mass, 2);

						const double cos_min = (parametrization.gamma * a * rho - 1.0) / (parametrization.gamma * b * reduced_rho);
						const double cos = cos_min + (1.0 - cos_min) * c;

						return (1.0 - cos_min) * DecayFunctions::differential_decay_function_integrand(rho, reduced_rho, cos, a, b, h0, parametrization, resonance);
					}, {0.0, 0.0}, {1.0, 1.0});

					const auto result = integrator.integrate();
					const double value = result.value;

					#pragma omp critical
					{
						grid[i * yE_count + j] = value;

						calculated_values++;

						const int base_precision = 5;
						const int z_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(z)));
						const int yE_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(yE)));

						std::cout << std::fixed << std::setprecision(base_precision);
						std::cout << "[grid] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
						std::cout << " / " << total_count;
						std::cout << ": " << value;
						std::cout << std::setprecision(z_precision) << " (z = " << z;
						std::cout << std::setprecision(yE_precision) << ", yE = " << yE << ")";
						std::cout << "\r" << std::flush;
					}
				}
			}
		}

		for (std::size_t i = 0; i < z_count; i++) {
			for (std::size_t j = 0; j < yE_count; j++) {
				file << grid[i * yE_count + j] << IO::endl;
			}
		}

		file.close();

		std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
	}

	void generate_decay_grids(
		const std::filesystem::path grid_folder, const std::vector<double> &z_bins, const std::vector<double> &yE_bins,
		const std::vector<double> &E_min_bins, const std::vector<DecayParametrization> &parametrization_bins, 
		const std::vector<Particle> &resonance_bins, const Particle &target, const Particle &lepton
	) {
		const std::size_t member_count = E_min_bins.size() * parametrization_bins.size() * resonance_bins.size();
		std::size_t member_index = 0;

		for (std::size_t parametrization_index = 0; parametrization_index < parametrization_bins.size(); parametrization_index++) {
			const DecayParametrization &parametrization = parametrization_bins[parametrization_index];

			for (std::size_t resonance_index = 0; resonance_index < resonance_bins.size(); resonance_index++) {
				const Particle &resonance = resonance_bins[resonance_index];

				for (std::size_t Emin_index = 0; Emin_index < E_min_bins.size(); Emin_index++) {
					const double E_min = E_min_bins[Emin_index];

					member_index++;

					const std::string filename = std::format(
						"dg_E_{:.2f}_r_{:.2f}_{:.2f}_t_{:.2f}_l_{:.2f}_p_{:.2f}_{:.2f}_{:.2f}_{:.2f}.dat",
						E_min, resonance.mass, resonance.lifetime, target.mass, lepton.mass, parametrization.N, parametrization.alpha, parametrization.beta, parametrization.gamma
					);

					const std::filesystem::path grid_path = grid_folder / filename;

					IO::create_directory_tree(grid_path);
					std::ofstream file(grid_path);

					file << "E_min = " << E_min << ", resonance = " << resonance << ", target = " << target << ", lepton = " << lepton;
					file << ", N = " << parametrization.N << ", alpha = " << parametrization.alpha << ", beta = " << parametrization.beta << ", gamma = " << parametrization.gamma;
					file << IO::endl;

					for (const auto &bin : {z_bins, yE_bins}) {
						for (const double value : bin) {
							file << value << " ";
						}
						file << IO::endl;
					}
					file << "-----" << IO::endl;

					const std::size_t z_count = z_bins.size();
					const std::size_t yE_count = yE_bins.size();
					const std::size_t total_count = z_count * yE_count;

					unsigned int calculated_values = 0;

					std::streamsize original_precision = std::cout.precision();

					std::vector<double> grid(total_count);

					#pragma omp parallel if(parallelize) num_threads(number_of_threads)
					{
						#pragma omp for collapse(2) schedule(static)
						for (std::size_t i = 0; i < z_count; i++) {
							for (std::size_t j = 0; j < yE_count; j++) {
								const double z = z_bins[i];
								const double yE = yE_bins[j];

								Integrator integrator([&](double input[], size_t, void *) {
									const double rho = input[0];
									const double c = input[1];

									const double h0 = z * yE;
									if (h0 < resonance.mass) { return 0.0; }
									const double pp0 = rho * h0;
									if (pp0 < E_min) { return 0.0; }

									const double mu = lepton.mass / h0;
									const double reduced_rho = sqrt(std::pow(rho, 2) - std::pow(mu, 2));

									const double a = std::pow(h0, 2) / std::pow(resonance.mass, 2);
									const double b = h0 * sqrt(std::pow(h0, 2) - std::pow(resonance.mass, 2)) / std::pow(resonance.mass, 2);

									const double cos_min = (parametrization.gamma * a * rho - 1.0) / (parametrization.gamma * b * reduced_rho);
									const double cos = cos_min + (1.0 - cos_min) * c;

									return (1.0 - cos_min) * DecayFunctions::differential_decay_function_integrand(rho, reduced_rho, cos, a, b, h0, parametrization, resonance);
								}, {0.0, 0.0}, {1.0, 1.0}, nullptr, IntegrationMethod::CubaSuave);
								integrator.cuba.maximum_evaluations = 100'000'000;
								integrator.cuba.maximum_relative_error = 1e-4;

								const auto result = integrator.integrate();
								const double value = result.value;

								#pragma omp critical
								{
									grid[i * yE_count + j] = value;

									calculated_values++;

									const int base_precision = 5;
									const int z_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(z)));
									const int yE_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(yE)));

									std::cout << std::fixed << std::setprecision(base_precision);
									std::cout << "[grid] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
									std::cout << " / " << total_count;
									std::cout << " [" << "grid " << IO::leading_zeroes(member_index, Math::number_of_digits(member_count));
									std::cout << " / " << member_count << "]";
									std::cout << ": " << value;
									std::cout << std::setprecision(z_precision) << " (z = " << z;
									std::cout << std::setprecision(yE_precision) << ", yE = " << yE << ")";
									std::cout << "\r" << std::flush;
								}
								
							}
						}
					}

					for (std::size_t i = 0; i < z_count; i++) {
						for (std::size_t j = 0; j < yE_count; j++) {
							file << grid[i * yE_count + j] << IO::endl;
						}
					}

					file.close();

					std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
				}
			}
		}
	}
};


#endif