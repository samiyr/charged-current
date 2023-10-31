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

#include "mdspan/mdspan.hpp"

struct GridGenerator {

	const bool parallelize;
	const unsigned int number_of_threads;

	GridGenerator(const unsigned int number_of_threads = Utility::get_default_thread_count() / 2) : parallelize(number_of_threads > 1), number_of_threads(number_of_threads) {}

	void generate_decay_grid(
		const std::filesystem::path grid_folder, const std::vector<double> &x_bins, const std::vector<double> &z_bins, const std::vector<double> &Q2_bins,
		const double E_min, const DecayParametrization &parametrization, const Particle &resonance, const Particle &target, const Particle &lepton
	) {
		using namespace Kokkos;

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

		for (const auto &bin : {x_bins, z_bins, Q2_bins}) {
			for (const double value : bin) {
				file << value << " ";
			}
			file << IO::endl;
		}
		file << "-----" << IO::endl;

		const std::size_t x_count = x_bins.size();
		const std::size_t z_count = z_bins.size();
		const std::size_t Q2_count = Q2_bins.size();
		const std::size_t total_count = x_count * z_count * Q2_count;

		unsigned int calculated_values = 0;

		std::streamsize original_precision = std::cout.precision();

		std::vector<double> backing_vector(total_count);
		auto grid = mdspan(backing_vector.data(), x_count, z_count, Q2_count);

		#pragma omp parallel if(parallelize) num_threads(number_of_threads)
		{
			#pragma omp for collapse(3) schedule(static)
			for (std::size_t i = 0; i < x_count; i++) {
				const double x = x_bins[i];
				for (std::size_t j = 0; j < z_count; j++) {
					const double z = z_bins[j];
					for (std::size_t k = 0; k < Q2_count; k++) {
						const double Q2 = Q2_bins[k];

						Integrator integrator([&](double input[], size_t, void *) {
							const double rho = input[0];
							const double c = input[1];

							const double h0 = z * Q2 / (2.0 * x * target.mass);
							if (h0 < resonance.mass) { return 0.0; }
							const double pp0 = rho * h0;
							if (pp0 < E_min) { return 0.0; }

							const double mu = lepton.mass / h0;
							const double reduced_rho = sqrt(std::pow(rho, 2) - std::pow(mu, 2));

							const double a = std::pow(h0, 2) / std::pow(resonance.mass, 2);
							const double b = h0 * sqrt(std::pow(h0, 2) - std::pow(resonance.mass, 2)) / std::pow(resonance.mass, 2);

							const double cos_min = (parametrization.gamma * a * rho - 1.0) / (parametrization.gamma * b * reduced_rho);
							const double cos = cos_min + (1.0 - cos_min) * c;

							return (1.0 - cos_min) * DecayFunctions::differential_decay_function(cos, rho, z, x, Q2, E_min, parametrization, resonance, target, lepton);
						}, {0.0, 0.0}, {1.0, 1.0});

						const auto result = integrator.integrate();
						const double value = result.value;

						#pragma omp critical
						{
							grid[i, j, k] = value;

							calculated_values++;

							const int base_precision = 5;
							const int x_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(x)));
							const int z_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(z)));
							const int Q2_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(Q2)));

							std::cout << std::fixed << std::setprecision(base_precision);
							std::cout << "[grid] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
							std::cout << " / " << total_count;
							std::cout << ": " << value;
							std::cout << std::setprecision(x_precision) << " (x = " << x;
							std::cout << std::setprecision(z_precision) << ", z = " << z;
							std::cout << std::setprecision(Q2_precision) << ", Q2 = " << Q2 << ")";
							std::cout << "\r" << std::flush;
						}
					}
				}
			}
		}

		for (std::size_t i = 0; i < x_count; i++) {
			for (std::size_t j = 0; j < z_count; j++) {
				for (std::size_t k = 0; k < Q2_count; k++) {
					file << grid[i, j, k] << IO::endl;
				}
			}
		}

		file.close();

		std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
	}

	void generate_decay_grids(
		const std::filesystem::path grid_folder, const std::vector<double> &x_bins, const std::vector<double> &z_bins, const std::vector<double> &Q2_bins,
		const std::vector<double> &E_min_bins, const std::vector<DecayParametrization> &parametrization_bins, 
		const std::vector<Particle> &resonance_bins, const Particle &target, const Particle &lepton
	) {
		using namespace Kokkos;

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

					for (const auto &bin : {x_bins, z_bins, Q2_bins}) {
						for (const double value : bin) {
							file << value << " ";
						}
						file << IO::endl;
					}
					file << "-----" << IO::endl;

					const std::size_t x_count = x_bins.size();
					const std::size_t z_count = z_bins.size();
					const std::size_t Q2_count = Q2_bins.size();
					const std::size_t total_count = x_count * z_count * Q2_count;

					unsigned int calculated_values = 0;

					std::streamsize original_precision = std::cout.precision();

					std::vector<double> backing_vector(total_count);
					auto grid = mdspan(backing_vector.data(), x_count, z_count, Q2_count);

					#pragma omp parallel if(parallelize) num_threads(number_of_threads)
					{
						#pragma omp for collapse(3) schedule(static)
						for (std::size_t i = 0; i < x_count; i++) {
							const double x = x_bins[i];
							for (std::size_t j = 0; j < z_count; j++) {
								const double z = z_bins[j];
								for (std::size_t k = 0; k < Q2_count; k++) {
									const double Q2 = Q2_bins[k];

									Integrator integrator([&](double input[], size_t, void *) {
										const double rho = input[0];
										const double c = input[1];

										const double h0 = z * Q2 / (2.0 * x * target.mass);
										if (h0 < resonance.mass) { return 0.0; }
										const double pp0 = rho * h0;
										if (pp0 < E_min) { return 0.0; }

										const double mu = lepton.mass / h0;
										if (rho <= mu) { return 0.0; }
										const double reduced_rho = sqrt(std::pow(rho, 2) - std::pow(mu, 2));

										const double a = std::pow(h0, 2) / std::pow(resonance.mass, 2);
										const double b = h0 * sqrt(std::pow(h0, 2) - std::pow(resonance.mass, 2)) / std::pow(resonance.mass, 2);

										const double cos_min = (parametrization.gamma * a * rho - 1.0) / (parametrization.gamma * b * reduced_rho);
										const double cos = cos_min + (1.0 - cos_min) * c;

										const double differential_decay_value = DecayFunctions::differential_decay_function_integrand(rho, reduced_rho, cos, a, b, h0, parametrization, resonance);

										return (1.0 - cos_min) * differential_decay_value;
									}, {0.0, 0.0}, {1.0, 1.0}, nullptr, IntegrationMethod::CubaSuave);
									integrator.cuba.maximum_evaluations = 100'000'000;
									integrator.cuba.maximum_relative_error = 1e-4;

									const auto result = integrator.integrate();
									const double value = result.value;

									#pragma omp critical
									{
										grid[i, j, k] = value;

										calculated_values++;

										const int base_precision = 5;
										const int x_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(x)));
										const int z_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(z)));
										const int Q2_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(Q2)));

										std::cout << std::fixed << std::setprecision(base_precision);
										std::cout << "[grid] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
										std::cout << " / " << total_count;
										std::cout << " [" << "grid " << IO::leading_zeroes(member_index, Math::number_of_digits(member_count));
										std::cout << " / " << member_count << "]";
										std::cout << ": " << value;
										std::cout << std::setprecision(x_precision) << " (x = " << x;
										std::cout << std::setprecision(z_precision) << ", z = " << z;
										std::cout << std::setprecision(Q2_precision) << ", Q2 = " << Q2 << ")";
										std::cout << "\r" << std::flush;
									}
								}
							}
						}
					}

					for (std::size_t i = 0; i < x_count; i++) {
						for (std::size_t j = 0; j < z_count; j++) {
							for (std::size_t k = 0; k < Q2_count; k++) {
								file << grid[i, j, k] << IO::endl;
							}
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