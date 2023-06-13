#ifndef SIDIS_H
#define SIDIS_H

#include "SIDISComputation.cpp"
#include "Flavor.cpp"
#include "Process.cpp"
#include "TRFKinematics.cpp"
#include "DecayFunctions.cpp"
#include "FragmentationConfiguration.cpp"
#include "ScaleDependence.cpp"
#include <optional>

template <typename PDFInterface, typename FFInterface, typename DecayFunction = decltype(DecayFunctions::trivial), typename ScaleFunction = decltype(ScaleDependence::trivial)>
struct SIDIS {
	const FlavorVector active_flavors;
	const FlavorVector active_antiflavors;

	PDFInterface pdf;
	FragmentationConfiguration<FFInterface, DecayFunction> ff;

	bool parallelize = true;
	int number_of_threads = 8;

	const size_t points;
	double max_chi_squared_deviation = 0.2;
	double max_relative_error = 1e-5;
	unsigned int iter_max = 10;

	const Process process;

	double global_sqrt_s;
	bool momentum_fraction_mass_corrections = false;

	const std::optional<ScaleFunction> factorization_scale;
	const std::optional<ScaleFunction> fragmentation_scale;

	SIDIS (
		const FlavorVector _active_flavors, 
		const PDFInterface _pdf, 
		const FFInterface _ff, 
		const size_t _points, 
		const Process _process,
		const std::optional<ScaleFunction> _factorization_scale = std::nullopt,
		const std::optional<ScaleFunction> _fragmentation_scale = std::nullopt)
	: active_flavors(_active_flavors),
	pdf(_pdf),
	ff({_ff}, {TrivialDecay}),
	points(_points),
	process(_process),
	factorization_scale(_factorization_scale),
	fragmentation_scale(_fragmentation_scale) { }
	SIDIS (
		const FlavorVector _active_flavors, 
		const PDFInterface _pdf, 
		const FragmentationConfiguration<FFInterface, DecayFunction> _ff, 
		const size_t _points, 
		const Process _process,
		const std::optional<ScaleFunction> _factorization_scale = std::nullopt,
		const std::optional<ScaleFunction> _fragmentation_scale = std::nullopt)
	: active_flavors(_active_flavors), 
	pdf(_pdf),
	ff(_ff),
	points(_points),
	process(_process),
	factorization_scale(_factorization_scale),
	fragmentation_scale(_fragmentation_scale) {	}

	PerturbativeResult compute_structure_function(StructureFunction F, const double x, const double z, const double Q2) {
		SIDISComputation sidis(global_sqrt_s, active_flavors, pdf, ff, points, max_chi_squared_deviation, max_relative_error, iter_max, process, momentum_fraction_mass_corrections, factorization_scale, fragmentation_scale);
		return sidis.structure_function(F, x, z, Q2);
	}
	PerturbativeResult F2(const double x, const double z, const double Q2) {
		return compute_structure_function(StructureFunction::F2, x, z, Q2);
	}
	PerturbativeResult FL(const double x, const double z, const double Q2) {
		return compute_structure_function(StructureFunction::FL, x, z, Q2);
	}
	PerturbativeResult F1(const double x, const double z, const double Q2) {
		const PerturbativeResult f2 = compute_structure_function(StructureFunction::F2, x, z, Q2);
		const PerturbativeResult fL = compute_structure_function(StructureFunction::FL, x, z, Q2);
		return (f2 - fL) / (2 * x);
	}
	PerturbativeResult xF3(const double x, const double z, const double Q2) {
		return compute_structure_function(StructureFunction::xF3, x, z, Q2);
	}

	PerturbativeResult differential_cross_section(const double x, const double z, const double Q2, const bool use_direct = false) {
		SIDISComputation sidis(global_sqrt_s, active_flavors, pdf, ff, points, max_chi_squared_deviation, max_relative_error, iter_max, process, momentum_fraction_mass_corrections, factorization_scale, fragmentation_scale);
		const PerturbativeResult differential_cs = use_direct ? sidis.differential_cross_section_direct(x, z, Q2) : sidis.differential_cross_section_indirect(x, z, Q2);
		return differential_cs;
	}

	PerturbativeResult lepton_pair_cross_section(const TRFKinematics kinematics) {
		SIDISComputation sidis(global_sqrt_s, active_flavors, pdf, ff, points, max_chi_squared_deviation, max_relative_error, iter_max, process, momentum_fraction_mass_corrections, factorization_scale, fragmentation_scale);
		const PerturbativeResult cs = sidis.lepton_pair_cross_section(kinematics);
		return cs;
	}

	void differential_cross_section(const std::vector<double> x_bins, const std::vector<double> z_bins, const std::vector<double> Q2_bins, const std::string filename) {
		const size_t x_step_count = x_bins.size();
		const size_t z_step_count = z_bins.size();
		const size_t Q2_step_count = Q2_bins.size();

		int calculated_values = 0;

		std::ofstream file(filename);

		#pragma omp parallel for if(parallelize) num_threads(number_of_threads) collapse(3)
		for (size_t i = 0; i < x_step_count; i++) {
			for (size_t j = 0; j < Q2_step_count; j++) {
				for (size_t k = 0; k < z_step_count; k++) {
					const double x = x_bins[i];
					const double z = z_bins[k];
					const double Q2 = Q2_bins[j];
					
					const PerturbativeResult differential_cs = differential_cross_section(x, z, Q2);

					#pragma omp critical
					{
						file << x << ", " << z << ", " << Q2 << ", " << differential_cs.lo << ", " << differential_cs.nlo << std::endl;
						file.flush();

						calculated_values++;
						std::cout << "Calculated value " << calculated_values << " / " << x_step_count * z_step_count * Q2_step_count << " (x = " << x << ", z = " << z << ", Q2 = " << Q2 << ")" << std::endl;
					}
				}
			}
		}

		file.close();
	}
	PerturbativeResult lepton_pair_cross_section_Q2_sqrt_s(const double x, const double Q2, const double sqrt_s) {
		TRFKinematics kinematics = TRFKinematics::Q2_sqrt_s(x, Q2, sqrt_s, process.target.mass, process.projectile.mass);
		const PerturbativeResult cs_xQ2 = lepton_pair_cross_section(kinematics);
		return cs_xQ2;
	}
	PerturbativeResult lepton_pair_cross_section_y_E(const double x, const double y, const double E_beam) {
		TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
		const PerturbativeResult cs_xQ2 = lepton_pair_cross_section(kinematics);
		const double jacobian = (kinematics.s - std::pow(process.target.mass, 2) - std::pow(process.projectile.mass, 2)) * x;
		const PerturbativeResult cs_xy = cs_xQ2 * jacobian;
		return cs_xy;

	}
	void lepton_pair_cross_section(const std::vector<double> x_bins, const std::vector<double> y_bins, const std::vector<double> E_beam_bins, const std::string filename) {
		const size_t x_step_count = x_bins.size();
		const size_t y_step_count = y_bins.size();
		const size_t E_beam_step_count = E_beam_bins.size();

		int calculated_values = 0;

		std::ofstream file(filename);

		#pragma omp parallel for if(parallelize) num_threads(number_of_threads) collapse(3)
		for (size_t i = 0; i < x_step_count; i++) {
			for (size_t j = 0; j < E_beam_step_count; j++) {
				for (size_t k = 0; k < y_step_count; k++) {
					const double x = x_bins[i];
					const double y = y_bins[k];
					const double E_beam = E_beam_bins[j];
					
					TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
					const PerturbativeResult cs_xy = lepton_pair_cross_section_y_E(x, y, E_beam);

					#pragma omp critical
					{
						file << x << ", " << y << ", " << E_beam << ", " << cs_xy.lo << ", " << cs_xy.nlo << std::endl;
						file.flush();

						calculated_values++;
						std::cout << "Calculated value " << calculated_values << " / " << x_step_count * y_step_count * E_beam_step_count << ": " << cs_xy.lo << ", " << cs_xy.nlo << " (x = " << x << ", y = " << y << ", s = " << kinematics.s << ", E_beam = " << E_beam << ", Q2 = " << kinematics.Q2 << ")" << std::endl;
					}
				}
			}
		}

		file.close();
	}

};

#endif