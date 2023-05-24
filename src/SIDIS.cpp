#ifndef SIDIS_H
#define SIDIS_H

#include "SIDISComputation.cpp"
#include "Flavor.cpp"
#include "Process.cpp"
#include "TRFKinematics.cpp"
#include "DecayFunctions.cpp"

template <typename PDFInterface, typename FFInterface>
struct SIDIS {
	const FlavorVector active_flavors;
	const FlavorVector active_antiflavors;

	PDFInterface pdf;
	FFInterface ff;

	bool parallelize = true;
	int number_of_threads = 8;

	const size_t points;
	double max_chi_squared_deviation = 0.2;
	double max_relative_error = 1e-5;
	unsigned int iter_max = 10;

	const Process process;

	double global_sqrt_s;

	SIDIS (const FlavorVector _active_flavors, const PDFInterface _pdf, const FFInterface _ff, const size_t _points, const Process _process)
	: active_flavors(_active_flavors), 
	pdf(_pdf),
	ff(_ff),
	points(_points),
	process(_process)
	{ }

	PerturbativeResult compute_structure_function(StructureFunction F, const double x, const double z, const double Q2) {
		SIDISComputation sidis(global_sqrt_s, active_flavors, pdf, ff, points, max_chi_squared_deviation, max_relative_error, iter_max, process);
		return sidis.structure_function(F, x, z, Q2);
	}
	PerturbativeResult F2(const double x, const double z, const double Q2) {
		return compute_structure_function(StructureFunction::F2, x, z, Q2);
	}
	PerturbativeResult FL(const double x, const double z, const double Q2) {
		return compute_structure_function(StructureFunction::FL, x, z, Q2);
	}
	PerturbativeResult xF3(const double x, const double z, const double Q2) {
		return compute_structure_function(StructureFunction::xF3, x, z, Q2);
	}

	PerturbativeResult differential_cross_section(const double x, const double z, const double Q2, const bool use_direct = false) {
		SIDISComputation sidis(global_sqrt_s, active_flavors, pdf, ff, points, max_chi_squared_deviation, max_relative_error, iter_max, process);
		const PerturbativeResult differential_cs = use_direct ? sidis.differential_cross_section_direct(x, z, Q2) : sidis.differential_cross_section_indirect(x, z, Q2);
		return differential_cs;
	}

	template <typename DecayFunction>
	PerturbativeResult lepton_pair_cross_section(const TRFKinematics kinematics, const DecayParametrization parametrization, const DecayFunction decay_function) {
		SIDISComputation sidis(global_sqrt_s, active_flavors, pdf, ff, points, max_chi_squared_deviation, max_relative_error, iter_max, process);
		Decay decay(parametrization, decay_function);
		const PerturbativeResult cs = sidis.lepton_pair_cross_section(kinematics, decay);
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
	void lepton_pair_cross_section(const std::vector<double> x_bins, const std::vector<double> y_bins, const std::vector<double> E_beam_bins, const DecayParametrization parametrization, const std::string filename) {
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
					
					TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target_mass, process.projectile_mass);
					const PerturbativeResult cs_xQ2 = lepton_pair_cross_section(kinematics, parametrization, DecayFunctions::decay_function);
					const double jacobian = (kinematics.s - std::pow(process.target_mass, 2) - std::pow(process.projectile_mass, 2)) * x;
					const PerturbativeResult cs_xy = cs_xQ2 * jacobian;

					#pragma omp critical
					{
						file << x << ", " << y << ", " << E_beam << ", " << cs_xy.lo << ", " << cs_xy.nlo << std::endl;
						file.flush();

						calculated_values++;
						std::cout << "Calculated value " << calculated_values << " / " << x_step_count * y_step_count * E_beam_step_count << " (x = " << x << ", y = " << y << ", s = " << kinematics.s << ", E_beam = " << E_beam << ", Q2 = " << kinematics.Q2 << ")" << std::endl;
					}
				}
			}
		}

		file.close();
	}

};

#endif