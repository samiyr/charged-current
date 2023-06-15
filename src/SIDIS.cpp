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

template <
	typename PDFInterface, 
	typename FFInterface, 
	typename DecayFunction = decltype(DecayFunctions::trivial), 
	typename FactorizationScaleFunction = decltype(ScaleDependence::trivial), 
	typename FragmentationScaleFunction = decltype(ScaleDependence::trivial)
>
struct SIDIS {
	const FlavorVector active_flavors;
	const FlavorVector active_antiflavors;

	PDFInterface pdf;
	FragmentationConfiguration<FFInterface, DecayFunction> ff;

	bool parallelize = true;
	size_t number_of_threads = 8;

	const size_t points;
	double max_chi_squared_deviation = 0.2;
	double max_relative_error = 1e-5;
	unsigned int iter_max = 10;

	const Process process;

	double global_sqrt_s;
	bool momentum_fraction_mass_corrections = false;

	const std::optional<FactorizationScaleFunction> factorization_scale;
	const std::optional<FragmentationScaleFunction> fragmentation_scale;

	bool maintain_order_separation = true;
	bool combine_integrals = false;
	bool compute_differential_cross_section_directly = false;

	double up_mass = 0.0;
	double down_mass = 0.0;
	double charm_mass = 0.0;
	double strange_mass = 0.0;
	double top_mass = 0.0;
	double bottom_mass = 0.0;

	SIDIS (
		const FlavorVector _active_flavors, 
		const PDFInterface _pdf, 
		const FFInterface _ff, 
		const size_t _points, 
		const Process _process,
		const std::optional<FactorizationScaleFunction> _factorization_scale = std::nullopt,
		const std::optional<FragmentationScaleFunction> _fragmentation_scale = std::nullopt)
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
		const std::optional<FactorizationScaleFunction> _factorization_scale = std::nullopt,
		const std::optional<FragmentationScaleFunction> _fragmentation_scale = std::nullopt)
	: active_flavors(_active_flavors), 
	pdf(_pdf),
	ff(_ff),
	points(_points),
	process(_process),
	factorization_scale(_factorization_scale),
	fragmentation_scale(_fragmentation_scale) {	}

	private:
	auto construct_computation() {
		SIDISComputation sidis(
			global_sqrt_s, active_flavors, 
			{
				top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
			}, 
			pdf, ff, 
			points, max_chi_squared_deviation, max_relative_error, iter_max, 
			process, 
			momentum_fraction_mass_corrections, factorization_scale, fragmentation_scale
		);
		return sidis;
	}
	void output_run_info(std::ofstream &file, const std::string comment) {
		file << "#cross_section = ds/dxdy" << std::endl;
		file << "#active_flavors = ";
		for (const FlavorType flavor : active_flavors) {
			file << flavor << " ";
		}
		file << std::endl;
		file << "#active_antiflavors = ";
		for (const FlavorType flavor : active_antiflavors) {
			file << flavor << " ";
		}
		file << std::endl;
		
		file << "#pdf = " << pdf.set_name << std::endl;
		file << "#ff = ";
		for (const auto &frag : ff.interfaces) {
			file << frag.set_name << " ";
		}
		file << std::endl;

		file << "#parallelize = " << Utility::bool_to_string(parallelize) << std::endl;
		file << "#number_of_threads = " << number_of_threads << std::endl;
		file << "#points = " << points << std::endl;
		file << "#max_chi_squared_deviation = " << max_chi_squared_deviation << std::endl;
		file << "#max_relative_error = " << max_relative_error << std::endl;
		file << "#iter_max = " << iter_max << std::endl;
		file << "#up_mass = " << up_mass << std::endl;
		file << "#down_mass = " << down_mass << std::endl;
		file << "#charm_mass = " << charm_mass << std::endl;
		file << "#strange_mass = " << strange_mass << std::endl;
		file << "#top_mass = " << top_mass << std::endl;
		file << "#bottom_mass = " << bottom_mass << std::endl;

		if (!comment.empty()) {
			file << "#comment = " << comment << std::endl;
		}
	}

	public:
	PerturbativeResult compute_structure_function(StructureFunction F, const double z, const TRFKinematics kinematics) {
		SIDISComputation sidis = construct_computation();
		return sidis.compute_structure_function(F, z, kinematics, combine_integrals);
	}
	PerturbativeResult compute_structure_function(StructureFunction F, const double x, const double z, const double Q2) {
		TRFKinematics kinematics = TRFKinematics::Q2_sqrt_s(x, Q2, global_sqrt_s, process.target.mass, process.projectile.mass);
		return compute_structure_function(F, z, kinematics);
	}

	PerturbativeResult F2(const double z, const TRFKinematics kinematics) {
		return compute_structure_function(StructureFunction::F2, z, kinematics);
	}
	PerturbativeResult FL(const double z, const TRFKinematics kinematics) {
		return compute_structure_function(StructureFunction::FL, z, kinematics);
	}
	PerturbativeResult xF3(const double z, const TRFKinematics kinematics) {
		return compute_structure_function(StructureFunction::xF3, z, kinematics);
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

	PerturbativeResult differential_cross_section_xQ2(const double z, const TRFKinematics kinematics) {
		SIDISComputation sidis = construct_computation();
		if (compute_differential_cross_section_directly) {
			return combine_integrals 
					? sidis.differential_cross_section_xQ2_direct_combined(z, kinematics)
					: sidis.differential_cross_section_xQ2_direct_separated(z, kinematics);
		} else {
			return combine_integrals
					? sidis.differential_cross_section_xQ2_indirect_combined(z, kinematics)
					: sidis.differential_cross_section_xQ2_indirect_separated(z, kinematics);
		}
	}

	PerturbativeResult differential_cross_section_xy(const double z, const TRFKinematics kinematics) {
		const PerturbativeResult xQ2 = differential_cross_section_xQ2(z, kinematics);
		const double jacobian = SIDISFunctions::xy_jacobian(kinematics, process);
		return xQ2 * jacobian;
	}

	PerturbativeResult lepton_pair_cross_section_xQ2(const TRFKinematics kinematics) {
		SIDISComputation sidis = construct_computation();
		return combine_integrals ? sidis.lepton_pair_cross_section_xQ2_combined(kinematics) : sidis.lepton_pair_cross_section_xQ2_separated(kinematics);
	}
	PerturbativeResult lepton_pair_cross_section_xy(const TRFKinematics kinematics) {
		const PerturbativeResult xQ2 = lepton_pair_cross_section_xQ2(kinematics);
		const double jacobian = SIDISFunctions::xy_jacobian(kinematics, process);
		return xQ2 * jacobian;
	}

	void differential_cross_section_xQ2(const std::vector<double> x_bins, const std::vector<double> z_bins, const std::vector<double> Q2_bins, const std::string filename) {
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
					
					const PerturbativeResult differential_cs = differential_cross_section_xQ2(x, z, Q2);

					#pragma omp critical
					{
						file << x << ", " << z << ", " << Q2 << ", " << differential_cs.lo << ", " << differential_cs.nlo << std::endl;
						file.flush();

						calculated_values++;
						std::cout << "Calculated value " << calculated_values << " / " << x_step_count * z_step_count * Q2_step_count;
						std::cout << " (x = " << x << ", z = " << z << ", Q2 = " << Q2 << ")" << std::endl;
					}
				}
			}
		}

		file.close();
	}

	void lepton_pair_cross_section_xy(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::string filename, 
		const std::string comment = "") {

		const size_t x_step_count = x_bins.size();
		const size_t y_step_count = y_bins.size();
		const size_t E_beam_step_count = E_beam_bins.size();

		int calculated_values = 0;

		std::ofstream file(filename);

		output_run_info(file, comment);

		file << "x,y,E,LO,NLO,Q2,factorization_scale,fragmentation_scale" << std::endl;

		SIDISComputation sidis = construct_computation();
		#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(sidis)
		{
			#pragma omp for collapse(3)
			for (size_t i = 0; i < x_step_count; i++) {
				for (size_t j = 0; j < E_beam_step_count; j++) {
					for (size_t k = 0; k < y_step_count; k++) {
						const double x = x_bins[i];
						const double y = y_bins[k];
						const double E_beam = E_beam_bins[j];
						
						TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
						const PerturbativeResult cross_section_xQ2 = combine_integrals 
																		? sidis.lepton_pair_cross_section_xQ2_combined(kinematics) 
																		: sidis.lepton_pair_cross_section_xQ2_separated(kinematics);
						const double jacobian = SIDISFunctions::xy_jacobian(kinematics, process);
						const PerturbativeResult cross_section_xy = cross_section_xQ2 * jacobian;

						const double Q2 = kinematics.Q2;
						const double factorization_scale = sidis.compute_factorization_scale(kinematics);
						const double fragmentation_scale = sidis.compute_fragmentation_scale(kinematics);

						#pragma omp critical
						{
							file << x << ", " << y << ", " << E_beam << ", " << cross_section_xy.lo << ", " << cross_section_xy.nlo;
							file << ", " << Q2 << ", " << factorization_scale << ", " << fragmentation_scale << std::endl;
							file.flush();

							calculated_values++;
							std::cout << "Calculated value " << calculated_values << " / " << x_step_count * y_step_count * E_beam_step_count ;
							std::cout << ": " << cross_section_xy.lo << ", " << cross_section_xy.nlo;
							std::cout << " (x = " << x << ", y = " << y << ", s = " << kinematics.s << ", E_beam = " << E_beam << ", Q2 = " << kinematics.Q2 << ")";
							std::cout << std::endl;
						}
					}
				}
			}
		}
		file.close();
	}
};

#endif