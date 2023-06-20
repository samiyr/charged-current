#ifndef DIS_EVALUATION_H
#define DIS_EVALUATION_H

#include "PDFConcept.cpp"
#include "Row.cpp"
#include <fstream>
#include <string>
#include "DISComputation.cpp"
#include "ScaleDependence.cpp"

template <PDFConcept PDFInterface, typename FactorizationScaleFunction = decltype(ScaleDependence::trivial)>
struct DIS {
	const FlavorVector active_flavors;

	const PDFInterface pdf;

	bool parallelize = true;
	unsigned int number_of_threads = Utility::get_default_thread_count();

	const size_t points;
	double max_chi_squared_deviation = 0.2;
	double max_relative_error = 1e-5;
	unsigned int iter_max = 10;

	const Process process;

	double global_sqrt_s;
	bool momentum_fraction_mass_corrections = false;

	const std::optional<FactorizationScaleFunction> factorization_scale;

	bool compute_differential_cross_section_directly = false;

	double up_mass = 0.0;
	double down_mass = 0.0;
	double charm_mass = 0.0;
	double strange_mass = 0.0;
	double top_mass = 0.0;
	double bottom_mass = 0.0;

	bool use_modified_cross_section_prefactor = false;

	DIS(
		const FlavorVector _active_flavors,
		const PDFInterface _pdf,
		const size_t _points,
		const Process _process,
		const std::optional<FactorizationScaleFunction> _factorization_scale = std::nullopt
	) : active_flavors(_active_flavors),
	pdf(_pdf),
	points(_points),
	process(_process),
	factorization_scale(_factorization_scale) { }

	private:
	auto construct_computation() const {
		DISComputation dis(
			global_sqrt_s, active_flavors, 
			{
				top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
			}, 
			pdf,
			points, max_chi_squared_deviation, max_relative_error, iter_max, 
			process, 
			momentum_fraction_mass_corrections, factorization_scale,
			use_modified_cross_section_prefactor
		);
		return dis;
	}
	void output_run_info(std::ofstream &file, const std::string comment) {
		file << "#cross_section = ds/dxdy" << std::endl;
		file << "#active_flavors = ";
		for (const FlavorType flavor : active_flavors) {
			file << flavor << " ";
		}
		file << std::endl;
		
		file << "#pdf = " << pdf.set_name << std::endl;

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
	PerturbativeQuantity compute_structure_function(const StructureFunction F, const TRFKinematics kinematics) const {
		DISComputation dis = construct_computation();
		return dis.compute_structure_function(F, kinematics);
	}
	PerturbativeQuantity compute_structure_function(const StructureFunction F, const double x, const double Q2) const {
		TRFKinematics kinematics = TRFKinematics::Q2_sqrt_s(x, Q2, global_sqrt_s, process.target.mass, process.projectile.mass);
		return compute_structure_function(F, kinematics);
	}

	PerturbativeQuantity F2(const TRFKinematics kinematics) const {
		return compute_structure_function(StructureFunction::F2, kinematics);
	}
	PerturbativeQuantity FL(const TRFKinematics kinematics) const {
		return compute_structure_function(StructureFunction::FL, kinematics);
	}
	PerturbativeQuantity xF3(const TRFKinematics kinematics) const {
		return compute_structure_function(StructureFunction::xF3, kinematics);
	}

	PerturbativeQuantity F2(const double x, const double Q2) const {
		return compute_structure_function(StructureFunction::F2, x, Q2);
	}
	PerturbativeQuantity FL(const double x, const double Q2) const {
		return compute_structure_function(StructureFunction::FL, x, Q2);
	}
	PerturbativeQuantity xF3(const double x, const double Q2) const {
		return compute_structure_function(StructureFunction::xF3, x, Q2);
	}

	PerturbativeQuantity differential_cross_section_xQ2(const TRFKinematics kinematics) const {
		DISComputation dis = construct_computation();
		if (compute_differential_cross_section_directly) {
			return dis.differential_cross_section_xQ2_direct(kinematics);
		} else {
			return dis.differential_cross_section_xQ2_indirect(kinematics);
		}
	}

	PerturbativeQuantity differential_cross_section_xy(const TRFKinematics kinematics) const {
		const PerturbativeQuantity xQ2 = differential_cross_section_xQ2(kinematics);
		const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
		return xQ2 * jacobian;
	}

	void differential_cross_section_xy(
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

		file << "x,y,E,LO,NLO,Q2,factorization_scale" << std::endl;

		DISComputation dis = construct_computation();
		#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(dis)
		{
			#pragma omp for collapse(3)
			for (size_t i = 0; i < x_step_count; i++) {
				for (size_t j = 0; j < E_beam_step_count; j++) {
					for (size_t k = 0; k < y_step_count; k++) {
						const double x = x_bins[i];
						const double y = y_bins[k];
						const double E_beam = E_beam_bins[j];
						
						TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
						const PerturbativeQuantity cross_section_xQ2 = compute_differential_cross_section_directly 
																		? dis.differential_cross_section_xQ2_direct(kinematics) 
																		: dis.differential_cross_section_xQ2_indirect(kinematics);
						const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
						const PerturbativeQuantity cross_section_xy = cross_section_xQ2 * jacobian;

						const double Q2 = kinematics.Q2;
						const double factorization_scale = dis.compute_factorization_scale(kinematics);

						#pragma omp critical
						{
							file << x << ", " << y << ", " << E_beam << ", " << cross_section_xy.lo << ", " << cross_section_xy.nlo;
							file << ", " << Q2 << ", " << factorization_scale << std::endl;
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