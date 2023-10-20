#ifndef DIS_EVALUATION_H
#define DIS_EVALUATION_H

#include <fstream>
#include <string>
#include <filesystem>
#include <chrono>
#include <format>

#include "PDF/PDFConcept.cpp"
#include "PDF/Interfaces/LHASetInterface.cpp"

#include "DIS/DISComputation.cpp"

#include "Common/ScaleDependence.cpp"

template <
	typename PDFInterface, 
	is_scale_dependence RenormalizationScale = decltype(ScaleDependence::trivial)::type,
	is_scale_dependence FactorizationScale = decltype(ScaleDependence::trivial)::type
> requires is_pdf_interface<PDFInterface> || is_instance<PDFInterface, LHASetInterface>
struct DIS {
	// Compile-time constant that tells whether the PDF type (PDFInterface) is a specialization of LHASetInterface,
	// which implies that PDF error sets are available.
	static constexpr bool has_pdf_error_sets = is_instance<PDFInterface, LHASetInterface>;

	const FlavorVector active_flavors;

	const PDFInterface pdf;

	bool parallelize = true;
	unsigned int number_of_threads = Utility::get_default_thread_count();

	IntegrationParameters integration_parameters = IntegrationParameters();

	const Process process;

	double global_sqrt_s;
	bool momentum_fraction_mass_corrections = false;

	const ScaleDependence::Function<RenormalizationScale> renormalization_scale;
	const ScaleDependence::Function<FactorizationScale> factorization_scale;

	double up_mass = 0.0;
	double down_mass = 0.0;
	double charm_mass = 0.0;
	double strange_mass = 0.0;
	double top_mass = 0.0;
	double bottom_mass = 0.0;

	bool use_modified_cross_section_prefactor = false;

	// Sets the minimum energy required for the muon coming out of the neutrino --> W + muon vertex, which introduces a maximum y value y_max = 1 - min / E_beam.
	// Enforced only in the integrated cross section.
	double primary_muon_min_energy = 0.0;

	// Sets the minimum hadronic energy (E_had = E_beam - E_primary muon), which introduces a minimum y value y_min = min / E_beam.
	// Enforced only in the integrated cross section.
	double hadronic_min_energy = 0.0;

	DIS(
		const FlavorVector _active_flavors,
		const PDFInterface _pdf,
		const Process _process,
		const ScaleDependence::Function<RenormalizationScale> _renormalization_scale = ScaleDependence::Function<RenormalizationScale>(),
		const ScaleDependence::Function<FactorizationScale> _factorization_scale = ScaleDependence::Function<FactorizationScale>()
	) : active_flavors(_active_flavors),
	pdf(_pdf),
	process(_process),
	renormalization_scale(_renormalization_scale),
	factorization_scale(_factorization_scale) { }

	private:
	auto construct_computation() const requires is_pdf_interface<PDFInterface> {
		DISComputation dis(
			global_sqrt_s, active_flavors, 
			{
				top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
			}, 
			pdf,
			integration_parameters,
			process, 
			momentum_fraction_mass_corrections, renormalization_scale, factorization_scale,
			use_modified_cross_section_prefactor,
			primary_muon_min_energy, hadronic_min_energy
		);
		return dis;
	}

	template <is_pdf_interface PDF>
	auto construct_computation(const PDF pdf_member) const requires has_pdf_error_sets {
		DISComputation sidis(
			global_sqrt_s, active_flavors, 
			{
				top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
			}, 
			pdf_member,
			integration_parameters,
			process, 
			momentum_fraction_mass_corrections, renormalization_scale, factorization_scale,
			use_modified_cross_section_prefactor,
			primary_muon_min_energy, hadronic_min_energy
		);
		return sidis;
	}

	void output_run_info(std::ofstream &file, const std::string comment) {
		file << "#timestamp = " << std::format("{:%d-%m-%Y %H:%M:%OS}", std::chrono::system_clock::now()) << IO::endl;
		file << "#cross_section = d^2s/dxdy" << IO::endl;
		file << "#active_flavors = ";
		for (const FlavorType flavor : active_flavors) {
			file << flavor << " ";
		}
		file << IO::endl;
		
		file << "#pdf = " << pdf.set_name << " [" << typeid(pdf).name() << "]" << IO::endl;
		file << "#process = " << process << IO::endl;

		file << "#parallelize = " << Conversion::bool_to_string(parallelize) << IO::endl;
		file << "#number_of_threads = " << number_of_threads << IO::endl;
		file << "#up_mass = " << up_mass << IO::endl;
		file << "#down_mass = " << down_mass << IO::endl;
		file << "#charm_mass = " << charm_mass << IO::endl;
		file << "#strange_mass = " << strange_mass << IO::endl;
		file << "#top_mass = " << top_mass << IO::endl;
		file << "#bottom_mass = " << bottom_mass << IO::endl;

		file << integration_parameters;

		if (!comment.empty()) {
			file << "#comment = " << comment << IO::endl;
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
	PerturbativeQuantity F3(const TRFKinematics kinematics) const {
		return compute_structure_function(StructureFunction::F3, kinematics);
	}

	PerturbativeQuantity F2(const double x, const double Q2) const {
		return compute_structure_function(StructureFunction::F2, x, Q2);
	}
	PerturbativeQuantity FL(const double x, const double Q2) const {
		return compute_structure_function(StructureFunction::FL, x, Q2);
	}
	PerturbativeQuantity F3(const double x, const double Q2) const {
		return compute_structure_function(StructureFunction::F3, x, Q2);
	}

	PerturbativeQuantity differential_cross_section_xQ2(const TRFKinematics kinematics) const {
		DISComputation dis = construct_computation();
		return dis.differential_cross_section_xQ2(kinematics);
	}

	PerturbativeQuantity differential_cross_section_xy(const TRFKinematics kinematics) const {
		const PerturbativeQuantity xQ2 = differential_cross_section_xQ2(kinematics);
		const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
		return xQ2 * jacobian;
	}

	PerturbativeQuantity integrated_cross_section(const double E_beam, const double Q2_min) const {
		DISComputation dis = construct_computation();

		const TRFKinematics placeholder_kinematics = TRFKinematics::y_E_beam(-1.0, -1.0, E_beam, process.target.mass, process.projectile.mass);
		const PerturbativeQuantity cross_section = dis.integrated_cross_section(placeholder_kinematics, Q2_min);

		return cross_section;
	}

	void differential_cross_section_xy(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path output, 
		const std::string comment = "") {

		const std::size_t x_step_count = x_bins.size();
		const std::size_t y_step_count = y_bins.size();
		const std::size_t E_beam_step_count = E_beam_bins.size();
		const std::size_t total_count = x_step_count * y_step_count * E_beam_step_count;

		int calculated_values = 0;

		IO::create_directory_tree(output);
		std::ofstream file(output);

		output_run_info(file, comment);

		file << "x,y,E,LO,NLO,NNLO,Q2,factorization_scale" << IO::endl;
		std::streamsize original_precision = std::cout.precision();

		DISComputation dis = construct_computation();
		#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(dis)
		{
			#pragma omp for collapse(3)
			for (std::size_t i = 0; i < x_step_count; i++) {
				for (std::size_t j = 0; j < E_beam_step_count; j++) {
					for (std::size_t k = 0; k < y_step_count; k++) {
						const double x = x_bins[i];
						const double y = y_bins[k];
						const double E_beam = E_beam_bins[j];
						
						TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
						const PerturbativeQuantity cross_section_xQ2 = dis.differential_cross_section_xQ2(kinematics);
						const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
						const PerturbativeQuantity cross_section_xy = cross_section_xQ2 * jacobian;

						const double Q2 = kinematics.Q2;
						const double factorization_scale = dis.factorization_scale_function(kinematics);

						#pragma omp critical
						{
							file << x << ", " << y << ", " << E_beam << ", " << cross_section_xy.lo << ", " << cross_section_xy.nlo;
							file << ", " << Q2 << ", " << factorization_scale << IO::endl;
							file.flush();

							calculated_values++;

							const int base_precision = 5;
							const int s_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.s)));
							const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));
							const int Q2_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.Q2)));

							std::cout << std::fixed << std::setprecision(base_precision);
							std::cout << "[DIS] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
							std::cout << " / " << total_count;
							std::cout << ": " << cross_section_xy;
							std::cout << " (x = " << x << ", y = " << y;
							std::cout << std::setprecision(s_precision) << ", s = " << kinematics.s;
							std::cout << std::setprecision(E_precision) << ", E_beam = " << E_beam;
							std::cout << std::setprecision(Q2_precision) << ", Q2 = " << kinematics.Q2 << ")";
							std::cout << "\r" << std::flush;
						}
					}
				}
			}
		}
		file.close();

		std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
	}

	void integrated_cross_section(
		const std::vector<double> E_beam_bins,
		const double Q2_min,
		const std::filesystem::path output, 
		const std::string comment = "") {

		const std::size_t E_beam_step_count = E_beam_bins.size();
		const std::size_t total_count = E_beam_step_count;

		int calculated_values = 0;

		IO::create_directory_tree(output);
		std::ofstream file(output);

		output_run_info(file, comment);

		file << "E,LO,NLO,NNLO" << IO::endl;
		std::streamsize original_precision = std::cout.precision();

		DISComputation dis = construct_computation();
		#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(dis)
		{
			#pragma omp for
			for (std::size_t i = 0; i < E_beam_step_count; i++) {
				const double E_beam = E_beam_bins[i];
				
				const TRFKinematics placeholder_kinematics = TRFKinematics::y_E_beam(-1.0, -1.0, E_beam, process.target.mass, process.projectile.mass);
				const PerturbativeQuantity cross_section = dis.integrated_cross_section(placeholder_kinematics, Q2_min);

				#pragma omp critical
				{
					file << E_beam << ", " << cross_section << IO::endl;
					file.flush();

					calculated_values++;

					const int base_precision = 5;
					const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));

					std::cout << std::fixed << std::setprecision(base_precision);
					std::cout << "[DIS] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
					std::cout << " / " << total_count;
					std::cout << ": " << cross_section;
					std::cout << std::setprecision(E_precision) << " (E_beam = " << E_beam << ")";
					std::cout << "\r" << std::flush;
				}
			}			
		}
		file.close();

		std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
	}

	void differential_cross_section_xy_error_sets(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path base_output, 
		const std::string comment = "") requires has_pdf_error_sets {

		const std::size_t member_count = pdf.size();
		const std::size_t x_step_count = x_bins.size();
		const std::size_t y_step_count = y_bins.size();
		const std::size_t E_beam_step_count = E_beam_bins.size();
		const std::size_t total_count = member_count * x_step_count * y_step_count * E_beam_step_count;

		int calculated_values = 0;

		std::streamsize original_precision = std::cout.precision();

		for (typename decltype(pdf)::size_type member_index = 0; member_index < member_count; member_index++) {
			const auto &pdf_member = pdf[member_index];
			DISComputation dis = construct_computation(pdf_member);

			const std::string path_trail = IO::leading_zeroes(pdf_member.set_member_number, 4);
			std::filesystem::path full_filename = base_output.stem();
			full_filename /= path_trail;
			full_filename.replace_extension(base_output.extension());
			std::filesystem::path output = base_output;
			output.replace_filename(full_filename);

			IO::create_directory_tree(output);
			std::ofstream file(output);

			output_run_info(file, comment);
			file << "#pdf_member = " << member_index << IO::endl;

			file << "x,y,E,LO,NLO,NNLO,Q2,renormalization_scale,factorization_scale" << IO::endl;

			#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(dis)
			{
				#pragma omp for collapse(3) schedule(guided)
				for (std::size_t i = 0; i < x_step_count; i++) {
					for (std::size_t j = 0; j < E_beam_step_count; j++) {
						for (std::size_t k = 0; k < y_step_count; k++) {
							const double x = x_bins[i];
							const double y = y_bins[k];
							const double E_beam = E_beam_bins[j];
							
							TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
							const PerturbativeQuantity cross_section_xQ2 = dis.differential_cross_section_xQ2(kinematics);
							const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
							const PerturbativeQuantity cross_section_xy = cross_section_xQ2 * jacobian;

							const double Q2 = kinematics.Q2;
							const double renormalization_scale = dis.renormalization_scale_function(kinematics);
							const double factorization_scale = dis.factorization_scale_function(kinematics);

							#pragma omp critical
							{
								file << x << ", " << y << ", " << E_beam << ", " << cross_section_xy.lo << ", " << cross_section_xy.nlo << ", " << cross_section_xy.nnlo;
								file << ", " << Q2 << ", " << renormalization_scale << ", " << factorization_scale << IO::endl;
								file.flush();

								calculated_values++;

								const int base_precision = 5;
								const int s_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.s)));
								const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));
								const int Q2_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.Q2)));

								std::cout << std::fixed << std::setprecision(base_precision);
								std::cout << "[DIS] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
								std::cout << " / " << total_count;
								std::cout << " [" << "pdf member " << IO::leading_zeroes(member_index + 1, Math::number_of_digits(member_count));
								std::cout << " / " << member_count << "]";
								std::cout << ": " << cross_section_xy;
								std::cout << " (x = " << x << ", y = " << y;
								std::cout << std::setprecision(s_precision) << ", s = " << kinematics.s;
								std::cout << std::setprecision(E_precision) << ", E_beam = " << E_beam;
								std::cout << std::setprecision(Q2_precision) << ", Q2 = " << kinematics.Q2 << ")";
								std::cout << "\r" << std::flush;
							}
						}
					}
				}
			}
			file.close();
		}
		std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
	}

	void integrated_cross_section_error_sets(
		const std::vector<double> E_beam_bins,
		const double Q2_min,
		const std::filesystem::path base_output, 
		const std::string comment = "") requires has_pdf_error_sets {

		const std::size_t member_count = pdf.size();
		const std::size_t E_beam_step_count = E_beam_bins.size();
		const std::size_t total_count = member_count * E_beam_step_count;

		int calculated_values = 0;

		std::streamsize original_precision = std::cout.precision();

		for (typename decltype(pdf)::size_type member_index = 0; member_index < member_count; member_index++) {
			const auto &pdf_member = pdf[member_index];
			DISComputation dis = construct_computation(pdf_member);

			const std::string path_trail = IO::leading_zeroes(pdf_member.set_member_number, 4);
			std::filesystem::path full_filename = base_output.stem();
			full_filename /= path_trail;
			full_filename.replace_extension(base_output.extension());
			std::filesystem::path output = base_output;
			output.replace_filename(full_filename);

			IO::create_directory_tree(output);
			std::ofstream file(output);

			output_run_info(file, comment);
			file << "#pdf_member = " << member_index << IO::endl;

			file << "E,LO,NLO,NNLO" << IO::endl;

			#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(dis)
			{
				#pragma omp for
				for (std::size_t i = 0; i < E_beam_step_count; i++) {
					const double E_beam = E_beam_bins[i];
					
					TRFKinematics placeholder_kinematics = TRFKinematics::y_E_beam(-1.0, -1.0, E_beam, process.target.mass, process.projectile.mass);
					const PerturbativeQuantity cross_section = dis.integrated_cross_section(placeholder_kinematics, Q2_min);

					#pragma omp critical
					{
						file << E_beam << ", " << cross_section << IO::endl;
						file.flush();

						calculated_values++;

						const int base_precision = 5;
						const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));

						std::cout << std::fixed << std::setprecision(base_precision);
						std::cout << "[DIS] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
						std::cout << " / " << total_count;
						std::cout << " [" << "pdf member " << IO::leading_zeroes(member_index + 1, Math::number_of_digits(member_count));
						std::cout << " / " << member_count << "]";
						std::cout << ": " << cross_section;
						std::cout << std::setprecision(E_precision) << " (E_beam = " << E_beam << ")";
						std::cout << "\r" << std::flush;
					}
				}
			}
			file.close();
		}
		std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
	}
};

#endif