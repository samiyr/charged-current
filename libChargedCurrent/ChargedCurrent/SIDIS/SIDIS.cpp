#ifndef SIDIS_H
#define SIDIS_H

#include <optional>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "SIDIS/SIDISComputation.cpp"

#include "Common/Flavor.cpp"
#include "Common/Process.cpp"
#include "Common/TRFKinematics.cpp"
#include "Common/ScaleDependence.cpp"
#include "Common/CommonFunctions.cpp"

#include "Decay/DecayFunctions.cpp"

#include "PDF/FragmentationConfiguration.cpp"
#include "PDF/PDFConcept.cpp"
#include "PDF/Interfaces/LHASetInterface.cpp"

template <
	typename PDFInterface,
	is_pdf_interface FFInterface, 
	is_decay_function DecayFunction = decltype(DecayFunctions::trivial), 
	is_scale_dependence RenormalizationScale = decltype(ScaleDependence::trivial)::type, 
	is_scale_dependence FactorizationScale = decltype(ScaleDependence::trivial)::type, 
	is_scale_dependence FragmentationScale = decltype(ScaleDependence::trivial)::type
> requires is_pdf_interface<PDFInterface> || is_instance<PDFInterface, LHASetInterface>
struct SIDIS {
	static constexpr bool has_pdf_error_sets = is_instance<PDFInterface, LHASetInterface>;

	const FlavorVector active_flavors;

	const PDFInterface pdf;
	const FragmentationConfiguration<FFInterface, DecayFunction> ff;

	// Parallelization with OpenMP over kinematical points (i.e. over different values of x, y and E). Enabled by default.
	bool parallelize = true;
	// If parallelization is enabled, controls the number of threads passed on to OpenMP. Defaults to half the number of cores on the machine.
	unsigned int number_of_threads = Utility::get_default_thread_count() / 2;

	bool use_aggressive_scaling = false;
	// Parameters passed on to all integrators.
	IntegrationParameters integration_parameters = IntegrationParameters();
	// The physical process considered, whether it's neutrino or antineutrino scattering and what are the target and projectile particles.
	const Process process;

	double global_sqrt_s;

	// A wrapped function that computes the renormalization scale in a given kinematical point.
	const ScaleDependence::Function<RenormalizationScale> renormalization_scale;
	const ScaleDependence::Function<FactorizationScale> factorization_scale;
	const ScaleDependence::Function<FragmentationScale> fragmentation_scale;

	bool maintain_order_separation = true;
	bool combine_integrals = false;
	bool compute_differential_cross_section_directly = false;

	double up_mass = 0.0;
	double down_mass = 0.0;
	double charm_mass = 0.0;
	double strange_mass = 0.0;
	double top_mass = 0.0;
	double bottom_mass = 0.0;

	bool use_modified_cross_section_prefactor = false;

	PerturbativeOrder order = PerturbativeOrder::NLO;
	bool use_nlp_nlo = false;

	ScaleVariation scale_variation = ScaleVariation::None;

	SIDIS (
		const FlavorVector _active_flavors, 
		const PDFInterface _pdf, 
		const FragmentationConfiguration<FFInterface, DecayFunction> _ff, 
		const Process _process,
		const ScaleDependence::Function<RenormalizationScale> _renormalization_scale = ScaleDependence::Function<RenormalizationScale>(),
		const ScaleDependence::Function<FactorizationScale> _factorization_scale = ScaleDependence::Function<FactorizationScale>(),
		const ScaleDependence::Function<FragmentationScale> _fragmentation_scale = ScaleDependence::Function<FragmentationScale>())
	: active_flavors(_active_flavors), 
	pdf(_pdf),
	ff(_ff),
	process(_process),
	renormalization_scale(_renormalization_scale),
	factorization_scale(_factorization_scale),
	fragmentation_scale(_fragmentation_scale) {	}

	SIDIS (
		const FlavorVector _active_flavors, 
		const PDFInterface _pdf, 
		const FFInterface _ff, 
		const Process _process,
		const ScaleDependence::Function<RenormalizationScale> _renormalization_scale = ScaleDependence::Function<RenormalizationScale>(),
		const ScaleDependence::Function<FactorizationScale> _factorization_scale = ScaleDependence::Function<FactorizationScale>(),
		const ScaleDependence::Function<FragmentationScale> _fragmentation_scale = ScaleDependence::Function<FragmentationScale>())
	: SIDIS(
		_active_flavors, 
		_pdf, 
		FragmentationConfiguration<FFInterface, DecayFunction>({_ff}), 
		_process, 
		_renormalization_scale, 
		_factorization_scale,
		_fragmentation_scale) { }

	private:
	auto construct_computation() const requires is_pdf_interface<PDFInterface> {
		SIDISComputation sidis(
			active_flavors, 
			{
				top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
			},
			pdf, ff,
			integration_parameters,
			process, 
			renormalization_scale, factorization_scale, fragmentation_scale,
			use_modified_cross_section_prefactor, order, use_nlp_nlo
		);
		return sidis;
	}
	template <is_pdf_interface PDF>
	auto construct_computation(const PDF pdf_member) const requires has_pdf_error_sets {
		SIDISComputation sidis(
			active_flavors, 
			{
				top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
			},
			pdf_member, ff,
			integration_parameters,
			process, 
			renormalization_scale, factorization_scale, fragmentation_scale,
			use_modified_cross_section_prefactor, order, use_nlp_nlo
		);
		return sidis;
	}
	auto construct_computation_scale_variation(const std::vector<double> scales) const requires is_pdf_interface<PDFInterface> {
		const auto renormalization = ScaleDependence::multiplicative(scales[0]);
		const auto factorization = ScaleDependence::multiplicative(scales[1]);
		const auto fragmentation = ScaleDependence::multiplicative(scales[2]);

		SIDISComputation sidis(
			active_flavors, 
			{
				top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
			},
			pdf, ff,
			integration_parameters,
			process, 
			renormalization, factorization, fragmentation,
			use_modified_cross_section_prefactor, order, use_nlp_nlo
		);
		return sidis;
	}
	void output_run_info(std::ofstream &file, const std::string comment) {
		file << "#cross_section = d^2s/dxdy" << IO::endl;
		file << "#active_flavors = ";
		for (const FlavorType flavor : active_flavors) {
			file << flavor << " ";
		}
		file << IO::endl;
		
		file << "#pdf = " << pdf.set_name << IO::endl;
		file << "#ff = ";
		for (const auto &frag : ff.interfaces) {
			file << frag.set_name << " ";
		}
		file << IO::endl;

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
	PerturbativeQuantity compute_structure_function(
		const StructureFunction F, const double z, const TRFKinematics kinematics
	) const requires is_pdf_interface<PDFInterface> {
		SIDISComputation sidis = construct_computation();
		return sidis.compute_structure_function(F, z, kinematics);
	}
	PerturbativeQuantity compute_structure_function(
		const StructureFunction F, const double x, const double z, const double Q2
	) const requires is_pdf_interface<PDFInterface> {
		TRFKinematics kinematics = TRFKinematics::Q2_sqrt_s(x, Q2, global_sqrt_s, process.target.mass, process.projectile.mass);
		return compute_structure_function(F, z, kinematics);
	}

	PerturbativeQuantity F2(const double z, const TRFKinematics kinematics) const requires is_pdf_interface<PDFInterface> {
		return compute_structure_function(StructureFunction::F2, z, kinematics);
	}
	PerturbativeQuantity FL(const double z, const TRFKinematics kinematics) const requires is_pdf_interface<PDFInterface> {
		return compute_structure_function(StructureFunction::FL, z, kinematics);
	}
	PerturbativeQuantity F3(const double z, const TRFKinematics kinematics) const requires is_pdf_interface<PDFInterface> {
		return compute_structure_function(StructureFunction::F3, z, kinematics);
	}

	PerturbativeQuantity F2(const double x, const double z, const double Q2) const requires is_pdf_interface<PDFInterface> {
		return compute_structure_function(StructureFunction::F2, x, z, Q2);
	}
	PerturbativeQuantity FL(const double x, const double z, const double Q2) const requires is_pdf_interface<PDFInterface> {
		return compute_structure_function(StructureFunction::FL, x, z, Q2);
	}
	PerturbativeQuantity F3(const double x, const double z, const double Q2) const requires is_pdf_interface<PDFInterface> {
		return compute_structure_function(StructureFunction::F3, x, z, Q2);
	}

	PerturbativeQuantity differential_cross_section_xQ2(const double z, const TRFKinematics kinematics) const requires is_pdf_interface<PDFInterface> {
		SIDISComputation sidis = construct_computation();
		return sidis.differential_cross_section_xQ2(z, kinematics);
	}

	PerturbativeQuantity differential_cross_section_xy(const double z, const TRFKinematics kinematics) const requires is_pdf_interface<PDFInterface> {
		const PerturbativeQuantity xQ2 = differential_cross_section_xQ2(z, kinematics);
		const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
		return xQ2 * jacobian;
	}

	PerturbativeQuantity lepton_pair_cross_section_xQ2(const TRFKinematics kinematics) const requires is_pdf_interface<PDFInterface>{
		SIDISComputation sidis = construct_computation();
		return sidis.lepton_pair_cross_section_xQ2(kinematics);
	}
	PerturbativeQuantity lepton_pair_cross_section_xy(const TRFKinematics kinematics) const requires is_pdf_interface<PDFInterface> {
		const PerturbativeQuantity xQ2 = lepton_pair_cross_section_xQ2(kinematics);
		const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
		return xQ2 * jacobian;
	}

	void differential_cross_section_xQ2(
		const std::vector<double> x_bins, 
		const std::vector<double> z_bins, 
		const std::vector<double> Q2_bins, 
		const std::filesystem::path output) requires is_pdf_interface<PDFInterface> {

		const std::size_t x_step_count = x_bins.size();
		const std::size_t z_step_count = z_bins.size();
		const std::size_t Q2_step_count = Q2_bins.size();

		int calculated_values = 0;

		IO::create_directory_tree(output);
		std::ofstream file(output);

		#pragma omp parallel for if(parallelize) num_threads(number_of_threads) collapse(3)
		for (std::size_t i = 0; i < x_step_count; i++) {
			for (std::size_t j = 0; j < Q2_step_count; j++) {
				for (std::size_t k = 0; k < z_step_count; k++) {
					const double x = x_bins[i];
					const double z = z_bins[k];
					const double Q2 = Q2_bins[j];
					
					const PerturbativeQuantity differential_cs = differential_cross_section_xQ2(x, z, Q2);

					#pragma omp critical
					{
						file << x << ", " << z << ", " << Q2 << ", " << differential_cs.lo << ", " << differential_cs.nlo << ", " << differential_cs.nnlo << IO::endl;
						file.flush();

						calculated_values++;
						std::cout << "Calculated value " << calculated_values << " / " << x_step_count * z_step_count * Q2_step_count;
						std::cout << " (x = " << x << ", z = " << z << ", Q2 = " << Q2 << ")" << IO::endl;
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
		const std::filesystem::path output, 
		const std::string comment = "") requires is_pdf_interface<PDFInterface> {

		if (scale_variation != ScaleVariation::None) { return lepton_pair_cross_section_xy_scale_variations(x_bins, y_bins, E_beam_bins, output, comment); }

		const std::size_t x_step_count = x_bins.size();
		const std::size_t y_step_count = y_bins.size();
		const std::size_t E_beam_step_count = E_beam_bins.size();

		int calculated_values = 0;

		IO::create_directory_tree(output);
		std::ofstream file(output);

		output_run_info(file, comment);

		file << "x,y,E,LO,NLO,NNLO,Q2,renormalization_scale,factorization_scale,fragmentation_scale" << IO::endl;

		SIDISComputation sidis = construct_computation();
		#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(sidis)
		{
			#pragma omp for collapse(3) schedule(guided)
			for (std::size_t i = 0; i < x_step_count; i++) {
				for (std::size_t j = 0; j < E_beam_step_count; j++) {
					for (std::size_t k = 0; k < y_step_count; k++) {
						const double x = x_bins[i];
						const double y = y_bins[k];
						const double E_beam = E_beam_bins[j];
						
						TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
						const PerturbativeQuantity cross_section_xQ2 = sidis.lepton_pair_cross_section_xQ2(kinematics);
						const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
						const PerturbativeQuantity cross_section_xy = cross_section_xQ2 * jacobian;

						const double Q2 = kinematics.Q2;
						const double renormalization_scale = sidis.renormalization_scale_function(kinematics);
						const double factorization_scale = sidis.factorization_scale_function(kinematics);
						const double fragmentation_scale = sidis.fragmentation_scale_function(kinematics);

						#pragma omp critical
						{
							file << x << ", " << y << ", " << E_beam << ", " << cross_section_xy.lo << ", " << cross_section_xy.nlo << ", " << cross_section_xy.nnlo;
							file << ", " << Q2 << ", " << renormalization_scale << ", " << factorization_scale << ", " << fragmentation_scale << IO::endl;
							file.flush();

							calculated_values++;
							std::cout << "Calculated value " << calculated_values << " / " << x_step_count * y_step_count * E_beam_step_count ;
							std::cout << ": " << cross_section_xy;
							std::cout << " (x = " << x << ", y = " << y << ", s = " << kinematics.s << ", E_beam = " << E_beam << ", Q2 = " << kinematics.Q2 << ")";
							std::cout << IO::endl;
						}
					}
				}
			}
		}
		file.close();
	}

	void lepton_pair_cross_section_xy_error_sets(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path base_output, 
		const std::string comment = "") requires has_pdf_error_sets {

		const std::size_t member_count = pdf.size();
		const std::size_t x_step_count = x_bins.size();
		const std::size_t y_step_count = y_bins.size();
		const std::size_t E_beam_step_count = E_beam_bins.size();

		int calculated_values = 0;

		for (typename decltype(pdf)::size_type member_index = 0; member_index < member_count; member_index++) {
			const auto &pdf_member = pdf[member_index];
			SIDISComputation sidis = construct_computation(pdf_member);

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

			file << "x,y,E,LO,NLO,NNLO,Q2,renormalization_scale,factorization_scale,fragmentation_scale" << IO::endl;

			#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(sidis)
			{
				#pragma omp for collapse(3) schedule(guided)
				for (std::size_t i = 0; i < x_step_count; i++) {
					for (std::size_t j = 0; j < E_beam_step_count; j++) {
						for (std::size_t k = 0; k < y_step_count; k++) {
							const double x = x_bins[i];
							const double y = y_bins[k];
							const double E_beam = E_beam_bins[j];
							
							TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
							const PerturbativeQuantity cross_section_xQ2 = sidis.lepton_pair_cross_section_xQ2(kinematics);
							const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
							const PerturbativeQuantity cross_section_xy = cross_section_xQ2 * jacobian;

							const double Q2 = kinematics.Q2;
							const double renormalization_scale = sidis.renormalization_scale_function(kinematics);
							const double factorization_scale = sidis.factorization_scale_function(kinematics);
							const double fragmentation_scale = sidis.fragmentation_scale_function(kinematics);

							#pragma omp critical
							{
								file << x << ", " << y << ", " << E_beam << ", " << cross_section_xy.lo << ", " << cross_section_xy.nlo << ", " << cross_section_xy.nnlo;
								file << ", " << Q2 << ", " << renormalization_scale << ", " << factorization_scale << ", " << fragmentation_scale << IO::endl;
								file.flush();

								calculated_values++;
								std::cout << "Calculated value " << calculated_values << " / " << member_count * x_step_count * y_step_count * E_beam_step_count;
								std::cout << " [" << "pdf member " << IO::leading_zeroes(member_index + 1, Math::number_of_digits(member_count));
								std::cout << " / " << member_count << "]";
								std::cout << ": " << cross_section_xy;
								std::cout << " (x = " << x << ", y = " << y << ", s = " << kinematics.s << ", E_beam = " << E_beam << ", Q2 = " << kinematics.Q2 << ")";
								std::cout << IO::endl;
							}
						}
					}
				}
			}
			file.close();
		}
	}

	private:
	void lepton_pair_cross_section_xy_scale_variations(
		const std::vector<double> x_bins, 
		const std::vector<double> y_bins, 
		const std::vector<double> E_beam_bins, 
		const std::filesystem::path base_output, 
		const std::string comment = "") requires is_pdf_interface<PDFInterface> {

		const std::size_t x_step_count = x_bins.size();
		const std::size_t y_step_count = y_bins.size();
		const std::size_t E_beam_step_count = E_beam_bins.size();

		int calculated_values = 0;

		const std::vector<double> scale_factors{0.25, 1.0, 4.0};
		std::vector<std::vector<double>> scales = Collections::tuples(scale_factors, 3);

		const double min_scale = *std::min_element(scale_factors.begin(), scale_factors.end());
		const double max_scale = *std::max_element(scale_factors.begin(), scale_factors.end());

		std::erase_if(scales, [=, this](const std::vector<double> tuple) {
			const double renormalization = tuple[0];
			const double factorization = tuple[1];
			const double fragmentation = tuple[2];

			const bool factorization_condition = min_scale <= factorization / renormalization && max_scale >= factorization / renormalization;
			const bool fragmentation_condition = min_scale <= fragmentation / renormalization && max_scale >= fragmentation / renormalization;

			switch (scale_variation) {
			case ScaleVariation::None:
				return !(renormalization == 1.0 && factorization == 1.0 && fragmentation == 1.0);			
			case ScaleVariation::All:
				return !(factorization_condition && fragmentation_condition);
			case ScaleVariation::RenormalizationFactorization:
				return !(factorization_condition && factorization == fragmentation);
			default: return true;
			}
		});

		const std::size_t scale_count = scales.size();

		if (use_aggressive_scaling) {
			std::vector<std::stringstream> output_streams(scale_count);

			#pragma omp parallel if(parallelize) num_threads(number_of_threads)
			{
				#pragma omp for collapse(4) schedule(guided)
				for (std::size_t scale_index = 0; scale_index < scale_count; scale_index++) {
					const std::vector<double> scale = scales[scale_index];
					SIDISComputation sidis = construct_computation_scale_variation(scale);

					for (std::size_t i = 0; i < x_step_count; i++) {
						for (std::size_t j = 0; j < E_beam_step_count; j++) {
							for (std::size_t k = 0; k < y_step_count; k++) {
								const double x = x_bins[i];
								const double y = y_bins[k];
								const double E_beam = E_beam_bins[j];
								
								TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
								const PerturbativeQuantity cross_section_xQ2 = sidis.lepton_pair_cross_section_xQ2(kinematics);
								const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
								const PerturbativeQuantity cross_section_xy = cross_section_xQ2 * jacobian;

								const double Q2 = kinematics.Q2;
								const double renormalization_scale = sidis.renormalization_scale_function(kinematics);
								const double factorization_scale = sidis.factorization_scale_function(kinematics);
								const double fragmentation_scale = sidis.fragmentation_scale_function(kinematics);

								#pragma omp critical
								{
									std::stringstream &ss = output_streams[scale_index];
									ss << x << ", " << y << ", " << E_beam << ", " << cross_section_xy.lo << ", " << cross_section_xy.nlo << ", " << cross_section_xy.nnlo;
									ss << ", " << Q2 << ", " << renormalization_scale << ", " << factorization_scale << ", " << fragmentation_scale << IO::endl;

									calculated_values++;
									std::cout << "Calculated value " << calculated_values << " / " << scale_count * x_step_count * y_step_count * E_beam_step_count;
									std::cout << " [" << "scale " << IO::leading_zeroes(scale_index + 1, Math::number_of_digits(scale_count)) << " / " << scale_count << "]";
									std::cout << ": " << cross_section_xy;
									std::cout << " (x = " << x << ", y = " << y << ", s = " << kinematics.s << ", E_beam = " << E_beam << ", Q2 = " << kinematics.Q2 << ")";
									std::cout << IO::endl;
								}
							}
						}
					}
				}
			}

			for (std::size_t scale_index = 0; scale_index < scale_count; scale_index++) {
				const std::vector<double> scale = scales[scale_index];
				const std::string path_trail = (scale[0] == 1.0 && scale[1] == 1.0 && scale[2] == 1.0) ? "base_scale" : "scale_" + std::to_string(scale_index);
				std::filesystem::path full_filename = base_output.stem();
				full_filename /= path_trail;
				full_filename.replace_extension(base_output.extension());
				std::filesystem::path output = base_output;
				output.replace_filename(full_filename);

				IO::create_directory_tree(output);
				std::ofstream file(output);

				output_run_info(file, comment);
				file << "#renormalization_scale = " << scale[0] << IO::endl;
				file << "#factorization_scale = " << scale[1] << IO::endl;
				file << "#fragmentation_scale = " << scale[2] << IO::endl;

				file << "x,y,E,LO,NLO,NNLO,Q2,renormalization_scale,factorization_scale,fragmentation_scale" << IO::endl;

				const std::stringstream &stream = output_streams[scale_index];
				file << stream.rdbuf();
				file.close();
			}
		} else {
			for (std::size_t scale_index = 0; scale_index < scale_count; scale_index++) {
				const std::vector<double> scale = scales[scale_index];
				SIDISComputation sidis = construct_computation_scale_variation(scale);

				const std::string path_trail = (scale[0] == 1.0 && scale[1] == 1.0 && scale[2] == 1.0) ? "base_scale" : "scale_" + std::to_string(scale_index);
				std::filesystem::path full_filename = base_output.stem();
				full_filename /= path_trail;
				full_filename.replace_extension(base_output.extension());
				std::filesystem::path output = base_output;
				output.replace_filename(full_filename);

				IO::create_directory_tree(output);
				std::ofstream file(output);

				output_run_info(file, comment);
				file << "#renormalization_scale = " << scale[0] << IO::endl;
				file << "#factorization_scale = " << scale[1] << IO::endl;
				file << "#fragmentation_scale = " << scale[2] << IO::endl;

				file << "x,y,E,LO,NLO,NNLO,Q2,renormalization_scale,factorization_scale,fragmentation_scale" << IO::endl;

				#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(sidis)
				{
					#pragma omp for collapse(3) schedule(guided)
					for (std::size_t i = 0; i < x_step_count; i++) {
						for (std::size_t j = 0; j < E_beam_step_count; j++) {
							for (std::size_t k = 0; k < y_step_count; k++) {
								const double x = x_bins[i];
								const double y = y_bins[k];
								const double E_beam = E_beam_bins[j];
								
								TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
								const PerturbativeQuantity cross_section_xQ2 = sidis.lepton_pair_cross_section_xQ2(kinematics);
								const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
								const PerturbativeQuantity cross_section_xy = cross_section_xQ2 * jacobian;

								const double Q2 = kinematics.Q2;
								const double renormalization_scale = sidis.renormalization_scale_function(kinematics);
								const double factorization_scale = sidis.factorization_scale_function(kinematics);
								const double fragmentation_scale = sidis.fragmentation_scale_function(kinematics);

								#pragma omp critical
								{
									file << x << ", " << y << ", " << E_beam << ", " << cross_section_xy.lo << ", " << cross_section_xy.nlo << ", " << cross_section_xy.nnlo;
									file << ", " << Q2 << ", " << renormalization_scale << ", " << factorization_scale << ", " << fragmentation_scale << IO::endl;
									file.flush();

									calculated_values++;
									std::cout << "Calculated value " << calculated_values << " / " << scale_count * x_step_count * y_step_count * E_beam_step_count;
									std::cout << " [" << "scale " << IO::leading_zeroes(scale_index + 1, Math::number_of_digits(scale_count)) << " / " << scale_count << "]";
									std::cout << ": " << cross_section_xy;
									std::cout << " (x = " << x << ", y = " << y << ", s = " << kinematics.s << ", E_beam = " << E_beam << ", Q2 = " << kinematics.Q2 << ")";
									std::cout << IO::endl;
								}
							}
						}
					}
				}
				file.close();
			}
		}
	}
};

#endif