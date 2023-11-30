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

struct DIS {
	const FlavorVector active_flavors;
	const Process process;
	const unsigned int number_of_threads;

	DIS(const FlavorVector &active_flavors, const Process &process, const unsigned int number_of_threads)
	: active_flavors(active_flavors), process(process), number_of_threads(number_of_threads) {}

	IntegrationParameters integration_parameters = IntegrationParameters();

	double up_mass = 0.0;
	double down_mass = 0.0;
	double strange_mass = 0.0;
	double top_mass = 0.0;
	double bottom_mass = 0.0;

	bool use_modified_cross_section_prefactor = true;

	bool parallelize() const {
		return number_of_threads > 1;
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* 								STRUCTURE FUNCTIONS								   */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	PerturbativeQuantity compute_structure_function(
		const StructureFunction F, const TRFKinematics kinematics,
		const double charm_mass, const double primary_muon_min_energy, const double hadronic_min_energy,
		const auto &pdf,
		const auto &renormalization_scale, const auto &factorization_scale
	) const {
		const DISComputation dis = construct_computation(
			charm_mass, primary_muon_min_energy, hadronic_min_energy,
			pdf,
			renormalization_scale, factorization_scale
		);
		return dis.compute_structure_function(F, kinematics);
	}

	PerturbativeQuantity F2(
		const TRFKinematics kinematics,
		const double charm_mass, const double primary_muon_min_energy, const double hadronic_min_energy,
		const auto &pdf,
		const auto &renormalization_scale, const auto &factorization_scale
	) const {
		return compute_structure_function(
			StructureFunction::F2, kinematics,
			charm_mass, primary_muon_min_energy, hadronic_min_energy,
			pdf,
			renormalization_scale, factorization_scale
		);
	}

	PerturbativeQuantity FL(
		const TRFKinematics kinematics,
		const double charm_mass, const double primary_muon_min_energy, const double hadronic_min_energy,
		const auto &pdf,
		const auto &renormalization_scale, const auto &factorization_scale
	) const {
		return compute_structure_function(
			StructureFunction::FL, kinematics,
			charm_mass, primary_muon_min_energy, hadronic_min_energy,
			pdf,
			renormalization_scale, factorization_scale
		);
	}

	PerturbativeQuantity F3(
		const TRFKinematics kinematics,
		const double charm_mass, const double primary_muon_min_energy, const double hadronic_min_energy,
		const auto &pdf,
		const auto &renormalization_scale, const auto &factorization_scale
	) const {
		return compute_structure_function(
			StructureFunction::F3, kinematics,
			charm_mass, primary_muon_min_energy, hadronic_min_energy,
			pdf,
			renormalization_scale, factorization_scale
		);
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* 							DIFFERENTIAL CROSS SECTIONS							   */
	/* 								INDIVIDUAL VALUES							  	   */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	PerturbativeQuantity differential_xQ2(
		const TRFKinematics &kinematics,
		const double charm_mass, const double primary_muon_min_energy, const double hadronic_min_energy,
		const auto &pdf,
		const auto &renormalization_scale, const auto &factorization_scale
	) const {
		const DISComputation dis = construct_computation(
			charm_mass, primary_muon_min_energy, hadronic_min_energy,
			pdf,
			renormalization_scale, factorization_scale
		);
		return dis.differential_cross_section_xQ2(kinematics);
	}

	PerturbativeQuantity differential_xy(
		const TRFKinematics &kinematics,
		const double charm_mass, const double primary_muon_min_energy, const double hadronic_min_energy,
		const auto &pdf,
		const auto &renormalization_scale, const auto &factorization_scale
	) const {
		const PerturbativeQuantity xQ2 = differential_xQ2(
			kinematics,
			charm_mass, primary_muon_min_energy, hadronic_min_energy,
			pdf,
			renormalization_scale, factorization_scale
		);
		const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
		return xQ2 * jacobian;
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* 							DIFFERENTIAL CROSS SECTIONS							   */
	/* 								MULTIPLE VALUES 							  	   */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	private:
	void differential_xy_base(
		const auto &dis,
		const std::vector<double> &xs,
		const std::vector<double> &ys,
		const std::vector<double> &Ebeams,
		const std::optional<std::string> variation,
		const std::size_t variation_index,
		const std::size_t variation_count,
		std::ofstream &file,
		const std::string comment
	) const {
		const std::size_t x_count = xs.size();
		const std::size_t y_count = ys.size();
		const std::size_t E_count = Ebeams.size();

		const std::size_t count = x_count * y_count * E_count * variation_count;

		std::size_t calculated_values = 0;
		
		output_run_info(file, dis, comment);

		file << "x,y,E,LO,NLO,NNLO,Q2,renormalization_scale,factorization_scale" << IO::endl;
		std::streamsize original_precision = std::cout.precision();

		#pragma omp parallel if(parallelize()) num_threads(number_of_threads) firstprivate(dis)
		{
			#pragma omp for collapse(3) schedule(guided)
			for (std::size_t i = 0; i < x_count; i++) {
				for (std::size_t j = 0; j < E_count; j++) {
					for (std::size_t k = 0; k < y_count; k++) {
						const double x = xs[i];
						const double y = ys[k];
						const double E_beam = Ebeams[j];
						
						TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
						const PerturbativeQuantity cross_section_xQ2 = dis.differential_cross_section_xQ2(kinematics);
						const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
						const PerturbativeQuantity cross_section_xy = cross_section_xQ2 * jacobian;

						const double Q2 = kinematics.Q2;
						const double renormalization_scale = dis.renormalization_scale_function(kinematics);
						const double factorization_scale = dis.factorization_scale_function(kinematics);

						#pragma omp critical
						{
							file << x << ", " << y << ", " << E_beam << ", " << cross_section_xy.lo << ", " << cross_section_xy.nlo;
							file << ", " << Q2 << ", " << renormalization_scale << ", " << factorization_scale << IO::endl;
							file.flush();

							calculated_values++;

							const int base_precision = 5;
							const int s_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.s)));
							const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));
							const int Q2_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.Q2)));

							std::cout << std::fixed << std::setprecision(base_precision);
							std::cout << "[DIS] " << IO::leading_zeroes(calculated_values + x_count * y_count * E_count * variation_index, Math::number_of_digits(count));
							std::cout << " / " << count;

							if (variation) {
								std::cout << " [" << *variation << " " << IO::leading_zeroes(variation_index + 1, Math::number_of_digits(variation_count));	
								std::cout << " / " << variation_count << "]";
							}

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

		if (variation_index == variation_count - 1) {
			std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
		}
	}

	public:
	void differential_xy(
		const std::vector<double> &xs,
		const std::vector<double> &ys,
		const std::vector<double> &Ebeams,
		const double charm_mass, const double primary_muon_min_energy, const double hadronic_min_energy,
		const auto &pdf,
		const auto &renormalization_scale, const auto &factorization_scale,
		const std::filesystem::path output,
		const std::string comment = ""
	) const {
		const DISComputation dis = construct_computation(
			charm_mass, primary_muon_min_energy, hadronic_min_energy,
			pdf,
			renormalization_scale, factorization_scale
		);

		IO::create_directory_tree(output);
		std::ofstream file(output);

		differential_xy_base(dis, xs, ys, Ebeams, std::nullopt, 0, 1, file, comment);
	}

	void differential_xy_errors(
		const std::vector<double> &xs,
		const std::vector<double> &ys,
		const std::vector<double> &Ebeams,
		const double charm_mass, const double primary_muon_min_energy, const double hadronic_min_energy,
		const auto &pdf,
		const auto &renormalization_scale, const auto &factorization_scale,
		const std::filesystem::path base_output,
		const std::optional<std::pair<std::size_t, std::size_t>> variation_range = std::nullopt,
		const std::string comment = ""
	) const {
		const bool custom_range = variation_range.has_value();
		const std::size_t variation_count = custom_range ? ((*variation_range).second - (*variation_range).first + 1) : pdf.size();
		const auto variation_start = static_cast<typename std::remove_reference_t<decltype(pdf)>::size_type>(custom_range ? (*variation_range).first : 0);
		const auto variation_end = static_cast<typename std::remove_reference_t<decltype(pdf)>::size_type>(custom_range ? (*variation_range).second + 1 : variation_count);

		for (typename std::remove_reference_t<decltype(pdf)>::size_type variation_index = variation_start; variation_index < variation_end; variation_index++) {
			const auto &pdf_member = pdf[variation_index];

			const DISComputation dis = construct_computation(
				charm_mass, primary_muon_min_energy, hadronic_min_energy,
				pdf_member,
				renormalization_scale, factorization_scale
			);

			const std::string path_trail = IO::leading_zeroes(pdf_member.set_member_number, 4);
			std::filesystem::path full_filename = base_output.stem();
			full_filename /= path_trail;
			full_filename.replace_extension(base_output.extension());
			std::filesystem::path output = base_output;
			output.replace_filename(full_filename);

			IO::create_directory_tree(output);
			std::ofstream file(output);
			file << "#pdf_member = " << variation_index << IO::endl;

			differential_xy_base(dis, xs, ys, Ebeams, "pdf set", variation_index - variation_start, variation_count, file, comment);
		}
	}

	void differential_xy_scales(
		const std::vector<double> &xs,
		const std::vector<double> &ys,
		const std::vector<double> &Ebeams,
		const double charm_mass, const double primary_muon_min_energy, const double hadronic_min_energy,
		const auto &pdf,
		const double min_renormalization_scale, const double min_factorization_scale,
		const std::filesystem::path base_output,
		const std::optional<std::pair<std::size_t, std::size_t>> variation_range = std::nullopt,
		const std::string comment = ""
	) const {

		const std::vector<std::vector<double>> scales = {
			{0.25, 0.25},
			{0.25, 1.0},
			{1.0, 0.25},
			{1.0, 1.0},
			{1.0, 4.0},
			{4.0, 1.0},
			{4.0, 4.0}
		};

		const bool custom_range = variation_range.has_value();
		const std::size_t variation_count = custom_range ? ((*variation_range).second - (*variation_range).first + 1) : scales.size();
		const std::size_t variation_start = custom_range ? (*variation_range).first : 0;
		const std::size_t variation_end = custom_range ? (*variation_range).second + 1 : variation_count;

		for (std::size_t variation_index = variation_start; variation_index < variation_end; variation_index++) {
			const std::vector<double> scale = scales[variation_index];

			const auto renormalization_scale = ScaleDependence::clamped_multiplicative(scale[0], min_renormalization_scale, std::numeric_limits<double>::max());
			const auto factorization_scale = ScaleDependence::clamped_multiplicative(scale[1], min_factorization_scale, std::numeric_limits<double>::max());

			const DISComputation dis = construct_computation(
				charm_mass, primary_muon_min_energy, hadronic_min_energy,
				pdf,
				renormalization_scale, factorization_scale
			);

			const std::string path_trail = (scale[0] == 1.0 && scale[1] == 1.0) ? "base_scale" : "scale_" + std::to_string(variation_index);
			std::filesystem::path full_filename = base_output.stem();
			full_filename /= path_trail;
			full_filename.replace_extension(base_output.extension());
			std::filesystem::path output = base_output;
			output.replace_filename(full_filename);

			IO::create_directory_tree(output);
			std::ofstream file(output);

			file << "#renormalization_scale = " << scale[0] << IO::endl;
			file << "#factorization_scale = " << scale[1] << IO::endl;

			differential_xy_base(dis, xs, ys, Ebeams, "scale", variation_index - variation_start, variation_count, file, comment);
		}
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* 							INTEGRATED CROSS SECTIONS							   */
	/* 								 MULTIPLE VALUES 							  	   */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	private:

	void x_integrated_base(
		const auto &dis,
		const std::vector<double> &Ebeams,
		const std::vector<double> &Q2s,
		const std::optional<std::string> variation,
		const std::size_t variation_index,
		const std::size_t variation_count,
		std::ofstream &file,
		const std::string comment
	) const {
		const std::size_t E_count = Ebeams.size();
		const std::size_t Q2_count = Q2s.size();

		const std::size_t count = E_count * Q2_count * variation_count;

		std::size_t calculated_values = 0;
		
		output_run_info(file, dis, comment);

		file << "E,Q2,LO,NLO,NNLO" << IO::endl;
		std::streamsize original_precision = std::cout.precision();

		#pragma omp parallel if(parallelize()) num_threads(number_of_threads) firstprivate(dis)
		{
			#pragma omp for collapse(2) schedule(guided)
			for (std::size_t i = 0; i < E_count; i++) {
				for (std::size_t j = 0; j < Q2_count; j++) {
					const double E_beam = Ebeams[i];
					const double Q2 = Q2s[j];
					
					const TRFKinematics placeholder_kinematics = TRFKinematics::y_E_beam(-1.0, -1.0, E_beam, process.target.mass, process.projectile.mass);
					const PerturbativeQuantity cross_section = dis.integrated_cross_section(placeholder_kinematics, Q2);

					#pragma omp critical
					{
						file << E_beam << ", " << Q2 << ", " << cross_section << IO::endl;
						file.flush();

						calculated_values++;

						const int base_precision = 5;
						const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));
						const int Q2_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(Q2)));

						std::cout << std::fixed << std::setprecision(base_precision);
						std::cout << "[DIS] " << IO::leading_zeroes(calculated_values + E_count * Q2_count * variation_index, Math::number_of_digits(count));
						std::cout << " / " << count;

						if (variation) {
							std::cout << " [" << *variation << " " << IO::leading_zeroes(variation_index + 1, Math::number_of_digits(variation_count));
							std::cout << " / " << variation_count << "]";
						}

						std::cout << ": " << cross_section;
						std::cout << std::setprecision(E_precision) << " (E_beam = " << E_beam;
						std::cout << std::setprecision(Q2_precision) << " , Q2 = " << Q2 << ")";
						std::cout << "\r" << std::flush;
					}
				}
			}			
		}
		file.close();

		if (variation_index == variation_count - 1) {
			std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
		}
	}
	void integrated_base(
		const auto &dis,
		const std::vector<double> &Ebeams,
		const double Q2_min,
		const std::optional<std::string> variation,
		const std::size_t variation_index,
		const std::size_t variation_count,
		std::ofstream &file,
		const std::string comment
	) const {
		const std::size_t E_count = Ebeams.size();

		const std::size_t count = E_count * variation_count;

		std::size_t calculated_values = 0;
		
		output_run_info(file, dis, comment);

		file << "E,LO,NLO,NNLO" << IO::endl;
		std::streamsize original_precision = std::cout.precision();

		#pragma omp parallel if(parallelize()) num_threads(number_of_threads) firstprivate(dis)
		{
			#pragma omp for
			for (std::size_t i = 0; i < E_count; i++) {
				const double E_beam = Ebeams[i];
				
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
					std::cout << "[DIS] " << IO::leading_zeroes(calculated_values + E_count * variation_index, Math::number_of_digits(count));
					std::cout << " / " << count;

					if (variation) {
						std::cout << " [" << *variation << " " << IO::leading_zeroes(variation_index + 1, Math::number_of_digits(variation_count));
						std::cout << " / " << variation_count << "]";
					}

					std::cout << ": " << cross_section;
					std::cout << std::setprecision(E_precision) << " (E_beam = " << E_beam << ")";
					std::cout << "\r" << std::flush;
				}
			}			
		}
		file.close();

		if (variation_index == variation_count - 1) {
			std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
		}
	}

	public:
	void x_integrated(
		const std::vector<double> &Ebeams,
		const std::vector<double> &Q2s,
		const double charm_mass, const double primary_muon_min_energy, const double hadronic_min_energy,
		const auto &pdf,
		const auto &renormalization_scale, const auto &factorization_scale,
		const std::filesystem::path output,
		const std::string comment = ""
	) const {
		const DISComputation dis = construct_computation(
			charm_mass, primary_muon_min_energy, hadronic_min_energy,
			pdf,
			renormalization_scale, factorization_scale
		);

		IO::create_directory_tree(output);
		std::ofstream file(output);

		x_integrated_base(dis, Ebeams, Q2s, std::nullopt, 0, 1, file, comment);
	}
	void integrated(
		const std::vector<double> &Ebeams,
		const double Q2_min,
		const double charm_mass, const double primary_muon_min_energy, const double hadronic_min_energy,
		const auto &pdf,
		const auto &renormalization_scale, const auto &factorization_scale,
		const std::filesystem::path output,
		const std::string comment = ""
	) const {
		const DISComputation dis = construct_computation(
			charm_mass, primary_muon_min_energy, hadronic_min_energy,
			pdf,
			renormalization_scale, factorization_scale
		);

		IO::create_directory_tree(output);
		std::ofstream file(output);

		integrated_base(dis, Ebeams, Q2_min, std::nullopt, 0, 1, file, comment);
	}

	void integrated_errors(
		const std::vector<double> &Ebeams,
		const double Q2_min,
		const double charm_mass, const double primary_muon_min_energy, const double hadronic_min_energy,
		const auto &pdf,
		const auto &renormalization_scale, const auto &factorization_scale,
		const std::filesystem::path base_output,
		const std::optional<std::pair<std::size_t, std::size_t>> variation_range = std::nullopt,
		const std::string comment = ""
	) const {
		const bool custom_range = variation_range.has_value();
		const std::size_t variation_count = custom_range ? ((*variation_range).second - (*variation_range).first + 1) : pdf.size();
		const auto variation_start = static_cast<typename std::remove_reference_t<decltype(pdf)>::size_type>(custom_range ? (*variation_range).first : 0);
		const auto variation_end = static_cast<typename std::remove_reference_t<decltype(pdf)>::size_type>(custom_range ? (*variation_range).second + 1 : variation_count);

		for (typename std::remove_reference_t<decltype(pdf)>::size_type variation_index = variation_start; variation_index < variation_end; variation_index++) {
			const auto &pdf_member = pdf[variation_index];

			const DISComputation dis = construct_computation(
				charm_mass, primary_muon_min_energy, hadronic_min_energy,
				pdf_member,
				renormalization_scale, factorization_scale
			);

			const std::string path_trail = IO::leading_zeroes(pdf_member.set_member_number, 4);
			std::filesystem::path full_filename = base_output.stem();
			full_filename /= path_trail;
			full_filename.replace_extension(base_output.extension());
			std::filesystem::path output = base_output;
			output.replace_filename(full_filename);

			IO::create_directory_tree(output);
			std::ofstream file(output);
			file << "#pdf_member = " << variation_index << IO::endl;

			integrated_base(dis, Ebeams, Q2_min, "pdf set", variation_index - variation_start, variation_count, file, comment);
		}
	}

	void integrated_scales(
		const std::vector<double> &Ebeams,
		const double Q2_min,
		const double charm_mass, const double primary_muon_min_energy, const double hadronic_min_energy,
		const auto &pdf,
		const double min_renormalization_scale, const double min_factorization_scale,
		const std::filesystem::path base_output,
		const std::optional<std::pair<std::size_t, std::size_t>> variation_range = std::nullopt,
		const std::string comment = ""
	) const {

		const std::vector<std::vector<double>> scales = {
			{0.25, 0.25},
			{0.25, 1.0},
			{1.0, 0.25},
			{1.0, 1.0},
			{1.0, 4.0},
			{4.0, 1.0},
			{4.0, 4.0}
		};

		const bool custom_range = variation_range.has_value();
		const std::size_t variation_count = custom_range ? ((*variation_range).second - (*variation_range).first + 1) : scales.size();
		const std::size_t variation_start = custom_range ? (*variation_range).first : 0;
		const std::size_t variation_end = custom_range ? (*variation_range).second + 1 : variation_count;

		for (std::size_t variation_index = variation_start; variation_index < variation_end; variation_index++) {
			const std::vector<double> scale = scales[variation_index];

			const auto renormalization_scale = ScaleDependence::clamped_multiplicative(scale[0], min_renormalization_scale, std::numeric_limits<double>::max());
			const auto factorization_scale = ScaleDependence::clamped_multiplicative(scale[1], min_factorization_scale, std::numeric_limits<double>::max());

			const DISComputation dis = construct_computation(
				charm_mass, primary_muon_min_energy, hadronic_min_energy,
				pdf,
				renormalization_scale, factorization_scale
			);

			const std::string path_trail = (scale[0] == 1.0 && scale[1] == 1.0 && scale[2] == 1.0) ? "base_scale" : "scale_" + std::to_string(variation_index);
			std::filesystem::path full_filename = base_output.stem();
			full_filename /= path_trail;
			full_filename.replace_extension(base_output.extension());
			std::filesystem::path output = base_output;
			output.replace_filename(full_filename);

			IO::create_directory_tree(output);
			std::ofstream file(output);

			file << "#renormalization_scale = " << scale[0] << IO::endl;
			file << "#factorization_scale = " << scale[1] << IO::endl;

			integrated_base(dis, Ebeams, Q2_min, "scale", variation_index - variation_start, variation_count, file, comment);
		}
	}

	private:
	auto construct_computation(
		const double charm_mass, const double primary_muon_min_energy, const double hadronic_min_energy,
		const auto &pdf,
		const auto &renormalization_scale, const auto &factorization_scale) const {
		return DISComputation(
			active_flavors,
			{
				top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
			},
			pdf,
			integration_parameters,
			process,
			renormalization_scale, factorization_scale,
			use_modified_cross_section_prefactor,
			primary_muon_min_energy, hadronic_min_energy
		);
	}

	void output_run_info(std::ofstream &file, const auto &dis, const std::string comment) const {
		file << "#timestamp = " << std::format("{:%d-%m-%Y %H:%M:%OS}", std::chrono::system_clock::now()) << IO::endl;
		file << "#cross_section = d^2s/dxdy" << IO::endl;
		file << "#active_flavors = ";
		for (const FlavorType flavor : active_flavors) {
			file << flavor << " ";
		}
		file << IO::endl;
		
		file << "#pdf = " << dis.pdf1.set_name << " [" << typeid(dis.pdf1).name() << "]" << IO::endl;
		file << "#process = " << process << IO::endl;

		file << "#parallelize = " << Conversion::bool_to_string(parallelize()) << IO::endl;
		file << "#number_of_threads = " << number_of_threads << IO::endl;
		file << "#up_mass = " << up_mass << IO::endl;
		file << "#down_mass = " << down_mass << IO::endl;
		file << "#charm_mass = " << dis.flavors.flavor_masses[Flavor::Charm + 6] << IO::endl;
		file << "#strange_mass = " << strange_mass << IO::endl;
		file << "#top_mass = " << top_mass << IO::endl;
		file << "#bottom_mass = " << bottom_mass << IO::endl;

		file << integration_parameters;

		if (!comment.empty()) {
			file << "#comment = " << comment << IO::endl;
		}
	}
};

// template <
// 	typename PDFInterface, 
// 	is_scale_dependence RenormalizationScale = decltype(ScaleDependence::trivial)::type,
// 	is_scale_dependence FactorizationScale = decltype(ScaleDependence::trivial)::type
// > requires is_pdf_interface<PDFInterface> || is_instance<PDFInterface, LHASetInterface>
// struct DIS {
// 	// Compile-time constant that tells whether the PDF type (PDFInterface) is a specialization of LHASetInterface,
// 	// which implies that PDF error sets are available.
// 	static constexpr bool has_pdf_error_sets = is_instance<PDFInterface, LHASetInterface>;

// 	const FlavorVector active_flavors;

// 	const PDFInterface pdf;

// 	bool parallelize = true;
// 	unsigned int number_of_threads = Utility::get_default_thread_count();

// 	IntegrationParameters integration_parameters = IntegrationParameters();

// 	const Process process;

// 	double global_sqrt_s;
// 	bool momentum_fraction_mass_corrections = false;

// 	const ScaleDependence::Function<RenormalizationScale> renormalization_scale;
// 	const ScaleDependence::Function<FactorizationScale> factorization_scale;

// 	double up_mass = 0.0;
// 	double down_mass = 0.0;
// 	double charm_mass = 0.0;
// 	double strange_mass = 0.0;
// 	double top_mass = 0.0;
// 	double bottom_mass = 0.0;

// 	bool use_modified_cross_section_prefactor = false;

// 	// Sets the minimum energy required for the muon coming out of the neutrino --> W + muon vertex, which introduces a maximum y value y_max = 1 - min / E_beam.
// 	// Enforced only in the integrated cross section.
// 	double primary_muon_min_energy = 0.0;

// 	// Sets the minimum hadronic energy (E_had = E_beam - E_primary muon), which introduces a minimum y value y_min = min / E_beam.
// 	// Enforced only in the integrated cross section.
// 	double hadronic_min_energy = 0.0;

// 	// If enabled, factorization scale will be frozen to the smallest value determined by the PDF.
// 	// Requesting the value of the PDF at a smaller Q^2 value will then give the value with the smallest
// 	// allowed Q^2 instead. Additionally, scale logarithms will be frozen as well. In essence,
// 	// f(x, Q^2 < Q^2_min) = f(x, Q^2_min) and the same for scale logarithms.
// 	bool freeze_factorization_scale = false;

// 	bool scale_variation = false;

// 	DIS(
// 		const FlavorVector _active_flavors,
// 		const PDFInterface _pdf,
// 		const Process _process,
// 		const ScaleDependence::Function<RenormalizationScale> _renormalization_scale = ScaleDependence::Function<RenormalizationScale>(),
// 		const ScaleDependence::Function<FactorizationScale> _factorization_scale = ScaleDependence::Function<FactorizationScale>()
// 	) : active_flavors(_active_flavors),
// 	pdf(_pdf),
// 	process(_process),
// 	renormalization_scale(_renormalization_scale),
// 	factorization_scale(_factorization_scale) { }

// 	private:
// 	auto construct_computation() const requires is_pdf_interface<PDFInterface> {
// 		DISComputation dis(
// 			global_sqrt_s, active_flavors, 
// 			{
// 				top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
// 			}, 
// 			pdf,
// 			integration_parameters,
// 			process, 
// 			momentum_fraction_mass_corrections, renormalization_scale, factorization_scale,
// 			use_modified_cross_section_prefactor,
// 			primary_muon_min_energy, hadronic_min_energy
// 		);
// 		return dis;
// 	}

// 	template <is_pdf_interface PDF>
// 	auto construct_computation(const PDF pdf_member) const requires has_pdf_error_sets {
// 		DISComputation dis(
// 			global_sqrt_s, active_flavors, 
// 			{
// 				top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
// 			}, 
// 			pdf_member,
// 			integration_parameters,
// 			process, 
// 			momentum_fraction_mass_corrections, renormalization_scale, factorization_scale,
// 			use_modified_cross_section_prefactor,
// 			primary_muon_min_energy, hadronic_min_energy
// 		);
// 		return dis;
// 	}

// 	auto construct_computation_scale_variation(const std::vector<double> scales) const requires is_pdf_interface<PDFInterface> {
// 		const double pdf_Q2_min = freeze_factorization_scale ? pdf.Q2_min() : std::numeric_limits<double>::min();
// 		const double pdf_Q2_max = freeze_factorization_scale ? pdf.Q2_max() : std::numeric_limits<double>::max();

// 		const auto renormalization = ScaleDependence::multiplicative(scales[0]);
// 		const auto factorization = ScaleDependence::clamped_multiplicative(scales[1], pdf_Q2_min, pdf_Q2_max);

// 		DISComputation dis(
// 			global_sqrt_s, active_flavors, 
// 			{
// 				top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
// 			}, 
// 			pdf,
// 			integration_parameters,
// 			process, 
// 			momentum_fraction_mass_corrections, renormalization, factorization,
// 			use_modified_cross_section_prefactor,
// 			primary_muon_min_energy, hadronic_min_energy
// 		);
// 		return dis;
// 	}

// 	void output_run_info(std::ofstream &file, const std::string comment) {
// 		file << "#timestamp = " << std::format("{:%d-%m-%Y %H:%M:%OS}", std::chrono::system_clock::now()) << IO::endl;
// 		file << "#cross_section = d^2s/dxdy" << IO::endl;
// 		file << "#active_flavors = ";
// 		for (const FlavorType flavor : active_flavors) {
// 			file << flavor << " ";
// 		}
// 		file << IO::endl;
		
// 		file << "#pdf = " << pdf.set_name << " [" << typeid(pdf).name() << "]" << IO::endl;
// 		file << "#process = " << process << IO::endl;

// 		file << "#parallelize = " << Conversion::bool_to_string(parallelize) << IO::endl;
// 		file << "#number_of_threads = " << number_of_threads << IO::endl;
// 		file << "#up_mass = " << up_mass << IO::endl;
// 		file << "#down_mass = " << down_mass << IO::endl;
// 		file << "#charm_mass = " << charm_mass << IO::endl;
// 		file << "#strange_mass = " << strange_mass << IO::endl;
// 		file << "#top_mass = " << top_mass << IO::endl;
// 		file << "#bottom_mass = " << bottom_mass << IO::endl;

// 		file << integration_parameters;

// 		if (!comment.empty()) {
// 			file << "#comment = " << comment << IO::endl;
// 		}
// 	}

// 	public:
// 	PerturbativeQuantity compute_structure_function(const StructureFunction F, const TRFKinematics kinematics) const {
// 		DISComputation dis = construct_computation();
// 		return dis.compute_structure_function(F, kinematics);
// 	}
// 	PerturbativeQuantity compute_structure_function(const StructureFunction F, const double x, const double Q2) const {
// 		TRFKinematics kinematics = TRFKinematics::Q2_sqrt_s(x, Q2, global_sqrt_s, process.target.mass, process.projectile.mass);
// 		return compute_structure_function(F, kinematics);
// 	}

// 	PerturbativeQuantity F2(const TRFKinematics kinematics) const {
// 		return compute_structure_function(StructureFunction::F2, kinematics);
// 	}
// 	PerturbativeQuantity FL(const TRFKinematics kinematics) const {
// 		return compute_structure_function(StructureFunction::FL, kinematics);
// 	}
// 	PerturbativeQuantity F3(const TRFKinematics kinematics) const {
// 		return compute_structure_function(StructureFunction::F3, kinematics);
// 	}

// 	PerturbativeQuantity F2(const double x, const double Q2) const {
// 		return compute_structure_function(StructureFunction::F2, x, Q2);
// 	}
// 	PerturbativeQuantity FL(const double x, const double Q2) const {
// 		return compute_structure_function(StructureFunction::FL, x, Q2);
// 	}
// 	PerturbativeQuantity F3(const double x, const double Q2) const {
// 		return compute_structure_function(StructureFunction::F3, x, Q2);
// 	}

// 	PerturbativeQuantity differential_cross_section_xQ2(const TRFKinematics kinematics) const {
// 		DISComputation dis = construct_computation();
// 		return dis.differential_cross_section_xQ2(kinematics);
// 	}

// 	PerturbativeQuantity differential_cross_section_xy(const TRFKinematics kinematics) const {
// 		const PerturbativeQuantity xQ2 = differential_cross_section_xQ2(kinematics);
// 		const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
// 		return xQ2 * jacobian;
// 	}

// 	PerturbativeQuantity integrated_cross_section(const double E_beam, const double Q2_min) const {
// 		DISComputation dis = construct_computation();

// 		const TRFKinematics placeholder_kinematics = TRFKinematics::y_E_beam(-1.0, -1.0, E_beam, process.target.mass, process.projectile.mass);
// 		const PerturbativeQuantity cross_section = dis.integrated_cross_section(placeholder_kinematics, Q2_min);

// 		return cross_section;
// 	}

// 	void differential_cross_section_xy(
// 		const std::vector<double> x_bins, 
// 		const std::vector<double> y_bins, 
// 		const std::vector<double> E_beam_bins, 
// 		const std::filesystem::path output, 
// 		const std::string comment = "") {

// 		const std::size_t x_step_count = x_bins.size();
// 		const std::size_t y_step_count = y_bins.size();
// 		const std::size_t E_beam_step_count = E_beam_bins.size();
// 		const std::size_t total_count = x_step_count * y_step_count * E_beam_step_count;

// 		int calculated_values = 0;

// 		IO::create_directory_tree(output);
// 		std::ofstream file(output);

// 		output_run_info(file, comment);

// 		file << "x,y,E,LO,NLO,NNLO,Q2,factorization_scale" << IO::endl;
// 		std::streamsize original_precision = std::cout.precision();

// 		DISComputation dis = construct_computation();
// 		#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(dis)
// 		{
// 			#pragma omp for collapse(3)
// 			for (std::size_t i = 0; i < x_step_count; i++) {
// 				for (std::size_t j = 0; j < E_beam_step_count; j++) {
// 					for (std::size_t k = 0; k < y_step_count; k++) {
// 						const double x = x_bins[i];
// 						const double y = y_bins[k];
// 						const double E_beam = E_beam_bins[j];
						
// 						TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
// 						const PerturbativeQuantity cross_section_xQ2 = dis.differential_cross_section_xQ2(kinematics);
// 						const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
// 						const PerturbativeQuantity cross_section_xy = cross_section_xQ2 * jacobian;

// 						const double Q2 = kinematics.Q2;
// 						const double factorization_scale = dis.factorization_scale_function(kinematics);

// 						#pragma omp critical
// 						{
// 							file << x << ", " << y << ", " << E_beam << ", " << cross_section_xy.lo << ", " << cross_section_xy.nlo;
// 							file << ", " << Q2 << ", " << factorization_scale << IO::endl;
// 							file.flush();

// 							calculated_values++;

// 							const int base_precision = 5;
// 							const int s_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.s)));
// 							const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));
// 							const int Q2_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.Q2)));

// 							std::cout << std::fixed << std::setprecision(base_precision);
// 							std::cout << "[DIS] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
// 							std::cout << " / " << total_count;
// 							std::cout << ": " << cross_section_xy;
// 							std::cout << " (x = " << x << ", y = " << y;
// 							std::cout << std::setprecision(s_precision) << ", s = " << kinematics.s;
// 							std::cout << std::setprecision(E_precision) << ", E_beam = " << E_beam;
// 							std::cout << std::setprecision(Q2_precision) << ", Q2 = " << kinematics.Q2 << ")";
// 							std::cout << "\r" << std::flush;
// 						}
// 					}
// 				}
// 			}
// 		}
// 		file.close();

// 		std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
// 	}

// 	void integrated_cross_section(
// 		const std::vector<double> E_beam_bins,
// 		const double Q2_min,
// 		const std::filesystem::path output, 
// 		const std::string comment = "") {

// 		if (scale_variation) { return integrated_cross_section_scale_variations(E_beam_bins, Q2_min, output, comment); }

// 		const std::size_t E_beam_step_count = E_beam_bins.size();
// 		const std::size_t total_count = E_beam_step_count;

// 		int calculated_values = 0;

// 		IO::create_directory_tree(output);
// 		std::ofstream file(output);

// 		output_run_info(file, comment);

// 		file << "E,LO,NLO,NNLO" << IO::endl;
// 		std::streamsize original_precision = std::cout.precision();

// 		DISComputation dis = construct_computation();
// 		#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(dis)
// 		{
// 			#pragma omp for
// 			for (std::size_t i = 0; i < E_beam_step_count; i++) {
// 				const double E_beam = E_beam_bins[i];
				
// 				const TRFKinematics placeholder_kinematics = TRFKinematics::y_E_beam(-1.0, -1.0, E_beam, process.target.mass, process.projectile.mass);
// 				const PerturbativeQuantity cross_section = dis.integrated_cross_section(placeholder_kinematics, Q2_min);

// 				#pragma omp critical
// 				{
// 					file << E_beam << ", " << cross_section << IO::endl;
// 					file.flush();

// 					calculated_values++;

// 					const int base_precision = 5;
// 					const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));

// 					std::cout << std::fixed << std::setprecision(base_precision);
// 					std::cout << "[DIS] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
// 					std::cout << " / " << total_count;
// 					std::cout << ": " << cross_section;
// 					std::cout << std::setprecision(E_precision) << " (E_beam = " << E_beam << ")";
// 					std::cout << "\r" << std::flush;
// 				}
// 			}			
// 		}
// 		file.close();

// 		std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
// 	}

// 	void differential_cross_section_xy_error_sets(
// 		const std::vector<double> x_bins, 
// 		const std::vector<double> y_bins, 
// 		const std::vector<double> E_beam_bins, 
// 		const std::filesystem::path base_output, 
// 		const std::string comment = "") requires has_pdf_error_sets {

// 		const std::size_t member_count = pdf.size();
// 		const std::size_t x_step_count = x_bins.size();
// 		const std::size_t y_step_count = y_bins.size();
// 		const std::size_t E_beam_step_count = E_beam_bins.size();
// 		const std::size_t total_count = member_count * x_step_count * y_step_count * E_beam_step_count;

// 		int calculated_values = 0;

// 		std::streamsize original_precision = std::cout.precision();

// 		for (typename decltype(pdf)::size_type member_index = 0; member_index < member_count; member_index++) {
// 			const auto &pdf_member = pdf[member_index];
// 			DISComputation dis = construct_computation(pdf_member);

// 			const std::string path_trail = IO::leading_zeroes(pdf_member.set_member_number, 4);
// 			std::filesystem::path full_filename = base_output.stem();
// 			full_filename /= path_trail;
// 			full_filename.replace_extension(base_output.extension());
// 			std::filesystem::path output = base_output;
// 			output.replace_filename(full_filename);

// 			IO::create_directory_tree(output);
// 			std::ofstream file(output);

// 			output_run_info(file, comment);
// 			file << "#pdf_member = " << member_index << IO::endl;

// 			file << "x,y,E,LO,NLO,NNLO,Q2,renormalization_scale,factorization_scale" << IO::endl;

// 			#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(dis)
// 			{
// 				#pragma omp for collapse(3) schedule(guided)
// 				for (std::size_t i = 0; i < x_step_count; i++) {
// 					for (std::size_t j = 0; j < E_beam_step_count; j++) {
// 						for (std::size_t k = 0; k < y_step_count; k++) {
// 							const double x = x_bins[i];
// 							const double y = y_bins[k];
// 							const double E_beam = E_beam_bins[j];
							
// 							TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
// 							const PerturbativeQuantity cross_section_xQ2 = dis.differential_cross_section_xQ2(kinematics);
// 							const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
// 							const PerturbativeQuantity cross_section_xy = cross_section_xQ2 * jacobian;

// 							const double Q2 = kinematics.Q2;
// 							const double renormalization_scale = dis.renormalization_scale_function(kinematics);
// 							const double factorization_scale = dis.factorization_scale_function(kinematics);

// 							#pragma omp critical
// 							{
// 								file << x << ", " << y << ", " << E_beam << ", " << cross_section_xy.lo << ", " << cross_section_xy.nlo << ", " << cross_section_xy.nnlo;
// 								file << ", " << Q2 << ", " << renormalization_scale << ", " << factorization_scale << IO::endl;
// 								file.flush();

// 								calculated_values++;

// 								const int base_precision = 5;
// 								const int s_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.s)));
// 								const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));
// 								const int Q2_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.Q2)));

// 								std::cout << std::fixed << std::setprecision(base_precision);
// 								std::cout << "[DIS] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
// 								std::cout << " / " << total_count;
// 								std::cout << " [" << "pdf member " << IO::leading_zeroes(member_index + 1, Math::number_of_digits(member_count));
// 								std::cout << " / " << member_count << "]";
// 								std::cout << ": " << cross_section_xy;
// 								std::cout << " (x = " << x << ", y = " << y;
// 								std::cout << std::setprecision(s_precision) << ", s = " << kinematics.s;
// 								std::cout << std::setprecision(E_precision) << ", E_beam = " << E_beam;
// 								std::cout << std::setprecision(Q2_precision) << ", Q2 = " << kinematics.Q2 << ")";
// 								std::cout << "\r" << std::flush;
// 							}
// 						}
// 					}
// 				}
// 			}
// 			file.close();
// 		}
// 		std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
// 	}

// 	void integrated_cross_section_error_sets(
// 		const std::vector<double> E_beam_bins,
// 		const double Q2_min,
// 		const std::filesystem::path base_output, 
// 		const std::string comment = "") requires has_pdf_error_sets {

// 		const std::size_t member_count = pdf.size();
// 		const std::size_t E_beam_step_count = E_beam_bins.size();
// 		const std::size_t total_count = member_count * E_beam_step_count;

// 		int calculated_values = 0;

// 		std::streamsize original_precision = std::cout.precision();

// 		for (typename decltype(pdf)::size_type member_index = 0; member_index < member_count; member_index++) {
// 			const auto &pdf_member = pdf[member_index];
// 			DISComputation dis = construct_computation(pdf_member);

// 			const std::string path_trail = IO::leading_zeroes(pdf_member.set_member_number, 4);
// 			std::filesystem::path full_filename = base_output.stem();
// 			full_filename /= path_trail;
// 			full_filename.replace_extension(base_output.extension());
// 			std::filesystem::path output = base_output;
// 			output.replace_filename(full_filename);

// 			IO::create_directory_tree(output);
// 			std::ofstream file(output);

// 			output_run_info(file, comment);
// 			file << "#pdf_member = " << member_index << IO::endl;

// 			file << "E,LO,NLO,NNLO" << IO::endl;

// 			#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(dis)
// 			{
// 				#pragma omp for
// 				for (std::size_t i = 0; i < E_beam_step_count; i++) {
// 					const double E_beam = E_beam_bins[i];
					
// 					TRFKinematics placeholder_kinematics = TRFKinematics::y_E_beam(-1.0, -1.0, E_beam, process.target.mass, process.projectile.mass);
// 					const PerturbativeQuantity cross_section = dis.integrated_cross_section(placeholder_kinematics, Q2_min);

// 					#pragma omp critical
// 					{
// 						file << E_beam << ", " << cross_section << IO::endl;
// 						file.flush();

// 						calculated_values++;

// 						const int base_precision = 5;
// 						const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));

// 						std::cout << std::fixed << std::setprecision(base_precision);
// 						std::cout << "[DIS] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
// 						std::cout << " / " << total_count;
// 						std::cout << " [" << "pdf member " << IO::leading_zeroes(member_index + 1, Math::number_of_digits(member_count));
// 						std::cout << " / " << member_count << "]";
// 						std::cout << ": " << cross_section;
// 						std::cout << std::setprecision(E_precision) << " (E_beam = " << E_beam << ")";
// 						std::cout << "\r" << std::flush;
// 					}
// 				}
// 			}
// 			file.close();
// 		}
// 		std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
// 	}

// 	private:
// 	void integrated_cross_section_scale_variations(
// 		const std::vector<double> E_beam_bins,
// 		const double Q2_min,
// 		const std::filesystem::path base_output, 
// 		const std::string comment = "") {

// 		const std::vector<std::vector<double>> scales = {
// 			{0.25, 0.25},
// 			{0.25, 1.0},
// 			{1.0, 0.25},
// 			{1.0, 1.0},
// 			{1.0, 4.0},
// 			{4.0, 1.0},
// 			{4.0, 4.0}
// 		};

// 		const std::size_t scale_count = scales.size();

// 		const std::size_t E_beam_step_count = E_beam_bins.size();

// 		const std::size_t total_count = scale_count * E_beam_step_count;

// 		int calculated_values = 0;

// 		std::streamsize original_precision = std::cout.precision();

// 		for (std::size_t scale_index = 0; scale_index < scale_count; scale_index++) {
// 			const std::vector<double> scale = scales[scale_index];
// 			DISComputation dis = construct_computation_scale_variation(scale);

// 			const std::string path_trail = (scale[0] == 1.0 && scale[1] == 1.0) ? "base_scale" : "scale_" + std::to_string(scale_index);
// 			std::filesystem::path full_filename = base_output.stem();
// 			full_filename /= path_trail;
// 			full_filename.replace_extension(base_output.extension());
// 			std::filesystem::path output = base_output;
// 			output.replace_filename(full_filename);

// 			IO::create_directory_tree(output);
// 			std::ofstream file(output);

// 			output_run_info(file, comment);
// 			file << "#renormalization_scale = " << scale[0] << IO::endl;
// 			file << "#factorization_scale = " << scale[1] << IO::endl;

// 			file << "E,LO,NLO,NNLO" << IO::endl;

// 			#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(dis)
// 			{
// 				#pragma omp for
// 				for (std::size_t i = 0; i < E_beam_step_count; i++) {
// 					const double E_beam = E_beam_bins[i];
					
// 					TRFKinematics placeholder_kinematics = TRFKinematics::y_E_beam(-1.0, -1.0, E_beam, process.target.mass, process.projectile.mass);
// 					const PerturbativeQuantity cross_section = dis.integrated_cross_section(placeholder_kinematics, Q2_min);

// 					#pragma omp critical
// 					{
// 						file << E_beam << ", " << cross_section << IO::endl;
// 						file.flush();

// 						calculated_values++;

// 						const int base_precision = 5;
// 						const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));

// 						std::cout << std::fixed << std::setprecision(base_precision);
// 						std::cout << "[DIS] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
// 						std::cout << " / " << total_count;
// 						std::cout << " [" << "scale " << IO::leading_zeroes(scale_index + 1, Math::number_of_digits(scale_count));
// 						std::cout << " / " << scale_count << "]";
// 						std::cout << ": " << cross_section;
// 						std::cout << std::setprecision(E_precision) << " (E_beam = " << E_beam << ")";
// 						std::cout << "\r" << std::flush;
// 					}
// 				}
// 			}
// 			file.close();
// 		}
// 		std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
// 	}
// };

#endif