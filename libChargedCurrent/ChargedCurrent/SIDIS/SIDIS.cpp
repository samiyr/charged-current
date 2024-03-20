#ifndef SIDIS_H
#define SIDIS_H

#include <optional>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <chrono>
#include <format>

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

#include "Threading/ThreadResult.cpp"
#include "Threading/shared_multi_future.cpp"

struct SIDIS {
	const FlavorVector active_flavors;
	const Process process;
	const unsigned int number_of_threads;

	SIDIS(const FlavorVector &active_flavors, const Process &process, const unsigned int number_of_threads)
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
		const StructureFunction F, const double z, const TRFKinematics kinematics,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_interface auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale
	) const {
		const SIDISComputation sidis = construct_computation(
			order, use_nlp_nlo, 
			charm_mass, primary_muon_min_energy, 
			pdf, ff,
			renormalization_scale, factorization_scale, fragmentation_scale
		);
		return sidis.compute_structure_function(F, z, kinematics);
	}

	PerturbativeQuantity F2(
		const double z, const TRFKinematics kinematics,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_interface auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale
	) const {
		return compute_structure_function(
			StructureFunction::F2, z, kinematics,
			order, use_nlp_nlo,
			charm_mass, primary_muon_min_energy,
			pdf, ff,
			renormalization_scale, factorization_scale, fragmentation_scale
		);
	}

	PerturbativeQuantity FL(
		const double z, const TRFKinematics kinematics,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_interface auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale
	) const {
		return compute_structure_function(
			StructureFunction::FL, z, kinematics,
			order, use_nlp_nlo,
			charm_mass, primary_muon_min_energy,
			pdf, ff,
			renormalization_scale, factorization_scale, fragmentation_scale
		);
	}

	PerturbativeQuantity F3(
		const double z, const TRFKinematics kinematics,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_interface auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale
	) const {
		return compute_structure_function(
			StructureFunction::F3, z, kinematics,
			order, use_nlp_nlo,
			charm_mass, primary_muon_min_energy,
			pdf, ff,
			renormalization_scale, factorization_scale, fragmentation_scale
		);
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* 							DIFFERENTIAL CROSS SECTIONS							   */
	/* 								INDIVIDUAL VALUES							  	   */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	PerturbativeQuantity differential_xQ2(
		const double z, const TRFKinematics &kinematics,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_interface auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale
	) const {
		const SIDISComputation sidis = construct_computation(
			order, use_nlp_nlo, 
			charm_mass, primary_muon_min_energy, 
			pdf, ff,
			renormalization_scale, factorization_scale, fragmentation_scale
		);
		return sidis.differential_cross_section_xQ2(z, kinematics);
	}

	PerturbativeQuantity differential_xy(
		const double z, const TRFKinematics &kinematics,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_interface auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale
	) const {
		const PerturbativeQuantity xQ2 = differential_xQ2(
			z, kinematics,
			order, use_nlp_nlo, 
			charm_mass, primary_muon_min_energy, 
			pdf, ff,
			renormalization_scale, factorization_scale, fragmentation_scale
		);
		const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
		return xQ2 * jacobian;
	}

	PerturbativeQuantity lepton_pair_zxQ2(
		const TRFKinematics &kinematics, const double z,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_interface auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale
	) const {
		const SIDISComputation sidis = construct_computation(
			order, use_nlp_nlo, 
			charm_mass, primary_muon_min_energy, 
			pdf, ff,
			renormalization_scale, factorization_scale, fragmentation_scale
		);
		return sidis.lepton_pair_cross_section_xQ2(kinematics, z);
	}

	PerturbativeQuantity lepton_pair_zxy(
		const TRFKinematics &kinematics, const double z,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_interface auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale
	) const {
		const PerturbativeQuantity xQ2 = lepton_pair_zxQ2(
			kinematics, z,
			order, use_nlp_nlo, 
			charm_mass, primary_muon_min_energy, 
			pdf, ff,
			renormalization_scale, factorization_scale, fragmentation_scale
		);
		const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
		return xQ2 * jacobian;
	}

	PerturbativeQuantity lepton_pair_xQ2(
		const TRFKinematics &kinematics,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_interface auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale
	) const {
		const SIDISComputation sidis = construct_computation(
			order, use_nlp_nlo, 
			charm_mass, primary_muon_min_energy, 
			pdf, ff,
			renormalization_scale, factorization_scale, fragmentation_scale
		);
		return sidis.lepton_pair_cross_section_xQ2(kinematics);
	}

	PerturbativeQuantity lepton_pair_xy(
		const TRFKinematics &kinematics,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_interface auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale
	) const {
		const PerturbativeQuantity xQ2 = lepton_pair_xQ2(
			kinematics,
			order, use_nlp_nlo, 
			charm_mass, primary_muon_min_energy, 
			pdf, ff,
			renormalization_scale, factorization_scale, fragmentation_scale
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
	void lepton_pair_zxy_base(
		const auto &sidis,
		const std::vector<double> &zs,
		const std::vector<double> &xs,
		const std::vector<double> &ys,
		const std::vector<double> &Ebeams,
		const std::optional<std::string> variation,
		const std::size_t variation_index,
		const std::size_t variation_count,
		std::ofstream &file,
		const std::string comment
	) const {
		const std::size_t z_count = zs.size();
		const std::size_t x_count = xs.size();
		const std::size_t y_count = ys.size();
		const std::size_t E_count = Ebeams.size();

		const std::size_t count = z_count * x_count * y_count * E_count * variation_count;

		std::size_t calculated_values = 0;
		
		output_run_info(file, sidis, comment);

		file << "z,x,y,E,LO,NLO,NNLO,Q2,renormalization_scale,factorization_scale,fragmentation_scale" << IO::endl;
		std::streamsize original_precision = std::cout.precision();

		#pragma omp parallel if(parallelize()) num_threads(number_of_threads) firstprivate(sidis)
		{
			#pragma omp for collapse(4) schedule(guided)
			for (std::size_t i = 0; i < x_count; i++) {
				for (std::size_t j = 0; j < E_count; j++) {
					for (std::size_t k = 0; k < y_count; k++) {
						for (std::size_t l = 0; l < z_count; l++) {
							const double x = xs[i];
							const double y = ys[k];
							const double E_beam = Ebeams[j];
							const double z = zs[l];
							
							TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
							const PerturbativeQuantity cross_section_zxQ2 = sidis.lepton_pair_cross_section_zxQ2(kinematics, z);
							const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
							const PerturbativeQuantity cross_section_zxy = cross_section_zxQ2 * jacobian;

							const double Q2 = kinematics.Q2;
							const double renormalization_scale = sidis.renormalization_scale_function(kinematics);
							const double factorization_scale = sidis.factorization_scale_function(kinematics);
							const double fragmentation_scale = sidis.fragmentation_scale_function(kinematics);

							#pragma omp critical
							{
								file << z << "," << x << ", " << y << ", " << E_beam << ", ";
								file << cross_section_zxy.lo << ", " << cross_section_zxy.nlo << ", " << cross_section_zxy.nnlo;
								file << ", " << Q2 << ", " << renormalization_scale << ", " << factorization_scale << ", " << fragmentation_scale << IO::endl;
								file.flush();

								calculated_values++;

								const int base_precision = 5;
								const int s_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.s)));
								const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));
								const int Q2_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.Q2)));

								std::cout << std::fixed << std::setprecision(base_precision);
								std::cout << "[SIDIS] " << IO::leading_zeroes(
									calculated_values + x_count * y_count * E_count * variation_index, Math::number_of_digits(count)
								);
								std::cout << " / " << count;

								if (variation) {
									std::cout << " [" << *variation << " " << IO::leading_zeroes(variation_index + 1, Math::number_of_digits(variation_count));	
									std::cout << " / " << variation_count << "]";
								}

								std::cout << ": " << cross_section_zxy;
								std::cout << " (z = " << z << ", x = " << x << ", y = " << y;
								std::cout << std::setprecision(s_precision) << ", s = " << kinematics.s;
								std::cout << std::setprecision(E_precision) << ", E_beam = " << E_beam;
								std::cout << std::setprecision(Q2_precision) << ", Q2 = " << kinematics.Q2 << ")";
								std::cout << "\r" << std::flush;
							}
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
	ThreadResult lepton_pair_xy_base(
		BS::thread_pool &thread_pool,
		const auto &sidis,
		const std::vector<double> &xs,
		const std::vector<double> &ys,
		const std::vector<double> &E_beams,
		const std::optional<std::string> variation,
		const std::size_t variation_index,
		const std::size_t variation_count,
		const std::filesystem::path &output,
		std::stringstream &header,
		const std::string comment
	) const {	
		output_run_info(header, sidis, comment);
		header << "x,y,E,LO,NLO,NNLO,Q2,renormalization_scale,factorization_scale,fragmentation_scale" << IO::endl;

		shared_multi_future<std::string> futures;

		for (const double x : xs) {
			for (const double y : ys) {
				for (const double E_beam : E_beams) {
					futures.push_back(thread_pool.submit_task([=] {
						std::stringstream stream;
						TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);

						const PerturbativeQuantity cross_section_xQ2 = sidis.lepton_pair_cross_section_xQ2(kinematics);
						const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
						const PerturbativeQuantity cross_section_xy = cross_section_xQ2 * jacobian;

						const double Q2 = kinematics.Q2;
						const double renormalization_scale = sidis.renormalization_scale_function(kinematics);
						const double factorization_scale = sidis.factorization_scale_function(kinematics);
						const double fragmentation_scale = sidis.fragmentation_scale_function(kinematics);

						stream << kinematics.x << ", " << kinematics.y << ", " << kinematics.E_beam << ", " << cross_section_xy.lo << ", " << cross_section_xy.nlo << ", " << cross_section_xy.nnlo;
						stream << ", " << Q2 << ", " << renormalization_scale << ", " << factorization_scale << ", " << fragmentation_scale << IO::endl;

						return stream.str();
					}).share());
				}
			}
		}

		return ThreadResult { futures, output, header.str() };
	}

	public:
	void lepton_pair_zxy(
		const std::vector<double> &zs,
		const std::vector<double> &xs,
		const std::vector<double> &ys,
		const std::vector<double> &Ebeams,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_interface auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale,
		const std::filesystem::path output,
		const std::string comment = ""
	) const {
		const SIDISComputation sidis = construct_computation(
			order, use_nlp_nlo,
			charm_mass, primary_muon_min_energy,
			pdf, ff,
			renormalization_scale, factorization_scale, fragmentation_scale
		);

		IO::create_directory_tree(output);
		std::ofstream file(output);

		lepton_pair_zxy_base(sidis, zs, xs, ys, Ebeams, std::nullopt, 0, 1, file, comment);
	}
	ThreadResult lepton_pair_xy(
		BS::thread_pool &thread_pool,
		const std::vector<double> &xs,
		const std::vector<double> &ys,
		const std::vector<double> &E_beams,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_interface auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale,
		const std::filesystem::path output,
		const std::string comment = ""
	) const {
		const SIDISComputation sidis = construct_computation(
			order, use_nlp_nlo,
			charm_mass, primary_muon_min_energy,
			pdf, ff,
			renormalization_scale, factorization_scale, fragmentation_scale
		);

		std::stringstream header;

		return lepton_pair_xy_base(thread_pool, sidis, xs, ys, E_beams, std::nullopt, 0, 1, output, header, comment);
	}

	std::vector<ThreadResult> lepton_pair_xy_errors(
		BS::thread_pool &thread_pool,
		const std::vector<double> &xs,
		const std::vector<double> &ys,
		const std::vector<double> &E_beams,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale,
		const std::filesystem::path base_output,
		const std::optional<std::pair<std::size_t, std::size_t>> variation_range = std::nullopt,
		const std::string comment = ""
	) const {
		const bool custom_range = variation_range.has_value();
		const std::size_t variation_count = custom_range ? ((*variation_range).second - (*variation_range).first + 1) : pdf.size();
		const auto variation_start = static_cast<typename std::remove_reference_t<decltype(pdf)>::size_type>(custom_range ? (*variation_range).first : 0);
		const auto variation_end = static_cast<typename std::remove_reference_t<decltype(pdf)>::size_type>(custom_range ? (*variation_range).second + 1 : variation_count);

		std::vector<ThreadResult> results;

		for (typename std::remove_reference_t<decltype(pdf)>::size_type variation_index = variation_start; variation_index < variation_end; variation_index++) {
			const auto &pdf_member = pdf[variation_index];

			const SIDISComputation sidis = construct_computation(
				order, use_nlp_nlo,
				charm_mass, primary_muon_min_energy,
				pdf_member, ff,
				renormalization_scale, factorization_scale, fragmentation_scale
			);

			const std::string path_trail = IO::leading_zeroes(pdf_member.set_member_number, 4);
			std::filesystem::path full_filename = base_output.stem();
			full_filename /= path_trail;
			full_filename.replace_extension(base_output.extension());
			std::filesystem::path output = base_output;
			output.replace_filename(full_filename);

			std::stringstream header;
			header << "#pdf_member = " << variation_index << IO::endl;

			results.push_back(lepton_pair_xy_base(thread_pool, sidis, xs, ys, E_beams, "pdf set", variation_index - variation_start, variation_count, output, header, comment));
		}

		return results;
	}

	std::vector<ThreadResult> lepton_pair_xy_scales(
		BS::thread_pool &thread_pool,
		const std::vector<double> &xs,
		const std::vector<double> &ys,
		const std::vector<double> &E_beams,
		const ScaleVariation scale_variation,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_interface auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization, const is_scale_dependence auto &factorization, const is_scale_dependence auto &fragmentation,
		const std::filesystem::path base_output,
		const std::optional<std::pair<std::size_t, std::size_t>> variation_range = std::nullopt,
		const std::string comment = ""
	) const {
		const std::vector<double> scale_factors{0.25, 1.0, 4.0};
		std::vector<std::vector<double>> scales = Collections::tuples(scale_factors, 3);

		const double min_scale = *std::min_element(scale_factors.begin(), scale_factors.end());
		const double max_scale = *std::max_element(scale_factors.begin(), scale_factors.end());

		std::erase_if(scales, [=](const std::vector<double> tuple) {
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

		const bool custom_range = variation_range.has_value();
		const std::size_t variation_count = custom_range ? ((*variation_range).second - (*variation_range).first + 1) : scales.size();
		const std::size_t variation_start = custom_range ? (*variation_range).first : 0;
		const std::size_t variation_end = custom_range ? (*variation_range).second + 1 : variation_count;

		std::vector<ThreadResult> results;

		for (std::size_t variation_index = variation_start; variation_index < variation_end; variation_index++) {
			const std::vector<double> scale = scales[variation_index];

			const auto renormalization_scale = scale[0] * renormalization;
			const auto factorization_scale = scale[1] * factorization;
			const auto fragmentation_scale = scale[2] * fragmentation;

			const SIDISComputation sidis = construct_computation(
				order, use_nlp_nlo,
				charm_mass, primary_muon_min_energy,
				pdf, ff,
				renormalization_scale, factorization_scale, fragmentation_scale
			);

			const std::string path_trail = (scale[0] == 1.0 && scale[1] == 1.0 && scale[2] == 1.0) ? "base_scale" : "scale_" + std::to_string(variation_index);
			std::filesystem::path full_filename = base_output.stem();
			full_filename /= path_trail;
			full_filename.replace_extension(base_output.extension());
			std::filesystem::path output = base_output;
			output.replace_filename(full_filename);

			std::stringstream header;

			header << "#renormalization_scale = " << scale[0] << IO::endl;
			header << "#factorization_scale = " << scale[1] << IO::endl;
			header << "#fragmentation_scale = " << scale[2] << IO::endl;

			results.push_back(lepton_pair_xy_base(thread_pool, sidis, xs, ys, E_beams, "scale", variation_index - variation_start, variation_count, output, header, comment));
		}

		return results;
	}

	std::vector<ThreadResult> lepton_pair_xy_decays(
		BS::thread_pool &thread_pool,
		const std::vector<double> &xs,
		const std::vector<double> &ys,
		const std::vector<double> &E_beams,
		const auto &decay_variations,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_interface auto &pdf,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale,
		const std::filesystem::path base_output,
		const std::optional<std::pair<std::size_t, std::size_t>> variation_range = std::nullopt,
		const std::string comment = ""
	) const {
		const bool custom_range = variation_range.has_value();
		const std::size_t variation_count = custom_range ? ((*variation_range).second - (*variation_range).first + 1) : decay_variations.size();
		const std::size_t variation_start = custom_range ? (*variation_range).first : 0;
		const std::size_t variation_end = custom_range ? (*variation_range).second + 1 : variation_count;

		std::vector<ThreadResult> results;

		for (std::size_t variation_index = variation_start; variation_index < variation_end; variation_index++) {
			const FragmentationConfiguration decay_variation = decay_variations[variation_index];

			const SIDISComputation sidis = construct_computation(
				order, use_nlp_nlo,
				charm_mass, primary_muon_min_energy,
				pdf, decay_variation,
				renormalization_scale, factorization_scale, fragmentation_scale
			);

			const std::string path_trail = IO::leading_zeroes(variation_index, 4);
			std::filesystem::path full_filename = base_output.stem();
			full_filename /= path_trail;
			full_filename.replace_extension(base_output.extension());
			std::filesystem::path output = base_output;
			output.replace_filename(full_filename);

			std::stringstream header;

			header << "#N = " << decay_variation.decays.front().parametrization.N << IO::endl;
			header << "#alpha = " << decay_variation.decays.front().parametrization.alpha << IO::endl;
			header << "#beta = " << decay_variation.decays.front().parametrization.beta << IO::endl;
			header << "#gamma = " << decay_variation.decays.front().parametrization.gamma << IO::endl;

			results.push_back(lepton_pair_xy_base(thread_pool, sidis, xs, ys, E_beams, "decay", variation_index - variation_start, variation_count, output, header, comment));
		}

		return results;
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* 							INTEGRATED CROSS SECTIONS							   */
	/* 								 MULTIPLE VALUES 							  	   */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

	private:

	void x_integrated_lepton_pair_base(
		const auto &sidis,
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
		
		output_run_info(file, sidis, comment);

		file << "E,Q2,LO,NLO,NNLO" << IO::endl;
		std::streamsize original_precision = std::cout.precision();

		#pragma omp parallel if(parallelize()) num_threads(number_of_threads) firstprivate(sidis)
		{
			#pragma omp for collapse(2) schedule(guided)
			for (std::size_t i = 0; i < E_count; i++) {
				for (std::size_t j = 0; j < Q2_count; j++) {
					const double E_beam = Ebeams[i];
					const double Q2 = Q2s[j];
					
					const TRFKinematics placeholder_kinematics = TRFKinematics::y_E_beam(-1.0, -1.0, E_beam, process.target.mass, process.projectile.mass);
					const PerturbativeQuantity cross_section = sidis.x_integrated_lepton_pair_cross_section(placeholder_kinematics, Q2);

					#pragma omp critical
					{
						file << E_beam << ", " << Q2 << ", " << cross_section << IO::endl;
						file.flush();

						calculated_values++;

						const int base_precision = 5;
						const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));
						const int Q2_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(Q2)));

						std::cout << std::fixed << std::setprecision(base_precision);
						std::cout << "[SIDIS] " << IO::leading_zeroes(calculated_values + E_count * Q2_count * variation_index, Math::number_of_digits(count));
						std::cout << " / " << count;

						if (variation) {
							std::cout << " [" << *variation << " " << IO::leading_zeroes(variation_index + 1, Math::number_of_digits(variation_count));
							std::cout << " / " << variation_count << "]";
						}

						std::cout << ": " << cross_section;
						std::cout << std::setprecision(E_precision) << " (E_beam = " << E_beam;
						std::cout << std::setprecision(Q2_precision) << ", Q2 = " << Q2 << ")";
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
	void integrated_lepton_pair_base(
		const auto &sidis,
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
		
		output_run_info(file, sidis, comment);

		file << "E,LO,NLO,NNLO" << IO::endl;
		std::streamsize original_precision = std::cout.precision();

		#pragma omp parallel if(parallelize()) num_threads(number_of_threads) firstprivate(sidis)
		{
			#pragma omp for
			for (std::size_t i = 0; i < E_count; i++) {
				const double E_beam = Ebeams[i];
				
				const TRFKinematics placeholder_kinematics = TRFKinematics::y_E_beam(-1.0, -1.0, E_beam, process.target.mass, process.projectile.mass);
				const PerturbativeQuantity cross_section = sidis.integrated_lepton_pair_cross_section(placeholder_kinematics, Q2_min);

				#pragma omp critical
				{
					file << E_beam << ", " << cross_section << IO::endl;
					file.flush();

					calculated_values++;

					const int base_precision = 5;
					const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));

					std::cout << std::fixed << std::setprecision(base_precision);
					std::cout << "[SIDIS] " << IO::leading_zeroes(calculated_values + E_count * variation_index, Math::number_of_digits(count));
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
	void x_integrated_lepton_pair(
		const std::vector<double> &Ebeams,
		const std::vector<double> &Q2s,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_interface auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale,
		const std::filesystem::path output,
		const std::string comment = ""
	) const {
		const SIDISComputation sidis = construct_computation(
			order, use_nlp_nlo,
			charm_mass, primary_muon_min_energy,
			pdf, ff,
			renormalization_scale, factorization_scale, fragmentation_scale
		);

		IO::create_directory_tree(output);
		std::ofstream file(output);

		x_integrated_lepton_pair_base(sidis, Ebeams, Q2s, std::nullopt, 0, 1, file, comment);
	}
	void integrated_lepton_pair(
		const std::vector<double> &Ebeams,
		const double Q2_min,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_interface auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale,
		const std::filesystem::path output,
		const std::string comment = ""
	) const {
		const SIDISComputation sidis = construct_computation(
			order, use_nlp_nlo,
			charm_mass, primary_muon_min_energy,
			pdf, ff,
			renormalization_scale, factorization_scale, fragmentation_scale
		);

		IO::create_directory_tree(output);
		std::ofstream file(output);

		integrated_lepton_pair_base(sidis, Ebeams, Q2_min, std::nullopt, 0, 1, file, comment);
	}

	void integrated_lepton_pair_errors(
		const std::vector<double> &Ebeams,
		const double Q2_min,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_set_interface auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale,
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

			const SIDISComputation sidis = construct_computation(
				order, use_nlp_nlo,
				charm_mass, primary_muon_min_energy,
				pdf_member, ff,
				renormalization_scale, factorization_scale, fragmentation_scale
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

			integrated_lepton_pair_base(sidis, Ebeams, Q2_min, "pdf set", variation_index - variation_start, variation_count, file, comment);
		}
	}

	void integrated_lepton_pair_scales(
		const std::vector<double> &Ebeams,
		const double Q2_min,
		const ScaleVariation scale_variation,
		const PerturbativeOrder order, const bool use_nlp_nlo,
		const double charm_mass, const double primary_muon_min_energy,
		const is_pdf_interface auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization, const is_scale_dependence auto &factorization, const is_scale_dependence auto &fragmentation,
		const std::filesystem::path base_output,
		const std::optional<std::pair<std::size_t, std::size_t>> variation_range = std::nullopt,
		const std::string comment = ""
	) const {
		const std::vector<double> scale_factors{0.25, 1.0, 4.0};
		std::vector<std::vector<double>> scales = Collections::tuples(scale_factors, 3);

		const double min_scale = *std::min_element(scale_factors.begin(), scale_factors.end());
		const double max_scale = *std::max_element(scale_factors.begin(), scale_factors.end());

		std::erase_if(scales, [=](const std::vector<double> tuple) {
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

		const bool custom_range = variation_range.has_value();
		const std::size_t variation_count = custom_range ? ((*variation_range).second - (*variation_range).first + 1) : scales.size();
		const std::size_t variation_start = custom_range ? (*variation_range).first : 0;
		const std::size_t variation_end = custom_range ? (*variation_range).second + 1 : variation_count;

		for (std::size_t variation_index = variation_start; variation_index < variation_end; variation_index++) {
			const std::vector<double> scale = scales[variation_index];

			const auto renormalization_scale = scale[0] * renormalization;
			const auto factorization_scale = scale[1] * factorization;
			const auto fragmentation_scale = scale[2] * fragmentation;

			const SIDISComputation sidis = construct_computation(
				order, use_nlp_nlo,
				charm_mass, primary_muon_min_energy,
				pdf, ff,
				renormalization_scale, factorization_scale, fragmentation_scale
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
			file << "#fragmentation_scale = " << scale[2] << IO::endl;

			integrated_lepton_pair_base(sidis, Ebeams, Q2_min, "scale", variation_index - variation_start, variation_count, file, comment);
		}
	}

	private:
	auto construct_computation(
		const PerturbativeOrder order,
		const bool use_nlp_nlo,
		const double charm_mass,
		const double primary_muon_min_energy,
		const auto &pdf, const auto &ff,
		const is_scale_dependence auto &renormalization_scale, const is_scale_dependence auto &factorization_scale, const is_scale_dependence auto &fragmentation_scale) const {
		return SIDISComputation(
			active_flavors,
			{
				top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
			},
			pdf, ff,
			integration_parameters,
			process,
			renormalization_scale, factorization_scale, fragmentation_scale,
			use_modified_cross_section_prefactor, order, use_nlp_nlo,
			primary_muon_min_energy
		);
	}

	void output_run_info(std::iostream &file, const auto &computation, const std::string comment) const {
		file << "#timestamp = " << std::format("{:%d-%m-%Y %H:%M:%OS}", std::chrono::system_clock::now()) << IO::endl;
		file << "#active_flavors = ";
		for (const FlavorType flavor : computation.flavors.active_flavors) {
			file << flavor << " ";
		}
		file << IO::endl;
		
		file << "#pdf = " << computation.pdf1.set_name << " [" << typeid(computation.pdf1).name() << "]" << IO::endl;
		file << "#ff = ";
		for (const auto &[frag, decay] : std::views::zip(computation.ff1.interfaces, computation.ff1.decays)) {
			const DecayParametrization p = decay.parametrization;
			file << frag.set_name << " [N = " << p.N << " | alpha = " << p.alpha << " | beta = " << p.beta << " | gamma = " << p.gamma << "] ";
		}
		file << IO::endl;

		file << "#process = " << computation.process << IO::endl;
		file << "#use_nlp_nlo = " << Conversion::bool_to_string(computation.use_nlp_nlo) << IO::endl;
		file << "#parallelize = " << Conversion::bool_to_string(parallelize()) << IO::endl;
		file << "#number_of_threads = " << number_of_threads << IO::endl;
		file << "#up_mass = " << up_mass << IO::endl;
		file << "#down_mass = " << down_mass << IO::endl;
		file << "#charm_mass = " << computation.flavors.flavor_masses[Flavor::Charm + 6] << IO::endl;
		file << "#strange_mass = " << strange_mass << IO::endl;
		file << "#top_mass = " << top_mass << IO::endl;
		file << "#bottom_mass = " << bottom_mass << IO::endl;
		file << computation.integration_parameters;

		if (!comment.empty()) {
			file << "#comment = " << comment << IO::endl;
		}
	}
};

// /// @brief 
// /// @tparam PDFInterface The type used for PDFs. Must be either a type that satisfies the is_pdf_interface concept or a specialization of LHASetInterface.
// /// @tparam FFInterface The type used for fragmentation functions. Must be a type that satisfies the is_pdf_interface concept.
// /// @tparam DecayFunction The type for the decay function, most likely a lambda. Must satisfy the is_decay_function concept. 
// /// Defaults to the type of DecayFunctions::trivial.
// /// @tparam RenormalizationScale The type for the renormalization scale function, most likely a lambda. Must satisfy the is_scale_dependence concept.
// /// Defaults to the type of ScaleDependence::trivial::type.
// /// @tparam FactorizationScale The type for the factorization scale function, most likely a lambda. Must satisfy the is_scale_dependence concept.
// /// Defaults to the type of ScaleDependence::trivial::type.
// /// @tparam FragmentationScale The type for the fragmentation scale function, most likely a lambda. Must satisfy the is_scale_dependence concept.
// /// Defaults to the type of ScaleDependence::trivial::type.
// template <
// 	typename PDFInterface,
// 	is_pdf_interface FFInterface, 
// 	is_decay_function DecayFunction = decltype(DecayFunctions::trivial), 
// 	is_scale_dependence RenormalizationScale = decltype(ScaleDependence::trivial)::type, 
// 	is_scale_dependence FactorizationScale = decltype(ScaleDependence::trivial)::type, 
// 	is_scale_dependence FragmentationScale = decltype(ScaleDependence::trivial)::type
// > requires is_pdf_interface<PDFInterface> || is_instance<PDFInterface, LHASetInterface>
// struct SIDIS {
// 	// Compile-time constant that tells whether the PDF type (PDFInterface) is a specialization of LHASetInterface,
// 	// which implies that PDF error sets are available.
// 	static constexpr bool has_pdf_error_sets = is_instance<PDFInterface, LHASetInterface>;

// 	const FlavorVector active_flavors;

// 	const PDFInterface pdf;
// 	const FragmentationConfiguration<FFInterface, DecayFunction> ff;

// 	// Parallelization with OpenMP over kinematical points (i.e. over different values of x, y and E). Enabled by default.
// 	bool parallelize = true;
// 	// If parallelization is enabled, controls the number of threads passed on to OpenMP. Defaults to half the number of cores on the machine.
// 	unsigned int number_of_threads = Utility::get_default_thread_count() / 2;
// 	// Parameters passed on to all integrators.
// 	IntegrationParameters integration_parameters = IntegrationParameters();
// 	// The physical process considered, whether it's neutrino or antineutrino scattering and what are the target and projectile particles.
// 	const Process process;
// 	// Global square root of Mandelstam s, used for some kinematics when not everything is specified.
// 	double global_sqrt_s;
// 	// A wrapped function that computes the renormalization scale in a given kinematical point.
// 	const ScaleDependence::Function<RenormalizationScale> renormalization_scale;
// 	// A wrapped function that computes the factorization scale in a given kinematical point.
// 	const ScaleDependence::Function<FactorizationScale> factorization_scale;
// 	// A wrapped function that computes the fragmentation scale in a given kinematical point.
// 	const ScaleDependence::Function<FragmentationScale> fragmentation_scale;

// 	// If enabled, factorization scale will be frozen to the smallest value determined by the PDF.
// 	// Requesting the value of the PDF at a smaller Q^2 value will then give the value with the smallest
// 	// allowed Q^2 instead. Additionally, scale logarithms will be frozen as well. In essence,
// 	// f(x, Q^2 < Q^2_min) = f(x, Q^2_min) and the same for scale logarithms.
// 	bool freeze_factorization_scale = false;
// 	// If enabled, fragmentation scale will be frozen to the smallest value determined by the FF.
// 	// Requesting the value of the FF at a smaller Q^2 value will then give the value with the smallest
// 	// allowed Q^2 instead. Additionally, scale logarithms will be frozen as well. In essence,
// 	// f(x, Q^2 < Q^2_min) = f(x, Q^2_min) and the same for scale logarithms.
// 	bool freeze_fragmentation_scale = true;

// 	bool maintain_order_separation = true;
// 	bool combine_integrals = false;
// 	bool compute_differential_cross_section_directly = false;

// 	double up_mass = 0.0;
// 	double down_mass = 0.0;
// 	double charm_mass = 0.0;
// 	double strange_mass = 0.0;
// 	double top_mass = 0.0;
// 	double bottom_mass = 0.0;

// 	bool use_modified_cross_section_prefactor = false;

// 	PerturbativeOrder order = PerturbativeOrder::NLO;
// 	bool use_nlp_nlo = false;

// 	ScaleVariation scale_variation = ScaleVariation::None;

// 	std::vector<DecayParametrization> decay_variations{};
// 	bool decay_variation = false;

// 	// Sets the minimum energy required for the muon coming out of the neutrino --> W + muon vertex, which introduces a maximum y value y_max = 1 - min / E_beam.
// 	// Enforced only in the integrated lepton-pair cross section.
// 	double primary_muon_min_energy = 0.0;

// 	SIDIS (
// 		const FlavorVector _active_flavors, 
// 		const PDFInterface _pdf, 
// 		const FragmentationConfiguration<FFInterface, DecayFunction> _ff, 
// 		const Process _process,
// 		const ScaleDependence::Function<RenormalizationScale> _renormalization_scale = ScaleDependence::Function<RenormalizationScale>(),
// 		const ScaleDependence::Function<FactorizationScale> _factorization_scale = ScaleDependence::Function<FactorizationScale>(),
// 		const ScaleDependence::Function<FragmentationScale> _fragmentation_scale = ScaleDependence::Function<FragmentationScale>())
// 	: active_flavors(_active_flavors), 
// 	pdf(_pdf),
// 	ff(_ff),
// 	process(_process),
// 	renormalization_scale(_renormalization_scale),
// 	factorization_scale(_factorization_scale),
// 	fragmentation_scale(_fragmentation_scale) {	}

// 	SIDIS (
// 		const FlavorVector _active_flavors, 
// 		const PDFInterface _pdf, 
// 		const FFInterface _ff, 
// 		const Process _process,
// 		const ScaleDependence::Function<RenormalizationScale> _renormalization_scale = ScaleDependence::Function<RenormalizationScale>(),
// 		const ScaleDependence::Function<FactorizationScale> _factorization_scale = ScaleDependence::Function<FactorizationScale>(),
// 		const ScaleDependence::Function<FragmentationScale> _fragmentation_scale = ScaleDependence::Function<FragmentationScale>())
// 	: SIDIS(
// 		_active_flavors, 
// 		_pdf, 
// 		FragmentationConfiguration<FFInterface, DecayFunction>({_ff}), 
// 		_process, 
// 		_renormalization_scale, 
// 		_factorization_scale,
// 		_fragmentation_scale) { }

// 	private:
// 	auto construct_computation() const requires is_pdf_interface<PDFInterface> {
// 		SIDISComputation sidis(
// 			active_flavors, 
// 			{
// 				top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
// 			},
// 			pdf, ff,
// 			integration_parameters,
// 			process, 
// 			renormalization_scale, factorization_scale, fragmentation_scale,
// 			use_modified_cross_section_prefactor, order, use_nlp_nlo,
// 			primary_muon_min_energy
// 		);
// 		return sidis;
// 	}
// 	template <is_pdf_interface PDF>
// 	auto construct_computation(const PDF pdf_member) const requires has_pdf_error_sets {
// 		SIDISComputation sidis(
// 			active_flavors, 
// 			{
// 				top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
// 			},
// 			pdf_member, ff,
// 			integration_parameters,
// 			process, 
// 			renormalization_scale, factorization_scale, fragmentation_scale,
// 			use_modified_cross_section_prefactor, order, use_nlp_nlo,
// 			primary_muon_min_energy
// 		);
// 		return sidis;
// 	}
// 	auto construct_computation_scale_variation(const std::vector<double> scales) const requires is_pdf_interface<PDFInterface> {
// 		const double pdf_Q2_min = freeze_factorization_scale ? pdf.Q2_min() : std::numeric_limits<double>::min();
// 		const double pdf_Q2_max = freeze_factorization_scale ? pdf.Q2_max() : std::numeric_limits<double>::max();

// 		const double ff_Q2_min = freeze_fragmentation_scale ? ff.Q2_min() : std::numeric_limits<double>::min();
// 		const double ff_Q2_max = freeze_fragmentation_scale ? ff.Q2_max() : std::numeric_limits<double>::max();

// 		const auto renormalization = ScaleDependence::multiplicative(scales[0]);
// 		const auto factorization = ScaleDependence::clamped_multiplicative(scales[1], pdf_Q2_min, pdf_Q2_max);
// 		const auto fragmentation = ScaleDependence::clamped_multiplicative(scales[2], ff_Q2_min, ff_Q2_max);

// 		SIDISComputation sidis(
// 			active_flavors, 
// 			{
// 				top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
// 			},
// 			pdf, ff,
// 			integration_parameters,
// 			process, 
// 			renormalization, factorization, fragmentation,
// 			use_modified_cross_section_prefactor, order, use_nlp_nlo,
// 			primary_muon_min_energy
// 		);
// 		return sidis;
// 	}
// 	auto construct_computation_ff_variation(const FragmentationConfiguration<FFInterface, DecayFunction> ff_variation) const requires is_pdf_interface<PDFInterface> {
// 		SIDISComputation sidis(
// 			active_flavors, 
// 			{
// 				top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
// 			},
// 			pdf, ff_variation,
// 			integration_parameters,
// 			process, 
// 			renormalization_scale, factorization_scale, fragmentation_scale,
// 			use_modified_cross_section_prefactor, order, use_nlp_nlo,
// 			primary_muon_min_energy
// 		);
// 		return sidis;
// 	}
// 	void output_run_info(std::ofstream &file, const auto &sidis, const std::string comment) {
// 		file << "#timestamp = " << std::format("{:%d-%m-%Y %H:%M:%OS}", std::chrono::system_clock::now()) << IO::endl;
// 		file << "#cross_section = d^2s/dxdy" << IO::endl;
// 		file << "#active_flavors = ";
// 		for (const FlavorType flavor : active_flavors) {
// 			file << flavor << " ";
// 		}
// 		file << IO::endl;
		
// 		file << "#pdf = " << pdf.set_name << " [" << typeid(pdf).name() << "]" << IO::endl;
// 		file << "#ff = ";
// 		for (const auto &[frag, decay] : std::views::zip(sidis.ff1.interfaces, sidis.ff1.decays)) {
// 			const DecayParametrization p = decay.parametrization;
// 			file << frag.set_name << " [N = " << p.N << " | alpha = " << p.alpha << " | beta = " << p.beta << " | gamma = " << p.gamma << "] ";
// 		}
// 		file << IO::endl;

// 		file << "#process = " << process << IO::endl;
// 		file << "#use_nlp_nlo = " << Conversion::bool_to_string(use_nlp_nlo) << IO::endl;
// 		file << "#parallelize = " << Conversion::bool_to_string(parallelize) << IO::endl;
// 		file << "#number_of_threads = " << number_of_threads << IO::endl;
// 		file << "#freeze_factorization_scale = " << freeze_factorization_scale << IO::endl;
// 		file << "#freeze_fragmentation_scale = " << freeze_fragmentation_scale << IO::endl;
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
// 	PerturbativeQuantity compute_structure_function(
// 		const StructureFunction F, const double z, const TRFKinematics kinematics
// 	) const requires is_pdf_interface<PDFInterface> {
// 		SIDISComputation sidis = construct_computation();
// 		return sidis.compute_structure_function(F, z, kinematics);
// 	}
// 	PerturbativeQuantity compute_structure_function(
// 		const StructureFunction F, const double x, const double z, const double Q2
// 	) const requires is_pdf_interface<PDFInterface> {
// 		TRFKinematics kinematics = TRFKinematics::Q2_sqrt_s(x, Q2, global_sqrt_s, process.target.mass, process.projectile.mass);
// 		return compute_structure_function(F, z, kinematics);
// 	}

// 	PerturbativeQuantity F2(const double z, const TRFKinematics kinematics) const requires is_pdf_interface<PDFInterface> {
// 		return compute_structure_function(StructureFunction::F2, z, kinematics);
// 	}
// 	PerturbativeQuantity FL(const double z, const TRFKinematics kinematics) const requires is_pdf_interface<PDFInterface> {
// 		return compute_structure_function(StructureFunction::FL, z, kinematics);
// 	}
// 	PerturbativeQuantity F3(const double z, const TRFKinematics kinematics) const requires is_pdf_interface<PDFInterface> {
// 		return compute_structure_function(StructureFunction::F3, z, kinematics);
// 	}

// 	PerturbativeQuantity F2(const double x, const double z, const double Q2) const requires is_pdf_interface<PDFInterface> {
// 		return compute_structure_function(StructureFunction::F2, x, z, Q2);
// 	}
// 	PerturbativeQuantity FL(const double x, const double z, const double Q2) const requires is_pdf_interface<PDFInterface> {
// 		return compute_structure_function(StructureFunction::FL, x, z, Q2);
// 	}
// 	PerturbativeQuantity F3(const double x, const double z, const double Q2) const requires is_pdf_interface<PDFInterface> {
// 		return compute_structure_function(StructureFunction::F3, x, z, Q2);
// 	}

// 	PerturbativeQuantity differential_cross_section_xQ2(const double z, const TRFKinematics kinematics) const requires is_pdf_interface<PDFInterface> {
// 		SIDISComputation sidis = construct_computation();
// 		return sidis.differential_cross_section_xQ2(z, kinematics);
// 	}

// 	PerturbativeQuantity integrated_lepton_pair_cross_section(const double E_beam, const double Q2_min) const {
// 		SIDISComputation sidis = construct_computation();

// 		const TRFKinematics placeholder_kinematics = TRFKinematics::y_E_beam(-1.0, -1.0, E_beam, process.target.mass, process.projectile.mass);
// 		const PerturbativeQuantity cross_section = sidis.integrated_lepton_pair_cross_section(placeholder_kinematics, Q2_min);

// 		return cross_section;
// 	}

// 	PerturbativeQuantity differential_cross_section_xy(const double z, const TRFKinematics kinematics) const requires is_pdf_interface<PDFInterface> {
// 		const PerturbativeQuantity xQ2 = differential_cross_section_xQ2(z, kinematics);
// 		const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
// 		return xQ2 * jacobian;
// 	}

// 	PerturbativeQuantity lepton_pair_cross_section_xQ2(const TRFKinematics kinematics) const requires is_pdf_interface<PDFInterface>{
// 		SIDISComputation sidis = construct_computation();
// 		return sidis.lepton_pair_cross_section_xQ2(kinematics);
// 	}
// 	PerturbativeQuantity lepton_pair_cross_section_xy(const TRFKinematics kinematics) const requires is_pdf_interface<PDFInterface> {
// 		const PerturbativeQuantity xQ2 = lepton_pair_cross_section_xQ2(kinematics);
// 		const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
// 		return xQ2 * jacobian;
// 	}

// 	void differential_cross_section_xQ2(
// 		const std::vector<double> x_bins, 
// 		const std::vector<double> z_bins, 
// 		const std::vector<double> Q2_bins, 
// 		const std::filesystem::path output) requires is_pdf_interface<PDFInterface> {

// 		const std::size_t x_step_count = x_bins.size();
// 		const std::size_t z_step_count = z_bins.size();
// 		const std::size_t Q2_step_count = Q2_bins.size();

// 		int calculated_values = 0;

// 		IO::create_directory_tree(output);
// 		std::ofstream file(output);

// 		#pragma omp parallel for if(parallelize) num_threads(number_of_threads) collapse(3)
// 		for (std::size_t i = 0; i < x_step_count; i++) {
// 			for (std::size_t j = 0; j < Q2_step_count; j++) {
// 				for (std::size_t k = 0; k < z_step_count; k++) {
// 					const double x = x_bins[i];
// 					const double z = z_bins[k];
// 					const double Q2 = Q2_bins[j];
					
// 					const PerturbativeQuantity differential_cs = differential_cross_section_xQ2(x, z, Q2);

// 					#pragma omp critical
// 					{
// 						file << x << ", " << z << ", " << Q2 << ", " << differential_cs.lo << ", " << differential_cs.nlo << ", " << differential_cs.nnlo << IO::endl;
// 						file.flush();

// 						calculated_values++;
// 						std::cout << "Calculated value " << calculated_values << " / " << x_step_count * z_step_count * Q2_step_count;
// 						std::cout << " (x = " << x << ", z = " << z << ", Q2 = " << Q2 << ")" << IO::endl;
// 					}
// 				}
// 			}
// 		}

// 		file.close();
// 	}

// 	void lepton_pair_cross_section_xy(
// 		const std::vector<double> x_bins, 
// 		const std::vector<double> y_bins, 
// 		const std::vector<double> E_beam_bins, 
// 		const std::filesystem::path output, 
// 		const std::string comment = "") requires is_pdf_interface<PDFInterface> {

// 		if (scale_variation != ScaleVariation::None) { return lepton_pair_cross_section_xy_scale_variations(x_bins, y_bins, E_beam_bins, output, comment); }
// 		if (decay_variation) { return lepton_pair_cross_section_xy_ff_variations(x_bins, y_bins, E_beam_bins, output, comment); }

// 		const std::size_t x_step_count = x_bins.size();
// 		const std::size_t y_step_count = y_bins.size();
// 		const std::size_t E_beam_step_count = E_beam_bins.size();
// 		const std::size_t total_count = x_step_count * y_step_count * E_beam_step_count;

// 		int calculated_values = 0;

// 		IO::create_directory_tree(output);
// 		std::ofstream file(output);

// 		SIDISComputation sidis = construct_computation();

// 		output_run_info(file, sidis, comment);

// 		file << "x,y,E,LO,NLO,NNLO,Q2,renormalization_scale,factorization_scale,fragmentation_scale" << IO::endl;
// 		std::streamsize original_precision = std::cout.precision();

// 		#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(sidis)
// 		{
// 			#pragma omp for collapse(3) schedule(guided)
// 			for (std::size_t i = 0; i < x_step_count; i++) {
// 				for (std::size_t j = 0; j < E_beam_step_count; j++) {
// 					for (std::size_t k = 0; k < y_step_count; k++) {
// 						const double x = x_bins[i];
// 						const double y = y_bins[k];
// 						const double E_beam = E_beam_bins[j];
						
// 						TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
// 						const PerturbativeQuantity cross_section_xQ2 = sidis.lepton_pair_cross_section_xQ2(kinematics);
// 						const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
// 						const PerturbativeQuantity cross_section_xy = cross_section_xQ2 * jacobian;

// 						const double Q2 = kinematics.Q2;
// 						const double renormalization_scale = sidis.renormalization_scale_function(kinematics);
// 						const double factorization_scale = sidis.factorization_scale_function(kinematics);
// 						const double fragmentation_scale = sidis.fragmentation_scale_function(kinematics);

// 						#pragma omp critical
// 						{
// 							file << x << ", " << y << ", " << E_beam << ", " << cross_section_xy.lo << ", " << cross_section_xy.nlo << ", " << cross_section_xy.nnlo;
// 							file << ", " << Q2 << ", " << renormalization_scale << ", " << factorization_scale << ", " << fragmentation_scale << IO::endl;
// 							file.flush();

// 							calculated_values++;

// 							const int base_precision = 5;
// 							const int s_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.s)));
// 							const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));
// 							const int Q2_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.Q2)));

// 							std::cout << std::fixed << std::setprecision(base_precision);
// 							std::cout << "[SIDIS] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
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

// 	void lepton_pair_cross_section_xy_error_sets(
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
// 			SIDISComputation sidis = construct_computation(pdf_member);

// 			const std::string path_trail = IO::leading_zeroes(pdf_member.set_member_number, 4);
// 			std::filesystem::path full_filename = base_output.stem();
// 			full_filename /= path_trail;
// 			full_filename.replace_extension(base_output.extension());
// 			std::filesystem::path output = base_output;
// 			output.replace_filename(full_filename);

// 			IO::create_directory_tree(output);
// 			std::ofstream file(output);

// 			output_run_info(file, sidis, comment);
// 			file << "#pdf_member = " << member_index << IO::endl;

// 			file << "x,y,E,LO,NLO,NNLO,Q2,renormalization_scale,factorization_scale,fragmentation_scale" << IO::endl;

// 			#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(sidis)
// 			{
// 				#pragma omp for collapse(3) schedule(guided)
// 				for (std::size_t i = 0; i < x_step_count; i++) {
// 					for (std::size_t j = 0; j < E_beam_step_count; j++) {
// 						for (std::size_t k = 0; k < y_step_count; k++) {
// 							const double x = x_bins[i];
// 							const double y = y_bins[k];
// 							const double E_beam = E_beam_bins[j];
							
// 							TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
// 							const PerturbativeQuantity cross_section_xQ2 = sidis.lepton_pair_cross_section_xQ2(kinematics);
// 							const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
// 							const PerturbativeQuantity cross_section_xy = cross_section_xQ2 * jacobian;

// 							const double Q2 = kinematics.Q2;
// 							const double renormalization_scale = sidis.renormalization_scale_function(kinematics);
// 							const double factorization_scale = sidis.factorization_scale_function(kinematics);
// 							const double fragmentation_scale = sidis.fragmentation_scale_function(kinematics);

// 							#pragma omp critical
// 							{
// 								file << x << ", " << y << ", " << E_beam << ", " << cross_section_xy.lo << ", " << cross_section_xy.nlo << ", " << cross_section_xy.nnlo;
// 								file << ", " << Q2 << ", " << renormalization_scale << ", " << factorization_scale << ", " << fragmentation_scale << IO::endl;
// 								file.flush();

// 								calculated_values++;

// 								const int base_precision = 5;
// 								const int s_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.s)));
// 								const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));
// 								const int Q2_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.Q2)));

// 								std::cout << std::fixed << std::setprecision(base_precision);
// 								std::cout << "[SIDIS] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
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

// 	void integrated_lepton_pair_cross_section(
// 		const std::vector<double> E_beam_bins,
// 		const double Q2_min,
// 		const std::filesystem::path output, 
// 		const std::string comment = "") {

// 		if (scale_variation != ScaleVariation::None) { return integrated_lepton_pair_cross_section_scale_variations(E_beam_bins, Q2_min, output, comment); }

// 		const std::size_t E_beam_step_count = E_beam_bins.size();
// 		const std::size_t total_count = E_beam_step_count;

// 		int calculated_values = 0;

// 		IO::create_directory_tree(output);
// 		std::ofstream file(output);

// 		SIDISComputation sidis = construct_computation();

// 		output_run_info(file, sidis, comment);

// 		file << "E,LO,NLO,NNLO" << IO::endl;
// 		std::streamsize original_precision = std::cout.precision();

// 		#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(sidis)
// 		{
// 			#pragma omp for
// 			for (std::size_t i = 0; i < E_beam_step_count; i++) {
// 				const double E_beam = E_beam_bins[i];
				
// 				const TRFKinematics placeholder_kinematics = TRFKinematics::y_E_beam(-1.0, -1.0, E_beam, process.target.mass, process.projectile.mass);
// 				const PerturbativeQuantity cross_section = sidis.integrated_lepton_pair_cross_section(placeholder_kinematics, Q2_min);

// 				#pragma omp critical
// 				{
// 					file << E_beam << ", " << cross_section << IO::endl;
// 					file.flush();

// 					calculated_values++;

// 					const int base_precision = 5;
// 					const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));

// 					std::cout << std::fixed << std::setprecision(base_precision);
// 					std::cout << "[SIDIS] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
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

// 	void integrated_lepton_pair_cross_section_error_sets(
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
// 			SIDISComputation sidis = construct_computation(pdf_member);

// 			const std::string path_trail = IO::leading_zeroes(pdf_member.set_member_number, 4);
// 			std::filesystem::path full_filename = base_output.stem();
// 			full_filename /= path_trail;
// 			full_filename.replace_extension(base_output.extension());
// 			std::filesystem::path output = base_output;
// 			output.replace_filename(full_filename);

// 			IO::create_directory_tree(output);
// 			std::ofstream file(output);

// 			output_run_info(file, sidis, comment);
// 			file << "#pdf_member = " << member_index << IO::endl;

// 			file << "E,LO,NLO,NNLO" << IO::endl;

// 			#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(sidis)
// 			{
// 				#pragma omp for
// 				for (std::size_t i = 0; i < E_beam_step_count; i++) {
// 					const double E_beam = E_beam_bins[i];
					
// 					TRFKinematics placeholder_kinematics = TRFKinematics::y_E_beam(-1.0, -1.0, E_beam, process.target.mass, process.projectile.mass);
// 					const PerturbativeQuantity cross_section = sidis.integrated_lepton_pair_cross_section(placeholder_kinematics, Q2_min);

// 					#pragma omp critical
// 					{
// 						file << E_beam << ", " << cross_section << IO::endl;
// 						file.flush();

// 						calculated_values++;

// 						const int base_precision = 5;
// 						const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));

// 						std::cout << std::fixed << std::setprecision(base_precision);
// 						std::cout << "[SIDIS] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
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
// 	void lepton_pair_cross_section_xy_scale_variations(
// 		const std::vector<double> x_bins, 
// 		const std::vector<double> y_bins, 
// 		const std::vector<double> E_beam_bins, 
// 		const std::filesystem::path base_output, 
// 		const std::string comment = "") requires is_pdf_interface<PDFInterface> {

// 		const std::size_t x_step_count = x_bins.size();
// 		const std::size_t y_step_count = y_bins.size();
// 		const std::size_t E_beam_step_count = E_beam_bins.size();

// 		int calculated_values = 0;

// 		const std::vector<double> scale_factors{0.25, 1.0, 4.0};
// 		std::vector<std::vector<double>> scales = Collections::tuples(scale_factors, 3);

// 		const double min_scale = *std::min_element(scale_factors.begin(), scale_factors.end());
// 		const double max_scale = *std::max_element(scale_factors.begin(), scale_factors.end());

// 		std::erase_if(scales, [=, this](const std::vector<double> tuple) {
// 			const double renormalization = tuple[0];
// 			const double factorization = tuple[1];
// 			const double fragmentation = tuple[2];

// 			const bool factorization_condition = min_scale <= factorization / renormalization && max_scale >= factorization / renormalization;
// 			const bool fragmentation_condition = min_scale <= fragmentation / renormalization && max_scale >= fragmentation / renormalization;

// 			switch (scale_variation) {
// 			case ScaleVariation::None:
// 				return !(renormalization == 1.0 && factorization == 1.0 && fragmentation == 1.0);			
// 			case ScaleVariation::All:
// 				return !(factorization_condition && fragmentation_condition);
// 			case ScaleVariation::RenormalizationFactorization:
// 				return !(factorization_condition && factorization == fragmentation);
// 			default: return true;
// 			}
// 		});

// 		const std::size_t scale_count = scales.size();

// 		const std::size_t total_count = scale_count * x_step_count * y_step_count * E_beam_step_count;

// 		std::streamsize original_precision = std::cout.precision();

// 		for (std::size_t scale_index = 0; scale_index < scale_count; scale_index++) {
// 			const std::vector<double> scale = scales[scale_index];
// 			SIDISComputation sidis = construct_computation_scale_variation(scale);

// 			const std::string path_trail = (scale[0] == 1.0 && scale[1] == 1.0 && scale[2] == 1.0) ? "base_scale" : "scale_" + std::to_string(scale_index);
// 			std::filesystem::path full_filename = base_output.stem();
// 			full_filename /= path_trail;
// 			full_filename.replace_extension(base_output.extension());
// 			std::filesystem::path output = base_output;
// 			output.replace_filename(full_filename);

// 			IO::create_directory_tree(output);
// 			std::ofstream file(output);

// 			output_run_info(file, sidis, comment);
// 			file << "#renormalization_scale = " << scale[0] << IO::endl;
// 			file << "#factorization_scale = " << scale[1] << IO::endl;
// 			file << "#fragmentation_scale = " << scale[2] << IO::endl;

// 			file << "x,y,E,LO,NLO,NNLO,Q2,renormalization_scale,factorization_scale,fragmentation_scale" << IO::endl;

// 			#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(sidis)
// 			{
// 				#pragma omp for collapse(3) schedule(guided)
// 				for (std::size_t i = 0; i < x_step_count; i++) {
// 					for (std::size_t j = 0; j < E_beam_step_count; j++) {
// 						for (std::size_t k = 0; k < y_step_count; k++) {
// 							const double x = x_bins[i];
// 							const double y = y_bins[k];
// 							const double E_beam = E_beam_bins[j];
							
// 							TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
// 							const PerturbativeQuantity cross_section_xQ2 = sidis.lepton_pair_cross_section_xQ2(kinematics);
// 							const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
// 							const PerturbativeQuantity cross_section_xy = cross_section_xQ2 * jacobian;

// 							const double Q2 = kinematics.Q2;
// 							const double renormalization_scale = sidis.renormalization_scale_function(kinematics);
// 							const double factorization_scale = sidis.factorization_scale_function(kinematics);
// 							const double fragmentation_scale = sidis.fragmentation_scale_function(kinematics);

// 							#pragma omp critical
// 							{
// 								file << x << ", " << y << ", " << E_beam << ", " << cross_section_xy.lo << ", " << cross_section_xy.nlo << ", " << cross_section_xy.nnlo;
// 								file << ", " << Q2 << ", " << renormalization_scale << ", " << factorization_scale << ", " << fragmentation_scale << IO::endl;
// 								file.flush();

// 								calculated_values++;

// 								const int base_precision = 5;
// 								const int s_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.s)));
// 								const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));
// 								const int Q2_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.Q2)));

// 								std::cout << std::fixed << std::setprecision(base_precision);
// 								std::cout << "[SIDIS] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
// 								std::cout << " / " << total_count;
// 								std::cout << " [" << "scale " << IO::leading_zeroes(scale_index + 1, Math::number_of_digits(scale_count)) << " / " << scale_count << "]";
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

// 	void integrated_lepton_pair_cross_section_scale_variations(
// 		const std::vector<double> E_beam_bins,
// 		const double Q2_min,
// 		const std::filesystem::path base_output, 
// 		const std::string comment = "") requires is_pdf_interface<PDFInterface> {

// 		const std::size_t E_beam_step_count = E_beam_bins.size();

// 		int calculated_values = 0;

// 		const std::vector<double> scale_factors{0.25, 1.0, 4.0};
// 		std::vector<std::vector<double>> scales = Collections::tuples(scale_factors, 3);

// 		const double min_scale = *std::min_element(scale_factors.begin(), scale_factors.end());
// 		const double max_scale = *std::max_element(scale_factors.begin(), scale_factors.end());

// 		std::erase_if(scales, [=, this](const std::vector<double> tuple) {
// 			const double renormalization = tuple[0];
// 			const double factorization = tuple[1];
// 			const double fragmentation = tuple[2];

// 			const bool factorization_condition = min_scale <= factorization / renormalization && max_scale >= factorization / renormalization;
// 			const bool fragmentation_condition = min_scale <= fragmentation / renormalization && max_scale >= fragmentation / renormalization;

// 			switch (scale_variation) {
// 			case ScaleVariation::None:
// 				return !(renormalization == 1.0 && factorization == 1.0 && fragmentation == 1.0);			
// 			case ScaleVariation::All:
// 				return !(factorization_condition && fragmentation_condition);
// 			case ScaleVariation::RenormalizationFactorization:
// 				return !(factorization_condition && factorization == fragmentation);
// 			default: return true;
// 			}
// 		});

// 		const std::size_t scale_count = scales.size();

// 		const std::size_t total_count = scale_count * E_beam_step_count;

// 		std::streamsize original_precision = std::cout.precision();

// 		for (std::size_t scale_index = 0; scale_index < scale_count; scale_index++) {
// 			const std::vector<double> scale = scales[scale_index];
// 			SIDISComputation sidis = construct_computation_scale_variation(scale);

// 			const std::string path_trail = (scale[0] == 1.0 && scale[1] == 1.0 && scale[2] == 1.0) ? "base_scale" : "scale_" + std::to_string(scale_index);
// 			std::filesystem::path full_filename = base_output.stem();
// 			full_filename /= path_trail;
// 			full_filename.replace_extension(base_output.extension());
// 			std::filesystem::path output = base_output;
// 			output.replace_filename(full_filename);

// 			IO::create_directory_tree(output);
// 			std::ofstream file(output);

// 			output_run_info(file, sidis, comment);
// 			file << "#renormalization_scale = " << scale[0] << IO::endl;
// 			file << "#factorization_scale = " << scale[1] << IO::endl;
// 			file << "#fragmentation_scale = " << scale[2] << IO::endl;

// 			file << "E,LO,NLO,NNLO" << IO::endl;

// 			#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(sidis)
// 			{
// 				#pragma omp for
// 				for (std::size_t i = 0; i < E_beam_step_count; i++) {
// 					const double E_beam = E_beam_bins[i];
					
// 					TRFKinematics placeholder_kinematics = TRFKinematics::y_E_beam(-1.0, -1.0, E_beam, process.target.mass, process.projectile.mass);
// 					const PerturbativeQuantity cross_section = sidis.integrated_lepton_pair_cross_section(placeholder_kinematics, Q2_min);

// 					#pragma omp critical
// 					{
// 						file << E_beam << ", " << cross_section << IO::endl;
// 						file.flush();

// 						calculated_values++;

// 						const int base_precision = 5;
// 						const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));

// 						std::cout << std::fixed << std::setprecision(base_precision);
// 						std::cout << "[SIDIS] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
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

// 	void lepton_pair_cross_section_xy_ff_variations(
// 		const std::vector<double> x_bins, 
// 		const std::vector<double> y_bins, 
// 		const std::vector<double> E_beam_bins, 
// 		const std::filesystem::path base_output, 
// 		const std::string comment = "") requires is_pdf_interface<PDFInterface> {

// 		const std::size_t x_step_count = x_bins.size();
// 		const std::size_t y_step_count = y_bins.size();
// 		const std::size_t E_beam_step_count = E_beam_bins.size();
// 		const std::size_t variation_count = decay_variations.size();
// 		const std::size_t total_count = variation_count * x_step_count * y_step_count * E_beam_step_count;

// 		int calculated_values = 0;

// 		std::streamsize original_precision = std::cout.precision();

// 		for (std::size_t variation_index = 0; variation_index < variation_count; variation_index++) {
// 			const DecayParametrization parametrization = decay_variations[variation_index];
// 			FragmentationConfiguration ff_variation = ff;

// 			for (auto &decay : ff_variation.decays) {
// 				decay.parametrization = parametrization;
// 			}

// 			SIDISComputation sidis = construct_computation_ff_variation(ff_variation);

// 			const std::string path_trail = IO::leading_zeroes(variation_index, 4);
// 			std::filesystem::path full_filename = base_output.stem();
// 			full_filename /= path_trail;
// 			full_filename.replace_extension(base_output.extension());
// 			std::filesystem::path output = base_output;
// 			output.replace_filename(full_filename);

// 			IO::create_directory_tree(output);
// 			std::ofstream file(output);

// 			output_run_info(file, sidis, comment);
// 			file << "#N = " << ff_variation.decays.front().parametrization.N << IO::endl;
// 			file << "#alpha = " << ff_variation.decays.front().parametrization.alpha << IO::endl;
// 			file << "#beta = " << ff_variation.decays.front().parametrization.beta << IO::endl;
// 			file << "#gamma = " << ff_variation.decays.front().parametrization.gamma << IO::endl;

// 			file << "x,y,E,LO,NLO,NNLO,Q2,renormalization_scale,factorization_scale,fragmentation_scale" << IO::endl;

// 			#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(sidis)
// 			{
// 				#pragma omp for collapse(3) schedule(guided)
// 				for (std::size_t i = 0; i < x_step_count; i++) {
// 					for (std::size_t j = 0; j < E_beam_step_count; j++) {
// 						for (std::size_t k = 0; k < y_step_count; k++) {
// 							const double x = x_bins[i];
// 							const double y = y_bins[k];
// 							const double E_beam = E_beam_bins[j];
							
// 							TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
// 							const PerturbativeQuantity cross_section_xQ2 = sidis.lepton_pair_cross_section_xQ2(kinematics);
// 							const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
// 							const PerturbativeQuantity cross_section_xy = cross_section_xQ2 * jacobian;

// 							const double Q2 = kinematics.Q2;
// 							const double renormalization_scale = sidis.renormalization_scale_function(kinematics);
// 							const double factorization_scale = sidis.factorization_scale_function(kinematics);
// 							const double fragmentation_scale = sidis.fragmentation_scale_function(kinematics);

// 							#pragma omp critical
// 							{
// 								file << x << ", " << y << ", " << E_beam << ", " << cross_section_xy.lo << ", " << cross_section_xy.nlo << ", " << cross_section_xy.nnlo;
// 								file << ", " << Q2 << ", " << renormalization_scale << ", " << factorization_scale << ", " << fragmentation_scale << IO::endl;
// 								file.flush();

// 								calculated_values++;

// 								const int base_precision = 5;
// 								const int s_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.s)));
// 								const int E_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(E_beam)));
// 								const int Q2_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(kinematics.Q2)));

// 								std::cout << std::fixed << std::setprecision(base_precision);
// 								std::cout << "[SIDIS] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
// 								std::cout << " / " << total_count;
// 								std::cout << " [" << "variation ";
// 								std::cout << IO::leading_zeroes(variation_index + 1, Math::number_of_digits(variation_count)) << " / " << variation_count;
// 								std::cout << "]";
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
// };

#endif