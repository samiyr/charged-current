#ifndef SIDIS_H
#define SIDIS_H

#include "SIDIS/SIDISComputation.cpp"
#include "Common/Flavor.cpp"
#include "Common/Process.cpp"
#include "Common/TRFKinematics.cpp"
#include "Decay/DecayFunctions.cpp"
#include "PDF/FragmentationConfiguration.cpp"
#include "Common/ScaleDependence.cpp"
#include <optional>
#include "Common/CommonFunctions.cpp"
#include "PDF/PDFConcept.cpp"
#include "Legacy/FilePath.cpp"

template <
	PDFConcept PDFInterface, 
	PDFConcept FFInterface, 
	DecayFunctions::Concept DecayFunction = decltype(DecayFunctions::trivial), 
	ScaleDependence::Concept RenormalizationScale = decltype(ScaleDependence::trivial)::type, 
	ScaleDependence::Concept FactorizationScale = decltype(ScaleDependence::trivial)::type, 
	ScaleDependence::Concept FragmentationScale = decltype(ScaleDependence::trivial)::type
>
struct SIDIS {
	const FlavorVector active_flavors;

	const PDFInterface pdf;
	const FragmentationConfiguration<FFInterface, DecayFunction> ff;

	bool parallelize = true;
	unsigned int number_of_threads = Utility::get_default_thread_count();

	IntegrationParameters integration_parameters = IntegrationParameters();

	const Process process;

	double global_sqrt_s;
	bool momentum_fraction_mass_corrections = false;

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
	: SIDIS(_active_flavors, _pdf, FragmentationConfiguration<FFInterface, DecayFunction>({_ff}), _process, _renormalization_scale, _factorization_scale, _fragmentation_scale) { }

	private:
	auto construct_computation() const {
		SIDISComputation sidis(
			global_sqrt_s, active_flavors, 
			{
				top_mass, bottom_mass, charm_mass, strange_mass, up_mass, down_mass, 0.0, down_mass, up_mass, strange_mass, charm_mass, bottom_mass, top_mass
			}, 
			pdf, ff,
			integration_parameters,
			process, 
			momentum_fraction_mass_corrections, renormalization_scale, factorization_scale, fragmentation_scale,
			use_modified_cross_section_prefactor
		);
		return sidis;
	}
	void output_run_info(std::ofstream &file, const std::string comment) {
		file << "#cross_section = ds/dxdy" << IO::endl;
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
	PerturbativeQuantity compute_structure_function(const StructureFunction F, const double z, const TRFKinematics kinematics) const {
		SIDISComputation sidis = construct_computation();
		return sidis.compute_structure_function(F, z, kinematics, combine_integrals);
	}
	PerturbativeQuantity compute_structure_function(const StructureFunction F, const double x, const double z, const double Q2) const {
		TRFKinematics kinematics = TRFKinematics::Q2_sqrt_s(x, Q2, global_sqrt_s, process.target.mass, process.projectile.mass);
		return compute_structure_function(F, z, kinematics);
	}

	PerturbativeQuantity F2(const double z, const TRFKinematics kinematics) const {
		return compute_structure_function(StructureFunction::F2, z, kinematics);
	}
	PerturbativeQuantity FL(const double z, const TRFKinematics kinematics) const {
		return compute_structure_function(StructureFunction::FL, z, kinematics);
	}
	PerturbativeQuantity xF3(const double z, const TRFKinematics kinematics) const {
		return compute_structure_function(StructureFunction::xF3, z, kinematics);
	}

	PerturbativeQuantity F2(const double x, const double z, const double Q2) const {
		return compute_structure_function(StructureFunction::F2, x, z, Q2);
	}
	PerturbativeQuantity FL(const double x, const double z, const double Q2) const {
		return compute_structure_function(StructureFunction::FL, x, z, Q2);
	}
	PerturbativeQuantity xF3(const double x, const double z, const double Q2) const {
		return compute_structure_function(StructureFunction::xF3, x, z, Q2);
	}

	PerturbativeQuantity differential_cross_section_xQ2(const double z, const TRFKinematics kinematics) const {
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

	PerturbativeQuantity differential_cross_section_xy(const double z, const TRFKinematics kinematics) const {
		const PerturbativeQuantity xQ2 = differential_cross_section_xQ2(z, kinematics);
		const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
		return xQ2 * jacobian;
	}

	PerturbativeQuantity lepton_pair_cross_section_xQ2(const TRFKinematics kinematics) const {
		SIDISComputation sidis = construct_computation();
		return combine_integrals ? sidis.lepton_pair_cross_section_xQ2_combined(kinematics) : sidis.lepton_pair_cross_section_xQ2_separated(kinematics);
	}
	PerturbativeQuantity lepton_pair_cross_section_xy(const TRFKinematics kinematics) const {
		const PerturbativeQuantity xQ2 = lepton_pair_cross_section_xQ2(kinematics);
		const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
		return xQ2 * jacobian;
	}

	void differential_cross_section_xQ2(const std::vector<double> x_bins, const std::vector<double> z_bins, const std::vector<double> Q2_bins, const FilePath output) {
		const size_t x_step_count = x_bins.size();
		const size_t z_step_count = z_bins.size();
		const size_t Q2_step_count = Q2_bins.size();

		int calculated_values = 0;

		std::ofstream file(output.path());

		#pragma omp parallel for if(parallelize) num_threads(number_of_threads) collapse(3)
		for (size_t i = 0; i < x_step_count; i++) {
			for (size_t j = 0; j < Q2_step_count; j++) {
				for (size_t k = 0; k < z_step_count; k++) {
					const double x = x_bins[i];
					const double z = z_bins[k];
					const double Q2 = Q2_bins[j];
					
					const PerturbativeQuantity differential_cs = differential_cross_section_xQ2(x, z, Q2);

					#pragma omp critical
					{
						file << x << ", " << z << ", " << Q2 << ", " << differential_cs.lo << ", " << differential_cs.nlo << IO::endl;
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
		const FilePath output, 
		const std::string comment = "") {

		const size_t x_step_count = x_bins.size();
		const size_t y_step_count = y_bins.size();
		const size_t E_beam_step_count = E_beam_bins.size();

		int calculated_values = 0;

		std::ofstream file(output.path());

		output_run_info(file, comment);

		file << "x,y,E,LO,NLO,Q2,renormalization_scale,factorization_scale,fragmentation_scale" << IO::endl;

		SIDISComputation sidis = construct_computation();
		#pragma omp parallel if(parallelize) num_threads(number_of_threads) firstprivate(sidis)
		{
			#pragma omp for collapse(3) schedule(guided)
			for (size_t i = 0; i < x_step_count; i++) {
				for (size_t j = 0; j < E_beam_step_count; j++) {
					for (size_t k = 0; k < y_step_count; k++) {
						const double x = x_bins[i];
						const double y = y_bins[k];
						const double E_beam = E_beam_bins[j];
						
						TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, process.target.mass, process.projectile.mass);
						const PerturbativeQuantity cross_section_xQ2 = combine_integrals 
																		? sidis.lepton_pair_cross_section_xQ2_combined(kinematics) 
																		: sidis.lepton_pair_cross_section_xQ2_separated(kinematics);
						const double jacobian = CommonFunctions::xy_jacobian(kinematics, process);
						const PerturbativeQuantity cross_section_xy = cross_section_xQ2 * jacobian;

						const double Q2 = kinematics.Q2;
						const double renormalization_scale = sidis.renormalization_scale_function(kinematics);
						const double factorization_scale = sidis.factorization_scale_function(kinematics);
						const double fragmentation_scale = sidis.fragmentation_scale_function(kinematics);

						#pragma omp critical
						{
							file << x << ", " << y << ", " << E_beam << ", " << cross_section_xy.lo << ", " << cross_section_xy.nlo;
							file << ", " << Q2 << ", " << renormalization_scale << ", " << factorization_scale << ", " << fragmentation_scale << IO::endl;
							file.flush();

							calculated_values++;
							std::cout << "Calculated value " << calculated_values << " / " << x_step_count * y_step_count * E_beam_step_count ;
							std::cout << ": " << cross_section_xy.lo << ", " << cross_section_xy.nlo;
							std::cout << " (x = " << x << ", y = " << y << ", s = " << kinematics.s << ", E_beam = " << E_beam << ", Q2 = " << kinematics.Q2 << ")";
							std::cout << IO::endl;
						}
					}
				}
			}
		}
		file.close();
	}
};

#endif