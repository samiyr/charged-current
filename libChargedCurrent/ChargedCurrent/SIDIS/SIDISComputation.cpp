#ifndef SIDIS_COMPUTATION_H
#define SIDIS_COMPUTATION_H

#include <vector>
#include <optional>

#include "Common/Flavor.cpp"
#include "Common/PerturbativeQuantity.cpp"
#include "Common/TRFKinematics.cpp"
#include "Common/StructureFunction.cpp"

#include "SIDIS/SIDISFunctions.cpp"

#include "Integration/Integrator.cpp"

#include "Decay/Decay.cpp"

#include "PDF/FragmentationConfiguration.cpp"
#include "PDF/PDFConcept.cpp"

#include "Common/ScaleDependence.cpp"

#include "PDF/Interfaces/LHAInterface.cpp"

template <
	is_pdf_interface PDFInterface, 
	is_pdf_interface FFInterface, 
	is_decay_function DecayFunction, 
	is_scale_dependence RenormalizationScale,
	is_scale_dependence FactorizationScale, 
	is_scale_dependence FragmentationScale
>
class SIDISComputation {
	public:
	FlavorInfo flavors;

	const PDFInterface pdf1;
	const FragmentationConfiguration<FFInterface, DecayFunction> ff1;

	const PDFInterface pdf2;
	const FragmentationConfiguration<FFInterface, DecayFunction> ff2;

	const IntegrationParameters integration_parameters;

	const Process process;

	const ScaleDependence::Function<RenormalizationScale> renormalization_scale_function;
	const ScaleDependence::Function<FactorizationScale> factorization_scale_function;
	const ScaleDependence::Function<FragmentationScale> fragmentation_scale_function;

	const bool use_modified_cross_section_prefactor;

	const PerturbativeOrder order;
	const bool use_nlp_nlo;

	const double primary_muon_min_energy;

	SIDISComputation (
		const FlavorVector _active_flavors, 
		const std::array<double, TOTAL_FLAVORS> _flavor_masses,
		const PDFInterface _pdf,
		const FragmentationConfiguration<FFInterface, DecayFunction> _ff,
		const IntegrationParameters _integration_parameters,
		const Process _process,
		const ScaleDependence::Function<RenormalizationScale> _renormalization_scale_function,
		const ScaleDependence::Function<FactorizationScale> _factorization_scale_function,
		const ScaleDependence::Function<FragmentationScale> _fragmentation_scale_function,
		const bool _use_modified_cross_section_prefactor,
		const PerturbativeOrder _order,
		const bool _use_nlp_nlo,
		const double primary_muon_min_energy
	) : flavors(_active_flavors, _flavor_masses),
	pdf1(_pdf), 
	ff1(_ff),
	pdf2(_pdf), 
	ff2(_ff),
	integration_parameters(_integration_parameters), 
	process(_process),
	renormalization_scale_function(_renormalization_scale_function),
	factorization_scale_function(_factorization_scale_function),
	fragmentation_scale_function(_fragmentation_scale_function),
	use_modified_cross_section_prefactor(_use_modified_cross_section_prefactor),
	order(_order),
	use_nlp_nlo(_use_nlp_nlo),
	primary_muon_min_energy(primary_muon_min_energy) { }

	double compute_alpha_s(const TRFKinematics &kinematics) const {
		const double renormalization_scale = renormalization_scale_function(kinematics);
		return pdf1.alpha_s(renormalization_scale);
	}
	
	PerturbativeQuantity structure_function(
		const double z, const TRFKinematics &kinematics, const auto &lo_integrand, const auto &nlo_integrand, const auto &nnlo_integrand, const int sign) const {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;

		double alpha_s = compute_alpha_s(kinematics);
		double nlo_coefficient = alpha_s / (2 * std::numbers::pi);
		double nnlo_coefficient = std::pow(nlo_coefficient, 2);

		const double renormalization_scale = renormalization_scale_function(kinematics);
		const double factorization_scale = factorization_scale_function(kinematics);
		const double fragmentation_scale = fragmentation_scale_function(kinematics);

		const double renormalization_scale_log = renormalization_scale == Q2 ? 0 : std::log(Q2 / renormalization_scale);
		const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);
		const double fragmentation_scale_log = fragmentation_scale == Q2 ? 0 : std::log(Q2 / fragmentation_scale);

		pdf1.evaluate(x, factorization_scale);
		ff1.evaluate(z, fragmentation_scale);

		SIDISFunctions::Parameters<PDFInterface, FFInterface, DecayFunction> params {
			pdf1, ff1, pdf2, ff2,
			flavors,
			nlo_coefficient,
			process, kinematics,
			z,
			renormalization_scale, factorization_scale, fragmentation_scale,
			renormalization_scale_log, factorization_scale_log, fragmentation_scale_log
		};

		const double lo = SIDISFunctions::evaluate<PDFInterface, FFInterface>({}, &params, lo_integrand, false, false, false, sign);

		double nlo = 0.0;
		if (order >= PerturbativeOrder::NLO) {
			Integrator nlo_integrator([nlo_integrand, sign](double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
				return SIDISFunctions::evaluate<PDFInterface, FFInterface>(input, params_in, nlo_integrand, true, true, false, sign);
			}, {x, z}, {1, 1}, integration_parameters, &params);
			const auto nlo_result = nlo_integrator.integrate();
			const double nlo_value = nlo_result.value;

			nlo = nlo_coefficient * nlo_value;
		}

		double nnlo = 0.0;
		if (order >= PerturbativeOrder::NNLO) {
			Integrator nnlo_integrator([nnlo_integrand, sign](double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
				return SIDISFunctions::evaluate<PDFInterface, FFInterface>(input, params_in, nnlo_integrand, true, true, false, sign);
			}, {x, z}, {1, 1}, integration_parameters, &params);
			const auto nnlo_result = nnlo_integrator.integrate();
			const double nnlo_value = nnlo_result.value;

			nnlo = nnlo_coefficient * nnlo_value;
		}

		return PerturbativeQuantity {lo, lo + nlo, lo + nlo + nnlo};
	}

	PerturbativeQuantity F2(const double z, const TRFKinematics &kinematics) const {
		return structure_function(
			z, kinematics, 
			SIDISFunctions::F2::LO::integrand, 
			use_nlp_nlo ? SIDISFunctions::F2::NLO_NLP::total_integrand : SIDISFunctions::F2::NLO::total_integrand, 
			SIDISFunctions::F2::NNLO_NLP::total_integrand, 1
		);
	}
	PerturbativeQuantity FL(const double z, const TRFKinematics &kinematics) const {
		return structure_function(
			z, kinematics, 
			SIDISFunctions::FL::LO::integrand, 
			use_nlp_nlo ? SIDISFunctions::F2::NLO_NLP::total_integrand : SIDISFunctions::F2::NLO::total_integrand, 
			SIDISFunctions::FL::NNLO_NLP::total_integrand, 1
		);
	}
	PerturbativeQuantity F3(const double z, const TRFKinematics &kinematics) const {
		return structure_function(
			z, kinematics, 
			SIDISFunctions::F3::LO::integrand, 
			use_nlp_nlo ? SIDISFunctions::F2::NLO_NLP::total_integrand : SIDISFunctions::F2::NLO::total_integrand, 
			SIDISFunctions::F3::NNLO_NLP::total_integrand, -1
		);
	}

	PerturbativeQuantity compute_structure_function(
		const StructureFunction structure_function, const double z, const TRFKinematics &kinematics) const {
		switch (structure_function) {
		case StructureFunction::F2: return F2(z, kinematics);
		case StructureFunction::FL: return FL(z, kinematics);
		case StructureFunction::F3: return F3(z, kinematics);
		default: return PerturbativeQuantity{0.0, 0.0};
		}
	}

	PerturbativeQuantity differential_cross_section_xQ2(const double z, const TRFKinematics &kinematics) const {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;

		double alpha_s = compute_alpha_s(kinematics);
		double nlo_coefficient = alpha_s / (2 * std::numbers::pi);
		double nnlo_coefficient = std::pow(nlo_coefficient, 2);

		const double renormalization_scale = renormalization_scale_function(kinematics);
		const double factorization_scale = factorization_scale_function(kinematics);
		const double fragmentation_scale = fragmentation_scale_function(kinematics);

		const double renormalization_scale_log = renormalization_scale == Q2 ? 0 : std::log(Q2 / renormalization_scale);
		const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);
		const double fragmentation_scale_log = fragmentation_scale == Q2 ? 0 : std::log(Q2 / fragmentation_scale);

		pdf1.evaluate(x, factorization_scale);
		ff1.evaluate(z, fragmentation_scale);

		SIDISFunctions::Parameters<PDFInterface, FFInterface, DecayFunction> params {
			pdf1, ff1, pdf2, ff2,
			flavors,
			nlo_coefficient,
			process, kinematics,
			z,
			renormalization_scale, factorization_scale, fragmentation_scale,
			renormalization_scale_log, factorization_scale_log, fragmentation_scale_log
		};

		const double lo = SIDISFunctions::cross_section<PDFInterface, FFInterface, DecayFunction>({}, &params, 
			SIDISFunctions::F2::LO::integrand, SIDISFunctions::FL::LO::integrand, SIDISFunctions::F3::LO::integrand,
			false, false, false
		);

		double nlo = 0.0;
		if (order >= PerturbativeOrder::NLO) {
			Integrator nlo_integrator([this](double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
				return SIDISFunctions::cross_section<PDFInterface, FFInterface, DecayFunction>(input, params_in, 
					use_nlp_nlo ? SIDISFunctions::F2::NLO_NLP::total_integrand : SIDISFunctions::F2::NLO::total_integrand, 
					use_nlp_nlo ? SIDISFunctions::FL::NLO_NLP::total_integrand : SIDISFunctions::FL::NLO::total_integrand, 
					use_nlp_nlo ? SIDISFunctions::F3::NLO_NLP::total_integrand : SIDISFunctions::F3::NLO::total_integrand,
					true, true, false
				);
			}, {x, z}, {1, 1}, integration_parameters, &params);
			const auto nlo_result = nlo_integrator.integrate();
			const double nlo_integral = nlo_result.value;

			nlo = nlo_coefficient * nlo_integral;
		}

		double nnlo = 0.0;
		if (order >= PerturbativeOrder::NNLO) {
			Integrator nnlo_integrator([](double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
				return SIDISFunctions::cross_section<PDFInterface, FFInterface, DecayFunction>(input, params_in, 
					SIDISFunctions::F2::NNLO_NLP::total_integrand, SIDISFunctions::FL::NNLO_NLP::total_integrand, SIDISFunctions::F3::NNLO_NLP::total_integrand,
					true, true, false
				);
			}, {x, z}, {1, 1}, integration_parameters, &params);
			const auto nnlo_result = nnlo_integrator.integrate();
			const double nnlo_integral = nnlo_result.value;

			nnlo = nnlo_coefficient * nnlo_integral;
		}


		const PerturbativeQuantity result = PerturbativeQuantity {lo, lo + nlo, lo + nlo + nnlo};
		const double prefactor = use_modified_cross_section_prefactor 
									? CommonFunctions::cross_section_modified_prefactor(kinematics) 
									: CommonFunctions::cross_section_prefactor(kinematics);
		return prefactor * result;
	}

	PerturbativeQuantity lepton_pair_cross_section_zxQ2(const TRFKinematics &kinematics, const double z) const {
		const double x = kinematics.x;
		const double Q2 = kinematics.Q2;

		std::vector<double> z_mins(ff1.decays.size());
		std::transform(ff1.decays.begin(), ff1.decays.end(), z_mins.begin(), [&kinematics](const auto &decay) {
			return SIDISFunctions::Helper::compute_z_min(kinematics, decay);
		});
		const double z_min = *std::min_element(z_mins.begin(), z_mins.end());

		if (z_min > 1.0 || z < z_min) { return PerturbativeQuantity{0.0, 0.0, 0.0}; }

		double alpha_s = compute_alpha_s(kinematics);
		const double nlo_coefficient = alpha_s / (2 * std::numbers::pi);
		const double nnlo_coefficient = std::pow(nlo_coefficient, 2);

		const double renormalization_scale = renormalization_scale_function(kinematics);
		const double factorization_scale = factorization_scale_function(kinematics);
		const double fragmentation_scale = fragmentation_scale_function(kinematics);

		const double renormalization_scale_log = renormalization_scale == Q2 ? 0 : std::log(Q2 / renormalization_scale);
		const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);
		const double fragmentation_scale_log = fragmentation_scale == Q2 ? 0 : std::log(Q2 / fragmentation_scale);

		pdf1.evaluate(x, factorization_scale);

		SIDISFunctions::Parameters<PDFInterface, FFInterface, DecayFunction> params {
			pdf1, ff1, pdf2, ff2,
			flavors,
			nlo_coefficient,
			process, kinematics,
			0.0,
			renormalization_scale, factorization_scale, fragmentation_scale,
			renormalization_scale_log, factorization_scale_log, fragmentation_scale_log
		};

		const std::vector<double> lo_input = {z};

		const double lo = SIDISFunctions::cross_section<PDFInterface, FFInterface, DecayFunction>(lo_input.data(), &params, 
			SIDISFunctions::F2::LO::integrand, SIDISFunctions::FL::LO::integrand, SIDISFunctions::F3::LO::integrand,
			false, false, true
		);

		double nlo = 0.0;
		if (order >= PerturbativeOrder::NLO) {
			Integrator nlo_integrator([&](double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
				double *total_input = new double[3];
				total_input[0] = input[0];
				total_input[1] = input[1];
				total_input[2] = z;

				const double value = SIDISFunctions::cross_section<PDFInterface, FFInterface, DecayFunction>(total_input, params_in, 
					use_nlp_nlo ? SIDISFunctions::F2::NLO_NLP::total_integrand : SIDISFunctions::F2::NLO::total_integrand, 
					use_nlp_nlo ? SIDISFunctions::FL::NLO_NLP::total_integrand : SIDISFunctions::FL::NLO::total_integrand, 
					use_nlp_nlo ? SIDISFunctions::F3::NLO_NLP::total_integrand : SIDISFunctions::F3::NLO::total_integrand,
					true, true, true
				);

				delete[] total_input;

				return value;
			}, {x, z_min}, {1.0, 1.0}, integration_parameters, &params);
			const auto nlo_result = nlo_integrator.integrate();
			const double nlo_integral = nlo_result.value;

			nlo = nlo_coefficient * nlo_integral;
		}

		double nnlo = 0.0;
		if (order >= PerturbativeOrder::NNLO) {
			Integrator nnlo_integrator([&](double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
				double *total_input = new double[3];
				total_input[0] = input[0];
				total_input[1] = input[1];
				total_input[2] = z;

				const double value = SIDISFunctions::cross_section<PDFInterface, FFInterface, DecayFunction>(total_input, params_in, 
					SIDISFunctions::F2::NNLO_NLP::total_integrand, SIDISFunctions::FL::NNLO_NLP::total_integrand, SIDISFunctions::F3::NNLO_NLP::total_integrand,
					true, true, true
				);

				delete[] total_input;

				return value;
			}, {x, z_min}, {1.0, 1.0}, integration_parameters, &params);
			const auto nnlo_result = nnlo_integrator.integrate();
			const double nnlo_integral = nnlo_result.value;

			nnlo = nnlo_coefficient * nnlo_integral;
		}

		const PerturbativeQuantity result = PerturbativeQuantity {lo, lo + nlo, lo + nlo + nnlo};
		const double prefactor = use_modified_cross_section_prefactor 
									? CommonFunctions::cross_section_modified_prefactor(kinematics) 
									: CommonFunctions::cross_section_prefactor(kinematics);
		return prefactor * result;
	}

	PerturbativeQuantity lepton_pair_cross_section_xQ2(const TRFKinematics &kinematics) const {
		const double x = kinematics.x;
		const double Q2 = kinematics.Q2;

		std::vector<double> z_mins(ff1.decays.size());
		std::transform(ff1.decays.begin(), ff1.decays.end(), z_mins.begin(), [&kinematics](const auto &decay) {
			return SIDISFunctions::Helper::compute_z_min(kinematics, decay);
		});
		const double z_min = *std::min_element(z_mins.begin(), z_mins.end());

		if (z_min > 1.0) { return PerturbativeQuantity{0.0, 0.0, 0.0}; }

		double alpha_s = compute_alpha_s(kinematics);
		const double nlo_coefficient = alpha_s / (2 * std::numbers::pi);
		const double nnlo_coefficient = std::pow(nlo_coefficient, 2);

		const double renormalization_scale = renormalization_scale_function(kinematics);
		const double factorization_scale = factorization_scale_function(kinematics);
		const double fragmentation_scale = fragmentation_scale_function(kinematics);

		const double renormalization_scale_log = renormalization_scale == Q2 ? 0 : std::log(Q2 / renormalization_scale);
		const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);
		const double fragmentation_scale_log = fragmentation_scale == Q2 ? 0 : std::log(Q2 / fragmentation_scale);

		pdf1.evaluate(x, factorization_scale);

		SIDISFunctions::Parameters<PDFInterface, FFInterface, DecayFunction> params {
			pdf1, ff1, pdf2, ff2,
			flavors,
			nlo_coefficient,
			process, kinematics,
			0.0,
			renormalization_scale, factorization_scale, fragmentation_scale,
			renormalization_scale_log, factorization_scale_log, fragmentation_scale_log
		};

		Integrator lo_integrator([&](double input[], [[maybe_unused]] std::size_t dim, void* params_in) {
			return SIDISFunctions::cross_section<PDFInterface, FFInterface, DecayFunction>(input, params_in, 
				SIDISFunctions::F2::LO::integrand, SIDISFunctions::FL::LO::integrand, SIDISFunctions::F3::LO::integrand,
				false, false, true
			);
		}, {z_min}, {1}, integration_parameters, &params);
		const auto lo_result = lo_integrator.integrate();
		const double lo = lo_result.value;

		double nlo = 0.0;
		if (order >= PerturbativeOrder::NLO) {
			Integrator nlo_integrator([&](double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
				return SIDISFunctions::cross_section<PDFInterface, FFInterface, DecayFunction>(input, params_in, 
					use_nlp_nlo ? SIDISFunctions::F2::NLO_NLP::total_integrand : SIDISFunctions::F2::NLO::total_integrand, 
					use_nlp_nlo ? SIDISFunctions::FL::NLO_NLP::total_integrand : SIDISFunctions::FL::NLO::total_integrand, 
					use_nlp_nlo ? SIDISFunctions::F3::NLO_NLP::total_integrand : SIDISFunctions::F3::NLO::total_integrand,
					true, true, true
				);
			}, {x, z_min, z_min}, {1.0, 1.0, 1.0}, integration_parameters, &params);
			const auto nlo_result = nlo_integrator.integrate();
			const double nlo_integral = nlo_result.value;

			nlo = nlo_coefficient * nlo_integral;
		}

		double nnlo = 0.0;
		if (order >= PerturbativeOrder::NNLO) {
			Integrator nnlo_integrator([&](double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
				return SIDISFunctions::cross_section<PDFInterface, FFInterface, DecayFunction>(input, params_in, 
					SIDISFunctions::F2::NNLO_NLP::total_integrand, SIDISFunctions::FL::NNLO_NLP::total_integrand, SIDISFunctions::F3::NNLO_NLP::total_integrand,
					true, true, true
				);
			}, {x, z_min, z_min}, {1.0, 1.0, 1.0}, integration_parameters, &params);
			const auto nnlo_result = nnlo_integrator.integrate();
			const double nnlo_integral = nnlo_result.value;

			nnlo = nnlo_coefficient * nnlo_integral;
		}

		const PerturbativeQuantity result = PerturbativeQuantity {lo, lo + nlo, lo + nlo + nnlo};
		const double prefactor = use_modified_cross_section_prefactor 
									? CommonFunctions::cross_section_modified_prefactor(kinematics) 
									: CommonFunctions::cross_section_prefactor(kinematics);
		return prefactor * result;
	}

	PerturbativeQuantity x_integrated_lepton_pair_cross_section(const TRFKinematics &placeholder_kinematics, const double Q2) const {
		// input = [z, x]
		const auto lo_integrand = [&](double input[]) {
			const double x = input[1];

			const double target_mass = placeholder_kinematics.target_mass;
			const double E_beam = placeholder_kinematics.E_beam;
			const double y = Q2 / (2.0 * x * target_mass * E_beam);

			const double y_max = 1.0 - primary_muon_min_energy / E_beam;
			if (y > y_max) { return 0.0; }

			const TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, target_mass, placeholder_kinematics.projectile_mass);

			std::vector<double> z_mins(ff1.decays.size());
			std::transform(ff1.decays.begin(), ff1.decays.end(), z_mins.begin(), [&kinematics](const auto &decay) {
				return SIDISFunctions::Helper::compute_z_min(kinematics, decay);
			});
			const double z_min = *std::min_element(z_mins.begin(), z_mins.end());

			if (input[0] < z_min) { return 0.0; }

			const double renormalization_scale = renormalization_scale_function(kinematics);
			const double factorization_scale = factorization_scale_function(kinematics);
			const double fragmentation_scale = fragmentation_scale_function(kinematics);

			const double renormalization_scale_log = renormalization_scale == Q2 ? 0 : std::log(Q2 / renormalization_scale);
			const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);
			const double fragmentation_scale_log = fragmentation_scale == Q2 ? 0 : std::log(Q2 / fragmentation_scale);

			pdf1.evaluate(x, factorization_scale);

			SIDISFunctions::Parameters<PDFInterface, FFInterface, DecayFunction> params {
				pdf1, ff1, pdf2, ff2,
				flavors,
				0.0,
				process, kinematics,
				0.0,
				renormalization_scale, factorization_scale, fragmentation_scale,
				renormalization_scale_log, factorization_scale_log, fragmentation_scale_log
			};

			const double differential_cross_section = SIDISFunctions::cross_section<PDFInterface, FFInterface, DecayFunction>(input, &params, 
				SIDISFunctions::F2::LO::integrand, SIDISFunctions::FL::LO::integrand, SIDISFunctions::F3::LO::integrand,
				false, false, true
			);

			const double prefactor = use_modified_cross_section_prefactor 
							? CommonFunctions::cross_section_modified_prefactor(kinematics) 
							: CommonFunctions::cross_section_prefactor(kinematics);

			return prefactor * differential_cross_section;
		};

		// input = [xi, xip, z, x]
		const auto nlo_integrand = [&](double input[]) {
			const double x = input[3];

			const double target_mass = placeholder_kinematics.target_mass;
			const double E_beam = placeholder_kinematics.E_beam;
			const double y = Q2 / (2.0 * x * target_mass * E_beam);

			const double y_max = 1.0 - primary_muon_min_energy / E_beam;
			if (y > y_max) { return 0.0; }

			const TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, target_mass, placeholder_kinematics.projectile_mass);

			std::vector<double> z_mins(ff1.decays.size());
			std::transform(ff1.decays.begin(), ff1.decays.end(), z_mins.begin(), [&kinematics](const auto &decay) {
				return SIDISFunctions::Helper::compute_z_min(kinematics, decay);
			});
			const double z_min = *std::min_element(z_mins.begin(), z_mins.end());

			if (input[1] < z_min || input[2] < z_min) { return 0.0; }

			double alpha_s = compute_alpha_s(kinematics);
			const double nlo_coefficient = alpha_s / (2 * std::numbers::pi);

			const double renormalization_scale = renormalization_scale_function(kinematics);
			const double factorization_scale = factorization_scale_function(kinematics);
			const double fragmentation_scale = fragmentation_scale_function(kinematics);

			const double renormalization_scale_log = renormalization_scale == Q2 ? 0 : std::log(Q2 / renormalization_scale);
			const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);
			const double fragmentation_scale_log = fragmentation_scale == Q2 ? 0 : std::log(Q2 / fragmentation_scale);

			pdf1.evaluate(x, factorization_scale);

			SIDISFunctions::Parameters<PDFInterface, FFInterface, DecayFunction> params {
				pdf1, ff1, pdf2, ff2,
				flavors,
				0.0,
				process, kinematics,
				0.0,
				renormalization_scale, factorization_scale, fragmentation_scale,
				renormalization_scale_log, factorization_scale_log, fragmentation_scale_log
			};

			const double differential_cross_section = SIDISFunctions::cross_section<PDFInterface, FFInterface, DecayFunction>(input, &params, 
				use_nlp_nlo ? SIDISFunctions::F2::NLO_NLP::total_integrand : SIDISFunctions::F2::NLO::total_integrand, 
				use_nlp_nlo ? SIDISFunctions::FL::NLO_NLP::total_integrand : SIDISFunctions::FL::NLO::total_integrand, 
				use_nlp_nlo ? SIDISFunctions::F3::NLO_NLP::total_integrand : SIDISFunctions::F3::NLO::total_integrand,
				true, true, true
			);

			const double prefactor = use_modified_cross_section_prefactor 
							? CommonFunctions::cross_section_modified_prefactor(kinematics) 
							: CommonFunctions::cross_section_prefactor(kinematics);

			return prefactor * nlo_coefficient * differential_cross_section;
		};

		// input = [xi, xip, z, x]
		const auto nnlo_integrand = [&](double input[]) {
			const double x = input[3];

			const double target_mass = placeholder_kinematics.target_mass;
			const double E_beam = placeholder_kinematics.E_beam;
			const double y = Q2 / (2.0 * x * target_mass * E_beam);

			const double y_max = 1.0 - primary_muon_min_energy / E_beam;
			if (y > y_max) { return 0.0; }

			const TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, target_mass, placeholder_kinematics.projectile_mass);

			std::vector<double> z_mins(ff1.decays.size());
			std::transform(ff1.decays.begin(), ff1.decays.end(), z_mins.begin(), [&kinematics](const auto &decay) {
				return SIDISFunctions::Helper::compute_z_min(kinematics, decay);
			});
			const double z_min = *std::min_element(z_mins.begin(), z_mins.end());

			if (input[1] < z_min || input[2] < z_min) { return 0.0; }

			double alpha_s = compute_alpha_s(kinematics);
			const double nnlo_coefficient = std::pow(alpha_s / (2 * std::numbers::pi), 2);

			const double renormalization_scale = renormalization_scale_function(kinematics);
			const double factorization_scale = factorization_scale_function(kinematics);
			const double fragmentation_scale = fragmentation_scale_function(kinematics);

			const double renormalization_scale_log = renormalization_scale == Q2 ? 0 : std::log(Q2 / renormalization_scale);
			const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);
			const double fragmentation_scale_log = fragmentation_scale == Q2 ? 0 : std::log(Q2 / fragmentation_scale);

			pdf1.evaluate(x, factorization_scale);

			SIDISFunctions::Parameters<PDFInterface, FFInterface, DecayFunction> params {
				pdf1, ff1, pdf2, ff2,
				flavors,
				0.0,
				process, kinematics,
				0.0,
				renormalization_scale, factorization_scale, fragmentation_scale,
				renormalization_scale_log, factorization_scale_log, fragmentation_scale_log
			};

			const double differential_cross_section = SIDISFunctions::cross_section<PDFInterface, FFInterface, DecayFunction>(input, &params, 
				SIDISFunctions::F2::NNLO_NLP::total_integrand, SIDISFunctions::FL::NNLO_NLP::total_integrand, SIDISFunctions::F3::NNLO_NLP::total_integrand,
				true, true, true
			);

			const double prefactor = use_modified_cross_section_prefactor 
							? CommonFunctions::cross_section_modified_prefactor(kinematics) 
							: CommonFunctions::cross_section_prefactor(kinematics);

			return prefactor * nnlo_coefficient * differential_cross_section;
		};

		const double target_mass = placeholder_kinematics.target_mass;
		const double E_beam = placeholder_kinematics.E_beam;

		const double x_min = Q2 / (2.0 * target_mass * E_beam);

		Integrator lo_integrator([&](double input[], std::size_t dim, [[maybe_unused]] void *params_in) {
			const double x = input[0];

			double *scaled_input = new double[dim];
			scaled_input[1] = x; // x
			scaled_input[0] = input[1]; // z

			const double result = lo_integrand(scaled_input);

			delete[] scaled_input;
			return result;
		}, {x_min /* x */, 0.0 /* z */}, {1.0, 1.0}, integration_parameters, nullptr);

		const auto lo_result = lo_integrator.integrate();
		const double lo = lo_result.value;

		double nlo = 0.0;
		if (order >= PerturbativeOrder::NLO) {
			Integrator nlo_integrator([&](double input[], std::size_t dim, [[maybe_unused]] void *params_in) {
				const double x = input[0];

				double *scaled_input = new double[dim];
				scaled_input[0] = x + (1.0 - x) * input[1]; // xi
				scaled_input[1] = input[2]; // xip
				scaled_input[2] = input[3]; // z
				scaled_input[3] = x; // x

				const double result = (1.0 - x) * nlo_integrand(scaled_input);

				delete[] scaled_input;
				return result;
			}, {x_min /* x */, 0.0 /* scaled xi */, 0.0 /* xip */, 0.0 /* z */}, {1.0, 1.0, 1.0, 1.0}, integration_parameters, nullptr);

			const auto nlo_result = nlo_integrator.integrate();
			nlo = nlo_result.value;
		}

		double nnlo = 0.0;
		if (order >= PerturbativeOrder::NNLO) {
			Integrator nnlo_integrator([&](double input[], std::size_t dim, [[maybe_unused]] void *params_in) {
				const double x = input[0];

				double *scaled_input = new double[dim];
				scaled_input[0] = x + (1.0 - x) * input[1]; // xi
				scaled_input[1] = input[2]; // xip
				scaled_input[2] = input[3]; // z
				scaled_input[3] = x; // x

				const double result = (1.0 - x) * nnlo_integrand(scaled_input);

				delete[] scaled_input;
				return result;
			}, {x_min /* x */, 0.0 /* scaled xi */, 0.0 /* xip */, 0.0 /* z */}, {1.0, 1.0, 1.0, 1.0}, integration_parameters, nullptr);
			const auto nnlo_result = nnlo_integrator.integrate();
			nnlo = nnlo_result.value;
		}

		const PerturbativeQuantity result = PerturbativeQuantity {lo, lo + nlo, lo + nlo + nnlo};
		return result;
	}
	PerturbativeQuantity integrated_lepton_pair_cross_section(const TRFKinematics &placeholder_kinematics, const double Q2_min) const {
		// input = [z, x, Q2]
		const auto lo_integrand = [&](double input[]) {
			const double x = input[1];
			const double Q2 = input[2];

			const double target_mass = placeholder_kinematics.target_mass;
			const double E_beam = placeholder_kinematics.E_beam;
			const double y = Q2 / (2.0 * x * target_mass * E_beam);

			const double y_max = 1.0 - primary_muon_min_energy / E_beam;
			if (y > y_max) { return 0.0; }

			const TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, target_mass, placeholder_kinematics.projectile_mass);

			std::vector<double> z_mins(ff1.decays.size());
			std::transform(ff1.decays.begin(), ff1.decays.end(), z_mins.begin(), [&kinematics](const auto &decay) {
				return SIDISFunctions::Helper::compute_z_min(kinematics, decay);
			});
			const double z_min = *std::min_element(z_mins.begin(), z_mins.end());

			if (input[0] < z_min) { return 0.0; }

			const double renormalization_scale = renormalization_scale_function(kinematics);
			const double factorization_scale = factorization_scale_function(kinematics);
			const double fragmentation_scale = fragmentation_scale_function(kinematics);

			const double renormalization_scale_log = renormalization_scale == Q2 ? 0 : std::log(Q2 / renormalization_scale);
			const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);
			const double fragmentation_scale_log = fragmentation_scale == Q2 ? 0 : std::log(Q2 / fragmentation_scale);

			pdf1.evaluate(x, factorization_scale);

			SIDISFunctions::Parameters<PDFInterface, FFInterface, DecayFunction> params {
				pdf1, ff1, pdf2, ff2,
				flavors,
				0.0,
				process, kinematics,
				0.0,
				renormalization_scale, factorization_scale, fragmentation_scale,
				renormalization_scale_log, factorization_scale_log, fragmentation_scale_log
			};

			const double differential_cross_section = SIDISFunctions::cross_section<PDFInterface, FFInterface, DecayFunction>(input, &params, 
				SIDISFunctions::F2::LO::integrand, SIDISFunctions::FL::LO::integrand, SIDISFunctions::F3::LO::integrand,
				false, false, true
			);

			const double prefactor = use_modified_cross_section_prefactor 
							? CommonFunctions::cross_section_modified_prefactor(kinematics) 
							: CommonFunctions::cross_section_prefactor(kinematics);

			return prefactor * differential_cross_section;
		};

		// input = [xi, xip, z, x, Q2]
		const auto nlo_integrand = [&](double input[]) {
			const double x = input[3];
			const double Q2 = input[4];

			const double target_mass = placeholder_kinematics.target_mass;
			const double E_beam = placeholder_kinematics.E_beam;
			const double y = Q2 / (2.0 * x * target_mass * E_beam);

			const double y_max = 1.0 - primary_muon_min_energy / E_beam;
			if (y > y_max) { return 0.0; }

			const TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, target_mass, placeholder_kinematics.projectile_mass);

			std::vector<double> z_mins(ff1.decays.size());
			std::transform(ff1.decays.begin(), ff1.decays.end(), z_mins.begin(), [&kinematics](const auto &decay) {
				return SIDISFunctions::Helper::compute_z_min(kinematics, decay);
			});
			const double z_min = *std::min_element(z_mins.begin(), z_mins.end());

			if (input[1] < z_min || input[2] < z_min) { return 0.0; }

			double alpha_s = compute_alpha_s(kinematics);
			const double nlo_coefficient = alpha_s / (2 * std::numbers::pi);

			const double renormalization_scale = renormalization_scale_function(kinematics);
			const double factorization_scale = factorization_scale_function(kinematics);
			const double fragmentation_scale = fragmentation_scale_function(kinematics);

			const double renormalization_scale_log = renormalization_scale == Q2 ? 0 : std::log(Q2 / renormalization_scale);
			const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);
			const double fragmentation_scale_log = fragmentation_scale == Q2 ? 0 : std::log(Q2 / fragmentation_scale);

			pdf1.evaluate(x, factorization_scale);

			SIDISFunctions::Parameters<PDFInterface, FFInterface, DecayFunction> params {
				pdf1, ff1, pdf2, ff2,
				flavors,
				0.0,
				process, kinematics,
				0.0,
				renormalization_scale, factorization_scale, fragmentation_scale,
				renormalization_scale_log, factorization_scale_log, fragmentation_scale_log
			};

			const double differential_cross_section = SIDISFunctions::cross_section<PDFInterface, FFInterface, DecayFunction>(input, &params, 
				use_nlp_nlo ? SIDISFunctions::F2::NLO_NLP::total_integrand : SIDISFunctions::F2::NLO::total_integrand, 
				use_nlp_nlo ? SIDISFunctions::FL::NLO_NLP::total_integrand : SIDISFunctions::FL::NLO::total_integrand, 
				use_nlp_nlo ? SIDISFunctions::F3::NLO_NLP::total_integrand : SIDISFunctions::F3::NLO::total_integrand,
				true, true, true
			);

			const double prefactor = use_modified_cross_section_prefactor 
							? CommonFunctions::cross_section_modified_prefactor(kinematics) 
							: CommonFunctions::cross_section_prefactor(kinematics);

			return prefactor * nlo_coefficient * differential_cross_section;
		};

		// input = [xi, xip, z, x, Q2]
		const auto nnlo_integrand = [&](double input[]) {
			const double x = input[3];
			const double Q2 = input[4];

			const double target_mass = placeholder_kinematics.target_mass;
			const double E_beam = placeholder_kinematics.E_beam;
			const double y = Q2 / (2.0 * x * target_mass * E_beam);

			const double y_max = 1.0 - primary_muon_min_energy / E_beam;
			if (y > y_max) { return 0.0; }

			const TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, target_mass, placeholder_kinematics.projectile_mass);

			std::vector<double> z_mins(ff1.decays.size());
			std::transform(ff1.decays.begin(), ff1.decays.end(), z_mins.begin(), [&kinematics](const auto &decay) {
				return SIDISFunctions::Helper::compute_z_min(kinematics, decay);
			});
			const double z_min = *std::min_element(z_mins.begin(), z_mins.end());

			if (input[1] < z_min || input[2] < z_min) { return 0.0; }

			double alpha_s = compute_alpha_s(kinematics);
			const double nnlo_coefficient = std::pow(alpha_s / (2 * std::numbers::pi), 2);

			const double renormalization_scale = renormalization_scale_function(kinematics);
			const double factorization_scale = factorization_scale_function(kinematics);
			const double fragmentation_scale = fragmentation_scale_function(kinematics);

			const double renormalization_scale_log = renormalization_scale == Q2 ? 0 : std::log(Q2 / renormalization_scale);
			const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);
			const double fragmentation_scale_log = fragmentation_scale == Q2 ? 0 : std::log(Q2 / fragmentation_scale);

			pdf1.evaluate(x, factorization_scale);

			SIDISFunctions::Parameters<PDFInterface, FFInterface, DecayFunction> params {
				pdf1, ff1, pdf2, ff2,
				flavors,
				0.0,
				process, kinematics,
				0.0,
				renormalization_scale, factorization_scale, fragmentation_scale,
				renormalization_scale_log, factorization_scale_log, fragmentation_scale_log
			};

			const double differential_cross_section = SIDISFunctions::cross_section<PDFInterface, FFInterface, DecayFunction>(input, &params, 
				SIDISFunctions::F2::NNLO_NLP::total_integrand, SIDISFunctions::FL::NNLO_NLP::total_integrand, SIDISFunctions::F3::NNLO_NLP::total_integrand,
				true, true, true
			);

			const double prefactor = use_modified_cross_section_prefactor 
							? CommonFunctions::cross_section_modified_prefactor(kinematics) 
							: CommonFunctions::cross_section_prefactor(kinematics);

			return prefactor * nnlo_coefficient * differential_cross_section;
		};

		const double target_mass = placeholder_kinematics.target_mass;
		const double E_beam = placeholder_kinematics.E_beam;

		const double x_min = Q2_min / (2.0 * target_mass * E_beam);

		Integrator lo_integrator([&](double input[], std::size_t dim, [[maybe_unused]] void *params_in) {
			const double x = input[0];
			const double Q2_max = 2.0 * x * target_mass * E_beam;

			double *scaled_input = new double[dim];
			scaled_input[1] = x; // x
			scaled_input[2] = Q2_min + (Q2_max - Q2_min) * input[1]; // Q^2
			scaled_input[0] = input[2]; // z

			const double result = (Q2_max - Q2_min) * lo_integrand(scaled_input);

			delete[] scaled_input;
			return result;
		}, {x_min /* x */, 0.0 /* scaled Q^2 */, 0.0 /* z */}, {1.0, 1.0, 1.0}, integration_parameters, nullptr);

		const auto lo_result = lo_integrator.integrate();
		const double lo = lo_result.value;

		double nlo = 0.0;
		if (order >= PerturbativeOrder::NLO) {
			Integrator nlo_integrator([&](double input[], std::size_t dim, [[maybe_unused]] void *params_in) {
				const double x = input[0];
				const double Q2_max = 2.0 * x * target_mass * E_beam;

				double *scaled_input = new double[dim];
				scaled_input[0] = x + (1.0 - x) * input[2]; // xi
				scaled_input[1] = input[3]; // xip
				scaled_input[2] = input[4]; // z
				scaled_input[3] = x; // x
				scaled_input[4] = Q2_min + (Q2_max - Q2_min) * input[1]; // Q^2

				const double result = (Q2_max - Q2_min) * (1.0 - x) * nlo_integrand(scaled_input);

				delete[] scaled_input;
				return result;
			}, {x_min /* x */, 0.0 /* scaled Q^2 */, 0.0 /* scaled xi */, 0.0 /* xip */, 0.0 /* z */}, {1.0, 1.0, 1.0, 1.0, 1.0}, integration_parameters, nullptr);

			const auto nlo_result = nlo_integrator.integrate();
			nlo = nlo_result.value;
		}

		double nnlo = 0.0;
		if (order >= PerturbativeOrder::NNLO) {
			Integrator nnlo_integrator([&](double input[], std::size_t dim, [[maybe_unused]] void *params_in) {
				const double x = input[0];
				const double Q2_max = 2.0 * x * target_mass * E_beam;

				double *scaled_input = new double[dim];
				scaled_input[0] = x + (1.0 - x) * input[2]; // xi
				scaled_input[1] = input[3]; // xip
				scaled_input[2] = input[4]; // z
				scaled_input[3] = x; // x
				scaled_input[4] = Q2_min + (Q2_max - Q2_min) * input[1]; // Q^2

				const double result = (Q2_max - Q2_min) * (1.0 - x) * nnlo_integrand(scaled_input);

				delete[] scaled_input;
				return result;
			}, {x_min /* x */, 0.0 /* scaled Q^2 */, 0.0 /* scaled xi */, 0.0 /* xip */, 0.0 /* z */}, {1.0, 1.0, 1.0, 1.0, 1.0}, integration_parameters, nullptr);
			const auto nnlo_result = nnlo_integrator.integrate();
			nnlo = nnlo_result.value;
		}

		const PerturbativeQuantity result = PerturbativeQuantity {lo, lo + nlo, lo + nlo + nnlo};
		return result;
	}
};

#endif