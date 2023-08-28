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
		const bool _use_nlp_nlo
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
	use_nlp_nlo(_use_nlp_nlo) { }

	double compute_alpha_s(const TRFKinematics &kinematics) const {
		const double renormalization_scale = renormalization_scale_function(kinematics);
		return pdf1.alpha_s(renormalization_scale);
	}
	
	template <typename LO, typename NLO, typename NNLO>
	PerturbativeQuantity structure_function(
		const double z, const TRFKinematics &kinematics, const LO lo_integrand, const NLO nlo_integrand, const NNLO nnlo_integrand, const int sign) const {
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

		const double lo = SIDISFunctions::construct<PDFInterface, FFInterface>({}, &params, lo_integrand, false, false, false, sign);

		double nlo = 0.0;
		if (order >= PerturbativeOrder::NLO) {
			Integrator nlo_integrator([nlo_integrand, sign](double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
				return SIDISFunctions::construct<PDFInterface, FFInterface>(input, params_in, nlo_integrand, true, true, false, sign);
			}, {x, z}, {1, 1}, integration_parameters, &params);
			const auto nlo_result = nlo_integrator.integrate();
			const double nlo_value = nlo_result.value;

			nlo = nlo_coefficient * nlo_value;
		}

		double nnlo = 0.0;
		if (order >= PerturbativeOrder::NNLO) {
			Integrator nnlo_integrator([nnlo_integrand, sign](double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
				return SIDISFunctions::construct<PDFInterface, FFInterface>(input, params_in, nnlo_integrand, true, true, false, sign);
			}, {x, z}, {1, 1}, integration_parameters, &params);
			const auto nnlo_result = nnlo_integrator.integrate();
			const double nnlo_value = nnlo_result.value;

			nnlo = nnlo_coefficient * nnlo_value;
		}

		return PerturbativeQuantity {lo, lo + nlo, lo + nlo + nnlo};
	}

	PerturbativeQuantity F2(const double z, const TRFKinematics &kinematics) const {
		return structure_function(
			z, kinematics, SIDISFunctions::F2::LO::integrand, SIDISFunctions::F2::NLO::total_integrand, SIDISFunctions::F2::NNLO_NLP::total_integrand, 1
		);
	}
	PerturbativeQuantity FL(const double z, const TRFKinematics &kinematics) const {
		return structure_function(
			z, kinematics, SIDISFunctions::FL::LO::integrand, SIDISFunctions::FL::NLO::total_integrand, SIDISFunctions::FL::NNLO_NLP::total_integrand, 1
		);
	}
	PerturbativeQuantity F3(const double z, const TRFKinematics &kinematics) const {
		return structure_function(
			z, kinematics, SIDISFunctions::F3::LO::integrand, SIDISFunctions::F3::NLO::total_integrand, SIDISFunctions::F3::NNLO_NLP::total_integrand, -1
		);
	}

	PerturbativeQuantity compute_structure_function(
		const StructureFunction structure_function, const double z, const TRFKinematics &kinematics) const {
		switch (structure_function) {
		case StructureFunction::F2: return F2(z, kinematics);
		case StructureFunction::FL: return FL(z, kinematics);
		case StructureFunction::F3: return F3(z, kinematics);
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

	PerturbativeQuantity lepton_pair_cross_section_xQ2(const TRFKinematics &kinematics) const {
		const double x = kinematics.x;
		const double Q2 = kinematics.Q2;

		std::vector<double> z_mins(ff1.decays.size());
		std::transform(ff1.decays.begin(), ff1.decays.end(), z_mins.begin(), [&kinematics](const auto &decay) {
			return SIDISFunctions::Helper::compute_z_min(kinematics, decay);
		});
		const double z_min = *std::min_element(z_mins.begin(), z_mins.end());

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
};

#endif