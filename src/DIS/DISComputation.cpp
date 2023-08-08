#ifndef DIS_COMPUTATION_H
#define DIS_COMPUTATION_H

#include "Integration/Integrator.cpp"
#include "Common/Constants.cpp"
#include "Utility/Utility.cpp"
#include "Common/PerturbativeQuantity.cpp"
#include "Common/Process.cpp"
#include "Common/StructureFunction.cpp"
#include "DIS/DISFunctions.cpp"
#include "Common/CommonFunctions.cpp"
#include <optional>
#include "PDF/PDFConcept.cpp"
#include "Common/ScaleDependence.cpp"

template <
	PDFConcept PDFInterface, 
	ScaleDependence::Concept RenormalizationScale,
	ScaleDependence::Concept FactorizationScale
>
class DISComputation {
	public:
	double sqrt_s;
	double s;

	FlavorInfo flavors;

	const PDFInterface pdf1;
	const PDFInterface pdf2;

	const IntegrationParameters integration_parameters;

	const Process process;
	const bool momentum_fraction_mass_corrections;

	const ScaleDependence::Function<RenormalizationScale> renormalization_scale_function;
	const ScaleDependence::Function<FactorizationScale> factorization_scale_function;
	const bool use_modified_cross_section_prefactor;

	DISComputation (
		const double _sqrt_s, 
		const FlavorVector _active_flavors,
		const std::array<double, TOTAL_FLAVORS> _flavor_masses, 
		const PDFInterface _pdf,
		const IntegrationParameters _integration_parameters,
		const Process _process,
		const bool _momentum_fraction_mass_corrections,
		const ScaleDependence::Function<RenormalizationScale> _renormalization_scale_function,
		const ScaleDependence::Function<FactorizationScale> _factorization_scale_function,
		const bool _use_modified_cross_section_prefactor
	) : sqrt_s(_sqrt_s), 
	s(_sqrt_s * _sqrt_s), 
	flavors(_active_flavors, _flavor_masses), 
	pdf1(_pdf), 
	pdf2(_pdf), 
	integration_parameters(_integration_parameters), 
	process(_process),
	momentum_fraction_mass_corrections(_momentum_fraction_mass_corrections),
	renormalization_scale_function(_renormalization_scale_function),
	factorization_scale_function(_factorization_scale_function),
	use_modified_cross_section_prefactor(_use_modified_cross_section_prefactor) { }

	double compute_alpha_s(const TRFKinematics &kinematics) const {
		const double renormalization_scale = renormalization_scale_function(kinematics);
		return pdf1.alpha_s(renormalization_scale);
	}

	PerturbativeQuantity F2(const TRFKinematics &kinematics) const {	
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;

		double alpha_s = compute_alpha_s(kinematics);
		double nlo_coefficient = alpha_s / (2 * std::numbers::pi);

		const double factorization_scale = factorization_scale_function(kinematics);
		const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);

		pdf1.evaluate(x, factorization_scale);

		DISFunctions::Parameters<PDFInterface> params {
			pdf1, pdf2,
			flavors,
			nlo_coefficient,
			process, kinematics,
			factorization_scale,
			factorization_scale_log
		};
		
		const double lo = DISFunctions::Evaluation::construct<PDFInterface>({}, &params, DISFunctions::F2::LO::integrand, false, 1);

		Integrator nlo_integrator([](double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
			return DISFunctions::Evaluation::construct<PDFInterface>(input, params_in, DISFunctions::F2::NLO::integrand, true, 1);
		}, {x}, {1.0}, integration_parameters, &params);
		auto nlo_result = nlo_integrator.integrate();
		const double nlo_value = nlo_result.value;

		const double nlo_delta = DISFunctions::Evaluation::construct<PDFInterface>({}, &params, DISFunctions::F2::NLO::delta, false, 1);

		const double nlo = nlo_coefficient * (nlo_value + nlo_delta);

		return PerturbativeQuantity {lo, lo + nlo};
	}

	PerturbativeQuantity FL(const TRFKinematics &kinematics) const {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;

		double alpha_s = compute_alpha_s(kinematics);
		double nlo_coefficient = alpha_s / (2 * std::numbers::pi);

		const double factorization_scale = factorization_scale_function(kinematics);
		const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);

		DISFunctions::Parameters<PDFInterface> params {
			pdf1, pdf2,
			flavors,
			nlo_coefficient,
			process, kinematics,
			factorization_scale,
			factorization_scale_log
		};

		// TODO: LO etc mass correction integrals

		Integrator nlo_integrator([](double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
			return DISFunctions::Evaluation::construct<PDFInterface>(input, params_in, DISFunctions::FL::NLO::integrand, true, 1);
		}, {x}, {1.0}, integration_parameters, &params);
		auto nlo_result = nlo_integrator.integrate();
		const double nlo_value = nlo_result.value;

		const double nlo = nlo_coefficient * nlo_value;

		return PerturbativeQuantity {0.0, nlo};
	}
	PerturbativeQuantity F3(const TRFKinematics &kinematics) const {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;

		double alpha_s = compute_alpha_s(kinematics);
		double nlo_coefficient = alpha_s / (2 * std::numbers::pi);
		
		const double factorization_scale = factorization_scale_function(kinematics);
		const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);

		pdf1.evaluate(x, factorization_scale);

		DISFunctions::Parameters<PDFInterface> params {
			pdf1, pdf2,
			flavors,
			nlo_coefficient,
			process, kinematics,
			factorization_scale,
			factorization_scale_log
		};

		const double lo = DISFunctions::Evaluation::construct<PDFInterface>({}, &params, DISFunctions::F3::LO::integrand, false, -1);

		Integrator nlo_integrator([](double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
			return DISFunctions::Evaluation::construct<PDFInterface>(input, params_in, DISFunctions::F3::NLO::integrand, true, -1);
		}, {x}, {1.0}, integration_parameters, &params);
		auto nlo_result = nlo_integrator.integrate();
		const double nlo_value = nlo_result.value;

		const double nlo_delta = DISFunctions::Evaluation::construct<PDFInterface>({}, &params, DISFunctions::F3::NLO::delta, false, -1);

		const double nlo = nlo_coefficient * (nlo_value + nlo_delta);

		return PerturbativeQuantity {lo, lo + nlo};
	}
	constexpr PerturbativeQuantity compute_structure_function(const StructureFunction F, const TRFKinematics &kinematics) const {
		switch (F) {
		case StructureFunction::F2: return F2(kinematics);
		case StructureFunction::FL: return FL(kinematics);
		case StructureFunction::F3: return F3(kinematics);
		}
	}

	PerturbativeQuantity differential_cross_section_xQ2_direct(const TRFKinematics &kinematics) const {
		if (!kinematics.is_valid()) { return PerturbativeQuantity {0.0, 0.0}; }

		const double prefactor = use_modified_cross_section_prefactor 
									? CommonFunctions::cross_section_modified_prefactor(kinematics) 
									: CommonFunctions::cross_section_prefactor(kinematics);

		const PerturbativeQuantity f2 = F2(kinematics);
		const PerturbativeQuantity fL = FL(kinematics);
		const PerturbativeQuantity f3 = F3(kinematics);

		const PerturbativeQuantity result = prefactor * CommonFunctions::make_cross_section_variable(kinematics, process, f2, fL, f3);

		return result;
	}

	PerturbativeQuantity differential_cross_section_xQ2_indirect(const TRFKinematics &kinematics) const {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;

		double alpha_s = compute_alpha_s(kinematics);
		double nlo_coefficient = alpha_s / (2 * std::numbers::pi);

		const double factorization_scale = factorization_scale_function(kinematics);
		const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);

		pdf1.evaluate(x, factorization_scale);

		DISFunctions::Parameters<PDFInterface> params {
			pdf1, pdf2,
			flavors,
			nlo_coefficient,
			process, kinematics,
			factorization_scale,
			factorization_scale_log
		};

		const double lo = DISFunctions::Evaluation::cross_section<PDFInterface>({}, &params, 
			DISFunctions::F2::LO::integrand, DISFunctions::FL::LO::integrand, DISFunctions::F3::LO::integrand,
			false
		);

		Integrator nlo_integrator([](double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
			return DISFunctions::Evaluation::cross_section<PDFInterface>(input, params_in, 
				DISFunctions::F2::NLO::integrand, DISFunctions::FL::NLO::integrand, DISFunctions::F3::NLO::integrand,
				true
			);
		}, {x}, {1.0}, integration_parameters, &params);
		const auto nlo_result = nlo_integrator.integrate();
		const double nlo_integral = nlo_result.value;

		const double delta_contribution = DISFunctions::Evaluation::cross_section<PDFInterface>({}, &params,
			DISFunctions::F2::NLO::delta, DISFunctions::FL::NLO::delta, DISFunctions::F3::NLO::delta, 
			false
		);

		const double nlo = nlo_coefficient * (nlo_integral + delta_contribution);

		const PerturbativeQuantity result = PerturbativeQuantity {lo, lo + nlo};
		const double prefactor = use_modified_cross_section_prefactor 
									? CommonFunctions::cross_section_modified_prefactor(kinematics) 
									: CommonFunctions::cross_section_prefactor(kinematics);
		return prefactor * result;
	}
};

#endif