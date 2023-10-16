#ifndef DIS_COMPUTATION_H
#define DIS_COMPUTATION_H

#include <optional>

#include "Common/Constants.cpp"
#include "Common/PerturbativeQuantity.cpp"
#include "Common/Process.cpp"
#include "Common/StructureFunction.cpp"
#include "Common/CommonFunctions.cpp"
#include "Common/ScaleDependence.cpp"

#include "Integration/Integrator.cpp"

#include "Utility/Utility.cpp"

#include "DIS/DISFunctions.cpp"

#include "PDF/PDFConcept.cpp"

template <
	is_pdf_interface PDFInterface, 
	is_scale_dependence RenormalizationScale,
	is_scale_dependence FactorizationScale
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

	template <typename LO, typename NLO, typename Delta>
	PerturbativeQuantity structure_function(
		const TRFKinematics &kinematics, 
		const LO lo_integrand, const NLO nlo_integrand, const Delta delta_integrand,
		const int sign) const {
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
		
		const double lo = DISFunctions::construct<PDFInterface>({}, &params, lo_integrand, false, sign);

		Integrator nlo_integrator([nlo_integrand, sign](double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
			return DISFunctions::construct<PDFInterface>(input, params_in, nlo_integrand, true, sign);
		}, {x}, {1.0}, integration_parameters, &params);
		auto nlo_result = nlo_integrator.integrate();
		const double nlo_value = nlo_result.value;

		const double nlo_delta = DISFunctions::construct<PDFInterface>({}, &params, delta_integrand, false, sign);

		const double nlo = nlo_coefficient * (nlo_value + nlo_delta);

		return PerturbativeQuantity {lo, lo + nlo};
	}

	PerturbativeQuantity F2(const TRFKinematics &kinematics) const {
		return structure_function(
			kinematics, DISFunctions::F2::LO::integrand, DISFunctions::F2::NLO::integrand, DISFunctions::F2::NLO::delta, 1
		);
	}

	PerturbativeQuantity FL(const TRFKinematics &kinematics) const {
		return structure_function(
			kinematics, DISFunctions::FL::LO::integrand, DISFunctions::FL::NLO::integrand, DISFunctions::FL::NLO::delta, 1
		);
	}
	PerturbativeQuantity F3(const TRFKinematics &kinematics) const {
		return structure_function(
			kinematics, DISFunctions::F3::LO::integrand, DISFunctions::F3::NLO::integrand, DISFunctions::F3::NLO::delta, -1
		);
	}
	constexpr PerturbativeQuantity compute_structure_function(const StructureFunction F, const TRFKinematics &kinematics) const {
		switch (F) {
		case StructureFunction::F2: return F2(kinematics);
		case StructureFunction::FL: return FL(kinematics);
		case StructureFunction::F3: return F3(kinematics);
		}
	}

	PerturbativeQuantity differential_cross_section_xQ2(const TRFKinematics &kinematics) const {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;

		const double alpha_s = compute_alpha_s(kinematics);
		const double nlo_coefficient = alpha_s / (2 * std::numbers::pi);

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

		const double lo = DISFunctions::cross_section<PDFInterface>({}, &params, 
			DISFunctions::F2::LO::integrand, DISFunctions::FL::LO::integrand, DISFunctions::F3::LO::integrand,
			false
		);

		Integrator nlo_integrator([](double input[], [[maybe_unused]] std::size_t dim, void *params_in) {
			return DISFunctions::cross_section<PDFInterface>(input, params_in, 
				DISFunctions::F2::NLO::integrand, DISFunctions::FL::NLO::integrand, DISFunctions::F3::NLO::integrand,
				true
			);
		}, {x}, {1.0}, integration_parameters, &params);
		const auto nlo_result = nlo_integrator.integrate();
		const double nlo_integral = nlo_result.value;

		const double delta_contribution = DISFunctions::cross_section<PDFInterface>({}, &params,
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

	PerturbativeQuantity integrated_cross_section(const TRFKinematics &placeholder_kinematics, const double Q2_min) const {
		const auto lo_integrand = [&](double input[]) {
			const double x = input[0];
			const double Q2 = input[1];

			const double target_mass = placeholder_kinematics.target_mass;
			const double E_beam = placeholder_kinematics.E_beam;
			const double y = Q2 / (2.0 * x * target_mass * E_beam);
			const TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, target_mass, placeholder_kinematics.projectile_mass);

			const double factorization_scale = factorization_scale_function(kinematics);
			const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);

			pdf1.evaluate(x, factorization_scale);

			DISFunctions::Parameters<PDFInterface> params {
				pdf1, pdf2,
				flavors,
				0.0,
				process, kinematics,
				factorization_scale,
				factorization_scale_log
			};

			const double prefactor = use_modified_cross_section_prefactor 
										? CommonFunctions::cross_section_modified_prefactor(kinematics) 
										: CommonFunctions::cross_section_prefactor(kinematics);

			const double differential_cross_section = DISFunctions::cross_section<PDFInterface>({}, &params, 
				DISFunctions::F2::LO::integrand, DISFunctions::FL::LO::integrand, DISFunctions::F3::LO::integrand,
				false
			);
			return prefactor * differential_cross_section;
		};

		const auto nlo_integrand = [&](double input[]) {
			const double x = input[1];
			const double Q2 = input[2];

			const double target_mass = placeholder_kinematics.target_mass;
			const double E_beam = placeholder_kinematics.E_beam;
			const double y = Q2 / (2.0 * x * target_mass * E_beam);
			const TRFKinematics kinematics = TRFKinematics::y_E_beam(x, y, E_beam, target_mass, placeholder_kinematics.projectile_mass);

			const double alpha_s = compute_alpha_s(kinematics);
			const double nlo_coefficient = alpha_s / (2 * std::numbers::pi);

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

			const double differential_cross_section = DISFunctions::cross_section<PDFInterface>(input, &params, 
				DISFunctions::F2::NLO::integrand, DISFunctions::FL::NLO::integrand, DISFunctions::F3::NLO::integrand,
				true
			) + DISFunctions::cross_section<PDFInterface>({}, &params,
				DISFunctions::F2::NLO::delta, DISFunctions::FL::NLO::delta, DISFunctions::F3::NLO::delta, 
				false
			) / (1.0 - x);

			const double prefactor = use_modified_cross_section_prefactor 
										? CommonFunctions::cross_section_modified_prefactor(kinematics) 
										: CommonFunctions::cross_section_prefactor(kinematics);

			return prefactor * differential_cross_section;
		};

		const double target_mass = placeholder_kinematics.target_mass;
		const double E_beam = placeholder_kinematics.E_beam;

		const double x_min = Q2_min / (2.0 * target_mass * E_beam);

		Integrator lo_integrator([&](double input[], std::size_t dim, [[maybe_unused]] void *params_in) {
			const double x = input[0];
			const double Q2_max = 2.0 * x * target_mass * E_beam;

			double *scaled_input = new double[dim];
			scaled_input[0] = x; // x
			scaled_input[1] = Q2_min + (Q2_max - Q2_min) * input[1]; // Q^2

			const double result = (Q2_max - Q2_min) * lo_integrand(scaled_input);

			delete[] scaled_input;
			return result;
		}, {x_min, 0.0}, {1.0 - 1e-9, 1.0 - 1e-9}, integration_parameters, nullptr);


		Integrator nlo_integrator([&](double input[], std::size_t dim, [[maybe_unused]] void *params_in) {
			const double x = input[0];
			const double Q2_max = 2.0 * x * target_mass * E_beam;

			double *scaled_input = new double[dim];
			scaled_input[0] = x + (1.0 - x) * input[2]; // xi
			scaled_input[1] = x; // x
			scaled_input[2] = Q2_min + (Q2_max - Q2_min) * input[1]; // Q^2

			const double result = (Q2_max - Q2_min) * (1.0 - x) * nlo_integrand(scaled_input);

			delete[] scaled_input;
			return result;
		}, {x_min, 0.0, 0.0}, {1.0 - 1e-9, 1.0 - 1e-9, 1.0 - 1e-9}, integration_parameters, nullptr);

		const auto nlo_result = nlo_integrator.integrate();
		const double nlo = nlo_result.value;

		const auto lo_result = lo_integrator.integrate();
		const double lo = lo_result.value;

		const PerturbativeQuantity result = PerturbativeQuantity {lo, lo + nlo};
		return result;
	}
};

#endif