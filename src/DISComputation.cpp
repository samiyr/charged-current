#ifndef DIS_COMPUTATION_H
#define DIS_COMPUTATION_H

#include "Integrator.cpp"
#include "Constants.cpp"
#include "Utility.cpp"
#include "PerturbativeQuantity.cpp"
#include "Process.cpp"
#include "StructureFunction.cpp"
#include "DISFunctions.cpp"
#include "CommonFunctions.cpp"
#include <optional>

template <typename PDFInterface, typename FactorizationScaleFunction>
class DISComputation {
	public:
	double sqrt_s;
	double s;

	FlavorInfo flavors;

	const PDFInterface pdf1;
	const PDFInterface pdf2;

	const size_t points;
	const double max_chi_squared_deviation;
	const double max_relative_error;
	const unsigned int iter_max;

	const Process process;
	const bool momentum_fraction_mass_corrections;

	const std::optional<FactorizationScaleFunction> factorization_scale_function;

	DISComputation (
		const double _sqrt_s, 
		const FlavorVector _active_flavors,
		const std::array<double, TOTAL_FLAVORS> _flavor_masses, 
		const PDFInterface _pdf,
		const size_t _points,
		const double _max_chi_squared_deviation, 
		const double _max_relative_error,
		const unsigned int _iter_max,
		const Process _process,
		const bool _momentum_fraction_mass_corrections,
		const std::optional<FactorizationScaleFunction> _factorization_scale_function
	) : sqrt_s(_sqrt_s), 
	s(_sqrt_s * _sqrt_s), 
	flavors(_active_flavors, _flavor_masses), 
	pdf1(_pdf), 
	pdf2(_pdf), 
	points(_points), 
	max_chi_squared_deviation(_max_chi_squared_deviation), 
	max_relative_error(_max_relative_error),
	iter_max(_iter_max),
	process(_process),
	momentum_fraction_mass_corrections(_momentum_fraction_mass_corrections),
	factorization_scale_function(_factorization_scale_function) { }

	constexpr bool nontrivial_factorization_scale() const { return factorization_scale_function.has_value(); }
	constexpr double compute_factorization_scale(const TRFKinematics &kinematics) const {
		if (nontrivial_factorization_scale()) {
			return (*factorization_scale_function)(kinematics);
		}
		return kinematics.Q2;
	}

	PerturbativeQuantity F2(const TRFKinematics &kinematics) const {	
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		const double factorization_scale = compute_factorization_scale(kinematics);
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
		
		const double lo = x * DISFunctions::Evaluation::construct<PDFInterface>({}, &params, DISFunctions::Integrands::F2x_lo_integrand, false, 1);

		Integrator nlo_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return DISFunctions::Evaluation::construct<PDFInterface>(input, params_in, DISFunctions::Integrands::F2x_nlo_integrand, true, 1);
		}, {x}, {1.0}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		auto nlo_result = nlo_integrator.integrate();
		const double nlo_value = nlo_result.value;

		const double nlo_delta = DISFunctions::Evaluation::construct<PDFInterface>({}, &params, DISFunctions::Integrands::F2x_delta_integrand, false, 1);

		const double nlo = nlo_coefficient * x * (nlo_value + nlo_delta);

		return PerturbativeQuantity {lo, lo + nlo};
	}

	PerturbativeQuantity FL(const TRFKinematics &kinematics) const {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		const double factorization_scale = compute_factorization_scale(kinematics);
		const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);

		DISFunctions::Parameters<PDFInterface> params {
			pdf1, pdf2,
			flavors,
			nlo_coefficient,
			process, kinematics,
			factorization_scale,
			factorization_scale_log
		};

		Integrator nlo_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return DISFunctions::Evaluation::construct<PDFInterface>(input, params_in, DISFunctions::Integrands::FLx_nlo_integrand, true, 1);
		}, {x}, {1.0}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		auto nlo_result = nlo_integrator.integrate();
		const double nlo_value = nlo_result.value;

		const double nlo = nlo_coefficient * x * nlo_value;

		return PerturbativeQuantity {0.0, nlo};
	}
	PerturbativeQuantity xF3(const TRFKinematics &kinematics) const {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);
		
		const double factorization_scale = compute_factorization_scale(kinematics);
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

		const double lo = x * DISFunctions::Evaluation::construct<PDFInterface>({}, &params, DISFunctions::Integrands::F3_lo_integrand, false, -1);

		Integrator nlo_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return DISFunctions::Evaluation::construct<PDFInterface>(input, params_in, DISFunctions::Integrands::F3_nlo_integrand, true, -1);
		}, {x}, {1.0}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		auto nlo_result = nlo_integrator.integrate();
		const double nlo_value = nlo_result.value;

		const double nlo_delta = DISFunctions::Evaluation::construct<PDFInterface>({}, &params, DISFunctions::Integrands::F3x_delta_integrand, false, -1);

		const double nlo = nlo_coefficient * x * (nlo_value + nlo_delta);

		return PerturbativeQuantity {lo, lo + nlo};
	}
	constexpr PerturbativeQuantity compute_structure_function(const StructureFunction F, const TRFKinematics &kinematics) const {
		switch (F) {
		case StructureFunction::F2: return F2(kinematics);
		case StructureFunction::FL: return FL(kinematics);
		case StructureFunction::xF3: return xF3(kinematics);
		}
	}

	PerturbativeQuantity differential_cross_section_xQ2_direct(const TRFKinematics &kinematics) const {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;
		const double y = kinematics.y;

		if (!kinematics.is_valid()) { return PerturbativeQuantity {0.0, 0.0}; }

		const double prefactor = CommonFunctions::cross_section_prefactor(Q2);

		const PerturbativeQuantity f2 = F2(kinematics);
		const PerturbativeQuantity fL = FL(kinematics);
		const PerturbativeQuantity xf3 = xF3(kinematics);

		const double term1 = 1 - y + 0.5 * y * y;
		const double term2 = - 0.5 * y * y;
		const double term3 = y * (1 - 0.5 * y);

		const PerturbativeQuantity result = prefactor * (term1 * f2 + term2 * fL + double(process.W_sign()) * term3 * xf3) / x;
		return result;
	}

	PerturbativeQuantity differential_cross_section_xQ2_indirect(const TRFKinematics &kinematics) const {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		const double factorization_scale = compute_factorization_scale(kinematics);
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
			DISFunctions::Integrands::F2x_lo_integrand, DISFunctions::Integrands::FLx_lo_integrand, DISFunctions::Integrands::F3_lo_integrand,
			false
		);

		Integrator nlo_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return DISFunctions::Evaluation::cross_section<PDFInterface>(input, params_in, 
				DISFunctions::Integrands::F2x_nlo_integrand, DISFunctions::Integrands::FLx_nlo_integrand, DISFunctions::Integrands::F3_nlo_integrand,
				true
			);
		}, {x}, {1.0}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto nlo_result = nlo_integrator.integrate();
		const double nlo_integral = nlo_result.value;

		const double delta_contribution = DISFunctions::Evaluation::cross_section<PDFInterface>({}, &params,
			DISFunctions::Integrands::F2x_delta_integrand, DISFunctions::Integrands::FLx_delta_integrand, DISFunctions::Integrands::F3x_delta_integrand, 
			false
		);

		const double nlo = nlo_coefficient * (nlo_integral + delta_contribution);

		const PerturbativeQuantity result = PerturbativeQuantity {lo, lo + nlo};
		const double prefactor = CommonFunctions::cross_section_prefactor(Q2);
		return prefactor * result;
	}
};

#endif