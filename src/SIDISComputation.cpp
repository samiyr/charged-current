#ifndef SIDIS_COMPUTATION_H
#define SIDIS_COMPUTATION_H

#include <vector>
#include "Flavor.cpp"
#include "PerturbativeQuantity.cpp"
#include "PDFCommon.cpp"
#include "SIDISFunctions.cpp"
#include "Integrator.cpp"
#include "Decay.cpp"
#include "TRFKinematics.cpp"
#include "FragmentationConfiguration.cpp"
#include <optional>

template <typename PDFInterface, typename FFInterface, typename DecayFunction, typename FactorizationScaleFunction, typename FragmentationScaleFunction>
class SIDISComputation/* : Utility::Traced<SIDISComputation<PDFInterface, FFInterface, DecayFunction, FactorizationScaleFunction, FragmentationScaleFunction>> */ {
	public:
	double sqrt_s;
	double s;
	
	FlavorInfo flavors;

	PDFInterface pdf1;
	FragmentationConfiguration<FFInterface, DecayFunction> ff1;

	PDFInterface pdf2;
	FragmentationConfiguration<FFInterface, DecayFunction> ff2;

	const size_t points;
	const double max_chi_squared_deviation;
	const double max_relative_error;
	const unsigned int iter_max;

	const Process process;
	const bool momentum_fraction_mass_corrections;

	const std::optional<FactorizationScaleFunction> factorization_scale_function;
	const std::optional<FragmentationScaleFunction> fragmentation_scale_function;

	SIDISComputation (
		const double _sqrt_s, 
		const FlavorVector _active_flavors, 
		const std::array<double, TOTAL_FLAVORS> _flavor_masses,
		const PDFInterface _pdf,
		const FragmentationConfiguration<FFInterface, DecayFunction> _ff,
		const size_t _points,
		const double _max_chi_squared_deviation, 
		const double _max_relative_error,
		const unsigned int _iter_max,
		const Process _process,
		const bool _momentum_fraction_mass_corrections,
		const std::optional<FactorizationScaleFunction> _factorization_scale_function,
		const std::optional<FragmentationScaleFunction> _fragmentation_scale_function
	) : sqrt_s(_sqrt_s),
	s(_sqrt_s * _sqrt_s), 
	flavors(_active_flavors, _flavor_masses),
	pdf1(_pdf), 
	ff1(_ff),
	pdf2(_pdf), 
	ff2(_ff),
	points(_points), 
	max_chi_squared_deviation(_max_chi_squared_deviation), 
	max_relative_error(_max_relative_error),
	iter_max(_iter_max),
	process(_process),
	momentum_fraction_mass_corrections(_momentum_fraction_mass_corrections),
	factorization_scale_function(_factorization_scale_function),
	fragmentation_scale_function(_fragmentation_scale_function) { }

	constexpr bool nontrivial_factorization_scale() { return factorization_scale_function.has_value(); }
	constexpr bool nontrivial_fragmentation_scale() { return fragmentation_scale_function.has_value(); }
	constexpr double compute_factorization_scale(const TRFKinematics &kinematics) {
		if (nontrivial_factorization_scale()) {
			return (*factorization_scale_function)(kinematics);
		}
		return kinematics.Q2;
	}
	constexpr double compute_fragmentation_scale(const TRFKinematics &kinematics) {
		if (fragmentation_scale_function) {
			return (*fragmentation_scale_function)(kinematics);
		}
		return kinematics.Q2;
	}

	PerturbativeQuantity F2_combined(const double z, const TRFKinematics &kinematics) {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		const double factorization_scale = compute_factorization_scale(kinematics);
		const double fragmentation_scale = compute_fragmentation_scale(kinematics);

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
			factorization_scale, fragmentation_scale,
			factorization_scale_log, fragmentation_scale_log
		};

		const double lo = x * SIDISFunctions::Evaluation::construct<PDFInterface, FFInterface>({}, &params, SIDISFunctions::Integrands::F2x_lo_integrand, false, false, false, 1);

		Integrator nlo_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::construct<PDFInterface, FFInterface>(input, params_in, SIDISFunctions::Integrands::F2x_nlo_integrand, true, true, false, 1);
		}, {x, z}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		auto nlo_result = nlo_integrator.integrate();
		const double nlo_value = nlo_result.value;

		const double nlo = nlo_coefficient * x * nlo_value;

		return PerturbativeQuantity {lo, lo + nlo};
	}
	PerturbativeQuantity FL_combined(const double z, const TRFKinematics &kinematics) {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		const double factorization_scale = compute_factorization_scale(kinematics);
		const double fragmentation_scale = compute_fragmentation_scale(kinematics);

		const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);
		const double fragmentation_scale_log = fragmentation_scale == Q2 ? 0 : std::log(Q2 / fragmentation_scale);

		SIDISFunctions::Parameters<PDFInterface, FFInterface, DecayFunction> params {
			pdf1, ff1, pdf2, ff2,
			flavors,
			nlo_coefficient,
			process, kinematics,
			z,
			factorization_scale, fragmentation_scale,
			factorization_scale_log, fragmentation_scale_log
		};

		Integrator xi_xip_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::construct<PDFInterface, FFInterface>(input, params_in, SIDISFunctions::Integrands::FLx_xi_xip_integrand, true, true, false, 1);
		}, {x, z}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		auto xi_xip_result = xi_xip_integrator.integrate();
		const double xi_xip_integral = xi_xip_result.value;

		const double nlo = nlo_coefficient * x * xi_xip_integral;

		return PerturbativeQuantity {0.0, nlo};
	}
	PerturbativeQuantity xF3_combined(const double z, const TRFKinematics &kinematics) {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);
		
		const double factorization_scale = compute_factorization_scale(kinematics);
		const double fragmentation_scale = compute_fragmentation_scale(kinematics);

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
			factorization_scale, fragmentation_scale,
			factorization_scale_log, fragmentation_scale_log
		};

		const double lo = x * SIDISFunctions::Evaluation::construct<PDFInterface, FFInterface>({}, &params, SIDISFunctions::Integrands::F3_lo_integrand, false, false, false, -1);

		Integrator nlo_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::construct<PDFInterface, FFInterface>(input, params_in, SIDISFunctions::Integrands::F3_nlo_integrand, true, true, false, -1);
		}, {x, z}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		auto nlo_result = nlo_integrator.integrate();
		const double nlo_value = nlo_result.value;

		const double nlo = nlo_coefficient * x * nlo_value;

		return PerturbativeQuantity {lo, lo + nlo};
	}

	PerturbativeQuantity F2_separated(const double z, const TRFKinematics &kinematics) {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		const double factorization_scale = compute_factorization_scale(kinematics);
		const double fragmentation_scale = compute_fragmentation_scale(kinematics);

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
			factorization_scale, fragmentation_scale,
			factorization_scale_log, fragmentation_scale_log
		};

		const double lo = x * SIDISFunctions::Evaluation::construct<PDFInterface, FFInterface>({}, &params, SIDISFunctions::Integrands::F2x_lo_integrand, false, false, false, 1);

		Integrator xi_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::construct<PDFInterface, FFInterface>(input, params_in, SIDISFunctions::Integrands::F2x_xi_integrand, true, false, false, 1);
		}, {x}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		auto xi_result = xi_integrator.integrate();
		const double xi_integral = xi_result.value;
		
		Integrator xip_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::construct<PDFInterface, FFInterface>(input, params_in, SIDISFunctions::Integrands::F2x_xip_integrand, false, true, false, 1);
		}, {z}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		auto xip_result = xip_integrator.integrate();
		const double xip_integral = xip_result.value;

		Integrator xi_xip_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::construct<PDFInterface, FFInterface>(input, params_in, SIDISFunctions::Integrands::F2x_xi_xip_integrand, true, true, false, 1);
		}, {x, z}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		auto xi_xip_result = xi_xip_integrator.integrate();
		const double xi_xip_integral = xi_xip_result.value;

		const double nlo1 = SIDISFunctions::Evaluation::construct<PDFInterface, FFInterface>({}, &params, SIDISFunctions::Integrands::F2x_delta_integrand, false, false, false, 1);
		const double nlo2 = xi_integral + xip_integral + xi_xip_integral;

		const double nlo = nlo_coefficient * x * (nlo1 + nlo2);

		return PerturbativeQuantity {lo, lo + nlo};
	}
	PerturbativeQuantity FL_separated(const double z, const TRFKinematics &kinematics) {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		const double factorization_scale = compute_factorization_scale(kinematics);
		const double fragmentation_scale = compute_fragmentation_scale(kinematics);

		const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);
		const double fragmentation_scale_log = fragmentation_scale == Q2 ? 0 : std::log(Q2 / fragmentation_scale);

		SIDISFunctions::Parameters<PDFInterface, FFInterface, DecayFunction> params {
			pdf1, ff1, pdf2, ff2,
			flavors,
			nlo_coefficient,
			process, kinematics,
			z,
			factorization_scale, fragmentation_scale,
			factorization_scale_log, fragmentation_scale_log
		};

		Integrator xi_xip_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::construct<PDFInterface, FFInterface>(input, params_in, SIDISFunctions::Integrands::FLx_xi_xip_integrand, true, true, false, 1);
		}, {x, z}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		auto xi_xip_result = xi_xip_integrator.integrate();
		const double xi_xip_integral = xi_xip_result.value;

		const double nlo = nlo_coefficient * x * xi_xip_integral;

		return PerturbativeQuantity {0, nlo};
	}
	PerturbativeQuantity xF3_separated(const double z, const TRFKinematics &kinematics) {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);
		
		const double factorization_scale = compute_factorization_scale(kinematics);
		const double fragmentation_scale = compute_fragmentation_scale(kinematics);

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
			factorization_scale, fragmentation_scale,
			factorization_scale_log, fragmentation_scale_log
		};

		const double lo = x * SIDISFunctions::Evaluation::construct<PDFInterface, FFInterface>({}, &params, SIDISFunctions::Integrands::F3_lo_integrand, false, false, false, -1);

		Integrator xi_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::construct<PDFInterface, FFInterface>(input, params_in, SIDISFunctions::Integrands::F3_xi_integrand, true, false, false, -1);
		}, {x}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		auto xi_result = xi_integrator.integrate();
		const double xi_integral = xi_result.value;
		
		Integrator xip_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::construct<PDFInterface, FFInterface>(input, params_in, SIDISFunctions::Integrands::F3_xip_integrand, false, true, false, -1);
		}, {z}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		auto xip_result = xip_integrator.integrate();
		const double xip_integral = xip_result.value;

		Integrator xi_xip_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::construct<PDFInterface, FFInterface>(input, params_in, SIDISFunctions::Integrands::F3_xi_xip_integrand, true, true, false, -1);
		}, {x, z}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		auto xi_xip_result = xi_xip_integrator.integrate();
		const double xi_xip_integral = xi_xip_result.value;

		const double nlo1 = SIDISFunctions::Evaluation::construct<PDFInterface, FFInterface>({}, &params, SIDISFunctions::Integrands::F3_delta_integrand, false, false, false, -1);
		const double nlo2 = xi_integral + xip_integral + xi_xip_integral;

		const double nlo = nlo_coefficient * x * (nlo1 + nlo2);

		return PerturbativeQuantity {lo, lo + nlo};
	}

	PerturbativeQuantity compute_structure_function(const StructureFunction structure_function, const double z, const TRFKinematics &kinematics, const bool combine_integrals) {
		switch (structure_function) {
		case StructureFunction::F2: return combine_integrals ? F2_combined(z, kinematics) : F2_separated(z, kinematics);
		case StructureFunction::FL: return combine_integrals ? FL_combined(z, kinematics) : FL_separated(z, kinematics);
		case StructureFunction::xF3: return combine_integrals ? xF3_combined(z, kinematics) : xF3_separated(z, kinematics);
		}
	}

	PerturbativeQuantity differential_cross_section_xQ2_direct(const TRFKinematics &kinematics, const PerturbativeQuantity f2, const PerturbativeQuantity fL, const PerturbativeQuantity xf3) {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;
		const double y = kinematics.y;
		const double prefactor = CommonFunctions::cross_section_prefactor(Q2);

		if (!kinematics.is_valid()) { return PerturbativeQuantity {0.0, 0.0}; }

		// const std::optional<double> y_opt = CommonFunctions::compute_y(x, Q2, s, process.target.mass, process.projectile.mass);
		// if (!y_opt.has_value()) {
		// 	return {0, 0};
		// }
		// const double y = *y_opt;

		const double term1 = 1 - y + 0.5 * y * y;
		const double term2 = - 0.5 * y * y;
		const double term3 = y * (1 - 0.5 * y);

		const PerturbativeQuantity result = (term1 * f2 + term2 * fL + double(process.W_sign()) * term3 * xf3) / x;
		return prefactor * result;
	}

	PerturbativeQuantity differential_cross_section_xQ2_direct_combined(const double z, const TRFKinematics &kinematics) {
		return differential_cross_section_xQ2_direct(kinematics, F2_combined(z, kinematics), FL_combined(z, kinematics), xF3_combined(z, kinematics));
	}
	PerturbativeQuantity differential_cross_section_xQ2_direct_separated(const double z, const TRFKinematics &kinematics) {
		return differential_cross_section_xQ2_direct(kinematics, F2_separated(z, kinematics), FL_separated(z, kinematics), xF3_separated(z, kinematics));
	}

	PerturbativeQuantity differential_cross_section_xQ2_indirect_separated(const double z, const TRFKinematics &kinematics) {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		const double factorization_scale = compute_factorization_scale(kinematics);
		const double fragmentation_scale = compute_fragmentation_scale(kinematics);

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
			factorization_scale, fragmentation_scale,
			factorization_scale_log, fragmentation_scale_log
		};

		const double lo = SIDISFunctions::Evaluation::cross_section<PDFInterface, FFInterface, DecayFunction>({}, &params, 
			SIDISFunctions::Integrands::F2x_lo_integrand, SIDISFunctions::Integrands::FLx_lo_integrand, SIDISFunctions::Integrands::F3_lo_integrand,
			false, false, false
		);

		Integrator xi_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::cross_section<PDFInterface, FFInterface, DecayFunction>(input, params_in, 
				SIDISFunctions::Integrands::F2x_xi_integrand, SIDISFunctions::Integrands::FLx_xi_integrand, SIDISFunctions::Integrands::F3_xi_integrand,
				true, false, false
			);
		}, {x}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto xi_result = xi_integrator.integrate();
		const double xi_integral = xi_result.value;
		
		Integrator xip_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::cross_section<PDFInterface, FFInterface, DecayFunction>(input, params_in, 
				SIDISFunctions::Integrands::F2x_xip_integrand, SIDISFunctions::Integrands::FLx_xip_integrand, SIDISFunctions::Integrands::F3_xip_integrand,
				false, true, false
			);
		}, {z}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto xip_result = xip_integrator.integrate();
		const double xip_integral = xip_result.value;

		Integrator xi_xip_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::cross_section<PDFInterface, FFInterface, DecayFunction>(input, params_in, 
				SIDISFunctions::Integrands::F2x_xi_xip_integrand, SIDISFunctions::Integrands::FLx_xi_xip_integrand, SIDISFunctions::Integrands::F3_xi_xip_integrand,
				true, true, false
			);
		}, {x, z}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto xi_xip_result = xi_xip_integrator.integrate();
		const double xi_xip_integral = xi_xip_result.value;

		const double nlo1 = SIDISFunctions::Evaluation::cross_section<PDFInterface, FFInterface, DecayFunction>({}, &params,
			SIDISFunctions::Integrands::F2x_delta_integrand, SIDISFunctions::Integrands::FLx_delta_integrand, SIDISFunctions::Integrands::F3_delta_integrand,
			false, false, false
		);
		const double nlo2 = xi_integral + xip_integral + xi_xip_integral;

		const double nlo = nlo_coefficient * (nlo1 + nlo2);

		const PerturbativeQuantity result = PerturbativeQuantity {lo, lo + nlo};
		const double prefactor = CommonFunctions::cross_section_prefactor(Q2);
		return prefactor * result;
	}

	PerturbativeQuantity differential_cross_section_xQ2_indirect_combined(const double z, const TRFKinematics &kinematics) {
		const double Q2 = kinematics.Q2;
		const double x = kinematics.x;
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		const double factorization_scale = compute_factorization_scale(kinematics);
		const double fragmentation_scale = compute_fragmentation_scale(kinematics);

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
			factorization_scale, fragmentation_scale,
			factorization_scale_log, fragmentation_scale_log
		};

		const double lo = SIDISFunctions::Evaluation::cross_section<PDFInterface, FFInterface, DecayFunction>({}, &params, 
			SIDISFunctions::Integrands::F2x_lo_integrand, SIDISFunctions::Integrands::FLx_lo_integrand, SIDISFunctions::Integrands::F3_lo_integrand,
			false, false, false
		);

		Integrator nlo_integrator([](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::cross_section<PDFInterface, FFInterface, DecayFunction>(input, params_in, 
				SIDISFunctions::Integrands::F2x_nlo_integrand, SIDISFunctions::Integrands::FLx_nlo_integrand, SIDISFunctions::Integrands::F3_nlo_integrand,
				true, true, false
			);
		}, {x, z}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto nlo_result = nlo_integrator.integrate();
		const double nlo_integral = nlo_result.value;

		const double nlo = nlo_coefficient * nlo_integral;

		const PerturbativeQuantity result = PerturbativeQuantity {lo, lo + nlo};
		const double prefactor = CommonFunctions::cross_section_prefactor(Q2);
		return prefactor * result;
	}

	PerturbativeQuantity lepton_pair_cross_section_xQ2_separated(const TRFKinematics &kinematics) {
		const double x = kinematics.x;
		const double Q2 = kinematics.Q2;

		std::vector<double> z_mins;
		for (const auto &decay : ff1.decays) {
			z_mins.push_back(SIDISFunctions::compute_z_min(kinematics, decay));
		}
		const double z_min = *std::min_element(z_mins.begin(), z_mins.end());

		const double alpha_s = pdf1.alpha_s(Q2);
		const double nlo_coefficient = alpha_s / (2 * M_PI);

		const double factorization_scale = compute_factorization_scale(kinematics);
		const double fragmentation_scale = compute_fragmentation_scale(kinematics);

		const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);
		const double fragmentation_scale_log = fragmentation_scale == Q2 ? 0 : std::log(Q2 / fragmentation_scale);

		pdf1.evaluate(x, factorization_scale);

		SIDISFunctions::Parameters<PDFInterface, FFInterface, DecayFunction> params {
			pdf1, ff1, pdf2, ff2,
			flavors,
			nlo_coefficient,
			process, kinematics,
			0.0,
			factorization_scale, fragmentation_scale,
			factorization_scale_log, fragmentation_scale_log
		};

		Integrator lo_integrator([&](double input[], [[maybe_unused]] size_t dim, void* params_in) {
			return SIDISFunctions::Evaluation::cross_section<PDFInterface, FFInterface, DecayFunction>(input, params_in, 
				SIDISFunctions::Integrands::F2x_lo_integrand, SIDISFunctions::Integrands::FLx_lo_integrand, SIDISFunctions::Integrands::F3_lo_integrand,
				false, false, true
			);
		}, {z_min}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto lo_result = lo_integrator.integrate();
		const double lo = lo_result.value;

		Integrator xi_integrator([&](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::cross_section<PDFInterface, FFInterface, DecayFunction>(input, params_in, 
				SIDISFunctions::Integrands::F2x_xi_integrand, SIDISFunctions::Integrands::FLx_xi_integrand, SIDISFunctions::Integrands::F3_xi_integrand,
				true, false, true
			);
		}, {x, z_min}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto xi_result = xi_integrator.integrate();
		const double xi_integral = xi_result.value;
		
		Integrator xip_integrator([&](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::cross_section<PDFInterface, FFInterface, DecayFunction>(input, params_in, 
				SIDISFunctions::Integrands::F2x_xip_integrand, SIDISFunctions::Integrands::FLx_xip_integrand, SIDISFunctions::Integrands::F3_xip_integrand,
				false, true, true
			);
		}, {z_min, z_min}, {1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto xip_result = xip_integrator.integrate();
		const double xip_integral = xip_result.value;

		Integrator xi_xip_integrator([&](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::cross_section<PDFInterface, FFInterface, DecayFunction>(input, params_in, 
				SIDISFunctions::Integrands::F2x_xi_xip_integrand, SIDISFunctions::Integrands::FLx_xi_xip_integrand, SIDISFunctions::Integrands::F3_xi_xip_integrand,
				true, true, true
			);
		}, {x, z_min, z_min}, {1, 1, 1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto xi_xip_result = xi_xip_integrator.integrate();
		const double xi_xip_integral = xi_xip_result.value;

		// Could combine with z_integrator, but then lose separation between LO and NLO
		Integrator delta_integrator([&](double input[], [[maybe_unused]] size_t dim, void* params_in) {
			return SIDISFunctions::Evaluation::cross_section<PDFInterface, FFInterface, DecayFunction>(input, params_in, 
				SIDISFunctions::Integrands::F2x_delta_integrand, SIDISFunctions::Integrands::FLx_delta_integrand, SIDISFunctions::Integrands::F3_delta_integrand,
				false, false, true
			);
		}, {z_min}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto delta_result = delta_integrator.integrate();
		const double delta_integral = delta_result.value;

		const double nlo = nlo_coefficient * (xi_integral + xip_integral + xi_xip_integral + delta_integral);

		const PerturbativeQuantity result = PerturbativeQuantity {lo, lo + nlo};
		const double prefactor = CommonFunctions::cross_section_prefactor(Q2);
		return prefactor * result;
	}
	PerturbativeQuantity lepton_pair_cross_section_xQ2_combined(const TRFKinematics &kinematics) {
		const double x = kinematics.x;
		const double Q2 = kinematics.Q2;

		std::vector<double> z_mins;
		for (const auto &decay : ff1.decays) {
			z_mins.push_back(SIDISFunctions::compute_z_min(kinematics, decay));
		}
		const double z_min = *std::min_element(z_mins.begin(), z_mins.end());

		const double alpha_s = pdf1.alpha_s(Q2);
		const double nlo_coefficient = alpha_s / (2 * M_PI);

		const double factorization_scale = compute_factorization_scale(kinematics);
		const double fragmentation_scale = compute_fragmentation_scale(kinematics);

		const double factorization_scale_log = factorization_scale == Q2 ? 0 : std::log(Q2 / factorization_scale);
		const double fragmentation_scale_log = fragmentation_scale == Q2 ? 0 : std::log(Q2 / fragmentation_scale);

		pdf1.evaluate(x, factorization_scale);

		SIDISFunctions::Parameters<PDFInterface, FFInterface, DecayFunction> params {
			pdf1, ff1, pdf2, ff2,
			flavors,
			nlo_coefficient,
			process, kinematics,
			0.0,
			factorization_scale, fragmentation_scale,
			factorization_scale_log, fragmentation_scale_log
		};

		Integrator lo_integrator([&](double input[], [[maybe_unused]] size_t dim, void* params_in) {
			return SIDISFunctions::Evaluation::cross_section<PDFInterface, FFInterface, DecayFunction>(input, params_in, 
				SIDISFunctions::Integrands::F2x_lo_integrand, SIDISFunctions::Integrands::FLx_lo_integrand, SIDISFunctions::Integrands::F3_lo_integrand,
				false, false, true
			);
		}, {z_min}, {1}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto lo_result = lo_integrator.integrate();
		const double lo = lo_result.value;

		Integrator nlo_integrator([&](double input[], [[maybe_unused]] size_t dim, void *params_in) {
			return SIDISFunctions::Evaluation::cross_section<PDFInterface, FFInterface, DecayFunction>(input, params_in, 
				SIDISFunctions::Integrands::F2x_nlo_integrand, SIDISFunctions::Integrands::FLx_nlo_integrand, SIDISFunctions::Integrands::F3_nlo_integrand,
				true, true, true
			) / (1.0 - input[2]);
		}, {x, z_min, z_min}, {1, 1, 1.0 - 1e-3}, points, &params, max_chi_squared_deviation, max_relative_error, iter_max);
		const auto nlo_result = nlo_integrator.integrate();
		const double nlo_integral = nlo_result.value / (1.0 - x);

		const double nlo = nlo_coefficient * nlo_integral;

		const PerturbativeQuantity result = PerturbativeQuantity {lo, lo + nlo};
		const double prefactor = CommonFunctions::cross_section_prefactor(Q2);
		return prefactor * result;
	}
};

#endif