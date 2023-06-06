#ifndef SIDIS_FUNCTIONS_H
#define SIDIS_FUNCTIONS_H

#include "Constants.cpp"
#include <cmath>
#include "Flavor.cpp"
#include "PDFCommon.cpp"
#include "Decay.cpp"
#include "FragmentationConfiguration.cpp"
#include "TRFKinematics.cpp"
#include "DecayFunctions.cpp"

namespace SIDISFunctions {
	template <typename PDFInterface, typename FFInterface, typename DecayFunction>
	struct Parameters {
		PDFInterface &pdf1;
		FragmentationConfiguration<FFInterface, DecayFunction> &ff1;
		PDFInterface &pdf2;
		FragmentationConfiguration<FFInterface, DecayFunction> &ff2;

		const FlavorInfo &flavors;

		const double nlo_coefficient;

		const Process &process;
		const TRFKinematics &kinematics;

		const double z;

		// Could some caching here
		// const double xq;
		// const double zq;
	};

	constexpr double delta_contribution(const double x, const double z) {
		return 2 * Constants::C_F * (std::pow(std::log(1 - x) + std::log(1 - z), 2) - 8) / z;
	}

	template <typename DecayFunction>
	constexpr double compute_z_min(const TRFKinematics kinematics, const Decay<DecayFunction> &decay) {
		return std::max({
			decay.lepton_momentum_min / (kinematics.y * kinematics.E_beam), 
			2 * kinematics.x * decay.resonance.mass * decay.hadron.mass / kinematics.Q2,
			decay.z_min_cutoff
		});
	}

	namespace Evaluation {
		template <typename PDFInterface, typename FFInterface, typename DecayFunction, typename Signature>
		constexpr static double construct(const double x, const double z, const double xi, const double xip, const Parameters<PDFInterface, FFInterface, DecayFunction> &params, const Signature integrand, const int sign, FFInterface &ff1, FFInterface &ff2, const FlavorType flavor1, const FlavorType flavor2, const FlavorType antiflavor1, const FlavorType antiflavor2) {			
			PDFInterface &pdf1 = params.pdf1;
			PDFInterface &pdf2 = params.pdf2;

			const double xg_hat = pdf2.xg();
			const double zg_hat = ff2.xg();

			const double xq = pdf1.xf(flavor1);
			const double zq = ff1.xf(flavor2);

			const double xq_hat = pdf2.xf(flavor1);
			const double zq_hat = ff2.xf(flavor2);

			const double anti_xq = pdf1.xf(antiflavor1);
			const double anti_zq = ff1.xf(antiflavor2);

			const double anti_xq_hat = pdf2.xf(antiflavor1);
			const double anti_zq_hat = ff2.xf(antiflavor2);

			const double xq_zq = xq * zq + sign * anti_xq * anti_zq;
			const double xq_hat_zq = xq_hat * zq + sign * anti_xq_hat * anti_zq;
			const double xq_zq_hat = xq * zq_hat + sign * anti_xq * anti_zq_hat;
			const double xq_hat_zq_hat = xq_hat * zq_hat + sign * anti_xq_hat * anti_zq_hat;
			
			const double xq_zg_hat = xq * zg_hat + sign * anti_xq * zg_hat;
			const double xg_hat_zq = xg_hat * zq + sign * xg_hat * anti_zq;
			const double xq_hat_zg_hat = xq_hat * zg_hat + sign * anti_xq_hat * zg_hat;
			const double xg_hat_zq_hat = xg_hat * zq_hat + sign * xg_hat * anti_zq_hat;

			const double value = integrand(xi, xip, x, z, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);

			return value;
		}
		template <typename PDFInterface, typename FFInterface, typename DecayFunction, typename Signature>
		constexpr static double construct(const double input[], const Parameters<PDFInterface, FFInterface, DecayFunction> &params, const Signature signature, const bool xi_int, const bool xip_int, const bool z_int, const int sign, FFInterface &ff1, FFInterface &ff2, Decay<DecayFunction> &decay, const double z_min) {
			const size_t xi_index = 0;
			const size_t xip_index = size_t(xi_int);
			const size_t z_index = size_t(xi_int) + size_t(xip_int);

			const double xi = xi_int ? input[xi_index] : 1.0;
			const double xip = xip_int ? input[xip_index] : 1.0;

			const TRFKinematics &kinematics = params.kinematics;

			const double x = kinematics.x;
			const double z = z_int ? input[z_index] : params.z;
			const double Q2 = kinematics.Q2;

			const double z_hat = z / xip;

			if (xi_int) {
				if (std::abs(xi - 1) < 1e-15) { return 0; }
			}
			if (xip_int) {
				if (std::abs(xip - 1) < 1e-15) { return 0; }
				ff2.evaluate(z_hat, Q2);
			}

			if (z_int) {
				if (xip < z) { return 0.0; }
				ff1.evaluate(z, Q2);
			}

			const FlavorInfo &flavors = params.flavors;

			const Process &process = params.process;
			const bool positive_W = process.positive_W();

			const FlavorVector &flavors1 = positive_W ? flavors.lower_flavors : flavors.upper_flavors;
			const FlavorVector &flavors2 = positive_W ? flavors.upper_flavors : flavors.lower_flavors;

			double sum = 0.0;
			for (const FlavorType flavor1 : flavors1) {
				const FlavorType antiflavor2 = Flavor::conjugate_flavor(flavor1);
				for (const FlavorType flavor2 : flavors2) {
					const FlavorType antiflavor1 = Flavor::conjugate_flavor(flavor2);

					const double V_ckm = CKM::squared(flavor1, flavor2);
					const double total_value = construct(x, z, xi, xip, params, signature, sign, ff1, ff2, flavor1, flavor2, antiflavor1, antiflavor2);
					const double summand = V_ckm * total_value;

					sum += summand;
				}
			}

			const double decay_value = decay(x, z, Q2, z_min);
			sum *= decay_value;

			return sum;
		}

		template <typename PDFInterface, typename FFInterface, typename DecayFunction = decltype(DecayFunctions::trivial), typename Signature>
		constexpr static double construct(const double input[], void *params_in, const Signature integrand, const bool xi_int, const bool xip_int, const bool z_int, const int sign) {
			const Parameters<PDFInterface, FFInterface, DecayFunction> &params = *static_cast<Parameters<PDFInterface, FFInterface, DecayFunction> *>(params_in);

			if (xi_int) {
				const double xi = input[0];

				const TRFKinematics &kinematics = params.kinematics;

				const double x = kinematics.x;
				const double Q2 = kinematics.Q2;

				const double x_hat = x / xi;

				params.pdf2.evaluate(x_hat, Q2);
			}

			FragmentationConfiguration<FFInterface, DecayFunction> &ffs1 = params.ff1;
			FragmentationConfiguration<FFInterface, DecayFunction> &ffs2 = params.ff2;

			double sum = 0.0;

			for (size_t ff_i = 0; ff_i < ffs1.size(); ff_i++) {
				FFInterface &ff1 = ffs1[ff_i];
				FFInterface &ff2 = ffs2[ff_i];
				auto &decay = ffs1.decays[ff_i];

				const double z_min = SIDISFunctions::compute_z_min(params.kinematics, decay);
				const double value = construct<PDFInterface, FFInterface, DecayFunction>(input, params, integrand, xi_int, xip_int, z_int, sign, ff1, ff2, decay, z_min);
				const double summand = value;

				sum += summand;
			}

			return sum;
		}

		template <typename PDFInterface, typename FFInterface, typename DecayFunction, typename Signature>
		constexpr static double cross_section(const double input[], void *params_in, const Signature F2, const Signature FL, const Signature xF3, const bool xi_int, const bool xip_int, const bool z_int) {
			const struct Parameters<PDFInterface, FFInterface, DecayFunction> &params = *static_cast<Parameters<PDFInterface, FFInterface, DecayFunction> *>(params_in);
			
			const double f2 = construct<PDFInterface, FFInterface, DecayFunction>(input, params_in, F2, xi_int, xip_int, z_int, 1);
			const double fL = construct<PDFInterface, FFInterface, DecayFunction>(input, params_in, FL, xi_int, xip_int, z_int, 1);
			const double xf3 = construct<PDFInterface, FFInterface, DecayFunction>(input, params_in, xF3, xi_int, xip_int, z_int, -1);

			const double cs = CommonFunctions::make_cross_section_variable(params.kinematics.x, params.kinematics.Q2, params.kinematics.s, params.process, f2, fL, xf3);

			return cs;
		}

		// template <typename PDFInterface, typename FFInterface, typename Signature, typename DecayFunction>
		// constexpr static double evaluate_gsl_sidis_integrand(double input[], [[maybe_unused]] size_t dim, void *params_in, Signature integrand, bool xi_int, bool xip_int, bool quark_minus) {
		// 	const struct Parameters<PDFInterface, FFInterface, DecayFunction> *params = static_cast<Parameters<PDFInterface, FFInterface, DecayFunction> *>(params_in);

		// 	const bool xi_xip_int = xi_int && xip_int;

		// 	const size_t xi_index = 0;
		// 	const size_t xip_index = xi_int ? 1 : 0;

		// 	const double xi = xi_int ? input[xi_index] : 1.0;
		// 	const double xip = xip_int ? input[xip_index] : 1.0;

		// 	const double x = params->x;
		// 	const double z = params->z;
		// 	const double Q2 = params->Q2;

		// 	const double x_hat = x / xi;
		// 	const double z_hat = z / xip;
		// 	const double xq_zq = params->xq_zq;

		// 	PDFInterface &pdf1 = params->pdf1;
		// 	FragmentationConfiguration<FFInterface, DecayFunction> &ff1 = params->ff1;
		// 	PDFInterface &pdf2 = params->pdf2;
		// 	FragmentationConfiguration<FFInterface, DecayFunction> &ff2 = params->ff2;

		// 	const Process &process = params->process;

		// 	if (xi_int && std::abs(xi - 1) < 1e-15) { return 0; }
		// 	if (xip_int && std::abs(xip - 1) < 1e-15) { return 0; }

		// 	if (xi_int) { pdf2.evaluate(x_hat, Q2); }
		// 	if (xip_int) { ff2.evaluate(z_hat, Q2); }

		// 	const FlavorInfo &flavors = params->flavors;

		// 	const double xq_hat_zq = xi_int ? PDFCommon::xq_zq_sum(pdf2, ff1, flavors, quark_minus, process) : 0.0;
		// 	const double xq_zq_hat = xip_int ? PDFCommon::xq_zq_sum(pdf1, ff2, flavors, quark_minus, process) : 0.0;
		// 	const double xq_hat_zq_hat = xi_xip_int ? PDFCommon::xq_zq_sum(pdf2, ff2, flavors, quark_minus, process) : 0.0;

		// 	const double xq_zg_hat = xip_int ? PDFCommon::xq_zg_sum(pdf1, ff2, flavors, quark_minus, process) : 0.0;
		// 	const double xg_hat_zq = xi_int ? PDFCommon::xg_zq_sum(pdf2, ff1, flavors, quark_minus, process) : 0.0;
		// 	const double xq_hat_zg_hat = xi_xip_int ? PDFCommon::xq_zg_sum(pdf2, ff2, flavors, quark_minus, process) : 0.0;
		// 	const double xg_hat_zq_hat = xi_xip_int ? PDFCommon::xg_zq_sum(pdf2, ff2, flavors, quark_minus, process) : 0.0;

		// 	const double integrand_value = integrand(
		// 		xi, xip, x, z, 
		// 		xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, 
		// 		xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat
		// 	);
		// 	return integrand_value;
		// }

		// template <typename PDFInterface, typename FFInterface, typename Signature, typename DecayFunction>
		// constexpr static double evaluate_gsl_sidis_decay_cross_section_integrand_ff(const double input[], [[maybe_unused]] const size_t dim, void *params_in, const Signature F2, const Signature FL, const Signature xF3, const Decay<DecayFunction> decay, const double z_min, const bool xi_int, const bool xip_int, const bool z_int, FFInterface &ff1, FFInterface &ff2) {
		// 	const struct Parameters<PDFInterface, FFInterface, DecayFunction> *params = static_cast<Parameters<PDFInterface, FFInterface, DecayFunction> *>(params_in);

		// 	const bool xi_xip_int = xi_int && xip_int;

		// 	const size_t xi_index = 0;
		// 	const size_t xip_index = size_t(xi_int);
		// 	const size_t z_index = size_t(xi_int) + size_t(xip_int);

		// 	const double xi = xi_int ? input[xi_index] : 1.0;
		// 	const double xip = xip_int ? input[xip_index] : 1.0;

		// 	const double x = params->x;
		// 	const double z = z_int ? input[z_index] : params->z;
		// 	const double s = params->s;
		// 	const double Q2 = params->Q2;

		// 	PDFInterface &pdf1 = params->pdf1;
		// 	PDFInterface &pdf2 = params->pdf2;

		// 	const Process &process = params->process;

		// 	const FlavorInfo &flavors = params->flavors;

		// 	const double x_hat = x / xi;
		// 	const double z_hat = z / xip;

		// 	if (z_int) {
		// 		if (xip < z) { return 0.0; }
		// 		ff1.evaluate(z, Q2);
		// 	}
		// 	const double xq_zq = z_int ? PDFCommon::xq_zq_sum(pdf1, ff1, flavors, false, process) : params->xq_zq;
		// 	const double xq_zq_3 = z_int ? PDFCommon::xq_zq_sum(pdf1, ff1, flavors, true, process) : params->xq_zq;

		// 	if (xi_int && std::abs(xi - 1) < 1e-15) { return 0; }
		// 	if (xip_int && std::abs(xip - 1) < 1e-15) { return 0; }

		// 	if (xi_int) { pdf2.evaluate(x_hat, Q2); }
		// 	if (xip_int) { ff2.evaluate(z_hat, Q2); }

		// 	const double xq_hat_zq = xi_int ? PDFCommon::xq_zq_sum(pdf2, ff1, flavors, false, process) : 0.0;
		// 	const double xq_zq_hat = xip_int ? PDFCommon::xq_zq_sum(pdf1, ff2, flavors, false, process) : 0.0;
		// 	const double xq_hat_zq_hat = xi_xip_int ? PDFCommon::xq_zq_sum(pdf2, ff2, flavors, false, process) : 0.0;

		// 	const double xq_zg_hat = xip_int ? PDFCommon::xq_zg_sum(pdf1, ff2, flavors, false, process) : 0.0;
		// 	const double xg_hat_zq = xi_int ? PDFCommon::xg_zq_sum(pdf2, ff1, flavors, false, process) : 0.0;
		// 	const double xq_hat_zg_hat = xi_xip_int ? PDFCommon::xq_zg_sum(pdf2, ff2, flavors, false, process) : 0.0;
		// 	const double xg_hat_zq_hat = xi_xip_int ? PDFCommon::xg_zq_sum(pdf2, ff2, flavors, false, process) : 0.0;


		// 	const double xq_hat_zq_3 = xi_int ? PDFCommon::xq_zq_sum(pdf2, ff1, flavors, true, process) : 0.0;
		// 	const double xq_zq_hat_3 = xip_int ? PDFCommon::xq_zq_sum(pdf1, ff2, flavors, true, process) : 0.0;
		// 	const double xq_hat_zq_hat_3 = xi_xip_int ? PDFCommon::xq_zq_sum(pdf2, ff2, flavors, true, process) : 0.0;

		// 	const double xq_zg_hat_3 = xip_int ? PDFCommon::xq_zg_sum(pdf1, ff2, flavors, true, process) : 0.0;
		// 	const double xg_hat_zq_3 = xi_int ? PDFCommon::xg_zq_sum(pdf2, ff1, flavors, true, process) : 0.0;
		// 	const double xq_hat_zg_hat_3 = xi_xip_int ? PDFCommon::xq_zg_sum(pdf2, ff2, flavors, true, process) : 0.0;
		// 	const double xg_hat_zq_hat_3 = xi_xip_int ? PDFCommon::xg_zq_sum(pdf2, ff2, flavors, true, process) : 0.0;

		// 	const double f2 = F2(
		// 		xi, xip, x, z, 
		// 		xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, 
		// 		xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat
		// 	);
		// 	const double fL = FL(
		// 		xi, xip, x, z, 
		// 		xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, 
		// 		xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat
		// 	);
		// 	const double xf3 = xF3(
		// 		xi, xip, x, z, 
		// 		xq_zq_3, xq_hat_zq_3, xq_zq_hat_3, xq_hat_zq_hat_3, 
		// 		xq_zg_hat_3, xg_hat_zq_3, xq_hat_zg_hat_3, xg_hat_zq_hat_3
		// 	);

		// 	const double integrand_value = CommonFunctions::make_cross_section_variable(x, Q2, s, process, f2, fL, xf3);
		// 	const double decay_function_value = decay(x, z, Q2, z_min);
		// 	const double final_value = decay_function_value * integrand_value;

		// 	return final_value;
		// }
		// template <typename PDFInterface, typename FFInterface, typename Signature, typename DecayFunction>
		// constexpr static double evaluate_gsl_sidis_decay_cross_section_integrand(const double input[], [[maybe_unused]] const size_t dim, void *params_in, const Signature F2, const Signature FL, const Signature xF3, const Decay<DecayFunction> decay, const double z_min, const bool xi_int, const bool xip_int, const bool z_int) {
		// 	const struct Parameters<PDFInterface, FFInterface, DecayFunction> *params = static_cast<Parameters<PDFInterface, FFInterface, DecayFunction> *>(params_in);
		// 	FragmentationConfiguration<FFInterface, DecayFunction> &ffs1 = params->ff1;
		// 	FragmentationConfiguration<FFInterface, DecayFunction> &ffs2 = params->ff2;

		// 	double sum = 0.0;
		// 	for (size_t i = 0; i < ffs1.size(); i++) {
		// 		FFInterface &ff1 = ffs1[i];
		// 		FFInterface &ff2 = ffs2[i];

		// 		const double value = evaluate_gsl_sidis_decay_cross_section_integrand_ff<PDFInterface, FFInterface>(input, dim, params_in, F2, FL, xF3, decay, z_min, xi_int, xip_int, z_int, ff1, ff2);
		// 		const double coefficient = ffs1.coefficients[i];

		// 		sum += coefficient * value;
		// 	}

		// 	return sum;
		// }
		// template <typename PDFInterface, typename FFInterface, typename Signature>
		// static double evaluate_gsl_sidis_cross_section_integrand(double input[], size_t dim, void *params_in, Signature F2, Signature FL, Signature xF3, bool xi_int, bool xip_int, bool z_int) {
		// 	const Decay decay = Decay(DecayFunctions::trivial);
		// 	return evaluate_gsl_sidis_decay_cross_section_integrand<PDFInterface, FFInterface>(input, dim, params_in, F2, FL, xF3, decay, 0.0, xi_int, xip_int, z_int);
		// }
	}

	namespace Integrands {
		static constexpr double F2_lo_integrand(
			[[maybe_unused]] const double xi, 
			[[maybe_unused]] const double xip, 
			[[maybe_unused]] const double x, 
			const double z, 
			const double xq_zq,
			[[maybe_unused]] const double xq_hat_zq, 
			[[maybe_unused]] const double xq_zq_hat,
			[[maybe_unused]] const double xq_hat_zq_hat, 
			[[maybe_unused]] const double xq_zg_hat,
			[[maybe_unused]] const double xg_hat_zq,
			[[maybe_unused]] const double xq_hat_zg_hat,
			[[maybe_unused]] const double xg_hat_zq_hat) {

			return 2 * xq_zq / z;
		}
		static constexpr double FL_lo_integrand(
			[[maybe_unused]] const double xi, 
			[[maybe_unused]] const double xip, 
			[[maybe_unused]] const double x, 
			[[maybe_unused]] const double z, 
			[[maybe_unused]] const double xq_zq,
			[[maybe_unused]] const double xq_hat_zq, 
			[[maybe_unused]] const double xq_zq_hat,
			[[maybe_unused]] const double xq_hat_zq_hat, 
			[[maybe_unused]] const double xq_zg_hat,
			[[maybe_unused]] const double xg_hat_zq,
			[[maybe_unused]] const double xq_hat_zg_hat,
			[[maybe_unused]] const double xg_hat_zq_hat) {

			return 0.0;
		}
		static constexpr double F3_lo_integrand(
			[[maybe_unused]] const double xi, 
			[[maybe_unused]] const double xip, 
			[[maybe_unused]] const double x, 
			const double z, 
			[[maybe_unused]] const double xq_zq,
			[[maybe_unused]] const double xq_hat_zq, 
			[[maybe_unused]] const double xq_zq_hat,
			[[maybe_unused]] const double xq_hat_zq_hat, 
			[[maybe_unused]] const double xq_zg_hat,
			[[maybe_unused]] const double xg_hat_zq,
			[[maybe_unused]] const double xq_hat_zg_hat,
			[[maybe_unused]] const double xg_hat_zq_hat) {
			return 2 * xq_zq / z;
		}
		static constexpr double F2_delta_integrand(
			[[maybe_unused]] const double xi, 
			[[maybe_unused]] const double xip, 
			const double x, 
			const double z, 
			const double xq_zq,
			[[maybe_unused]] const double xq_hat_zq, 
			[[maybe_unused]] const double xq_zq_hat,
			[[maybe_unused]] const double xq_hat_zq_hat, 
			[[maybe_unused]] const double xq_zg_hat,
			[[maybe_unused]] const double xg_hat_zq,
			[[maybe_unused]] const double xq_hat_zg_hat,
			[[maybe_unused]] const double xg_hat_zq_hat) {

			return SIDISFunctions::delta_contribution(x, z) * xq_zq;
		}
		static constexpr double FL_delta_integrand(
			[[maybe_unused]] const double xi, 
			[[maybe_unused]] const double xip, 
			[[maybe_unused]] const double x, 
			[[maybe_unused]] const double z, 
			[[maybe_unused]] const double xq_zq,
			[[maybe_unused]] const double xq_hat_zq, 
			[[maybe_unused]] const double xq_zq_hat,
			[[maybe_unused]] const double xq_hat_zq_hat, 
			[[maybe_unused]] const double xq_zg_hat,
			[[maybe_unused]] const double xg_hat_zq,
			[[maybe_unused]] const double xq_hat_zg_hat,
			[[maybe_unused]] const double xg_hat_zq_hat) {

			return 0.0;
		}
		static constexpr double F3_delta_integrand(
			[[maybe_unused]] const double xi, 
			[[maybe_unused]] const double xip, 
			const double x, 
			const double z, 
			const double xq_zq,
			[[maybe_unused]] const double xq_hat_zq, 
			[[maybe_unused]] const double xq_zq_hat,
			[[maybe_unused]] const double xq_hat_zq_hat, 
			[[maybe_unused]] const double xq_zg_hat,
			[[maybe_unused]] const double xg_hat_zq,
			[[maybe_unused]] const double xq_hat_zg_hat,
			[[maybe_unused]] const double xg_hat_zq_hat) {

			return SIDISFunctions::delta_contribution(x, z) * xq_zq;
		}
		static constexpr double F2_xi_integrand(
			const double xi, 
			[[maybe_unused]] const double xip, 
			[[maybe_unused]] const double x, 
			const double z, 
			const double xq_zq,
			const double xq_hat_zq, 
			[[maybe_unused]] const double xq_zq_hat,
			[[maybe_unused]] const double xq_hat_zq_hat, 
			[[maybe_unused]] const double xq_zg_hat,
			const double xg_hat_zq,
			[[maybe_unused]] const double xq_hat_zg_hat,
			[[maybe_unused]] const double xg_hat_zq_hat) {

			const double log_term = std::log((1 - xi) / xi);
			const double logm1 = Utility::logm1(z);

			const double term1 = (1 - xi) * (1 + log_term + logm1) - (2 * xi * std::log(xi)) / (1 - xi);
			const double term2 = 2 * (Utility::logm1(xi) + logm1);
			const double term3 = (xi * xq_hat_zq - xq_zq)  / (1 - xi);

			const double quark_contribution = Constants::C_F * (xq_hat_zq * term1 + term2 * term3);

			const double term4 = 1 - (std::pow(xi, 2) + std::pow(1 - xi, 2)) * (1 - log_term);
			const double term5 = logm1 * (1 - 2 * xi * (1 - xi));

			const double gluon_contribution = Constants::T_R * xg_hat_zq * (term4 + term5);

			return 2 * (quark_contribution + gluon_contribution) / z;
		}
		static constexpr double F3_xi_integrand(
			const double xi, 
			const double xip, 
			const double x, 
			const double z, 
			const double xq_zq,
			const double xq_hat_zq, 
			const double xq_zq_hat,
			const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			const double xg_hat_zq,
			const double xq_hat_zg_hat,
			const double xg_hat_zq_hat) {
			return F2_xi_integrand(xi, xip, x, z, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
		}

		static constexpr double F2_xip_integrand(
			[[maybe_unused]] const double xi, 
			const double xip, 
			const double x, 
			const double z, 
			const double xq_zq,
			[[maybe_unused]] const double xq_hat_zq, 
			const double xq_zq_hat,
			[[maybe_unused]] const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			[[maybe_unused]] const double xg_hat_zq,
			[[maybe_unused]] const double xq_hat_zg_hat,
			[[maybe_unused]] const double xg_hat_zq_hat) {

			const double log_term = std::log(xip * (1 - xip));
			const double logm1 = Utility::logm1(x);
			
			const double term1 = (1 - xip) * (1 + log_term + logm1) + (2 * xip * std::log(xip)) / (1 - xip);
			const double term2 = 2 * (Utility::logm1(xip) + logm1);
			const double term3 = (xip * xq_zq_hat - xq_zq)  / (1 - xip);

			const double quark_contribution = Constants::C_F * (xq_zq_hat * term1 + term2 * term3);

			const double term4 = xip + log_term * (1 + std::pow(1 - xip, 2)) / xip;
			const double term5 = logm1 * (xip + 2 * (1 - xip) / xip);

			const double gluon_contribution = Constants::C_F * xq_zg_hat * (term4 + term5);

			return 2 * (quark_contribution + gluon_contribution) / z;
		}
		static constexpr double F3_xip_integrand(
			const double xi, 
			const double xip, 
			const double x, 
			const double z, 
			const double xq_zq,
			const double xq_hat_zq, 
			const double xq_zq_hat,
			const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			const double xg_hat_zq,
			const double xq_hat_zg_hat,
			const double xg_hat_zq_hat) {
			return F2_xip_integrand(xi, xip, x, z, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
		}

		static constexpr double F2_xi_xip_integrand(const double xi, 
		const double xip, 
			[[maybe_unused]] const double x, 
			const double z, 
		const double xq_zq,
		const double xq_hat_zq, 
		const double xq_zq_hat,
		const double xq_hat_zq_hat, 
		const double xq_zg_hat,
		const double xg_hat_zq,
		const double xq_hat_zg_hat,
		const double xg_hat_zq_hat) {
			const double term1 = (1 - xi) / (1 - xip) * (xq_hat_zq_hat - xq_hat_zq);
			const double term2 = (1 - xip) / (1 - xi) * (xq_hat_zq_hat - xq_zq_hat);
			const double term3 = 2 * (xi * xip * xq_hat_zq_hat - xi * xq_hat_zq - xip * xq_zq_hat + xq_zq) / ((1 - xi) * (1 - xip));
			const double term4 = 6 * xi * xip * xq_hat_zq_hat;

			const double quark_contribution = Constants::C_F * (term1 + term2 + term3 + term4);

			const double term5 = xq_hat_zg_hat * (6 * xi * (1 - xip) + (1 - xi) / xip);
			const double term6 = (xq_hat_zg_hat * (xip + 2 * xi * (1 - xip) / xip) - xq_zg_hat * (xip + 2 * (1 - xip) / xip)) / (1 - xi);
			const double gluon_contribution_1 = Constants::C_F * (term5 + term6);

			const double term7 = xg_hat_zq_hat * (12 * xi * (1 - xi) + (1 - xip) / xip);
			const double term8 = (xg_hat_zq_hat * (xip - 2 * xi * (1 - xi) / xip) - xg_hat_zq * (1 - 2 * xi * (1 - xi))) / (1 - xip);

			const double gluon_contribution_2 = Constants::T_R * (term7 + term8);

			return 2 * (quark_contribution + gluon_contribution_1 + gluon_contribution_2) / z;
		}
		static constexpr double FL_xi_integrand(
			[[maybe_unused]] const double xi, 
			[[maybe_unused]] const double xip, 
			[[maybe_unused]] const double x, 
			[[maybe_unused]] const double z, 
			[[maybe_unused]] const double xq_zq,
			[[maybe_unused]] const double xq_hat_zq, 
			[[maybe_unused]] const double xq_zq_hat,
			[[maybe_unused]] const double xq_hat_zq_hat, 
			[[maybe_unused]] const double xq_zg_hat,
			[[maybe_unused]] const double xg_hat_zq,
			[[maybe_unused]] const double xq_hat_zg_hat,
			[[maybe_unused]] const double xg_hat_zq_hat) {

			return 0.0;
		}
		static constexpr double FL_xip_integrand(
			[[maybe_unused]] const double xi, 
			[[maybe_unused]] const double xip, 
			[[maybe_unused]] const double x, 
			[[maybe_unused]] const double z, 
			[[maybe_unused]] const double xq_zq,
			[[maybe_unused]] const double xq_hat_zq, 
			[[maybe_unused]] const double xq_zq_hat,
			[[maybe_unused]] const double xq_hat_zq_hat, 
			[[maybe_unused]] const double xq_zg_hat,
			[[maybe_unused]] const double xg_hat_zq,
			[[maybe_unused]] const double xq_hat_zg_hat,
			[[maybe_unused]] const double xg_hat_zq_hat) {

			return 0.0;
		}
		static constexpr double FL_xi_xip_integrand(const double xi, 
		const double xip, 
			[[maybe_unused]] const double x, 
			const double z, 
		[[maybe_unused]] const double xq_zq,
		[[maybe_unused]] const double xq_hat_zq, 
		[[maybe_unused]] const double xq_zq_hat,
		const double xq_hat_zq_hat, 
		[[maybe_unused]] const double xq_zg_hat,
		[[maybe_unused]] const double xg_hat_zq,
		const double xq_hat_zg_hat,
		const double xg_hat_zq_hat) {
			const double term1 = xq_hat_zq_hat * Constants::C_F * 4 * xi * xip;
			const double term2 = xq_hat_zg_hat * Constants::C_F * 4 * xi * (1 - xip);
			const double term3 = xg_hat_zq_hat * Constants::T_R * 8 * xi * (1 - xi);
			
			const double value = term1 + term2 + term3;
			return 2 * value / z;
		}

		static constexpr double F3_xi_xip_integrand(const double xi, 
		const double xip, 
			const double x, 
			const double z, 
		const double xq_zq,
		const double xq_hat_zq, 
		const double xq_zq_hat,
		const double xq_hat_zq_hat, 
		const double xq_zg_hat,
		const double xg_hat_zq,
		const double xq_hat_zg_hat,
		const double xg_hat_zq_hat) {
			const double F2_value = F2_xi_xip_integrand(xi, xip, x, z, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
			const double term1 = Constants::C_F * xq_hat_zq_hat * (6 * xi * xip + 2 * (1 - xi - xip));
			const double term2 = Constants::C_F * xq_hat_zg_hat * (4 * xi * (1 - xip) + 2 * (1 - xi) * xip);
			const double term3 = Constants::T_R * xg_hat_zq_hat * (12 * xi * (1 - xi) + 2 * (1 - 2 * xi * (1 - xi)) / xip - 2);

			const double value = F2_value - 2 * (term1 + term2 + term3) / z;

			return value;
		}
	}

	// Signature make_cross_section_integrand(const double x, const double Q2, const double s, Process process, Signature F2, Signature FL, Signature xF3) {
		// const double prefactor = CommonFunctions::cross_section_prefactor(Q2);
		// const std::optional<double> y = CommonFunctions::compute_y(x, Q2, s);

		// const double term1 = y.has_value() ? 1 - *y + 0.5 * *y * *y : 0.0;
		// const double term2 = y.has_value() ? - 0.5 * *y * *y : 0.0;
		// const double term3 = y.has_value() ? *y * (1 - 0.5 * *y) : 0.0;

	// 	Signature integrand = [&](const double xi, 
	// 	const double xip, 
	// 	const double x, 
	// 	const double z, 
	// 	const double xq_zq,
	// 	const double xq_hat_zq, 
	// 	const double xq_zq_hat,
	// 	const double xq_hat_zq_hat, 
	// 	const double xq_zg_hat,
	// 	const double xg_hat_zq,
	// 	const double xq_hat_zg_hat,
	// 	const double xg_hat_zq_hat) {
	// 		const double f2 = F2(xi, xip, x, z, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
	// 		const double fL = FL(xi, xip, x, z, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
	// 		const double xf3 = xF3(xi, xip, x, z, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);

	// 		const double cs = prefactor * (term1 * f2 + term2 * fL + double(process.W_sign()) * term3 * xf3) / x;

	// 		return cs;
	// 	};

	// 	return integrand;
	// }
}

#endif