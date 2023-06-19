#ifndef SIDIS_FUNCTIONS_H
#define SIDIS_FUNCTIONS_H

#include "Constants.cpp"
#include <cmath>
#include "Flavor.cpp"
#include "Decay.cpp"
#include "FragmentationConfiguration.cpp"
#include "TRFKinematics.cpp"
#include "DecayFunctions.cpp"
#include "CKM.cpp"

namespace SIDISFunctions {
	template <typename PDFInterface, typename FFInterface, typename DecayFunction>
	struct Parameters {
		const PDFInterface &pdf1;
		const FragmentationConfiguration<FFInterface, DecayFunction> &ff1;
		const PDFInterface &pdf2;
		const FragmentationConfiguration<FFInterface, DecayFunction> &ff2;

		const FlavorInfo &flavors;

		const double nlo_coefficient;

		const Process &process;
		const TRFKinematics &kinematics;

		const double z;

		const double factorization_scale;
		const double fragmentation_scale;

		const double factorization_scale_log;
		const double fragmentation_scale_log;

		// Could some caching here
		// const double xq;
		// const double zq;
	};

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
		constexpr static double construct(const double x, const double z, const double xi, const double xip, const Parameters<PDFInterface, FFInterface, DecayFunction> &params, const Signature integrand, const int sign, const FFInterface &ff1, const FFInterface &ff2, const FlavorType flavor1, const FlavorType flavor2, const FlavorType antiflavor1, const FlavorType antiflavor2) {			
			const PDFInterface &pdf1 = params.pdf1;
			const PDFInterface &pdf2 = params.pdf2;

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

			const double value = integrand(xi, xip, x, z, params.factorization_scale_log, params.fragmentation_scale_log, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);

			return value;
		}
		template <typename PDFInterface, typename FFInterface, typename DecayFunction, typename Signature>
		constexpr static double construct(const double input[], const Parameters<PDFInterface, FFInterface, DecayFunction> &params, const Signature signature, const bool xi_int, const bool xip_int, const bool z_int, const int sign, const FFInterface &ff1, const FFInterface &ff2, const Decay<DecayFunction> &decay, const double z_min) {
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

			const double factorization_scale = params.factorization_scale;
			const double fragmentation_scale = params.fragmentation_scale;

			if (xi_int) {
				// if (xi < x) { return 0; }
				if (xi == 1.0) { return 0.0; }
			}
			if (z_int) {
				if (xip < z) { return 0.0; }
				ff1.evaluate(z, fragmentation_scale);
			}
			if (xip_int) {
				// if (xip < z) {
				// 	if (!z_int) {
				// 		std::cout << xip << std::endl; 
				// 	}
				// 	return 0; 
				// 	}
				if (xip == 1.0) { return 0.0; }
				ff2.evaluate(z_hat, fragmentation_scale);
			}

			const FlavorInfo &flavors = params.flavors;

			const Process &process = params.process;
			const bool positive_W = process.positive_W();

			const FlavorVector &flavors1 = positive_W ? flavors.lower_flavors : flavors.upper_flavors;
			const FlavorVector &flavors2 = positive_W ? flavors.upper_flavors : flavors.lower_flavors;

			const PDFInterface &pdf1 = params.pdf1;
			const PDFInterface &pdf2 = params.pdf2;

			double sum = 0.0;
			for (const FlavorType outgoing : flavors2) {
				const FlavorType anti_incoming = Flavor::conjugate_flavor(outgoing);
				for (const FlavorType incoming : flavors1) {
					const FlavorType anti_outgoing = Flavor::conjugate_flavor(incoming);

					const double x_mass = CommonFunctions::compute_momentum_fraction_mass_correction(x, Q2, flavors.mass(Flavor::Charm), 0.0);

					if (xi_int && xi < x_mass) { continue; }
					pdf1.evaluate(x_mass, factorization_scale);
					pdf2.evaluate(x_mass / xi, factorization_scale);

					const double V_ckm = CKM::squared(incoming, outgoing);
					const double total_value = construct(x_mass, z, xi, xip, params, signature, sign, ff1, ff2, incoming, outgoing, anti_incoming, anti_outgoing);
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
			const TRFKinematics &kinematics = params.kinematics;

			if (xi_int) {
				const double xi = input[0];
				const double x = kinematics.x;
				const double x_hat = x / xi;
				params.pdf2.evaluate(x_hat, params.factorization_scale);
			}

			const FragmentationConfiguration<FFInterface, DecayFunction> &ffs1 = params.ff1;
			const FragmentationConfiguration<FFInterface, DecayFunction> &ffs2 = params.ff2;

			double sum = 0.0;

			for (size_t ff_i = 0; ff_i < ffs1.size(); ff_i++) {
				const FFInterface &ff1 = ffs1[ff_i];
				const FFInterface &ff2 = ffs2[ff_i];
				const auto &decay = ffs1.decays[ff_i];

				const double z_min = SIDISFunctions::compute_z_min(kinematics, decay);
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
	}

	constexpr double delta_contribution(const double x, const double z) {
		return 2 * Constants::C_F * (std::pow(std::log(1 - x) + std::log(1 - z), 2) - 8) / (x * z);
	}

	namespace Integrands {
		static constexpr double F2x_lo_integrand(
			[[maybe_unused]] const double xi, 
			[[maybe_unused]] const double xip, 
			const double x, 
			const double z,
			[[maybe_unused]] const double factorization_scale_log,
			[[maybe_unused]] const double fragmentation_scale_log,
			const double xq_zq,
			[[maybe_unused]] const double xq_hat_zq, 
			[[maybe_unused]] const double xq_zq_hat,
			[[maybe_unused]] const double xq_hat_zq_hat, 
			[[maybe_unused]] const double xq_zg_hat,
			[[maybe_unused]] const double xg_hat_zq,
			[[maybe_unused]] const double xq_hat_zg_hat,
			[[maybe_unused]] const double xg_hat_zq_hat) {

			return 2 * xq_zq / (x * z);
		}
		static constexpr double FLx_lo_integrand(
			[[maybe_unused]] const double xi, 
			[[maybe_unused]] const double xip, 
			[[maybe_unused]] const double x, 
			[[maybe_unused]] const double z, 
			[[maybe_unused]] const double factorization_scale_log,
			[[maybe_unused]] const double fragmentation_scale_log,
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
			const double x, 
			const double z, 
			[[maybe_unused]] const double factorization_scale_log,
			[[maybe_unused]] const double fragmentation_scale_log,
			[[maybe_unused]] const double xq_zq,
			[[maybe_unused]] const double xq_hat_zq, 
			[[maybe_unused]] const double xq_zq_hat,
			[[maybe_unused]] const double xq_hat_zq_hat, 
			[[maybe_unused]] const double xq_zg_hat,
			[[maybe_unused]] const double xg_hat_zq,
			[[maybe_unused]] const double xq_hat_zg_hat,
			[[maybe_unused]] const double xg_hat_zq_hat) {
			return 2 * xq_zq / (x * z);
		}


		static constexpr double scale_delta_integrand(
			[[maybe_unused]] const double xi,
			[[maybe_unused]] const double xip,
			const double x, 
			const double z, 
			const double factorization_scale_log,
			const double fragmentation_scale_log,
			const double xq_zq,
			[[maybe_unused]] const double xq_hat_zq, 
			[[maybe_unused]] const double xq_zq_hat,
			[[maybe_unused]] const double xq_hat_zq_hat, 
			[[maybe_unused]] const double xq_zg_hat,
			[[maybe_unused]] const double xg_hat_zq,
			[[maybe_unused]] const double xq_hat_zg_hat,
			[[maybe_unused]] const double xg_hat_zq_hat) {
			return 2 * Constants::C_F * xq_zq * ((2 * std::log(1 - z) + 1.5) * fragmentation_scale_log + (2 * std::log(1 - x) + 1.5) * factorization_scale_log) / (x * z);
		}
		static constexpr double scale_xi_integrand(
			const double xi,
			[[maybe_unused]] const double xip,
			const double x, 
			const double z, 
			const double factorization_scale_log,
			[[maybe_unused]] const double fragmentation_scale_log,
			const double xq_zq,
			const double xq_hat_zq, 
			[[maybe_unused]] const double xq_zq_hat,
			[[maybe_unused]] const double xq_hat_zq_hat, 
			[[maybe_unused]] const double xq_zg_hat,
			const double xg_hat_zq,
			[[maybe_unused]] const double xq_hat_zg_hat,
			[[maybe_unused]] const double xg_hat_zq_hat) {
			const double term1 = Constants::C_F * ((1 + xi * xi) * xq_hat_zq - 2 * xq_zq) / (1 - xi);
			const double term2 = Constants::T_R * xg_hat_zq * (std::pow(xi, 2) + std::pow(1 - xi, 2));
			const double result = 2 * factorization_scale_log * (term1 + term2) / (x * z);
			return result;
		}
		static constexpr double scale_xip_integrand(
			[[maybe_unused]] const double xi,
			const double xip,
			const double x, 
			const double z, 
			[[maybe_unused]] const double factorization_scale_log,
			const double fragmentation_scale_log,
			const double xq_zq,
			[[maybe_unused]] const double xq_hat_zq, 
			const double xq_zq_hat,
			[[maybe_unused]] const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			[[maybe_unused]] const double xg_hat_zq,
			[[maybe_unused]] const double xq_hat_zg_hat,
			[[maybe_unused]] const double xg_hat_zq_hat) {
			const double term1 = Constants::C_F * ((1 + std::pow(xip, 2)) * xq_zq_hat - 2 * xq_zq) / (1 - xip);
			const double term2 = Constants::C_F * xq_zg_hat * (1 + std::pow(1 - xip, 2)) / xip;
			const double result = 2 * fragmentation_scale_log * (term1 + term2) / (x * z);
			return result;
		}
		static constexpr double F2x_delta_integrand(
			[[maybe_unused]] const double xi, 
			[[maybe_unused]] const double xip, 
			const double x, 
			const double z, 
			const double factorization_scale_log,
			const double fragmentation_scale_log,
			const double xq_zq,
			[[maybe_unused]] const double xq_hat_zq, 
			[[maybe_unused]] const double xq_zq_hat,
			[[maybe_unused]] const double xq_hat_zq_hat, 
			[[maybe_unused]] const double xq_zg_hat,
			[[maybe_unused]] const double xg_hat_zq,
			[[maybe_unused]] const double xq_hat_zg_hat,
			[[maybe_unused]] const double xg_hat_zq_hat) {

			const double term1 = SIDISFunctions::delta_contribution(x, z) * xq_zq;
			const double term2 = (factorization_scale_log == 0 && fragmentation_scale_log == 0) ? 0 : scale_delta_integrand(xi, xip, x, z, factorization_scale_log, fragmentation_scale_log, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
			return term1 + term2;
		}
		static constexpr double FLx_delta_integrand(
			[[maybe_unused]] const double xi, 
			[[maybe_unused]] const double xip, 
			[[maybe_unused]] const double x, 
			[[maybe_unused]] const double z, 
			[[maybe_unused]] const double fragmentation_scale_log,
			[[maybe_unused]] const double factorization_scale_log,
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
			const double xi, 
			const double xip, 
			const double x, 
			const double z, 
			const double factorization_scale_log,
			const double fragmentation_scale_log,
			const double xq_zq,
			const double xq_hat_zq, 
			const double xq_zq_hat,
			const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			const double xg_hat_zq,
			const double xq_hat_zg_hat,
			const double xg_hat_zq_hat) {

			return F2x_delta_integrand(xi, xip, x, z, factorization_scale_log, fragmentation_scale_log, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
		}
		static constexpr double F2x_xi_integrand(
			const double xi, 
			[[maybe_unused]] const double xip, 
			const double x, 
			const double z, 
			const double factorization_scale_log,
			const double fragmentation_scale_log,
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

			const double coefficient_contribution = 2 * (quark_contribution + gluon_contribution) / (x * z);
			const double scale_contribution = (factorization_scale_log == 0 && fragmentation_scale_log == 0) ? 0 : scale_xi_integrand(xi, xip, x, z, factorization_scale_log, fragmentation_scale_log, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);

			return coefficient_contribution + scale_contribution;
		}
		static constexpr double F3_xi_integrand(
			const double xi, 
			const double xip, 
			const double x, 
			const double z, 
			const double factorization_scale_log,
			const double fragmentation_scale_log,
			const double xq_zq,
			const double xq_hat_zq, 
			const double xq_zq_hat,
			const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			const double xg_hat_zq,
			const double xq_hat_zg_hat,
			const double xg_hat_zq_hat) {
			return F2x_xi_integrand(xi, xip, x, z, factorization_scale_log, fragmentation_scale_log, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
		}

		static constexpr double F2x_xip_integrand(
			[[maybe_unused]] const double xi, 
			const double xip, 
			const double x, 
			const double z, 
			const double factorization_scale_log,
			const double fragmentation_scale_log,
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

			const double coefficient_contribution = 2 * (quark_contribution + gluon_contribution) / (x * z);
			const double scale_contribution = (factorization_scale_log == 0 && fragmentation_scale_log == 0) ? 0 : scale_xip_integrand(xi, xip, x, z, factorization_scale_log, fragmentation_scale_log, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);

			return coefficient_contribution + scale_contribution;
		}
		static constexpr double F3_xip_integrand(
			const double xi, 
			const double xip, 
			const double x, 
			const double z, 
			const double factorization_scale_log,
			const double fragmentation_scale_log,
			const double xq_zq,
			const double xq_hat_zq, 
			const double xq_zq_hat,
			const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			const double xg_hat_zq,
			const double xq_hat_zg_hat,
			const double xg_hat_zq_hat) {
			return F2x_xip_integrand(xi, xip, x, z, factorization_scale_log, fragmentation_scale_log, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
		}

		static constexpr double F2x_xi_xip_integrand(const double xi, 
			const double xip, 
			const double x, 
			const double z, 
			[[maybe_unused]] const double factorization_scale_log,
			[[maybe_unused]] const double fragmentation_scale_log,
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

			return 2 * (quark_contribution + gluon_contribution_1 + gluon_contribution_2) / (x * z);
		}
		static constexpr double FLx_xi_integrand(
			[[maybe_unused]] const double xi, 
			[[maybe_unused]] const double xip, 
			[[maybe_unused]] const double x, 
			[[maybe_unused]] const double z, 
			[[maybe_unused]] const double fragmentation_scale_log,
			[[maybe_unused]] const double factorization_scale_log,
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
		static constexpr double FLx_xip_integrand(
			[[maybe_unused]] const double xi, 
			[[maybe_unused]] const double xip, 
			[[maybe_unused]] const double x, 
			[[maybe_unused]] const double z, 
			[[maybe_unused]] const double factorization_scale_log,
			[[maybe_unused]] const double fragmentation_scale_log,
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
		static constexpr double FLx_xi_xip_integrand(const double xi, 
			const double xip, 
			const double x, 
			const double z, 
			[[maybe_unused]] const double factorization_scale_log,
			[[maybe_unused]] const double fragmentation_scale_log,
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
			return 2 * value / (x * z);
		}

		static constexpr double F3_xi_xip_integrand(const double xi, 
			const double xip, 
			const double x, 
			const double z, 
			const double factorization_scale_log,
			const double fragmentation_scale_log,
			const double xq_zq,
			const double xq_hat_zq, 
			const double xq_zq_hat,
			const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			const double xg_hat_zq,
			const double xq_hat_zg_hat,
			const double xg_hat_zq_hat) {
			const double F2_value = F2x_xi_xip_integrand(xi, xip, x, z, factorization_scale_log, fragmentation_scale_log, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
			const double term1 = Constants::C_F * xq_hat_zq_hat * (6 * xi * xip + 2 * (1 - xi - xip));
			const double term2 = Constants::C_F * xq_hat_zg_hat * (4 * xi * (1 - xip) + 2 * (1 - xi) * xip);
			const double term3 = Constants::T_R * xg_hat_zq_hat * (12 * xi * (1 - xi) + 2 * (1 - 2 * xi * (1 - xi)) / xip - 2);

			const double value = F2_value - 2 * (term1 + term2 + term3) / (x * z);

			return value;
		}

		template <typename Signature>
		static constexpr double make_nlo_integrand(
			const Signature delta_integrand_f,
			const Signature xi_integrand_f,
			const Signature xip_integrand_f,
			const Signature xi_xip_integrand_f,
			const double xi,
			const double xip, 
			const double x, 
			const double z, 
			const double factorization_scale_log,
			const double fragmentation_scale_log,
			const double xq_zq,
			const double xq_hat_zq, 
			const double xq_zq_hat,
			const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			const double xg_hat_zq,
			const double xq_hat_zg_hat,
			const double xg_hat_zq_hat) {
			
			const double delta_integrand = delta_integrand_f(xi, xip, x, z, factorization_scale_log, fragmentation_scale_log, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
			const double xi_integrand = xi_integrand_f(xi, xip, x, z, factorization_scale_log, fragmentation_scale_log, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
			const double xip_integrand = xip_integrand_f(xi, xip, x, z, factorization_scale_log, fragmentation_scale_log, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
			const double xi_xip_integrand = xi_xip_integrand_f(xi, xip, x, z, factorization_scale_log, fragmentation_scale_log, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);

			// const double integrand = Utility::kahan_sum({delta_integrand, (1.0 - x) * xi_integrand, (1.0 - z) * xip_integrand, (1.0 - x) * (1.0 - z) * xi_xip_integrand});
			// const double integrand = delta_integrand / ((1.0 - x) * (1.0 - z)) + xi_integrand / (1.0 - z) + xip_integrand / (1.0 - x) + xi_xip_integrand;
			const double integrand = delta_integrand + (1.0 - x) * xi_integrand + (1.0 - z) * xip_integrand + (1.0 - x) * (1.0 - z) * xi_xip_integrand;
			return integrand;
		}

		static constexpr double F2x_nlo_integrand(
			const double xi,
			const double xip, 
			const double x, 
			const double z, 
			const double factorization_scale_log,
			const double fragmentation_scale_log,
			const double xq_zq,
			const double xq_hat_zq, 
			const double xq_zq_hat,
			const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			const double xg_hat_zq,
			const double xq_hat_zg_hat,
			const double xg_hat_zq_hat) {
			
			return make_nlo_integrand(F2x_delta_integrand, F2x_xi_integrand, F2x_xip_integrand, F2x_xi_xip_integrand, xi, xip, x, z, factorization_scale_log, fragmentation_scale_log, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
		}
		static constexpr double FLx_nlo_integrand(
			const double xi,
			const double xip, 
			const double x, 
			const double z, 
			const double factorization_scale_log,
			const double fragmentation_scale_log,
			const double xq_zq,
			const double xq_hat_zq, 
			const double xq_zq_hat,
			const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			const double xg_hat_zq,
			const double xq_hat_zg_hat,
			const double xg_hat_zq_hat) {
			
			return make_nlo_integrand(FLx_delta_integrand, FLx_xi_integrand, FLx_xip_integrand, FLx_xi_xip_integrand, xi, xip, x, z, factorization_scale_log, fragmentation_scale_log, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
		}
		static constexpr double F3_nlo_integrand(
			const double xi,
			const double xip, 
			const double x, 
			const double z, 
			const double factorization_scale_log,
			const double fragmentation_scale_log,
			const double xq_zq,
			const double xq_hat_zq, 
			const double xq_zq_hat,
			const double xq_hat_zq_hat, 
			const double xq_zg_hat,
			const double xg_hat_zq,
			const double xq_hat_zg_hat,
			const double xg_hat_zq_hat) {
			
			return make_nlo_integrand(F3_delta_integrand, F3_xi_integrand, F3_xip_integrand, F3_xi_xip_integrand, xi, xip, x, z, factorization_scale_log, fragmentation_scale_log, xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat);
		}
	}
}

#endif