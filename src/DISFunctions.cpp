#ifndef DIS_FUNCTIONS_H
#define DIS_FUNCTIONS_H

#include "Constants.cpp"
#include <cmath>
#include "Flavor.cpp"
#include <functional>
#include "TRFKinematics.cpp"
#include "CommonFunctions.cpp"
#include "CKM.cpp"
#include "PDFConcept.cpp"

namespace DISFunctions {
	template <PDFConcept PDFInterface>
	struct Parameters {
		const PDFInterface &pdf1;
		const PDFInterface &pdf2;

		const FlavorInfo &flavors;

		const double nlo_coefficient;

		const Process &process;
		const TRFKinematics &kinematics;

		const double factorization_scale;
		const double factorization_scale_log;
	};

	static constexpr double delta_contribution(const double x) {
		const double log = std::log(1 - x);
		return -2.0 * Constants::C_F * (9.0 / 2.0 + (M_PI * M_PI) / 3.0 - log * log + 1.5 * log) / x;
	}

	namespace Evaluation {
		template <PDFConcept PDFInterface, typename Signature>
		constexpr static double construct(const double x, const double xi, const Parameters<PDFInterface> &params, const Signature integrand, const int sign, const FlavorType flavor1, const FlavorType antiflavor1) {			
			const PDFInterface &pdf1 = params.pdf1;
			const PDFInterface &pdf2 = params.pdf2;

			const double xg_hat = pdf2.xg();

			const double xq = pdf1.xf(flavor1);
			const double xq_hat = pdf2.xf(flavor1);

			const double anti_xq = pdf1.xf(antiflavor1);
			const double anti_xq_hat = pdf2.xf(antiflavor1);

			const double xq_total = xq + sign * anti_xq;
			const double xq_hat_total = xq_hat + sign * anti_xq_hat;

			const double value = integrand(xi, x, params.factorization_scale_log, xq_total, xq_hat_total, 2.0 * xg_hat);

			return value;
		}
		template <PDFConcept PDFInterface, typename Signature>
		constexpr static double construct(const double input[], const Parameters<PDFInterface> &params, const Signature signature, const bool xi_int, const int sign) {
			const double xi = xi_int ? input[0] : 1.0;

			const TRFKinematics &kinematics = params.kinematics;

			const double x = kinematics.x;
			const double Q2 = kinematics.Q2;

			const double factorization_scale = params.factorization_scale;

			if (xi_int) {
				// if (xi < x) { return 0; }
				if (xi == 1.0) { return 0.0; }
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
					const double x_mass = CommonFunctions::compute_momentum_fraction_mass_correction(x, Q2, flavors.mass(Flavor::Charm), 0.0);

					if (xi_int && xi < x_mass) { continue; }
					pdf1.evaluate(x_mass, factorization_scale);
					pdf2.evaluate(x_mass / xi, factorization_scale);

					const double V_ckm = CKM::squared(incoming, outgoing);
					const double total_value = construct(x_mass, xi, params, signature, sign, incoming, anti_incoming);
					const double summand = V_ckm * total_value;

					sum += summand;
				}
			}

			return sum;
		}

		template <PDFConcept PDFInterface, typename Signature>
		constexpr static double construct(const double input[], void *params_in, const Signature integrand, const bool xi_int, const int sign) {
			const Parameters<PDFInterface> &params = *static_cast<Parameters<PDFInterface> *>(params_in);
			const TRFKinematics &kinematics = params.kinematics;

			if (xi_int) {
				const double xi = input[0];
				const double x = kinematics.x;
				const double x_hat = x / xi;
				params.pdf2.evaluate(x_hat, params.factorization_scale);
			}

			const double value = construct<PDFInterface>(input, params, integrand, xi_int, sign);

			return value;
		}

		template <PDFConcept PDFInterface, typename Signature>
		constexpr static double cross_section(const double input[], void *params_in, const Signature F2, const Signature FL, const Signature xF3, const bool xi_int) {
			const struct Parameters<PDFInterface> &params = *static_cast<Parameters<PDFInterface> *>(params_in);
			
			const double f2 = construct<PDFInterface>(input, params_in, F2, xi_int, 1);
			const double fL = construct<PDFInterface>(input, params_in, FL, xi_int, 1);
			const double xf3 = construct<PDFInterface>(input, params_in, xF3, xi_int, -1);

			const double cs = CommonFunctions::make_cross_section_variable(params.kinematics.x, params.kinematics.Q2, params.kinematics.s, params.process, f2, fL, xf3);

			return cs;
		}
	}

	namespace Integrands {
		static constexpr double F2x_lo_integrand([[maybe_unused]] const double xi, const double x, [[maybe_unused]] const double factorization_scale_log, const double xq, [[maybe_unused]] const double xq_hat, [[maybe_unused]] const double xg_hat) {
			return 2 * xq / x;
		}
		static constexpr double F2x_delta_integrand([[maybe_unused]] const double xi, const double x, const double factorization_scale_log, const double xq, [[maybe_unused]] const double xq_hat, [[maybe_unused]] const double xg_hat) {
			const double term1 = DISFunctions::delta_contribution(x) * xq;
			return term1;
		}
		static constexpr double F2x_nlo_integrand(const double xi, const double x, const double factorization_scale_log, const double xq, const double xq_hat, const double xg_hat) {
			const double term1 = std::log(1 - xi) * ((1 + xi * xi) * xq_hat - 2 * xq) / (1 - xi);
			const double term2 = (xq_hat - xq) / (1 - xi);
			const double term3 = xq_hat * (- (1 + xi * xi) * std::log(xi) / (1 - xi) + 3 + 2 * xi);

			const double quark_contribution = Constants::C_F * (term1 - 1.5 * term2 + term3);

			const double term4 = (xi * xi + (1 - xi) * (1 - xi)) * std::log((1 - xi) / xi);
			const double term5 = -1 + 8 * xi * (1 - xi);

			const double gluon_contribution = 0.5 * xg_hat * (term4 + term5);

			const double total_contribution = 2 * (quark_contribution + gluon_contribution) / x;
			return total_contribution;
		}

		static constexpr double FLx_lo_integrand([[maybe_unused]] const double xi, [[maybe_unused]] const double x, [[maybe_unused]] const double factorization_scale_log, [[maybe_unused]] const double xq, [[maybe_unused]] const double xq_hat, [[maybe_unused]] const double xg_hat) {
			return 0.0;
		}
		static constexpr double FLx_delta_integrand([[maybe_unused]] const double xi, [[maybe_unused]] const double x, [[maybe_unused]] const double factorization_scale_log, [[maybe_unused]] const double xq, [[maybe_unused]] const double xq_hat, [[maybe_unused]] const double xg_hat) {
			return 0.0;
		}
		static constexpr double FLx_nlo_integrand(const double xi, const double x, [[maybe_unused]] const double factorization_scale_log, [[maybe_unused]] const double xq, const double xq_hat, const double xg_hat) {
			const double quark_contribution = Constants::C_F * xq_hat * 2 * xi;
			const double gluon_contribution = 2 * xi * (1 - xi) * xg_hat;

			const double total_contribution = 2 * (quark_contribution + gluon_contribution) / x;
			return total_contribution;
		}
		static constexpr double F3_lo_integrand([[maybe_unused]] const double xi, const double x, [[maybe_unused]] const double factorization_scale_log, const double xq, [[maybe_unused]] const double xq_hat, [[maybe_unused]] const double xg_hat) {
			return 2 * xq / x;
		}
		static constexpr double F3x_delta_integrand(const double xi, const double x, const double factorization_scale_log, const double xq, const double xq_hat, const double xg_hat) {
			return F2x_delta_integrand(xi, x, factorization_scale_log, xq, xq_hat, xg_hat);
		}
		static constexpr double F3_nlo_integrand(const double xi, const double x, const double factorization_scale_log, const double xq, const double xq_hat, [[maybe_unused]] const double xg_hat) {
			const double term1 = std::log(1 - xi) * ((1 + xi * xi) * xq_hat - 2 * xq) / (1 - xi);
			const double term2 = (xq_hat - xq) / (1 - xi);
			const double term3 = xq_hat * (- (1 + xi * xi) * std::log(xi) / (1 - xi) + 2 + xi);

			const double quark_contribution = 2 * Constants::C_F * (term1 - 1.5 * term2 + term3) / x;

			return quark_contribution;
		}
	}
}

#endif