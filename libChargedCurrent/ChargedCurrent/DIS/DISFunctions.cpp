#ifndef DIS_FUNCTIONS_H
#define DIS_FUNCTIONS_H

#include <cmath>
#include <functional>

#include "Common/Constants.cpp"
#include "Common/Flavor.cpp"
#include "Common/TRFKinematics.cpp"
#include "Common/CommonFunctions.cpp"
#include "Common/CKM.cpp"

#include "PDF/PDFConcept.cpp"

#include "DIS/Coefficients/F2.cpp"
#include "DIS/Coefficients/FL.cpp"
#include "DIS/Coefficients/F3.cpp"

namespace DISFunctions {
	template <is_pdf_interface PDFInterface>
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

	template <is_pdf_interface PDFInterface, typename Signature>
	constexpr static double evaluate_integrand(
		const double x, const double xi, 
		const Parameters<PDFInterface> &params, const Signature integrand, const int sign, 
		const FlavorType flavor1, const FlavorType antiflavor1,
		const double m2,
		const double Q2) {		

		const PDFInterface &pdf1 = params.pdf1;
		const PDFInterface &pdf2 = params.pdf2;

		const double xg_hat = pdf2.xg();

		const double xq = pdf1.xf(flavor1);
		const double xq_hat = pdf2.xf(flavor1);

		const double anti_xq = pdf1.xf(antiflavor1);
		const double anti_xq_hat = pdf2.xf(antiflavor1);

		const double xq_total = xq + sign * anti_xq;
		const double xq_hat_total = xq_hat + sign * anti_xq_hat;

		return integrand(
			xi, x, 
			params.factorization_scale_log, 
			xq_total, xq_hat_total, 2.0 * xg_hat, 
			m2, Q2
		);
	}
	template <is_pdf_interface PDFInterface, typename Signature>
	constexpr static double evaluate_integral(
		const double input[], const Parameters<PDFInterface> &params, const Signature signature, 
		const bool xi_int, const int sign) {

		const double xi = xi_int ? input[0] : 1.0;

		const TRFKinematics &kinematics = params.kinematics;

		const double x = kinematics.x;
		const double Q2 = kinematics.Q2;

		const double factorization_scale = params.factorization_scale;

		if (xi_int) {
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
				const double mass = flavors.mass(outgoing);
				const double x_mass = CommonFunctions::compute_momentum_fraction_mass_correction(x, Q2, mass, 0.0);

				if (xi_int && xi < x_mass) { continue; }
				if (x_mass >= 1.0) { continue; }
				pdf1.evaluate(x_mass, factorization_scale);
				pdf2.evaluate(x_mass / xi, factorization_scale);

				const double V_ckm = CKM::squared(incoming, outgoing);
				const double total_value = evaluate_integrand(
					x_mass, xi,
					params, signature, sign, 
					incoming, anti_incoming,
					std::pow(mass, 2), Q2
				);
				const double summand = V_ckm * total_value;

				sum += summand;
			}
		}

		return sum;
	}

	template <is_pdf_interface PDFInterface, typename Signature>
	constexpr static double evaluate(const double input[], void *params_in, const Signature integrand, const bool xi_int, const int sign) {
		const Parameters<PDFInterface> &params = *static_cast<Parameters<PDFInterface> *>(params_in);
		const TRFKinematics &kinematics = params.kinematics;

		if (xi_int) {
			const double xi = input[0];
			const double x = kinematics.x;
			const double x_hat = x / xi;
			params.pdf2.evaluate(x_hat, params.factorization_scale);
		}

		const double value = evaluate_integral<PDFInterface>(input, params, integrand, xi_int, sign);

		return value;
	}

	template <is_pdf_interface PDFInterface, typename Signature>
	constexpr static double cross_section(
		const double input[], void *params_in,
		const Signature F2, const Signature FL, const Signature F3, 
		const bool xi_int) {

		const struct Parameters<PDFInterface> &params = *static_cast<Parameters<PDFInterface> *>(params_in);
		
		const double f2 = evaluate<PDFInterface>(input, params_in, F2, xi_int, 1);
		const double fL = evaluate<PDFInterface>(input, params_in, FL, xi_int, 1);
		const double f3 = evaluate<PDFInterface>(input, params_in, F3, xi_int, -1);

		const double cs = CommonFunctions::make_cross_section_variable(params.kinematics, params.process, f2, fL, f3);

		return cs;
	}
}

#endif