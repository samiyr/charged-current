#ifndef SIDIS_FUNCTIONS_H
#define SIDIS_FUNCTIONS_H

#include <cmath>
#include <ranges>

#include "Common/Constants.cpp"
#include "Common/Flavor.cpp"
#include "Common/TRFKinematics.cpp"
#include "Common/CKM.cpp"
#include "Common/CommonFunctions.cpp"

#include "Decay/Decay.cpp"
#include "Decay/DecayFunctions.cpp"

#include "PDF/FragmentationConfiguration.cpp"
#include "PDF/PDFConcept.cpp"

#include "Utility/Math.cpp"

#include "SIDIS/Coefficients/F2.cpp"
#include "SIDIS/Coefficients/FL.cpp"
#include "SIDIS/Coefficients/F3.cpp"
#include "SIDIS/Coefficients/EvaluationParameters.cpp"

namespace SIDISFunctions {
	template <is_pdf_interface PDFInterface, is_pdf_interface FFInterface, is_decay_function DecayFunction>
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

		const double renormalization_scale;
		const double factorization_scale;
		const double fragmentation_scale;

		const double renormalization_scale_log;
		const double factorization_scale_log;
		const double fragmentation_scale_log;
	};
}

namespace SIDISFunctions {
	template <is_pdf_interface PDFInterface, is_pdf_interface FFInterface, is_decay_function DecayFunction, typename Signature>
	constexpr static double evaluate_integrand(
		const double x, 
		const double z, 
		const double xi, 
		const double xip, 
		const Parameters<PDFInterface, FFInterface, DecayFunction> &params, 
		const Signature integrand, 
		const FFInterface &ff1, 
		const FFInterface &ff2, 
		const FlavorType incoming, 
		const FlavorType outgoing, 
		const double m2,
		const double Q2) {

		const PDFInterface &pdf1 = params.pdf1;
		const PDFInterface &pdf2 = params.pdf2;

		const double pdf = pdf1.xf(incoming);
		const double pdf_hat = pdf2.xf(incoming);

		const double ff = ff1.xf(outgoing);
		const double ff_hat = ff2.xf(outgoing);

		return integrand({
			.xi = xi, .xip = xip, .x = x, .z = z,
			.pdf = pdf, .pdf_hat = pdf_hat,
			.ff = ff, .ff_hat = ff_hat,
			.renormalization_scale_log = params.renormalization_scale_log,
			.factorization_scale_log = params.factorization_scale_log, .fragmentation_scale_log = params.fragmentation_scale_log,
			.m2 = m2, .Q2 = Q2
		});
	}

	template <is_pdf_interface PDFInterface, is_pdf_interface FFInterface, is_decay_function DecayFunction>
	constexpr static double evaluate_integral_with_decay(
		const double input[], 
		const Parameters<PDFInterface, FFInterface, DecayFunction> &params, 
		const auto &&quark_to_quark, const auto &&quark_to_gluon, const auto &&gluon_to_quark,
		const bool xi_int, const bool xip_int, const bool z_int, const int sign, 
		const FFInterface &ff1, const FFInterface &ff2, 
		const Decay<DecayFunction> &decay, const double z_min) {

		const std::size_t xi_index = 0;
		const std::size_t xip_index = std::size_t(xi_int);
		const std::size_t z_index = std::size_t(xi_int) + std::size_t(xip_int);

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
			if (xi == 1.0) [[unlikely]] { return 0.0; }
		}
		if (z_int) {
			if (xip < z) { return 0.0; }
			ff1.evaluate(z, fragmentation_scale);
		}
		if (xip_int) {
			// if (xip < z) {
			// 	if (!z_int) {
			// 		std::cout << xip << IO::endl; 
			// 	}
			// 	return 0; 
			// 	}
			if (xip == 1.0) [[unlikely]] { return 0.0; }
			ff2.evaluate(z_hat, fragmentation_scale);
		}

		const FlavorInfo &flavors = params.flavors;

		const Process &process = params.process;
		const bool positive_W = process.positive_W();

		const FlavorVector &incoming_quark_flavors = positive_W ? flavors.lower_flavors : flavors.upper_flavors;
		const FlavorVector &outgoing_quark_flavors = positive_W ? flavors.upper_flavors : flavors.lower_flavors;

		const FlavorVector &incoming_antiquark_flavors = positive_W ? flavors.upper_antiflavors : flavors.lower_antiflavors;
		const FlavorVector &outgoing_antiquark_flavors = positive_W ? flavors.lower_antiflavors : flavors.upper_antiflavors;

		const PDFInterface &pdf1 = params.pdf1;
		const PDFInterface &pdf2 = params.pdf2;

		const auto evaluate_channel = [&](const FlavorType incoming, const FlavorType outgoing, const auto &&integrand) {
			const double incoming_mass = flavors.mass(incoming);
			const double outgoing_mass = flavors.mass(outgoing);

			const double mass = std::max(incoming_mass, outgoing_mass);
			const double chi = CommonFunctions::compute_momentum_fraction_mass_correction(x, Q2, mass, 0.0);

			if (xi_int && xi < chi) { return 0.0; }

			pdf1.evaluate(chi, factorization_scale);
			pdf2.evaluate(chi / xi, factorization_scale);

			return evaluate_integrand(
				chi, z, xi, xip,
				params, integrand, ff1, ff2,
				incoming, outgoing, std::pow(mass, 2), Q2
			);
		};

		double sum = 0.0;

		for (const FlavorType incoming : incoming_quark_flavors) {
			for (const FlavorType outgoing: outgoing_quark_flavors) {
				const double qq = evaluate_channel(incoming, outgoing, quark_to_quark);
				const double gq = evaluate_channel(incoming, Flavor::Gluon, quark_to_gluon);
				const double qg = evaluate_channel(Flavor::Gluon, outgoing, gluon_to_quark);

				sum += CKM::squared(incoming, outgoing) * (qq + gq + qg);
			}
		}

		for (const FlavorType incoming : incoming_antiquark_flavors) {
			for (const FlavorType outgoing: outgoing_antiquark_flavors) {
				const double qq = evaluate_channel(incoming, outgoing, quark_to_quark);
				const double gq = evaluate_channel(incoming, Flavor::Gluon, quark_to_gluon);
				const double qg = evaluate_channel(Flavor::Gluon, outgoing, gluon_to_quark);

				sum += static_cast<double>(sign) * CKM::squared(incoming, outgoing) * (qq + gq + qg);
			}
		}

		const double decay_value = decay(x, z, Q2, z_min);
		sum *= decay_value;

		return sum;
	}

	template <is_pdf_interface PDFInterface, is_pdf_interface FFInterface, is_decay_function DecayFunction = decltype(DecayFunctions::trivial)>
	constexpr static double evaluate(
		const double input[], void *params_in, const auto &&quark_to_quark, const auto &&quark_to_gluon, const auto &&gluon_to_quark, 
		const bool xi_int, const bool xip_int, const bool z_int, const int sign) {

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

		for (const auto &[ff1, ff2, decay] : std::views::zip(ffs1.interfaces, ffs2.interfaces, ffs1.decays)) {
			const double z_min = SIDISFunctions::Helper::compute_z_min(kinematics, decay);
			const double value = evaluate_integral_with_decay<PDFInterface, FFInterface, DecayFunction>(
				input, params,
				quark_to_quark, quark_to_gluon, gluon_to_quark, 
				xi_int, xip_int, z_int, sign, 
				ff1, ff2, decay, z_min
			);
			const double summand = value;

			sum += summand;
		}

		return sum;
	}

	template <is_pdf_interface PDFInterface, is_pdf_interface FFInterface, is_decay_function DecayFunction>
	constexpr static double cross_section(
		const double input[], void *params_in, 
		const auto &&F2_qq, const auto &&F2_gq, const auto &&F2_qg,
		const auto &&FL_qq, const auto &&FL_gq, const auto &&FL_qg,
		const auto &&F3_qq, const auto &&F3_gq, const auto &&F3_qg,
		const bool xi_int, const bool xip_int, const bool z_int) {
		const auto &params = *static_cast<Parameters<PDFInterface, FFInterface, DecayFunction> *>(params_in);
		
		const double f2 = evaluate<PDFInterface, FFInterface, DecayFunction>(input, params_in, F2_qq, F2_gq, F2_qg, xi_int, xip_int, z_int, 1);
		const double fL = evaluate<PDFInterface, FFInterface, DecayFunction>(input, params_in, FL_qq, FL_gq, FL_qg, xi_int, xip_int, z_int, 1);
		const double f3 = evaluate<PDFInterface, FFInterface, DecayFunction>(input, params_in, F3_qq, F3_gq, F3_qg, xi_int, xip_int, z_int, -1);

		const double cs = CommonFunctions::make_cross_section_variable(params.kinematics, params.process, f2, fL, f3);

		return cs;
	}
}

#endif