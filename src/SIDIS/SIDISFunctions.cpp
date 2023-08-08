#ifndef SIDIS_FUNCTIONS_H
#define SIDIS_FUNCTIONS_H

#include "Common/Constants.cpp"
#include <cmath>
#include "Common/Flavor.cpp"
#include "Decay/Decay.cpp"
#include "PDF/FragmentationConfiguration.cpp"
#include "Common/TRFKinematics.cpp"
#include "Decay/DecayFunctions.cpp"
#include "Common/CKM.cpp"
#include "PDF/PDFConcept.cpp"
#include <ranges>
#include "Utility/Math.cpp"

#include "SIDIS/Coefficients/F2.cpp"
#include "SIDIS/Coefficients/FL.cpp"
#include "SIDIS/Coefficients/F3.cpp"

namespace SIDISFunctions {
	template <PDFConcept PDFInterface, PDFConcept FFInterface, DecayFunctions::Concept DecayFunction>
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
	};
}

namespace SIDISFunctions {
	template <PDFConcept PDFInterface, PDFConcept FFInterface, DecayFunctions::Concept DecayFunction, typename Signature>
	constexpr static double construct(
		const double x, 
		const double z, 
		const double xi, 
		const double xip, 
		const Parameters<PDFInterface, FFInterface, DecayFunction> &params, 
		const Signature integrand, 
		const int sign, 
		const FFInterface &ff1, 
		const FFInterface &ff2, 
		const FlavorType flavor1, 
		const FlavorType flavor2, 
		const FlavorType antiflavor1, 
		const FlavorType antiflavor2, 
		const double log1mx, 
		const double log1mz, 
		const double logxi, 
		const double logxip, 
		const double log1mxi,
		const double log1mxip,
		const double m2,
		const double Q2) {

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

		const double value = integrand(
								xi, xip, x, z, 
								params.factorization_scale_log, params.fragmentation_scale_log, 
								xq_zq, xq_hat_zq, xq_zq_hat, xq_hat_zq_hat, xq_zg_hat, 
								xg_hat_zq, xq_hat_zg_hat, xg_hat_zq_hat, 
								log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
								m2, Q2);

		return value;
	}
	template <PDFConcept PDFInterface, PDFConcept FFInterface, DecayFunctions::Concept DecayFunction, typename Signature>
	constexpr static double construct(
		const double input[], 
		const Parameters<PDFInterface, FFInterface, DecayFunction> &params, 
		const Signature signature, 
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

		const FlavorVector &flavors1 = positive_W ? flavors.lower_flavors : flavors.upper_flavors;
		const FlavorVector &flavors2 = positive_W ? flavors.upper_flavors : flavors.lower_flavors;

		const PDFInterface &pdf1 = params.pdf1;
		const PDFInterface &pdf2 = params.pdf2;

		const double log1mz = Math::log1m(z);
		
		const double logxi = std::log(xi);
		const double logxip = std::log(xip);

		const double log1mxi = Math::log1m(xi);
		const double log1mxip = Math::log1m(xip);

		double sum = 0.0;
		for (const FlavorType outgoing : flavors2) {
			const FlavorType anti_incoming = Flavor::conjugate_flavor(outgoing);
			for (const FlavorType incoming : flavors1) {
				const FlavorType anti_outgoing = Flavor::conjugate_flavor(incoming);

				const double mass = flavors.mass(Flavor::Charm);
				const double x_mass = CommonFunctions::compute_momentum_fraction_mass_correction(x, Q2, mass, 0.0);

				if (xi_int && xi < x_mass) { continue; }
				pdf1.evaluate(x_mass, factorization_scale);
				pdf2.evaluate(x_mass / xi, factorization_scale);

				const double log1mx = Math::log1m(x_mass);

				const double V_ckm = CKM::squared(incoming, outgoing);
				const double total_value = construct(
											x_mass, z, xi, xip, 
											params, signature, sign, ff1, ff2, 
											incoming, outgoing, anti_incoming, anti_outgoing, 
											log1mx, log1mz, logxi, logxip, log1mxi, log1mxip,
											std::pow(mass, 2), Q2);
				const double summand = V_ckm * total_value;

				sum += summand;
			}
		}
		const double decay_value = decay(x, z, Q2, z_min);
		sum *= decay_value;

		return sum;
	}

	template <PDFConcept PDFInterface, PDFConcept FFInterface, DecayFunctions::Concept DecayFunction = decltype(DecayFunctions::trivial), typename Signature>
	constexpr static double construct(
		const double input[], void *params_in, const Signature integrand, 
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
			const double value = construct<PDFInterface, FFInterface, DecayFunction>(
									input, params, integrand, 
									xi_int, xip_int, z_int, sign, 
									ff1, ff2, decay, z_min);
			const double summand = value;

			sum += summand;
		}

		return sum;
	}

	template <PDFConcept PDFInterface, PDFConcept FFInterface, DecayFunctions::Concept DecayFunction, typename Signature>
	constexpr static double cross_section(
		const double input[], void *params_in, 
		const Signature F2, const Signature FL, const Signature F3, 
		const bool xi_int, const bool xip_int, const bool z_int) {
		const auto &params = *static_cast<Parameters<PDFInterface, FFInterface, DecayFunction> *>(params_in);
		
		const double f2 = construct<PDFInterface, FFInterface, DecayFunction>(input, params_in, F2, xi_int, xip_int, z_int, 1);
		const double fL = construct<PDFInterface, FFInterface, DecayFunction>(input, params_in, FL, xi_int, xip_int, z_int, 1);
		const double f3 = construct<PDFInterface, FFInterface, DecayFunction>(input, params_in, F3, xi_int, xip_int, z_int, -1);

		const double cs = CommonFunctions::make_cross_section_variable(params.kinematics, params.process, f2, fL, f3);

		return cs;
	}
}

#endif