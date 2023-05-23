#ifndef PDF_COMMON_H
#define PDF_COMMON_H

#include <string>
#include "Process.cpp"
#include "Flavor.cpp"
#include "CKM.cpp"

#define CHARGED_CURRENT 1

#if CHARGED_CURRENT == 1

namespace PDFCommon {
	template <typename PDFInterface>
	constexpr double xf_sum(const PDFInterface &pdf, const FlavorVector &flavors, const FlavorVector &reflected) {
		double sum = 0;
		for (const auto flavor : flavors) {
			const double ckm_contribution = CKM::sum_of_squares(reflected, flavor);
			sum += ckm_contribution * pdf.xf(flavor);
		}
		return sum;
	}
	template <typename PDFInterface>
	constexpr double xg_sum(const PDFInterface &pdf, const FlavorInfo &flavors) {
		double ckm_sum = 0;
		for (const auto flavor : flavors.upper_flavors) {
			const double ckm_contribution = CKM::sum_of_squares(flavors.lower_flavors, flavor);
			ckm_sum += ckm_contribution;
		}
		return 2 * ckm_sum * pdf.xg();
	}
	template <typename PDFInterface>
	constexpr double xq_sum(const PDFInterface &pdf,
	const FlavorInfo &flavors, 
	const bool quark_minus, 
	const Process process,
	const bool sum_first = true) {
		const bool positive_W = process.positive_W();
		
		// This condition follows by drawing a truth table: lower flavors are summed first iff positive_W and sum_first are both either true of false
		const bool lower_first = positive_W == sum_first;

		const double xq_sum = lower_first ? xf_sum(pdf, flavors.lower_flavors, flavors.upper_flavors) : xf_sum(pdf, flavors.upper_flavors, flavors.lower_flavors);
		const double anti_xq_sum = lower_first ? xf_sum(pdf, flavors.upper_antiflavors, flavors.lower_antiflavors) : xf_sum(pdf, flavors.lower_antiflavors, flavors.upper_antiflavors);
		const double xq = quark_minus ? xq_sum - anti_xq_sum : xq_sum + anti_xq_sum;
		return xq;
	}

	template <typename PDFInterface, typename FFInterface>
	constexpr static double xq_zq_sum(const PDFInterface &pdf, 
	const FFInterface &ff, 
	const FlavorInfo &flavors, 
	const bool quark_minus, 
	const Process process) {
		const bool positive_W = process.positive_W();

		const FlavorVector &flavors1 = positive_W ? flavors.lower_flavors : flavors.upper_flavors;
		const FlavorVector &flavors2 = positive_W ? flavors.upper_flavors : flavors.lower_flavors;

		double sum = 0.0;
		for (const FlavorType flavor1 : flavors1) {
			const FlavorType antiflavor2 = Flavor::conjugate_flavor(flavor1);
			for (const FlavorType flavor2 : flavors2) {
				const FlavorType antiflavor1 = Flavor::conjugate_flavor(flavor2);

				const double V_ckm = CKM::squared(flavor1, flavor2);

				const double xq = pdf.xf(flavor1);
				const double zD = ff.xf(flavor2);

				const double anti_xq = pdf.xf(antiflavor1);
				const double anti_zD = ff.xf(antiflavor2);

				const double contribution = quark_minus ? xq * zD - anti_xq * anti_zD : xq * zD + anti_xq * anti_zD;
				// std::cout << xq << std::endl;
				// std::cout << zD << std::endl;
				// std::cout << anti_xq << std::endl;
				// std::cout << anti_zD << std::endl << std::endl;
				sum += V_ckm * contribution;
			}
		}
		// std::cout << sum << std::endl;
		return sum;
	}

	template <typename PDFInterface, typename FFInterface>
	constexpr static double xq_zg_sum(const PDFInterface &pdf, 
	const FFInterface &ff, 
	const FlavorInfo &flavors, 
	const bool quark_minus, 
	const Process process) {
		const double xq_sum = PDFCommon::xq_sum(pdf, flavors, quark_minus, process, true);
		const double zg_sum = ff.xg();
		const double sum = xq_sum * zg_sum;

		return sum;
	}

	template <typename PDFInterface, typename FFInterface>
	constexpr static double xg_zq_sum(const PDFInterface &pdf, 
	const FFInterface &ff, 
	const FlavorInfo &flavors, 
	const bool quark_minus, 
	const Process process) {
		const double xg_sum = pdf.xg();
		const double zq_sum = PDFCommon::xq_sum(ff, flavors, quark_minus, process, false);
		const double sum = xg_sum * zq_sum;

		return sum;
	}
}

#else

namespace PDFCommon {
	template <typename PDFInterface>
	double xf_sum(const PDFInterface &pdf, const FlavorVector &flavors, const FlavorVector &reflected) {
		double sum = 0;
		for (const auto flavor : flavors) {
			sum += pdf.xf(flavor);
		}
		return sum;
	}
	template <typename PDFInterface>
	double xg_sum(const PDFInterface &pdf, const FlavorInfo &flavors) {
		double ckm_sum = 0;
		for (const auto flavor : flavors.active_flavors) {
			ckm_sum += 1.0;
		}
		return 2 * ckm_sum * pdf.xg();
	}
	template <typename PDFInterface>
	double xq_sum(const PDFInterface &pdf,
	const FlavorInfo &flavors, 
	const bool quark_minus, 
	const Process process,
	const bool sum_first = true) {
		const double xq_sum = xf_sum(pdf, flavors.active_flavors, flavors.active_flavors);
		const double anti_xq_sum = xf_sum(pdf, flavors.active_antiflavors, flavors.active_antiflavors);
		const double xq = quark_minus ? xq_sum - anti_xq_sum : xq_sum + anti_xq_sum;
		return xq;
	}

	template <typename PDFInterface, typename FFInterface>
	static double xq_zq_sum(const PDFInterface &pdf, 
	const FFInterface &ff, 
	const FlavorInfo &flavors, 
	const bool quark_minus, 
	const Process process) {
		const FlavorVector &flavors1 = flavors.active_flavors;
		const FlavorVector &flavors2 = flavors.active_flavors;

		double sum = 0.0;
		for (const FlavorType flavor1 : flavors1) {
			const FlavorType antiflavor2 = Flavor::conjugate_flavor(flavor1);
			for (const FlavorType flavor2 : flavors2) {
				const FlavorType antiflavor1 = Flavor::conjugate_flavor(flavor2);

				const double V_ckm = 1.0;

				const double xq = pdf.xf(flavor1);
				const double zD = ff.xf(flavor2);

				const double anti_xq = pdf.xf(antiflavor1);
				const double anti_zD = ff.xf(antiflavor2);

				const double contribution = quark_minus ? xq * zD - anti_xq * anti_zD : xq * zD + anti_xq * anti_zD;
				// std::cout << xq << std::endl;
				// std::cout << zD << std::endl;
				// std::cout << anti_xq << std::endl;
				// std::cout << anti_zD << std::endl << std::endl;
				sum += V_ckm * contribution;
			}
		}
		// std::cout << sum << std::endl;
		return sum;
	}

	template <typename PDFInterface, typename FFInterface>
	static double xq_zg_sum(const PDFInterface &pdf, 
	const FFInterface &ff, 
	const FlavorInfo &flavors, 
	const bool quark_minus, 
	const Process process) {
		const double xq_sum = PDFCommon::xq_sum(pdf, flavors, quark_minus, process, true);
		const double zg_sum = ff.xg();
		const double sum = xq_sum * zg_sum;

		return sum;
	}

	template <typename PDFInterface, typename FFInterface>
	static double xg_zq_sum(const PDFInterface &pdf, 
	const FFInterface &ff, 
	const FlavorInfo &flavors, 
	const bool quark_minus, 
	const Process process) {
		const double xg_sum = pdf.xg();
		const double zq_sum = PDFCommon::xq_sum(ff, flavors, quark_minus, process, false);
		const double sum = xg_sum * zq_sum;

		return sum;
	}
}

#endif

#endif