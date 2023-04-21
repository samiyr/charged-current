#ifndef PDFFF_H
#define PDFFF_H

#include "PDFInterface.cpp"
#include "FFInterface.cpp"
#include <string>
#include "Process.cpp"
#include "Flavor.cpp"

namespace PDF_FF {
	static double xq_zq_sum(const PDFInterface &pdf, const FFInterface &ff, const FlavorVector &upper_flavors, 
	const FlavorVector &lower_flavors, 
	const Process process) {
		const bool positive_W = process.positive_W();

		const FlavorVector &flavors1 = positive_W ? lower_flavors : upper_flavors;
		const FlavorVector &flavors2 = positive_W ? upper_flavors : lower_flavors;

		double sum = 0.0;
		for (size_t i = 0; i < flavors1.size(); i++) {
			for (size_t j = 0; j < flavors2.size(); j++) {
				const FlavorType flavor1 = flavors1[i];
				const FlavorType flavor2 = flavors2[j];

				const double V_ckm = CKM::squared(flavor1, flavor2);
				const double xq = pdf.xf(flavor1);
				const double zD = ff.xf(flavor2);

				sum += V_ckm * xq * zD;
			}
		}
		return sum;
	}

	static double xq_zg_sum(const PDFInterface &pdf, const FFInterface &ff, const FlavorVector &upper_flavors, 
	const FlavorVector &lower_flavors, 
	const Process process) {
		const bool positive_W = process.positive_W();

		const double xq_sum = positive_W ? pdf.xf_sum(lower_flavors, upper_flavors) : pdf.xf_sum(upper_flavors, lower_flavors);
		const double zDg = ff.xg();

		return xq_sum * zDg;
	}

	static double xg_zq_sum(const PDFInterface &pdf, const FFInterface &ff, const FlavorVector &upper_flavors, 
	const FlavorVector &lower_flavors, 
	const Process process) {
		const bool positive_W = process.positive_W();

		const double zDq_sum = positive_W ? ff.xf_sum(lower_flavors, upper_flavors) : ff.xf_sum(upper_flavors, lower_flavors);
		const double xg = pdf.xg();

		return xg * zDq_sum;
	}
}

#endif