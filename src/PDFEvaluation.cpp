#ifndef PDF_EVALUATION_H
#define PDF_EVALUATION_H

#include <string>
#include <stdexcept>
#include "LHAPDF/LHAPDF.h"
#include "Flavor.cpp"
#include "CKM.cpp"
#include "Process.cpp"

using namespace LHAPDF;

struct PDFEvaluation {
	std::string pdf_set_name;
	int set_member_number;

	std::vector<FlavorType> available_flavors;

	PDFEvaluation(std::string _pdf_set_name, int _set_member_number = 0)
	: pdf_set_name(_pdf_set_name), set_member_number(_set_member_number) {
		LHAPDF::Info &cfg = LHAPDF::getConfig();
		cfg.set_entry("Verbosity", 0);
		_pdf = mkPDF(pdf_set_name, set_member_number);
		available_flavors = _pdf->flavors();
	}

	// PDFEvaluation(const PDFEvaluation &p) {
	// 	pdf_set_name = p.pdf_set_name;
	// 	set_member_number = p.set_member_number;
	// 	available_flavors = p.available_flavors;
	// 	_pdf = p._pdf;
	// 	flavor_values = p.flavor_values;
	// 	std::cout << p.pdf_set_name << std::endl;
	// }

	void evaluate(const double x, const double Q2) {
		_pdf->xfxQ2(x, Q2, flavor_values);
	}
	double xf_evaluate(const FlavorType flavor, const double x, const double Q2) {
		return _pdf->xfxQ2(flavor, x, Q2);
	}
	double xf(const FlavorType flavor) {
		return flavor_values[flavor + 6];
	}
	double alpha_s(const double Q2) {
		return _pdf->alphasQ2(Q2);
	}
	double xf_sum(const std::vector<FlavorType> &flavors, const std::vector<FlavorType> &reflected) {
		double sum = 0;
		for (const auto flavor : flavors) {
			const double ckm_contribution = CKM::sum_of_squares(reflected, flavor);
			sum += ckm_contribution * xf(flavor);
		}
		return sum;
	}
	double xg() {
		return xf(Flavor::Gluon);
	}
	double xg_sum(const std::vector<FlavorType> &flavors, const std::vector<FlavorType> &reflected) {
		double ckm_sum = 0;
		for (const auto flavor : flavors) {
			const double ckm_contribution = CKM::sum_of_squares(reflected, flavor);
			ckm_sum += ckm_contribution;
		}
		return 2 * ckm_sum * xg();
	}
	double xq_sum(const FlavorVector &upper_flavors, 
	const FlavorVector &lower_flavors, 
	const FlavorVector &upper_antiflavors, 
	const FlavorVector &lower_antiflavors, 
	const bool quark_minus, 
	const Process process) {
		const bool positive_W = process.positive_W();
		const double xq_sum = positive_W ? xf_sum(lower_flavors, upper_flavors) : xf_sum(upper_flavors, lower_flavors);
		const double anti_xq_sum = positive_W ? xf_sum(upper_antiflavors, lower_antiflavors) : xf_sum(lower_antiflavors, upper_antiflavors);
		const double xq = quark_minus ? xq_sum - anti_xq_sum : xq_sum + anti_xq_sum;
		return xq;
	}
	private:
	PDF *_pdf;
	std::vector<double> flavor_values;
};

#endif