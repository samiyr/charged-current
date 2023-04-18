#ifndef PDF_EVALUATION_H
#define PDF_EVALUATION_H

#include <string>
#include <stdexcept>
#include "LHAPDF/LHAPDF.h"
#include "Flavor.cpp"
#include "CKM.cpp"

using namespace LHAPDF;

struct PDFEvaluation {
	std::string pdf_set_name;
	int set_member_number;

	std::vector<FlavorType> available_flavors;

	PDFEvaluation(std::string _pdf_set_name, int _set_member_number)
	: pdf_set_name(_pdf_set_name), set_member_number(_set_member_number) {
		_pdf = mkPDF(pdf_set_name, set_member_number);
		available_flavors = _pdf->flavors();
	}

	void evaluate(const double x, const double Q2) {
		_pdf->xfxQ2(x, Q2, flavor_values);
	}
	double xf_evaluate(const FlavorType flavor, const double x, const double Q2) {
		evaluate(x, Q2);
		return xf(flavor);
	}
	double xf(const FlavorType flavor) {
		return flavor_values[flavor];
	}
	double alpha_s(const double Q2) {
		return _pdf->alphasQ2(Q2);
	}
	double xf_sum(std::vector<FlavorType> flavors) {
		double sum = 0;
		for (const auto &flavor : flavors) {
			const double ckm_contribution = CKM::sum_of_squares(reflect_flavors(flavors), flavor);
			sum += ckm_contribution * xf(flavor);
		}
		return sum;
	}
	double xg() {
		return xf(Flavor::Gluon);
	}
	private:
	PDF *_pdf;
	std::map<FlavorType, double> flavor_values;
};

#endif